library(shiny)
library(shinyWidgets)
library(leaflet)
library(sf)
library(tigris)
library(leaflet.extras)
library(geojsonio)
library(jsonlite)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(zip)
library(esri2sf)
library(stringr)
library(shinycssloaders)
library(raster)
library(gstat)
library(stars)
library(data.table)
library(lubridate)
library(httr)

# Load Texas counties data and transform to WGS84
texas_counties <- counties(state = "TX", cb = TRUE)
texas_counties <- st_transform(texas_counties, 4326)

# Load emission factors and county lookup table
emission_factors <- read.csv("data/Emission_Factors.csv")
county_lookup <- read.csv("data/counties.csv")
hourly_factors_df <- read.csv("data/Hourly_Factors.csv")

# 1. Load the geometries from the tigris package for the year 2020
options(tigris_use_cache = TRUE)
tx_counties_geo <- counties("TX", year = 2020)
tx_tracts_geo <- tracts("TX", year = 2020)
tx_blkgps_geo <- block_groups("TX", year = 2020)

# 2. Load the mortality rate CSV data
county_mr_data <- read.csv("data/County_mr.csv")
tract_mr_data <- read.csv("data/Tract_mr.csv")
blkgp_mr_data <- read.csv("data/Blkgp_MR.csv")
county_mr_data$FIPS<- as.character(county_mr_data$FIPS)
tract_mr_data$Census.Tract <- as.character(tract_mr_data$Census.Tract)
blkgp_mr_data$Census.Block.Group <- as.character(blkgp_mr_data$Census.Block.Group)
# 3. Join the CSV data to the geometries
# Ensure GEOID columns are characters for a stable join
county_mr_sf <- left_join(tx_counties_geo, county_mr_data, by = c("GEOID" = "FIPS"))
tract_mr_sf <- left_join(tx_tracts_geo, tract_mr_data, by = c("GEOID" = "Census.Tract"))
blkgp_mr_sf <- left_join(tx_blkgps_geo, blkgp_mr_data, by = c("GEOID" = "Census.Block.Group"))

# Rename columns for consistency and convert to numeric
# Note: Using names() is safer than dplyr::rename in case of missing columns
names(county_mr_sf)[names(county_mr_sf) == "Population"] <- "Total_Pop"
names(county_mr_sf)[names(county_mr_sf) == "Mortality.Rate..Deaths.100000."] <- "Mortality_Rate"
names(tract_mr_sf)[names(tract_mr_sf) == "Population"] <- "Total_Pop"
names(tract_mr_sf)[names(tract_mr_sf) == "Mortality.Rate..Deaths.100000."] <- "Mortality_Rate"
names(blkgp_mr_sf)[names(blkgp_mr_sf) == "Population"] <- "Total_Pop"
names(blkgp_mr_sf)[names(blkgp_mr_sf) == "Mortality.Rate..Deaths.100000."] <- "Mortality_Rate"

# Convert columns to numeric, handling potential errors
county_mr_sf$Total_Pop <- as.numeric(county_mr_sf$Total_Pop)
county_mr_sf$Mortality_Rate <- as.numeric(county_mr_sf$Mortality_Rate) / 100000 # Convert to a rate per person
tract_mr_sf$Total_Pop <- as.numeric(tract_mr_sf$Total_Pop)
tract_mr_sf$Mortality_Rate <- as.numeric(tract_mr_sf$Mortality_Rate) / 100000
blkgp_mr_sf$Total_Pop <- as.numeric(blkgp_mr_sf$Total_Pop)
blkgp_mr_sf$Mortality_Rate <- as.numeric(blkgp_mr_sf$Mortality_Rate) / 100000

# Define Relative Risk (RR) values
rr_values <- list(
  "PM2.5" = list(
    "1.06 (Pope et al., 2002)" = 1.06,
    "1.08 (COMEAP, 2018)" = 1.08
  ),
  "PM10" = list(
    "1.04 (Pope et al., 2002)" = 1.04
  ),
  "NOx (as NO2)" = list(
    "1.0046 (Faustini et al., 2014)" = 1.0046
  ),
  "CO" = list(
    "Not Recommended" = NA
  )
)



# Function to read AERMOD POSTFILE
read_pos <- function(infile) {
  lines <- readLines(infile)
  data_lines <- lines[grepl("^\\s{2}\\d", lines)]
  data_text <- paste(data_lines, collapse = "\n")
  pos_data <- fread(data_text, header = FALSE)
  col_names <- c("X", "Y", "AVERAGE_CONC", "ZELEV", "ZHILL", "ZFLAG", "AVE", "GRP", "DATE", "NET_ID")
  setnames(pos_data, colnames(pos_data)[1:ncol(pos_data)], col_names[1:ncol(pos_data)])
  return(pos_data)
}


# Function to generate near-road receptors
generate_receptors <- function(project_area) {
  project_area <- st_transform(project_area, 32614)  # UTM Zone 14N
  large_polygon <- st_union(project_area)
  outer_boundary_polygon <- large_polygon
  distances_m <- c(5, 50, 100, 200, 300, 500) * 3.28084 * 0.3048  # Feet to meters
  buffers <- lapply(distances_m, function(dist) st_buffer(outer_boundary_polygon, dist = dist))
  bbox <- st_bbox(buffers[[length(buffers)]])
  bbox_polygon <- st_as_sfc(bbox, crs = st_crs(project_area))
  spacings_m <- c(25, 50, 75, 100) * 3.28084 * 0.3048  # Feet to meters
  receptors_all <- list()
  for (i in seq_along(spacings_m)) {
    grid_points <- st_make_grid(bbox_polygon, cellsize = spacings_m[i], what = "centers")
    receptors <- st_sf(geometry = grid_points, crs = st_crs(project_area))
    if (i < length(distances_m)) {
      receptors <- receptors[!st_intersects(receptors, buffers[[i]], sparse = FALSE)[,1], ]
      if (i < length(buffers)) {
        receptors <- receptors[st_intersects(receptors, buffers[[i + 1]], sparse = FALSE)[,1], ]
      }
    }
    receptors_all[[i]] <- receptors
  }
  receptors_sf <- do.call(rbind, receptors_all)
  return(receptors_sf)
}

get_met_url <- function(county) {
  county_lower <- tolower(gsub(" ", "", county))
  if(county_lower %in% c("harris", "travis")){
    if(county_lower == "harris") {
      return(list(default = "https://www.tceq.texas.gov/downloads/permitting/air/modeling/air-dispersion/counties-h-k/harris-iah5y.zip",
                  alternative = "https://www.tceq.texas.gov/downloads/permitting/air/modeling/air-dispersion/counties-h-k/harris-iah5y.zip"))
    }
    if(county_lower == "travis") {
      return(list(default = "https://www.tceq.texas.gov/downloads/permitting/air/modeling/air-dispersion/counties-s-z/travis-aus5y.zip",
                  alternative = "https://www.tceq.texas.gov/downloads/permitting/air/modeling/air-dispersion/counties-s-z/travis-aus5y.zip"))
    }
  } else {
    first <- substr(county_lower, 1, 1)
    if(first >= "a" && first <= "c"){
      folder <- "counties-a-c"
    } else if(first >= "d" && first <= "g"){
      folder <- "counties-d-g"
    } else if(first >= "h" && first <= "k"){
      folder <- "counties-h-k"
    } else if(first >= "l" && first <= "r"){
      folder <- "counties-l-r"
    } else {
      folder = "counties-s-z"
    }
    url <- paste0("https://www.tceq.texas.gov/downloads/permitting/air/modeling/air-dispersion/", 
                  folder, "/", county_lower, "5y.zip")
    return(url)
  }
}

# Function to package AERMOD input and MET files into a ZIP
package_met_files <- function(aermod_input_text, met_dir, met_year, met_roughness, county_name, zip_filename = "AERMOD_Package.zip") {
  county_data <- county_lookup[county_lookup$County_Name == county_name, ]
  if (nrow(county_data) == 0) {
    county_data <- county_lookup[county_lookup$County_Name == "Bexar", ]
    showNotification("County not found in counties.csv, using Bexar as fallback.", type = "warning")
  }
  surface_id <- county_data$Surface_ID
  upper_air_id <- county_data$Upper_Air_ID
  county_upper <- toupper(gsub(" ", "", county_name))
  file_prefix <- paste0(county_upper, "_", surface_id, upper_air_id)
  sfc_file <- paste0(file_prefix, met_year, met_roughness, ".SFC")
  pfl_file <- paste0(file_prefix, met_year, met_roughness, ".PFL")
  
  met_files <- list.files(met_dir, full.names = TRUE, recursive = TRUE)
  chosen_files <- met_files[grepl(sfc_file, met_files, fixed = TRUE) | grepl(pfl_file, met_files, fixed = TRUE)]
  
  if (length(chosen_files) < 2) {
    pattern <- paste0(met_year, met_roughness, "\\.(SFC|PFL)$")
    chosen_files <- met_files[grepl(pattern, met_files)]
  }
  
  if (length(chosen_files) < 2) {
    stop("Could not find both .SFC and .PFL files for the chosen MET options.")
  }
  
  aermod_exe <- "data/aermod.exe"
  if (!file.exists(aermod_exe)) {
    stop("aermod.exe not found in www folder.")
  }
  
  temp_dir <- tempfile()
  dir.create(temp_dir)
  
  batch_file <- file.path(temp_dir, "run_aermod.bat")
  writeLines(c("cd %~dp0", "aermod.exe"), con = batch_file)
  
  aermod_input_file <- file.path(temp_dir, "AERMOD.INP")
  writeLines(aermod_input_text, con = aermod_input_file)
  
  file.copy(aermod_exe, temp_dir)
  file.copy(chosen_files, temp_dir)
  
  zip_file <- file.path(tempdir(), zip_filename)
  zipr(zipfile = zip_file, files = list.files(temp_dir, full.names = TRUE), root = temp_dir)
  
  unlink(temp_dir, recursive = TRUE)
  return(zip_file)
}

# Function to package scenario files
package_scenario_files <- function(baseline_input, scenario_input, scenario_title, met_dir, met_year, met_roughness, county_name) {
  main_dir <- tempfile()
  dir.create(main_dir)
  baseline_dir <- file.path(main_dir, "Baseline")
  dir.create(baseline_dir)
  scenario_dir_name <- gsub("[^A-Za-z0-9_]", "_", scenario_title) # Sanitize title for folder name
  scenario_dir <- file.path(main_dir, scenario_dir_name)
  dir.create(scenario_dir)
  
  # Function to copy files to a directory
  prepare_run_dir <- function(target_dir, aermod_input_text) {
    county_data <- county_lookup[county_lookup$County_Name == county_name, ]
    surface_id <- county_data$Surface_ID
    upper_air_id <- county_data$Upper_Air_ID
    county_upper <- toupper(gsub(" ", "", county_name))
    file_prefix <- paste0(county_upper, "_", surface_id, upper_air_id)
    sfc_file <- paste0(file_prefix, met_year, met_roughness, ".SFC")
    pfl_file <- paste0(file_prefix, met_year, met_roughness, ".PFL")
    
    met_files <- list.files(met_dir, full.names = TRUE, recursive = TRUE)
    chosen_files <- met_files[grepl(sfc_file, met_files, fixed = TRUE) | grepl(pfl_file, met_files, fixed = TRUE)]
    
    if (length(chosen_files) < 2) {
      pattern <- paste0(met_year, met_roughness, "\\.(SFC|PFL)$")
      chosen_files <- met_files[grepl(pattern, met_files)]
    }
    
    aermod_exe <- "data/aermod.exe"
    batch_file_path <- file.path(target_dir, "run_aermod.bat")
    writeLines(c("cd %~dp0", "aermod.exe"), con = batch_file_path)
    
    aermod_input_file_path <- file.path(target_dir, "AERMOD.INP")
    writeLines(aermod_input_text, con = aermod_input_file_path)
    
    file.copy(aermod_exe, target_dir)
    file.copy(chosen_files, target_dir)
  }
  
  # Prepare both directories
  prepare_run_dir(baseline_dir, baseline_input)
  prepare_run_dir(scenario_dir, scenario_input)
  
  # Zip the main directory
  zip_file <- file.path(tempdir(), "AERMOD_Scenarios_Package.zip")
  zip::zipr(zipfile = zip_file, files = list.files(main_dir, full.names = TRUE), recurse = TRUE, compression_level = 9)
  
  unlink(main_dir, recursive = TRUE)
  return(zip_file)
}

# Helper function to get population data from EPA EnviroAtlas
get_population_raster <- function(template_raster) {
  # Project extent of template raster to Web Mercator (EPSG:3857)
  template_proj <- projectExtent(template_raster, crs = "+init=epsg:3857")
  bbox <- extent(template_proj)
  
  # Define the image size based on 30m resolution
  width <- round((bbox@xmax - bbox@xmin) / 30)
  height <- round((bbox@ymax - bbox@ymin) / 30)
  
  # Construct the URL for the EPA EnviroAtlas ImageServer
  base_url <- "https://enviroatlas.epa.gov/arcgis/rest/services/Rasters/Dasymetric_2020/ImageServer/exportImage"
  query <- list(
    bbox = paste(bbox@xmin, bbox@ymin, bbox@xmax, bbox@ymax, sep = ","),
    bboxSR = 3857,
    size = paste(width, height, sep = ","),
    imageSR = 3857,
    format = "tiff",
    pixelType = "F32",
    noDataInterpretation = "esriNoDataMatchAny",
    interpolation = "RSP_BilinearInterpolation",
    f = "image"
  )
  
  request_url <- httr::modify_url(base_url, query = query)
  temp_tiff <- tempfile(fileext = ".tif")
  
  tryCatch({
    response <- httr::GET(request_url, httr::write_disk(temp_tiff, overwrite = TRUE), httr::timeout(300))
    if (httr::status_code(response) != 200) {
      stop("Failed to download population data. Server responded with status: ", httr::status_code(response))
    }
    population_raster <- raster(temp_tiff)
    return(population_raster)
  }, error = function(e) {
    showNotification(paste("Error fetching population data:", e$message), type = "error", duration = 15)
    return(NULL)
  })
}

# Helper function to calculate weighted exposure
calculate_weighted_exposure <- function(conc_raster, pop_raster) {
  pop_raster_proj <- projectRaster(pop_raster, crs = crs(conc_raster), method = 'bilinear')
  conc_resampled <- resample(conc_raster, pop_raster_proj, method = 'bilinear')
  pop_raster_proj <- crop(pop_raster_proj, extent(conc_resampled))
  conc_resampled[is.na(conc_resampled)] <- 0
  pop_raster_proj[is.na(pop_raster_proj) | pop_raster_proj < 0] <- 0
  
  product_raster <- conc_resampled * pop_raster_proj
  sum_product <- cellStats(product_raster, 'sum', na.rm = TRUE)
  sum_pop <- cellStats(pop_raster_proj, 'sum', na.rm = TRUE)
  
  weighted_avg <- if (sum_pop > 0) sum_product / sum_pop else 0
  
  return(list(
    weighted_avg = weighted_avg,
    conc_resampled = conc_resampled,
    pop_resampled = pop_raster_proj
  ))
}


# User Interface (UI)
ui <- fluidPage(
  tabsetPanel(
    id = "maintabs",
    type = "hidden",
    # Landing Page
    tabPanel(
      title = "Landing Page",
      fluidRow(
        column(
          12,
          tags$head(
            tags$style(HTML("
              .hover-highlight {
                transition: opacity 0.3s;
                opacity: 0; 
                position: absolute;
                background-color: rgba(255,255,0,0.4);
                z-index: 9999;
                pointer-events: none;
              }
              .total-impact-box {
                border: 2px solid #007bff;
                border-radius: 5px;
                padding: 20px;
                text-align: center;
                background-color: #f0f8ff;
                margin-top: 20px;
              }
            "))
          ),
          div(
            style = "position: relative; display: inline-block;",
            img(src = "model_chain.png", width = "1800", usemap = "#modelMap"),
            tags$map(
              name = "modelMap",
              tags$area(shape = "rect", coords = "450,700,650,750", href = "#", alt = "Emission Modeling",
                        onmouseover = "document.getElementById('emissionOverlay').style.opacity=1;",
                        onmouseout = "document.getElementById('emissionOverlay').style.opacity=0;",
                        onclick = "Shiny.setInputValue('clicked_block','Emission Modeling');"),
              tags$area(shape = "rect", coords = "625,750,825,800", href = "#", alt = "Air Quality Modeling",
                        onmouseover = "document.getElementById('airQualityOverlay').style.opacity=1;",
                        onmouseout = "document.getElementById('airQualityOverlay').style.opacity=0;",
                        onclick = "Shiny.setInputValue('clicked_block','Air Quality Modeling');"),
              tags$area(shape = "rect", coords = "825,800,1025,900", href = "#", alt = "Exposure",
                        onmouseover = "document.getElementById('exposureOverlay').style.opacity=1;",
                        onmouseout = "document.getElementById('exposureOverlay').style.opacity=0;",
                        onclick = "Shiny.setInputValue('clicked_block','Exposure Assessment');"),
              tags$area(shape = "rect", coords = "1050,800,1200,950", href = "#", alt = "Health",
                        onmouseover = "document.getElementById('healthOverlay').style.opacity=1;",
                        onmouseout = "document.getElementById('healthOverlay').style.opacity=0;",
                        onclick = "Shiny.setInputValue('clicked_block','Health Impact Assessment');")
            ),
            div(id = "emissionOverlay", class = "hover-highlight", style = "top: 250px; left: 200px; width: 650px; height: 300px;"),
            div(id = "airQualityOverlay", class = "hover-highlight", style = "top: 250px; left: 200px; width: 950px; height: 300px;"),
            div(id = "exposureOverlay", class = "hover-highlight", style = "top: 250px; left: 200px; width: 1250px; height: 300px;"),
            div(id = "healthOverlay", class = "hover-highlight", style = "top: 250px; left: 200px; width: 1500px; height: 300px;")
          ),
          p("Click an area to proceed to the next page.")
        )
      )
    ),
    
    # Pipeline Selection
    tabPanel(
      title = "Pipeline Selection",
      fluidRow(
        column(
          12,
          h3("Customize your pipeline"),
          textOutput("userChoice", inline = TRUE),
          br(),
          uiOutput("analysisTypeUI"),
          br(),
          uiOutput("dataSourceUI"),
          br(),
          actionButton("goNext", "Run Model"),
          br(),
          actionButton("back_to_landing", "Back to Landing Page")
        )
      )
    ),
    
    # Emission Calculation
    tabPanel(
      title = "Emission Calculation",
      fluidRow(
        column(12, uiOutput("emCalcUI"))
      )
    ),
    
    # Dispersion Modeling
    tabPanel(title = "Dispersion Modeling",
             uiOutput("dispersionUI")
    ),
    # Visualize AERMOD Results
    tabPanel(title = "Visualize AERMOD Results",
             sidebarLayout(
               sidebarPanel(
                 h4("Upload AERMOD Output"),
                 fileInput("pst_file_base", "Upload Baseline POSTFILE (.pst)", accept = ".pst"),
                 textInput("scenario_name_base", "Baseline Scenario Name", value = "Baseline"),
                 hr(),
                 checkboxInput("compare_scenarios_viz", "Compare with Alternate Scenario?", FALSE),
                 conditionalPanel(
                   condition = "input.compare_scenarios_viz == true",
                   fileInput("pst_file_alt", "Upload Alternate Scenario POSTFILE (.pst)", accept = ".pst"),
                   textInput("scenario_name_alt", "Alternate Scenario Name", value = "Alternate Scenario")
                 ),
                 hr(),
                 h4("Visualization Options"),
                 selectInput("interpolation_res", "Interpolation Grid Resolution", choices = c("30m" = 30, "50m" = 50, "100m" = 100), selected = 30),
                 radioButtons("interpolation_method", "Interpolation Method:",
                              choices = c("Smooth (IDW)" = "idw", "Fast (Nearest Neighbor)" = "nn"),
                              selected = "idw"),
                 actionButton("visualize_results", "Generate Visualizations"),
                 br(),
                 conditionalPanel(
                   condition = "output.is_exposure_assessment_flow || output.is_hia_flow",
                   tagList(
                     br(),
                     actionButton("goto_exposure_tab_btn", "Proceed to Exposure Assessment", icon = icon("arrow-right"), class = "btn-info")
                   )
                 ),
                 br(),
                 actionButton("back_to_pipeline", "Back to Pipeline Selection")
               ),
               mainPanel(
                 uiOutput("results_visuals_ui")
               )
             )
    ),
    # Exposure Assessment Tab
    tabPanel(title = "Exposure Assessment",
             sidebarLayout(
               sidebarPanel(
                 h4("Exposure Assessment"),
                 p("This tool uses the AERMOD results from the previous tab and EPA's EnviroAtlas population data to estimate exposure."),
                 p("Click the button below to run the analysis. This may take a few moments as it downloads population data."),
                 selectInput("exposure_geo", "Visualize Exposure by:",
                             choices = c("Census Tract" = "Tract", "Census Block Group" = "Block Group")),
                 actionButton("run_exposure_assessment", "Calculate & Visualize Exposure"),
                 hr(),
                 conditionalPanel(
                   condition = "output.is_hia_flow",
                   actionButton("goto_hia_tab_btn", "Proceed to Health Impact Assessment", icon = icon("arrow-right"), class = "btn-success")
                 ),
                 actionButton("back_to_viz_results_btn", "Back to Results Upload")
               ),
               mainPanel(
                 withSpinner(uiOutput("exposure_results_display_ui"), type = 5)
               )
             )
    ),
    # Health Impact Assessment Tab
    tabPanel(title = "Health Impact Assessment",
             sidebarLayout(
               sidebarPanel(
                 h4("Health Impact Assessment Setup"),
                 p("Calculate the change in all-cause mortality based on the change in pollutant concentration."),
                 selectInput("hia_pollutant", "Select Pollutant:",
                             choices = names(rr_values)),
                 uiOutput("hia_rr_ui"),
                 selectInput("hia_mortality_geo", "Select Mortality Data:",
                             choices = c("Census Tract" = "Tract", "Census Block Group" = "Block Group", "County" = "County")),
                 actionButton("run_hia", "Run Health Impact Assessment", icon = icon("calculator")),
                 hr(),
                 actionButton("back_to_exposure_btn", "Back to Exposure Assessment")
               ),
               mainPanel(
                 withSpinner(uiOutput("hia_results_display_ui"), type = 5)
               )
             )
    )
  )
)
options(shiny.maxRequestSize=1000*1024^2)
# Server Logic
server <- function(input, output, session) {
  
  # Reactive values
  rv <- reactiveValues(
    stage = NULL,
    selected_county = NULL,
    hpms_links = NULL,
    selected_links = NULL,
    alternate_links = NULL, # For scenario comparison
    link_selection_type = NULL,
    visual_type = NULL,
    drawn_polygon = NULL,
    show_visuals = FALSE,
    show_emissions = FALSE,
    show_scenarios_ui = FALSE,
    scenario_mode = FALSE,
    scenario_title = NULL,
    link_polygons = NULL,
    link_polygons_alternate = NULL,
    near_road_receptors = NULL,
    near_road_receptors_alternate = NULL,
    gridded_receptors = NULL,
    gridded_receptors_alternate = NULL,
    met_dir = NULL,
    clicked_block = NULL,
    baseline_raster = NULL,
    alternate_raster = NULL,
    difference_raster = NULL,
    processed_results_data = NULL,
    is_hourly = FALSE,
    pop_raster = NULL,
    exposure_conc_base = NULL,
    exposure_conc_alt = NULL,
    resampled_conc_raster_base = NULL,
    resampled_conc_raster_alt = NULL,
    pop_raster_resampled = NULL,
    exposure_by_geo = NULL, # For exposure map
    hia_results = NULL, # For HIA map
    total_hia_impact = NULL # For HIA summary
  )
  
  # UI modal feedback wrapper
  with_modal <- function(expr, message = "Processing...") {
    showModal(modalDialog(message, footer = NULL))
    on.exit(removeModal(), add = TRUE)
    force(expr)
  }
  
  # Display chosen modeling block
  output$userChoice <- renderText({
    req(input$clicked_block)
    paste("You selected:", input$clicked_block)
  })
  
  # Switch to Pipeline Selection tab
  observeEvent(input$clicked_block, {
    rv$clicked_block <- input$clicked_block
    updateTabsetPanel(session, "maintabs", selected = "Pipeline Selection")
  })
  
  # Back to Landing Page from Pipeline Selection
  observeEvent(input$back_to_landing, {
    # Full reset
    for(name in names(rv)) {
      rv[[name]] <- NULL
    }
    session$sendCustomMessage(type = "resetInput", message = list(id = "clicked_block"))
    updateTabsetPanel(session, "maintabs", selected = "Landing Page")
  })
  
  # Back button navigations
  observeEvent(input$back_to_pipeline, { updateTabsetPanel(session, "maintabs", selected = "Pipeline Selection") })
  observeEvent(input$back_to_viz_results_btn, { updateTabsetPanel(session, "maintabs", selected = "Visualize AERMOD Results") })
  observeEvent(input$back_to_exposure_btn, { updateTabsetPanel(session, "maintabs", selected = "Exposure Assessment") })
  
  # Analysis Type Selection UI
  output$analysisTypeUI <- renderUI({
    req(input$clicked_block)
    if (input$clicked_block %in% c("Emission Modeling", "Air Quality Modeling")) {
      choices <- if(input$clicked_block == "Emission Modeling") c("Basic Calculation", "Scenario Comparison") else c("Single AERMOD Run", "Scenario Analysis (AERMOD)", "PM Hot Spot Analysis")
      radioGroupButtons(inputId = "analysisType", label = "Select Analysis Type:", choices = choices, justified = TRUE)
    } else if (input$clicked_block %in% c("Exposure Assessment", "Health Impact Assessment")) {
      tagList(
        h4("Instructions"),
        p("This workflow requires AERMOD results. Please use the button below to proceed to the 'Visualize AERMOD Results' tab to upload your data."),
        actionButton("goto_viz_for_exposure", "Go to AERMOD Results Upload")
      )
    }
  })
  
  observeEvent(input$goto_viz_for_exposure, { updateTabsetPanel(session, "maintabs", selected = "Visualize AERMOD Results") })
  observeEvent(input$analysisType, { rv$analysis <- input$analysisType })
  
  # Data Source Selection UI
  output$dataSourceUI <- renderUI({
    req(rv$analysis)
    if (input$clicked_block == "Emission Modeling") {
      tagList(
        h4("Select Traffic Data Source:"),
        prettyRadioButtons(inputId = "traffic_source", label = NULL,
                           choices = c("Built-in (HPMS - Texas Only)", "Custom"), selected = rv$traffic,
                           status = "primary", fill = TRUE),
        if (isTRUE(input$traffic_source == "Custom")) {
          tagList(fileInput("traffic_csv", "Upload Traffic Data CSV", accept = ".csv"),
                  uiOutput("trafficTemplateUI"))
        },
        br(),
        h4("Select Emission Data Source:"),
        prettyRadioButtons(inputId = "emission_source", label = NULL,
                           choices = c("Built-in (MOVES-ERLT - Texas Only)", "Custom"), selected = rv$emission,
                           status = "primary", fill = TRUE),
        if (isTRUE(input$emission_source == "Custom")) {
          tagList(fileInput("emission_csv", "Upload Emission Data CSV", accept = ".csv"),
                  uiOutput("emissionTemplateUI"))
        }
      )
    } else if (input$clicked_block == "Air Quality Modeling") {
      tagList(
        h4("Select Traffic Data Source:"),
        prettyRadioButtons(inputId = "traffic_source", label = NULL,
                           choices = if (rv$analysis == "PM Hot Spot Analysis") c("Custom") else c("Built-in (HPMS - Texas Only)", "Custom"),
                           selected = rv$traffic, status = "primary", fill = TRUE),
        if (isTRUE(input$traffic_source == "Custom")) {
          if (rv$analysis == "PM Hot Spot Analysis") {
            fileInput("traffic_shapefile", "Upload Traffic Shapefile (zipped)", accept = ".zip")
          } else {
            fileInput("traffic_csv", "Upload Traffic Data CSV", accept = ".csv")
          }
        },
        br(),
        h4("Select Emission Data Source:"),
        prettyRadioButtons(inputId = "emission_source", label = NULL,
                           choices = c("Built-in (MOVES-ERLT - Texas Only)", "Custom"), selected = rv$emission,
                           status = "primary", fill = TRUE),
        if (isTRUE(input$emission_source == "Custom")) {
          fileInput("emission_csv", "Upload Emission Data CSV", accept = ".csv")
        },
        br(),
        h4("Select Meteorological Data Source:"),
        prettyRadioButtons(inputId = "met_source", label = NULL,
                           choices = c("Built-in (TCEQ Preprocessed - Texas Only)", "Custom"), selected = rv$met,
                           status = "primary", fill = TRUE),
        if (isTRUE(input$met_source == "Custom")) {
          tagList(fileInput("met_sfc", "Upload MET .SFC File", accept = ".sfc"),
                  fileInput("met_pfl", "Upload MET .PFL File", accept = ".pfl"))
        },
        hr(),
        p("Already have AERMOD results?"),
        actionButton("visualize_aermod_results_btn", "Visualize Results")
      )
    }
  })
  
  observeEvent(input$visualize_aermod_results_btn, {
    updateTabsetPanel(session, "maintabs", selected = "Visualize AERMOD Results")
  })
  
  observeEvent(input$traffic_source, { rv$traffic <- input$traffic_source })
  observeEvent(input$emission_source, { rv$emission <- input$emission_source })
  observeEvent(input$met_source, { rv$met <- input$met_source })
  
  # Run Model
  observeEvent(input$goNext, {
    if (is.null(input$clicked_block) || input$clicked_block %in% c("Exposure Assessment", "Health Impact Assessment")) {
      # Do nothing if exposure or HIA is selected, as the flow is different
      return()
    }
    
    if (input$clicked_block == "Emission Modeling") {
      updateTabsetPanel(session, "maintabs", selected = "Emission Calculation")
    } else if (input$clicked_block == "Air Quality Modeling") {
      updateTabsetPanel(session, "maintabs", selected = "Emission Calculation")
    }
  })
  
  # Emission Calculation UI
  output$emCalcUI <- renderUI({
    if (is.null(rv$stage) || rv$stage == "select_county") {
      tagList(
        h3("Select a County"),
        withSpinner(leafletOutput("county_map"), type = 5),
        actionButton("next_county", "Next"),
        actionButton("back_to_landing_emission", "Back to Landing Page")
      )
    } else if (rv$stage == "display_links") {
      tagList(
        h3("HPMS Links in Selected County"),
        withSpinner(leafletOutput("links_map"), type = 5),
        radioButtons(inputId = "link_selection", label = "Select Links to Use:",
                     choices = c("Use all links in county" = "all", "Use links in selected box" = "selected"), selected = "all"),
        actionButton("visualize_traffic", "Next - Visualize Traffic"),
        conditionalPanel(
          condition = "output.showVisuals == true",
          selectInput(inputId = "visual_type", label = "Select Visual:",
                      choices = c("Map with Color-Coded Links" = "map_links", "AADT by Tract_ID" = "map_tract",
                                  "Column Chart" = "chart", "Hourly VMT Distribution" = "hourly_vmt")),
          uiOutput("visual_output"),
          actionButton("visualize_emissions", "Next - Visualize Emissions"),
          conditionalPanel(
            condition = "output.showEmissions == true",
            fluidRow(
              column(4, selectInput(inputId = "pollutant", label = "Select Pollutant:", choices = unique(emission_factors$Pollutant))),
              column(4, selectInput(inputId = "analysis_year", label = "Analysis Year:", choices = unique(emission_factors$Year))),
              column(4, selectInput(inputId = "emission_visual_type", label = "Select Emission Visual:",
                                    choices = c("Emission Map of All Links" = "emission_map_links",
                                                "Emissions by Tract" = "emission_map_tract",
                                                "Column Plot - Emissions by Road Type" = "emission_chart",
                                                "Hourly Variation of Emissions" = "hourly_emissions")))
            ),
            uiOutput("emission_visual_output"),
            uiOutput("proceed_dispersion_ui"),
            actionButton("show_scenarios_btn", "Compare Scenarios", icon = icon("balance-scale")),
            uiOutput("scenario_ui"), # New UI for scenarios
            actionButton("back_to_landing_emission", "Back to Landing Page")
          )
        )
      )
    }
  })
  
  # Scenario UI
  output$scenario_ui <- renderUI({
    req(rv$show_scenarios_ui)
    tagList(
      hr(),
      h4("Scenario Comparison"),
      radioButtons("scenario_type", "Select Scenario Type:",
                   choices = c("ZEV Adoption", "Operational Improvement"), inline = TRUE),
      conditionalPanel(
        condition = "input.scenario_type == 'ZEV Adoption'",
        selectInput("zev_adoption_rate", "ZEV Adoption Rate:",
                    choices = c("25%", "60%", "70%"))
      ),
      conditionalPanel(
        condition = "input.scenario_type == 'Operational Improvement'",
        selectInput("op_improvement_type", "Operational Improvement Type:",
                    choices = c("Improve restricted road (freeway) speed by 10%",
                                "Improve unrestricted roadway (arterials and others) speed by 10%"))
      ),
      actionButton("calculate_scenario", "Calculate Scenario Emissions"),
      withSpinner(plotOutput("scenario_comparison_plot"), type = 5)
    )
  })
  
  observeEvent(input$show_scenarios_btn, {
    rv$show_scenarios_ui <- !rv$show_scenarios_ui
  })
  
  # Conditional "Proceed to Dispersion Modeling" button
  output$proceed_dispersion_ui <- renderUI({
    if (input$clicked_block != "Emission Modeling") {
      actionButton("proceed_dispersion", "Proceed to Dispersion Modeling")
    }
  })
  # Back to Landing Page from Emission Calculation
  observeEvent(input$back_to_landing_emission, {
    # Full reset
    for(name in names(rv)) {
      rv[[name]] <- NULL
    }
    rv$show_visuals <- FALSE
    rv$show_emissions <- FALSE
    rv$show_scenarios_ui <- FALSE
    rv$scenario_mode <- FALSE
    session$sendCustomMessage(type = "resetInput", message = list(id = "clicked_block"))
    updateTabsetPanel(session, "maintabs", selected = "Landing Page")
  })
  
  output$showVisuals <- reactive({ rv$show_visuals })
  outputOptions(output, "showVisuals", suspendWhenHidden = FALSE)
  output$showEmissions <- reactive({ rv$show_emissions })
  outputOptions(output, "showEmissions", suspendWhenHidden = FALSE)
  
  # New reactive for conditional UI
  output$is_exposure_assessment_flow <- reactive({
    isTRUE(rv$clicked_block == "Exposure Assessment")
  })
  outputOptions(output, "is_exposure_assessment_flow", suspendWhenHidden = FALSE)
  
  # County map
  output$county_map <- renderLeaflet({
    leaflet(texas_counties) %>%
      addTiles() %>%
      addPolygons(layerId = ~GEOID, label = ~NAME,
                  highlightOptions = highlightOptions(color = "white", weight = 2, bringToFront = TRUE))
  })
  
  observeEvent(input$county_map_shape_click, {
    click <- input$county_map_shape_click
    if (!is.null(click)) {
      rv$selected_county <- click$id
      rv$highlighted_county <- click$id
      leafletProxy("county_map") %>%
        clearGroup("highlight") %>%
        addPolygons(data = texas_counties[texas_counties$GEOID == rv$highlighted_county, ],
                    color = "yellow", weight = 3, fillOpacity = 0.5, group = "highlight")
    }
  })
  
  # Load HPMS data from ArcGIS Online after county selection
  observeEvent(input$next_county, {
    if (!is.null(rv$selected_county)) {
      county_fips <- as.integer(rv$selected_county)
      with_modal({
        tryCatch({
          subset <- esri2sf(
            url = "https://services1.arcgis.com/qr14biwnHA6Vis6l/arcgis/rest/services/HPMS_Input/FeatureServer/0",
            where = paste0("CTY_FIPS = ", county_fips)
          )
          
          if (length(subset) == 0) {
            showModal(modalDialog(title = "No Data Found", paste("No HPMS data found for county FIPS:", county_fips), easyClose = TRUE))
            return()
          }
          
          hpms_subset <- st_transform(subset, 4326)
          hpms_subset$ID <- seq_len(nrow(hpms_subset))
          
          rv$hpms_links <- hpms_subset
          rv$stage <- "display_links"
        }, error = function(e) {
          showModal(modalDialog(title = "Error Loading Data", paste("Failed to load HPMS data:", conditionMessage(e)), easyClose = TRUE))
        })
      }, paste("Loading HPMS data for county FIPS:", county_fips))
    } else {
      showModal(modalDialog(title = "No County Selected", "Please select a county on the map before proceeding.", easyClose = TRUE))
    }
  })
  
  
  # Links map with drawing tools
  output$links_map <- renderLeaflet({
    tryCatch({
      req(rv$stage == "display_links", rv$selected_county, rv$hpms_links)
      filtered_links <- rv$hpms_links
      selected_geom <- texas_counties[texas_counties$GEOID == rv$selected_county, ]
      filtered_links <- sf::st_zm(filtered_links, drop = TRUE)
      selected_geom <- sf::st_zm(selected_geom, drop = TRUE)
      if (!all(sf::st_is_valid(selected_geom))) selected_geom <- sf::st_make_valid(selected_geom)
      if (!all(sf::st_is_valid(filtered_links))) filtered_links <- sf::st_make_valid(filtered_links)
      leaflet() %>%
        addTiles() %>%
        addPolygons(data = selected_geom, color = "blue", weight = 2, fillOpacity = 0.2, label = ~NAME) %>%
        addPolylines(data = filtered_links, color = "red", weight = 2, opacity = 0.7, group = "links") %>%
        addDrawToolbar(targetGroup = "draw", editOptions = editToolbarOptions(selectedPathOptions = selectedPathOptions()))
    }, error = function(e) {
      leaflet() %>% addTiles() %>% addPopups(lng = -98.5, lat = 31.0, popup = paste("Error rendering map:", conditionMessage(e)))
    })
  })
  
  observeEvent(input$links_map_draw_new_feature, {
    new_feature <- input$links_map_draw_new_feature
    feature_collection <- list(type = "FeatureCollection", features = list(new_feature))
    drawn_sf <- geojsonio::geojson_sf(jsonlite::toJSON(feature_collection, auto_unbox = TRUE))
    drawn_sf <- st_set_crs(drawn_sf, st_crs(rv$hpms_links))
    filtered_links <- rv$hpms_links
    rv$selected_links <- st_intersection(filtered_links, drawn_sf)
    leafletProxy("links_map") %>%
      clearGroup("selected_links") %>%
      addPolylines(data = rv$selected_links, color = "green", weight = 3, group = "selected_links")
  })
  
  observeEvent(input$link_selection, { rv$link_selection_type <- input$link_selection })
  
  observeEvent(input$visualize_traffic, {
    req(rv$selected_county, rv$hpms_links)
    filtered_links <- rv$hpms_links
    if (rv$link_selection_type == "all") {
      rv$selected_links <- filtered_links
    } else if (rv$link_selection_type == "selected" && !is.null(rv$selected_links)) {
      # Already filtered by polygon
    } else {
      showNotification("Please draw a polygon to select links.", type = "warning")
      return()
    }
    rv$show_visuals <- TRUE
  })
  
  observeEvent(input$visualize_emissions, { rv$show_emissions <- TRUE })
  
  observeEvent(input$proceed_dispersion, {
    updateTabsetPanel(session, "maintabs", selected = "Dispersion Modeling")
  })
  
  # New observer for the button to go to exposure assessment
  observeEvent(input$goto_exposure_tab_btn, {
    updateTabsetPanel(session, "maintabs", selected = "Exposure Assessment")
  })
  
  # Reusable function to calculate emissions (long format for plotting)
  calculate_hourly_emissions_fn <- function(links_df, pollutant, year) {
    req(links_df, pollutant, year)
    ef_filtered <- emission_factors %>%
      filter(Pollutant == pollutant, Year == year)
    
    link_emissions_df <- links_df %>%
      left_join(ef_filtered, by = c("DIST_NM" = "District", "Road_Type" = "Road_Type", "Average_Sp" = "Average_Speed")) %>%
      mutate(
        Daily_Emissions = AADT * Length_mil * Emission_Factor_g_mile,
        Length_m = Length_mil * 1609.34
      ) %>% filter(!is.na(Daily_Emissions))
    
    hourly_emissions_df <- link_emissions_df %>%
      left_join(hourly_factors_df, by = c("DIST_NM" = "District")) %>%
      mutate(
        Width_m = Lanes * 12 * 0.3048,
        Hourly_Emissions_g_hr = Daily_Emissions * Hourly_Fraction,
        Hourly_Emissions = Hourly_Emissions_g_hr / (3600 * Length_m * Width_m) # g/s.m²
      )
    
    return(hourly_emissions_df)
  }
  
  # Prepare Scenario Links
  observeEvent(input$calculate_scenario, {
    with_modal({
      req(rv$selected_links, input$scenario_type)
      
      alt_links <- rv$selected_links
      scenario_title <- ""
      
      if (input$scenario_type == "ZEV Adoption") {
        rate <- as.numeric(gsub("%", "", input$zev_adoption_rate)) / 100
        alt_links$AADT <- alt_links$AADT * (1 - rate)
        scenario_title <- paste0("ZEV Adoption ", input$zev_adoption_rate)
      } else if (input$scenario_type == "Operational Improvement") {
        if (input$op_improvement_type == "Improve restricted road (freeway) speed by 10%") {
          filter_cond <- str_detect(alt_links$Road_Type, "Restricted")
          scenario_title <- "Improve Freeway Speed by 10%"
        } else {
          filter_cond <- str_detect(alt_links$Road_Type, "Unrestricted")
          scenario_title <- "Improve Arterial Speed by 10%"
        }
        filter_cond[is.na(filter_cond)] <- FALSE
        alt_links$Average_Sp[filter_cond] <- alt_links$Average_Sp[filter_cond] * 1.10
        alt_links$Lanes[filter_cond] <- alt_links$Lanes[filter_cond] + 2
      }
      
      rv$alternate_links <- alt_links
      rv$scenario_title <- scenario_title
      rv$scenario_mode <- TRUE
    }, "Preparing scenario data...")
  })
  
  # Scenario Comparison Plot
  output$scenario_comparison_plot <- renderPlot({
    req(rv$scenario_mode, rv$selected_links, rv$alternate_links, input$pollutant, input$analysis_year)
    
    baseline_hourly <- calculate_hourly_emissions_fn(rv$selected_links, input$pollutant, input$analysis_year)
    alternate_hourly <- calculate_hourly_emissions_fn(rv$alternate_links, input$pollutant, input$analysis_year)
    
    req(baseline_hourly, alternate_hourly)
    
    baseline_data <- baseline_hourly %>% mutate(Scenario = "Baseline")
    alternate_data <- alternate_hourly %>% mutate(Scenario = rv$scenario_title)
    
    combined_data <- rbind(baseline_data, alternate_data)
    
    
    ggboxplot(combined_data, x = "Hour_ID", y = "Hourly_Emissions_g_hr", color = "black",fill = "Scenario",  scales   = "free_x",
              palette  = "jco",
              outlier.shape = 1,
              outlier.size  = 0.1) +
      stat_compare_means(aes(group = Scenario), method = "t.test", label = "p.signif") + 
      yscale("log10",.format= TRUE) +
      geom_vline(
        xintercept = seq(0.5, length(unique(combined_data$Hour_ID)) + 0.5, by = 1),
        linetype    = "solid",
        color       = "gray",
        size        = 0.5
      ) + 
      labs(title = "Hourly Emissions Comparison by Scenario", x = "Hour", y = "Hourly Emissions (g/h)") +
      theme_minimal()
  })
  
  # Calculate link-level daily emissions (for single run visuals)
  link_emissions <- reactive({
    req(rv$selected_links, input$pollutant, input$analysis_year)
    ef_filtered <- emission_factors %>%
      filter(Pollutant == input$pollutant, Year == input$analysis_year)
    rv$selected_links %>%
      left_join(ef_filtered, by = c("DIST_NM" = "District", "Road_Type" = "Road_Type", "Average_Sp" = "Average_Speed")) %>%
      mutate(Daily_Emissions = AADT * Length_mil * Emission_Factor_g_mile,
             Length_m = Length_mil * 1609.34)
  })
  # Calculate hourly emissions in g/s.m²
  hourly_emissions <- reactive({
    req(link_emissions(), rv$selected_links)
    hourly_factors <- read.csv("data/Hourly_Factors.csv")
    link_data <- link_emissions() %>%
      left_join(hourly_factors, by = c("DIST_NM" = "District")) %>%
      mutate(
        Width_m = Lanes * 12 * 0.3048,  # Assuming 12 ft per lane, converted to meters
        Hourly_Emissions_g_hr = Daily_Emissions * Hourly_Fraction,  # g/hr
        Hourly_Emissions = Hourly_Emissions_g_hr / (3600 * Length_m * Width_m)  # g/s.m²
      )
    link_data
  }) 
  # Reusable function to generate buffered polygons with hourly emissions for AERMOD
  generate_link_polygons <- function(links_sf, pollutant, year) {
    with_modal({
      req(links_sf, pollutant, year)
      
      ef_filtered <- emission_factors %>%
        filter(Pollutant == pollutant, Year == year)
      
      link_emissions_df <- links_sf %>%
        left_join(ef_filtered, by = c("DIST_NM" = "District", "Road_Type" = "Road_Type", "Average_Sp" = "Average_Speed")) %>%
        mutate(
          Daily_Emissions = AADT * Length_mil * Emission_Factor_g_mile,
          Length_m = Length_mil * 1609.34,
          Width_m = Lanes * 12 * 0.3048
        ) %>%
        filter(!is.na(Daily_Emissions))
      
      links_proj <- st_transform(link_emissions_df, 32614)
      buffer_dist <- links_proj$Width_m / 2
      polygons_proj <- st_buffer(links_proj, dist = buffer_dist, endCapStyle = "FLAT")
      polygons_proj <- st_make_valid(polygons_proj)
      polygons_proj <- polygons_proj[st_area(polygons_proj) > units::set_units(1, "m^2"), ]
      
      if (nrow(polygons_proj) == 0) {
        showNotification("No valid polygons generated after buffering.", type = "error")
        return(NULL)
      }
      
      polygons_proj$Area_m2 <- as.numeric(st_area(polygons_proj))
      
      hourly_factors_wide <- hourly_factors_df %>%
        dplyr::select(-Period) %>%
        pivot_wider(names_from = Hour_ID, values_from = Hourly_Fraction, names_prefix = "Hour_")
      
      polygons_with_hourly <- polygons_proj %>%
        left_join(hourly_factors_wide, by = c("DIST_NM" = "District"))
      
      for(h in 1:24){
        frac_col <- paste0("Hour_", h)
        aermod_col <- paste0("Hour_", h, "_g_s_m2")
        if(frac_col %in% names(polygons_with_hourly)) {
          polygons_with_hourly[[aermod_col]] <- (polygons_with_hourly$Daily_Emissions * polygons_with_hourly[[frac_col]]) / (3600 * polygons_with_hourly$Area_m2)
        }
      }
      
      return(st_transform(polygons_with_hourly, 4326))
    }, "Generating link polygons...")
  }
  
  # Observers for generating polygons in different modes
  observeEvent(input$generate_polygons, {
    rv$link_polygons <- generate_link_polygons(rv$selected_links, input$pollutant, input$analysis_year)
  })
  observeEvent(input$generate_polygons_base, {
    rv$link_polygons <- generate_link_polygons(rv$selected_links, input$pollutant, input$analysis_year)
  })
  observeEvent(input$generate_polygons_alt, {
    req(rv$alternate_links)
    rv$link_polygons_alternate <- generate_link_polygons(rv$alternate_links, input$pollutant, input$analysis_year)
  })
  
  # Generate near-road receptors
  observeEvent(input$generate_near_road_receptors, {
    req(rv$link_polygons)
    grid_points <- generate_receptors(rv$link_polygons)
    rv$near_road_receptors <- st_transform(st_sf(geometry = grid_points, crs = 32614), 4326)
    print("Generated near-road receptors")
  })
  
  # Generate gridded receptors
  observeEvent(input$generate_gridded_receptors, {
    req(rv$link_polygons, input$grid_resolution)
    resolution <- as.numeric(gsub("m", "", input$grid_resolution))
    bbox <- st_bbox(st_transform(rv$link_polygons, 32614))
    bbox_polygon <- st_as_sfc(bbox, crs = 32614)
    grid_points <- st_make_grid(bbox_polygon, cellsize = resolution, what = "centers")
    rv$gridded_receptors <- st_transform(st_sf(geometry = grid_points, crs = 32614), 4326)
    print(paste("Generated gridded receptors at", input$grid_resolution, "resolution"))
  })
  
  # Render polygon map
  output$polygon_map <- renderLeaflet({
    req(rv$link_polygons)
    pal <- colorNumeric("viridis", domain = rv$link_polygons$Daily_Emissions)
    leaflet(rv$link_polygons) %>%
      addTiles() %>%
      addPolygons(color = ~pal(Daily_Emissions), weight = 2,
                  label = ~paste("Emissions:", round(Daily_Emissions, 2), "g/day")) %>%
      leaflet::addLegend(pal = pal, values = rv$link_polygons$Daily_Emissions, title = "Emissions (g/day)", position = "bottomright")
  })
  
  # Update map with receptors
  observe({
    req(rv$link_polygons)
    map <- leafletProxy("polygon_map") %>% clearGroup("receptors")
    if (!is.null(rv$near_road_receptors)) {
      map <- map %>% addCircleMarkers(data = rv$near_road_receptors, radius = 2, color = "red", fillOpacity = 0.5,
                                      label = "Near-Road Receptor", group = "receptors")
    }
    if (!is.null(rv$gridded_receptors)) {
      map <- map %>% addCircleMarkers(data = rv$gridded_receptors, radius = 2, color = "green", fillOpacity = 0.5,
                                      label = "Gridded Receptor", group = "receptors")
    }
  })
  
  # Receptor table for AERMOD
  receptor_table <- reactive({
    req(rv$gridded_receptors)
    receptors_sf <- st_transform(rv$gridded_receptors, 32614)
    coords <- st_coordinates(receptors_sf)
    data.frame(ID = seq_len(nrow(coords)), X = coords[, "X"], Y = coords[, "Y"])
  })
  
  # Function to simplify & convert polygons to a table (limiting vertices to 20)
  convert_polygons_to_table <- function(polygons) {
    polygons <- st_transform(polygons, 32614)
    polygons <- st_make_valid(polygons)
    polygons <- suppressWarnings(st_cast(polygons, "POLYGON"))
    
    valid_polygons <- polygons[sapply(st_geometry(polygons), function(g) {
      coords <- st_coordinates(g)
      nrow(coords) >= 4 && st_is_valid(g)
    }), ]
    
    if (nrow(valid_polygons) == 0) stop("No valid polygons.")
    
    max_vertices <- 20
    coord_lists <- lapply(st_geometry(valid_polygons), function(geom) {
      coords <- st_coordinates(st_zm(geom, drop = TRUE))[, 1:2]
      n <- nrow(coords)
      if (n > max_vertices) {
        warning(sprintf("Polygon ID %s truncated from %d to %d vertices", valid_polygons$ID[1], n, max_vertices))
        coords <- coords[1:max_vertices, , drop = FALSE]
      }
      coords
    })
    
    coords_df <- do.call(rbind, lapply(coord_lists, function(coords) {
      coords_padded <- matrix(NA, nrow = 1, ncol = max_vertices * 2)
      coords_flat <- c(rbind(coords[, "X"], coords[, "Y"]))
      coords_padded[1, 1:length(coords_flat)] <- coords_flat
      coords_padded
    }))
    
    colnames(coords_df) <- c(rbind(paste0("Xv.", 1:max_vertices), paste0("Yv.", 1:max_vertices)))
    
    hourly_cols <- grep("_g_s_m2$", names(valid_polygons), value = TRUE)
    hourly_df <- st_drop_geometry(valid_polygons)[, hourly_cols, drop = FALSE]
    
    final_df <- data.frame(
      ID = valid_polygons$ID,
      Year = valid_polygons$Year,
      Length_m = valid_polygons$Length_mil,
      Width_m = valid_polygons$Width_m,
      Nvert = sapply(coord_lists, nrow),
      coords_df,
      hourly_df,
      stringsAsFactors = FALSE
    )
    return(final_df)
  }
  
  source_table <- reactive({
    req(rv$link_polygons)
    convert_polygons_to_table(rv$link_polygons)
  })
  
  # AERMOD input generation
  generate_aermod_text <- function(receptor_df, source_df, poll, title, averaging_time, county_name, met_year, met_roughness) {
    RE <- receptor_df
    SO <- source_df
    
    county_data <- county_lookup[county_lookup$County_Name == county_name, ]
    if (nrow(county_data) == 0) county_data <- county_lookup[county_lookup$County_Name == "Bexar", ]
    
    surface_wban <- county_data$Surface_Station_ID
    upper_wban <- county_data$Upper_Air_Station_ID
    elevation <- county_data$Surface_Station_Elevation_meters
    roughness <- county_data$Bowen_Ratio
    population <- case_when(
      county_name == "Bexar" ~ 2100000, county_name == "Harris" ~ 5000000,
      county_name == "Dallas" ~ 2800000, county_name == "Tarrant" ~ 2300000,
      county_name == "Travis" ~ 1400000, TRUE ~ 100000
    )
    
    surface_id <- county_data$Surface_ID
    upper_air_id <- county_data$Upper_Air_ID
    county_upper <- toupper(gsub(" ", "", county_name))
    file_prefix <- paste0(county_upper, "_", surface_id, upper_air_id)
    sfc_file <- paste0(file_prefix, met_year, met_roughness, ".SFC")
    pfl_file <- paste0(file_prefix, met_year, met_roughness, ".PFL")
    met_year_full <- paste0("20", met_year)
    
    co_lines <- c("CO STARTING", paste("CO TITLEONE", title), "CO MODELOPT FLAT CONC", paste("CO POLLUTID", poll),
                  paste("CO AVERTIME", averaging_time), "CO FLAGPOLE 1.8", sprintf("CO URBANOPT %d %.1f", population, roughness),
                  "CO RUNORNOT RUN", "CO FINISHED")
    
    source_starting <- c("SO STARTING", "SO ELEVUNIT METERS")
    source_location_lines <- sapply(1:nrow(SO), function(i) sprintf("SO LOCATION %s AREAPOLY %.1f %.1f", SO$ID[i], SO$Xv.1[i], SO$Yv.1[i]))
    source_param_lines <- sapply(1:nrow(SO), function(i) sprintf("SO SRCPARAM %s %.0f %.1f %.0f %.1f", SO$ID[i], 1, 2.6, SO$Nvert[i], 2.37))
    generate_so_areavert_line <- function(row) {
      nvert <- as.numeric(row["Nvert"])
      parts <- c(sprintf("SO AREAVERT %s", row["ID"]))
      for (j in 1:nvert) {
        parts <- c(parts, sprintf("%.1f %.1f", as.numeric(row[[paste0("Xv.", j)]]), as.numeric(row[[paste0("Yv.", j)]])))
      }
      paste(parts, collapse = " ")
    }
    so_areavert_lines <- apply(SO, 1, generate_so_areavert_line)
    source_emi_lines <- sapply(1:nrow(SO), function(i) {
      hourly_cols <- paste0("Hour_", 1:24, "_g_s_m2")
      hourly_rates <- as.numeric(SO[i, hourly_cols])
      hourly_rates_fmt <- sprintf("%.6e", hourly_rates)
      paste("SO EMISFACT", SO$ID[i], "HROFDY", paste(hourly_rates_fmt, collapse = " "))
    })
    source_finished <- c("SO URBANSRC ALL", "SO SRCGROUP ALL", "SO FINISHED")
    
    receptor_lines <- c("RE STARTING", "RE ELEVUNIT METERS",
                        sapply(1:nrow(RE), function(i) sprintf("RE DISCCART %.1f %.1f", RE$X[i], RE$Y[i])),
                        "RE FINISHED")
    
    met_lines <- c("ME STARTING", sprintf("ME SURFDATA %s %s", surface_wban, met_year_full),
                   sprintf("ME UAIRDATA %s %s", upper_wban, met_year_full), sprintf("ME SITEDATA 99999 %s TEXAS", met_year_full),
                   sprintf("ME PROFBASE %.1f METERS", elevation), sprintf("ME SURFFILE %s", sfc_file),
                   sprintf("ME PROFFILE %s", pfl_file), "ME FINISHED")
    
    ou_lines <- c("OU STARTING", paste("OU RECTABLE", averaging_time, "FIRST"), paste("OU MAXTABLE", averaging_time, "1"),
                  paste("OU DAYTABLE", averaging_time), paste("OU PLOTFILE", averaging_time, "ALL", "FIRST", paste0("ALL", averaging_time, "FIRST.plt"), "10000"),
                  paste("OU POSTFILE", averaging_time, "ALL", "PLOT", paste0("ALL", averaging_time, ".pst"), "10001"), "OU FINISHED")
    
    paste(c(co_lines, source_starting, source_location_lines, source_param_lines, so_areavert_lines,
            source_emi_lines, source_finished, receptor_lines, met_lines, ou_lines), collapse = "\n")
  }
  
  aermod_input_baseline <- reactive({
    req(receptor_table(), source_table(), input$pollutant, input$aermod_title, input$averaging_time, rv$selected_county, input$met_year, input$met_roughness)
    county_name <- texas_counties$NAME[texas_counties$GEOID == rv$selected_county]
    generate_aermod_text(receptor_table(), source_table(), input$pollutant, input$aermod_title, input$averaging_time, county_name, input$met_year, input$met_roughness)
  })
  
  aermod_input_alternate <- reactive({
    req(rv$gridded_receptors_alternate, rv$link_polygons_alternate, input$pollutant, input$aermod_title, input$averaging_time, rv$selected_county, input$met_year, input$met_roughness)
    
    receptors_alt_sf <- st_transform(rv$gridded_receptors_alternate, 32614)
    coords_alt <- st_coordinates(receptors_alt_sf)
    receptor_df_alt <- data.frame(ID = seq_len(nrow(coords_alt)), X = coords_alt[, "X"], Y = coords_alt[, "Y"])
    
    source_df_alt <- convert_polygons_to_table(rv$link_polygons_alternate)
    
    county_name <- texas_counties$NAME[texas_counties$GEOID == rv$selected_county]
    title_alt <- paste(input$aermod_title, "-", rv$scenario_title)
    
    generate_aermod_text(receptor_df_alt, source_df_alt, input$pollutant, title_alt, input$averaging_time, county_name, input$met_year, input$met_roughness)
  })
  
  output$download_aermod <- downloadHandler(
    filename = function() {
      paste("AERMOD_Input_", input$pollutant, "_", input$averaging_time, ".txt", sep = "")
    },
    content = function(file) {
      writeLines(aermod_input_baseline(), file)
    }
  )
  
  output$met_county_text <- renderText({
    req(rv$selected_county)
    county_name <- texas_counties$NAME[texas_counties$GEOID == rv$selected_county]
    if(length(county_name) == 0) county_name <- "Unknown"
    county_name
  })
  
  output$met_option_ui <- renderUI({
    req(rv$selected_county)
    county_name <- texas_counties$NAME[texas_counties$GEOID == rv$selected_county]
    county_lower <- tolower(gsub(" ", "", county_name))
    if(county_lower %in% c("harris", "travis")) {
      radioButtons("met_option", "Select MET Option:",
                   choices = c("Default" = "default", "Alternative" = "alternative"), selected = "default")
    }
  })
  
  observeEvent(input$download_met_data, {
    req(rv$selected_county)
    with_modal({
      county_name <- texas_counties$NAME[texas_counties$GEOID == rv$selected_county]
      met_url_info <- get_met_url(county_name)
      met_url <- if(is.list(met_url_info)) met_url_info[[input$met_option %||% "default"]] else met_url_info
      
      temp_zip <- tempfile(fileext = ".zip")
      download.file(met_url, destfile = temp_zip, mode = "wb")
      
      met_temp_dir <- tempfile()
      dir.create(met_temp_dir)
      unzip(temp_zip, exdir = met_temp_dir, overwrite = TRUE, junkpaths = TRUE)
      rv$met_dir <- met_temp_dir
      met_files <- list.files(met_temp_dir, full.names = TRUE, recursive = TRUE)
      years_found <- unique(stringr::str_extract(basename(met_files), "[0-9]{2}"))
      years_found <- years_found[!is.na(years_found)]
      if(length(years_found)==0) years_found <- "20"
      updateSelectInput(session, "met_year", choices = sort(years_found), selected = ifelse("20" %in% years_found, "20", years_found[1]))
    }, "Downloading MET Data...")
  })
  
  observeEvent(input$package_and_download, {
    req(rv$met_dir, aermod_input_baseline(), rv$selected_county)
    with_modal({
      county_name <- texas_counties$NAME[texas_counties$GEOID == rv$selected_county]
      
      if (isTRUE(rv$scenario_mode)) {
        req(aermod_input_alternate(), rv$scenario_title)
        zip_file <- package_scenario_files(aermod_input_baseline(), aermod_input_alternate(), rv$scenario_title, rv$met_dir, input$met_year, input$met_roughness, county_name)
        filename <- "AERMOD_Scenarios_Package.zip"
      } else {
        zip_file <- package_met_files(aermod_input_baseline(), rv$met_dir, input$met_year, input$met_roughness, county_name)
        filename <- "AERMOD_Package.zip"
      }
      
      output$download_met_package_ui <- renderUI({
        downloadButton("download_met_package", "Download AERMOD Package")
      })
      output$download_met_package <- downloadHandler(
        filename = function() { filename },
        content = function(file) { file.copy(zip_file, file) },
        contentType = "application/zip"
      )
    }, "Packaging AERMOD Package...")
  })
  
  output$visual_output <- renderUI({
    req(rv$visual_type)
    switch(rv$visual_type,
           "map_links" = withSpinner(leafletOutput("map_links_output"), type = 5),
           "map_tract" = withSpinner(leafletOutput("map_tract_output"), type = 5),
           "chart" = withSpinner(plotOutput("chart_output"), type = 5),
           "hourly_vmt" = withSpinner(plotOutput("hourly_vmt_plot"), type = 5))
  })
  
  output$map_links_output <- renderLeaflet({
    req(rv$visual_type == "map_links", rv$selected_links)
    links <- rv$selected_links
    pal <- colorFactor(palette = c("orange", "yellow", "blue", "lightblue"),
                       domain = c("Rural Unrestricted Access", "Rural Restricted Access", 
                                  "Urban Unrestricted Access", "Urban Restricted Access"))
    weights <- ifelse(links$Road_Type %in% c("Rural Unrestricted Access", "Urban Unrestricted Access"), 2, 4)
    leaflet(links) %>%
      addTiles() %>%
      addPolylines(color = ~pal(Road_Type), weight = weights, opacity = 0.8, label = ~Road_Type)
  })
  
  output$map_tract_output <- renderLeaflet({
    req(rv$visual_type == "map_tract", rv$selected_links)
    tracts <- tracts(state = "TX", county = substr(rv$selected_county, 3, 5), cb = TRUE)
    tracts <- st_transform(tracts, 4326)
    aadt_by_tract <- rv$selected_links %>%
      st_drop_geometry() %>%
      group_by(TRACT_ID) %>%
      summarise(Total_VMT = sum(AADT * Length_mil, na.rm = TRUE)) %>%
      mutate(TRACT_ID = as.character(TRACT_ID)) %>%
      left_join(tracts, by = c("TRACT_ID" = "GEOID"))
    aadt_sf <- st_as_sf(aadt_by_tract)
    pal <- colorBin("YlOrRd", domain = aadt_sf$Total_VMT, bins = 5)
    leaflet(aadt_sf) %>%
      addTiles() %>%
      addPolygons(fillColor = ~pal(Total_VMT), fillOpacity = 0.7, color = "white", weight = 1,
                  label = ~paste("Tract:", TRACT_ID, "VMT:", Total_VMT)) %>%
      addLegend(pal = pal, values = ~Total_VMT, title = "Total VMT", position = "bottomright")
  })
  
  output$chart_output <- renderPlot({
    req(rv$visual_type == "chart", rv$selected_links)
    chart_data <- rv$selected_links %>%
      st_drop_geometry() %>%
      group_by(Road_Type) %>%
      summarise(Total_VMT = sum(AADT * Length_mil, na.rm = TRUE))
    total_links <- nrow(rv$selected_links)
    ggplot(chart_data, aes(x = Road_Type, y = Total_VMT)) +
      geom_col(fill = "steelblue") +
      labs(title = "Total VMT by Road Type", x = "Road Type", y = "Total VMT",
           caption = paste("Total number of links:", total_links)) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$hourly_vmt_plot <- renderPlot({
    req(rv$visual_type == "hourly_vmt", rv$selected_links)
    selected_links <- rv$selected_links %>% mutate(VMT = AADT * Length_mil)
    vmt_by_type_district <- selected_links %>%
      st_drop_geometry() %>%
      group_by(Road_Type, DIST_NM) %>%
      summarise(total_VMT = sum(VMT, na.rm = TRUE))
    hourly_vmt <- vmt_by_type_district %>%
      left_join(hourly_factors_df, by = c("DIST_NM" = "District")) %>%
      mutate(hourly_VMT = total_VMT * Hourly_Fraction)
    total_hourly_vmt <- hourly_vmt %>%
      group_by(Road_Type, Hour_ID) %>%
      summarise(total_hourly_VMT = sum(hourly_VMT, na.rm = TRUE))
    ggplot(total_hourly_vmt, aes(x = Hour_ID, y = total_hourly_VMT, fill = Road_Type)) +
      geom_area() +
      labs(x = "Hour", y = "Hourly VMT", title = "Hourly VMT Distribution by Road Type") +
      scale_x_continuous(breaks = 1:24) +
      theme_minimal()
  })
  
  observeEvent(input$visual_type, { rv$visual_type <- input$visual_type })
  
  output$emission_visual_output <- renderUI({
    req(input$emission_visual_type)
    switch(input$emission_visual_type,
           "emission_map_links" = leafletOutput("emission_map_links"),
           "emission_map_tract" = leafletOutput("emission_map_tract"),
           "emission_chart" = plotOutput("emission_chart"),
           "hourly_emissions" = plotOutput("hourly_emissions"))
  })
  
  output$emission_map_links <- renderLeaflet({
    req(link_emissions())
    pal <- colorNumeric("viridis", domain = link_emissions()$Daily_Emissions)
    leaflet(link_emissions()) %>%
      addTiles() %>%
      addPolylines(color = ~pal(Daily_Emissions), weight = 2,
                   label = ~paste("Emissions:", round(Daily_Emissions, 2), "g/day")) %>%
      leaflet::addLegend(pal = pal, values = link_emissions()$Daily_Emissions, title = "Emissions (g/day)", position = "bottomright")
  })
  
  output$emission_map_tract <- renderLeaflet({
    req(link_emissions())
    emissions_by_tract <- link_emissions() %>%
      st_drop_geometry() %>%
      group_by(TRACT_ID) %>%
      summarise(Total_Emissions = sum(Daily_Emissions, na.rm = TRUE)) %>%
      mutate(TRACT_ID = as.character(TRACT_ID)) %>%
      left_join(tracts(state = "TX", county = substr(rv$selected_county, 3, 5), cb = TRUE), by = c("TRACT_ID" = "GEOID"))
    emissions_sf <- st_as_sf(emissions_by_tract)
    pal <- colorNumeric("viridis", domain = emissions_sf$Total_Emissions)
    leaflet(emissions_sf) %>%
      addTiles() %>%
      addPolygons(fillColor = ~pal(Total_Emissions), fillOpacity = 0.7, color = "white", weight = 1,
                  label = ~paste("Tract:", TRACT_ID, "Emissions:", round(Total_Emissions, 2), "g/day")) %>%
      leaflet::addLegend(pal = pal, values = emissions_sf$Total_Emissions, title = "Total Emissions (g/day)", position = "bottomright")
  })
  
  output$emission_chart <- renderPlot({
    req(link_emissions())
    chart_data <- link_emissions() %>%
      st_drop_geometry() %>%
      group_by(Road_Type) %>%
      summarise(Total_Emissions = sum(Daily_Emissions, na.rm = TRUE))
    ggplot(chart_data, aes(x = Road_Type, y = Total_Emissions)) +
      geom_col(fill = "steelblue") +
      labs(title = "Total Emissions by Road Type", x = "Road Type", y = "Total Emissions (g/day)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$hourly_emissions <- renderPlot({
    req(hourly_emissions())
    hourly_data <- hourly_emissions() %>%
      st_drop_geometry() %>%
      group_by(Road_Type, Hour_ID) %>%
      summarise(Total_Hourly_Emissions = sum(Hourly_Emissions_g_hr, na.rm = TRUE))
    ggplot(hourly_data, aes(x = Hour_ID, y = Total_Hourly_Emissions, fill = Road_Type)) +
      geom_area() +
      labs(title = "Hourly Variation of Emissions", x = "Hour", y = "Hourly Emissions (g/h)") +
      scale_x_continuous(breaks = 1:24) +
      theme_minimal()
  })
  
  # --- Dispersion Modeling UI and Logic for Scenarios ---
  output$dispersionUI <- renderUI({
    if (isTRUE(rv$scenario_mode)) {
      # Two-column layout for scenario comparison
      fluidRow(
        column(6, 
               h3("Baseline Scenario"),
               actionButton("generate_polygons_base", "Generate Link Polygons"),
               actionButton("generate_receptors_base", "Generate Receptors"),
               selectInput("grid_resolution_base", "Grid Resolution:", choices = c("250m", "500m", "1000m"), selected = "500m"),
               withSpinner(leafletOutput("polygon_map_base"), type = 5)
        ),
        column(6,
               h3(paste("Alternate Scenario:", rv$scenario_title)),
               actionButton("generate_polygons_alt", "Generate Link Polygons"),
               actionButton("generate_receptors_alt", "Generate Receptors"),
               selectInput("grid_resolution_alt", "Grid Resolution:", choices = c("250m", "500m", "1000m"), selected = "500m"),
               withSpinner(leafletOutput("polygon_map_alt"), type = 5)
        ),
        column(12, hr(), 
               h3("AERMOD Run Setup (Common for Both Scenarios)"),
               fluidRow(
                 column(6,
                        textInput("aermod_title", "AERMOD Title:", placeholder = "Enter base title for AERMOD run"),
                        selectInput("averaging_time", "Select Averaging Time:", choices = c("Hourly" = "1", "Daily" = "24"), selected = "1"),
                        h4("MET County: ", textOutput("met_county_text", inline = TRUE)),
                        uiOutput("met_option_ui"),
                        actionButton("download_met_data", "Download & Unpack MET ZIP")
                 ),
                 column(6,
                        selectInput("met_year", "Select MET Year", choices = NULL),
                        selectInput("met_roughness", "Select Surface Roughness", choices = c("L", "M", "H"), selected = "M"),
                        actionButton("package_and_download", "Package and Download AERMOD Scenarios"),
                        uiOutput("download_met_package_ui"),
                        br(),
                        actionButton("back_to_landing_dispersion", "Back to Landing Page")
                 )
               )
        )
      )
    } else {
      # Original single-run layout
      fluidRow(
        column(6,
               h3("Dispersion Modeling"),
               actionButton("generate_polygons", "Generate Link Polygons"),
               actionButton("generate_near_road_receptors", "Generate Near Road Receptors"),
               selectInput("grid_resolution", "Grid Resolution:", choices = c("250m", "500m", "1000m"), selected = "500m"),
               actionButton("generate_gridded_receptors", "Generate Gridded Receptors"),
               textInput("aermod_title", "AERMOD Title:", placeholder = "Enter title for AERMOD run"),
               selectInput("averaging_time", "Select Averaging Time:", choices = c("Hourly" = "1", "Daily" = "24"), selected = "1"),
               actionButton("generate_aermod_input", "Generate AERMOD Input"),
               downloadButton("download_aermod", "Download AERMOD Input"),
               withSpinner(leafletOutput("polygon_map"), type = 5)
        ),
        column(6,
               h3("MET Data for Dispersion Modeling"),
               h4("MET County: ", textOutput("met_county_text", inline = TRUE)),
               uiOutput("met_option_ui"),
               actionButton("download_met_data", "Download & Unpack MET ZIP"),
               br(), br(),
               selectInput("met_year", "Select MET Year", choices = NULL),
               selectInput("met_roughness", "Select Surface Roughness", choices = c("L", "M", "H"), selected = "M"),
               actionButton("package_and_download", "Package and Download AERMOD Package"),
               uiOutput("download_met_package_ui"),
               br(),
               actionButton("back_to_landing_dispersion", "Back to Landing Page")
        )
      )
    }
  })
  
  # Baseline Map and Logic
  output$polygon_map_base <- renderLeaflet({
    req(rv$link_polygons)
    pal <- colorNumeric("viridis", domain = rv$link_polygons$Daily_Emissions)
    leaflet(rv$link_polygons) %>% addTiles() %>%
      addPolygons(color = ~pal(Daily_Emissions), weight = 2, label = ~paste("Emissions:", round(Daily_Emissions, 2), "g/day")) %>%
      leaflet::addLegend(pal = pal, values = ~Daily_Emissions, title = "Baseline Emissions (g/day)", position = "bottomright")
  })
  observe({
    req(rv$link_polygons)
    proxy <- leafletProxy("polygon_map_base") %>% clearGroup("receptors")
    if (!is.null(rv$near_road_receptors)) proxy %>% addCircleMarkers(data = rv$near_road_receptors, radius = 2, color = "red", group = "receptors")
    if (!is.null(rv$gridded_receptors)) proxy %>% addCircleMarkers(data = rv$gridded_receptors, radius = 2, color = "green", group = "receptors")
  })
  observeEvent(input$generate_receptors_base, {
    req(rv$link_polygons)
    res <- as.numeric(gsub("m", "", input$grid_resolution_base))
    bbox <- st_bbox(st_transform(rv$link_polygons, 32614))
    grid <- st_make_grid(st_as_sfc(bbox), cellsize = res, what = "centers")
    rv$gridded_receptors <- st_transform(st_sf(geometry = grid, crs = 32614), 4326)
  })
  
  # Alternate Map and Logic
  output$polygon_map_alt <- renderLeaflet({
    req(rv$link_polygons_alternate)
    pal <- colorNumeric("viridis", domain = rv$link_polygons_alternate$Daily_Emissions)
    leaflet(rv$link_polygons_alternate) %>% addTiles() %>%
      addPolygons(color = ~pal(Daily_Emissions), weight = 2, label = ~paste("Emissions:", round(Daily_Emissions, 2), "g/day")) %>%
      leaflet::addLegend(pal = pal, values = ~Daily_Emissions, title = "Alternate Emissions (g/day)", position = "bottomright")
  })
  observe({
    req(rv$link_polygons_alternate)
    proxy <- leafletProxy("polygon_map_alt") %>% clearGroup("receptors_alt")
    if (!is.null(rv$near_road_receptors_alternate)) proxy %>% addCircleMarkers(data = rv$near_road_receptors_alternate, radius = 2, color = "red", group = "receptors_alt")
    if (!is.null(rv$gridded_receptors_alternate)) proxy %>% addCircleMarkers(data = rv$gridded_receptors_alternate, radius = 2, color = "green", group = "receptors_alt")
  })
  observeEvent(input$generate_receptors_alt, {
    req(rv$link_polygons_alternate)
    res <- as.numeric(gsub("m", "", input$grid_resolution_alt))
    bbox <- st_bbox(st_transform(rv$link_polygons_alternate, 32614))
    grid <- st_make_grid(st_as_sfc(bbox), cellsize = res, what = "centers")
    rv$gridded_receptors_alternate <- st_transform(st_sf(geometry = grid, crs = 32614), 4326)
  })
  
  
  # --- Navigation Buttons ---
  observeEvent(input$goto_exposure_tab_btn, { updateTabsetPanel(session, "maintabs", selected = "Exposure Assessment") })
  observeEvent(input$goto_hia_tab_btn, { updateTabsetPanel(session, "maintabs", selected = "Health Impact Assessment") })
  
  # Conditional outputs for flows
  output$is_exposure_assessment_flow <- reactive({ rv$clicked_block == "Exposure Assessment" })
  output$is_hia_flow <- reactive({ rv$clicked_block == "Health Impact Assessment" })
  outputOptions(output, "is_exposure_assessment_flow", suspendWhenHidden = FALSE)
  outputOptions(output, "is_hia_flow", suspendWhenHidden = FALSE)
  
  # --- Visualize AERMOD Results Logic ---
  
  observeEvent(input$visualize_results, {
    with_modal({
      req(input$pst_file_base, input$interpolation_method)
      
      # Helper function for processing and interpolation
      process_and_interpolate <- function(file_path, resolution, method) {
        aermod_data <- read_pos(file_path)
        avg_conc <- aermod_data[, .(AVG_CONC = mean(AVERAGE_CONC, na.rm = TRUE)), by = .(X, Y)]
        points_sf <- st_as_sf(avg_conc, coords = c("X", "Y"), crs = 32614)
        bbox <- st_bbox(points_sf)
        grid <- st_make_grid(st_as_sfc(bbox), cellsize = as.numeric(resolution), what = "centers")
        
        if (method == "idw") {
          # Smooth (IDW) interpolation
          idw_gstat <- gstat(formula = AVG_CONC ~ 1, data = points_sf)
          interp_grid <- predict(idw_gstat, newdata = grid)
        } else {
          # Fast (Nearest Neighbor) interpolation
          # Using krige with nmax=1 and no variogram model performs nearest neighbor interpolation
          interp_grid <- krige(formula = AVG_CONC ~ 1, locations = points_sf, newdata = grid, nmax = 1)
        }
        
        stars_grid <- st_rasterize(interp_grid["var1.pred"])
        raster_grid <- as(stars_grid, "Raster")
        return(raster_grid)
      }
      
      # Process baseline
      rv$baseline_raster <- process_and_interpolate(input$pst_file_base$datapath, input$interpolation_res, input$interpolation_method)
      
      # Process alternate scenario if requested
      if (input$compare_scenarios_viz && !is.null(input$pst_file_alt)) {
        rv$alternate_raster <- process_and_interpolate(input$pst_file_alt$datapath, input$interpolation_res, input$interpolation_method)
        rv$alternate_raster <- resample(rv$alternate_raster, rv$baseline_raster)
        
        # Calculate percentage difference, handling potential division by zero
        rv$difference_raster <- ((rv$baseline_raster - rv$alternate_raster))
        rv$difference_raster[is.infinite(rv$difference_raster)] <- NA # Handle cases where baseline is 0
        
      }
      
      # Process data for boxplots
      base_data <- read_pos(input$pst_file_base$datapath) %>%
        mutate(date_time = ymd_h(paste0("20", substr(DATE, 1, 6), substr(DATE, 7, 8))),
               day_of_week = wday(date_time, label = TRUE),
               month = month(date_time, label = TRUE),
               hour = hour(date_time),
               Scenario = input$scenario_name_base)
      
      rv$is_hourly <- grepl("ALL1.pst", input$pst_file_base$name, ignore.case = TRUE)
      
      if (input$compare_scenarios_viz && !is.null(input$pst_file_alt)) {
        alt_data <- read_pos(input$pst_file_alt$datapath) %>%
          mutate(date_time = ymd_h(paste0("20", substr(DATE, 1, 6), substr(DATE, 7, 8))),
                 day_of_week = wday(date_time, label = TRUE),
                 month = month(date_time, label = TRUE),
                 hour = hour(date_time),
                 Scenario = input$scenario_name_alt)
        rv$processed_results_data <- rbind(base_data, alt_data)
        rv$processed_results_data$Scenario <- factor(rv$processed_results_data$Scenario, 
                                                     levels = c(input$scenario_name_base, input$scenario_name_alt))
      } else {
        rv$processed_results_data <- base_data
      }
      
    }, "Processing and Interpolating Results...")
  })
  
  output$results_visuals_ui <- renderUI({
    req(rv$processed_results_data)
    
    if(input$compare_scenarios_viz && !is.null(rv$difference_raster)) {
      # Scenario comparison layout
      tagList(
        fluidRow(
          column(6, withSpinner(leafletOutput("baseline_map_viz"))),
          column(6, withSpinner(leafletOutput("alternate_map_viz")))
        ),
        hr(),
        h4("Percentage Difference in Concentration"),
        withSpinner(leafletOutput("difference_map_viz")),
        hr(),
        h4("Concentration Box Plots"),
        plotOutput("day_of_week_plot"),
        plotOutput("month_plot"),
        if(rv$is_hourly) { plotOutput("hour_plot") }
      )
    } else {
      # Single scenario layout
      tagList(
        withSpinner(leafletOutput("baseline_map_viz")),
        hr(),
        h4("Concentration Box Plots"),
        plotOutput("day_of_week_plot"),
        plotOutput("month_plot"),
        if(rv$is_hourly) { plotOutput("hour_plot") }
      )
    }
  })
  
  # Baseline map for visualization tab
  output$baseline_map_viz <- renderLeaflet({
    req(rv$baseline_raster)
    raster_wgs84 <- projectRaster(rv$baseline_raster, crs = "EPSG:4326", method = "bilinear")
    pal <- colorNumeric("viridis", domain = values(raster_wgs84), na.color = "transparent", reverse = TRUE)
    leaflet() %>% addTiles() %>%
      addRasterImage(raster_wgs84, colors = pal, opacity = 0.7) %>%
      addLegend(pal = pal, values = values(raster_wgs84), title = paste(input$scenario_name_base, "Concentration"))
  })
  
  # Alternate map for visualization tab
  output$alternate_map_viz <- renderLeaflet({
    req(rv$alternate_raster)
    raster_wgs84 <- projectRaster(rv$alternate_raster, crs = "EPSG:4326", method = "bilinear")
    # Use the same domain as the baseline for consistent scale
    domain_vals <- values(projectRaster(rv$baseline_raster, crs = "EPSG:4326", method = "bilinear"))
    pal <- colorNumeric("viridis", domain = domain_vals, na.color = "transparent", reverse = TRUE)
    leaflet() %>% addTiles() %>%
      addRasterImage(raster_wgs84, colors = pal, opacity = 0.7) %>%
      addLegend(pal = pal, values = domain_vals, title = paste(input$scenario_name_alt, "Concentration"))
  })
  
  # Difference map for visualization tab
  output$difference_map_viz <- renderLeaflet({
    req(rv$difference_raster)
    raster_wgs84 <- projectRaster(rv$difference_raster, crs = "EPSG:4326", method = "bilinear")
    max_abs_val <- max(abs(values(raster_wgs84)), na.rm = TRUE)
    pal <- colorNumeric(c("blue", "white", "red"), domain = c(-max_abs_val, max_abs_val), na.color = "transparent")
    leaflet() %>% addTiles() %>%
      addRasterImage(raster_wgs84, colors = pal, opacity = 0.7) %>%
      addLegend(pal = pal, values = values(raster_wgs84), title = "Difference in Concentration")
  })
  
  # Box Plots
  output$day_of_week_plot <- renderPlot({
    req(rv$processed_results_data)
    p <- ggboxplot(rv$processed_results_data, x = "day_of_week", y = "AVERAGE_CONC", color = "black", scales = "free_x", palette = "jco", outlier.shape = 1, outlier.size = 0.1)
    if(input$compare_scenarios_viz) {
      p <- ggboxplot(rv$processed_results_data, x = "day_of_week", y = "AVERAGE_CONC", color = "black", fill="Scenario", scales = "free_x", palette = "jco", outlier.shape = 1, outlier.size = 0.1)+stat_compare_means(aes(group = Scenario), method = "t.test", label = "p.signif")
    }
    p + yscale("log10", .format = TRUE) +
      geom_vline(xintercept = seq(0.5, length(unique(rv$processed_results_data$day_of_week)) + 0.5, by = 1), linetype = "solid", color = "gray", size = 0.5) +
      labs(title = "Concentration by Day of Week", x = "Day of Week", y = "Concentration") +
      theme_minimal()
  })
  
  output$month_plot <- renderPlot({
    req(rv$processed_results_data)
    p <- ggboxplot(rv$processed_results_data, x = "month", y = "AVERAGE_CONC", color = "black", scales = "free_x", palette = "jco", outlier.shape = 1, outlier.size = 0.1)
    if(input$compare_scenarios_viz) {
      p <- ggboxplot(rv$processed_results_data, x = "month", y = "AVERAGE_CONC", color = "black", fill="Scenario", scales = "free_x", palette = "jco", outlier.shape = 1, outlier.size = 0.1)+ stat_compare_means(aes(group = Scenario), method = "t.test", label = "p.signif")
    }
    p + yscale("log10", .format = TRUE) +
      geom_vline(xintercept = seq(0.5, length(unique(rv$processed_results_data$month)) + 0.5, by = 1), linetype = "solid", color = "gray", size = 0.5) +
      labs(title = "Concentration by Month", x = "Month", y = "Concentration") +
      theme_minimal()
  })
  
  output$hour_plot <- renderPlot({
    req(rv$processed_results_data, rv$is_hourly)
    p <- ggboxplot(rv$processed_results_data, x = "hour", y = "AVERAGE_CONC", color = "black", scales = "free_x", palette = "jco", outlier.shape = 1, outlier.size = 0.1)
    if(input$compare_scenarios_viz) {
      p <- ggboxplot(rv$processed_results_data, x = "hour", y = "AVERAGE_CONC", color = "black", fill="Scenario", scales = "free_x", palette = "jco", outlier.shape = 1, outlier.size = 0.1)+ stat_compare_means(aes(group = Scenario), method = "t.test", label = "p.signif")
    }
    p + yscale("log10", .format = TRUE) +
      geom_vline(xintercept = seq(0.5, length(unique(rv$processed_results_data$hour)) + 0.5, by = 1), linetype = "solid", color = "gray", size = 0.5) +
      labs(title = "Concentration by Hour of Day", x = "Hour", y = "Concentration") +
      theme_minimal()
  })
  
  
  # --- Exposure Assessment Logic ---
  
  observeEvent(input$run_exposure_assessment, {
    if (is.null(rv$baseline_raster)) {
      showNotification("Please generate visualizations on the 'Visualize AERMOD Results' tab first.", type = "error", duration = 10)
      return()
    }
    
    with_modal({
      rv$pop_raster <- get_population_raster(rv$baseline_raster)
      req(rv$pop_raster)
      
      # Overall weighted exposures
      results_base <- calculate_weighted_exposure(rv$baseline_raster, rv$pop_raster)
      rv$exposure_conc_base <- results_base$weighted_avg
      rv$resampled_conc_raster_base <- results_base$conc_resampled
      rv$pop_raster_resampled <- results_base$pop_resampled
      
      if (!is.null(rv$alternate_raster)) {
        results_alt <- calculate_weighted_exposure(rv$alternate_raster, rv$pop_raster)
        rv$exposure_conc_alt <- results_alt$weighted_avg
        rv$resampled_conc_raster_alt <- results_alt$conc_resampled
      }
      
      # Aggregated weighted exposures per geo unit
      geo_data <- switch(input$exposure_geo,
                         "Tract" = tract_mr_sf,
                         "Block Group" = blkgp_mr_sf)
      
      # Crop geo_data to study area
      raster_bbox <- st_bbox(rv$baseline_raster)
      bbox_sf <- st_as_sfc(raster_bbox, crs = crs(rv$baseline_raster))
      bbox_sf <- st_transform(bbox_sf, st_crs(geo_data))
      geo_data_cropped <- geo_data[st_intersects(geo_data, bbox_sf, sparse = FALSE)[,1], ]
      geo_data_proj <- st_transform(geo_data_cropped, crs(rv$baseline_raster))
      
      # Calculate weighted exposure for base
      exposure_base <- numeric(nrow(geo_data_proj))
      for(i in 1:nrow(geo_data_proj)){
        poly <- geo_data_proj[i,]
        conc_masked <- mask(crop(rv$resampled_conc_raster_base, extent(poly)), poly)
        pop_masked <- mask(crop(rv$pop_raster_resampled, extent(poly)), poly)
        sum_product <- cellStats(conc_masked * pop_masked, 'sum', na.rm = TRUE)
        sum_pop <- cellStats(pop_masked, 'sum', na.rm = TRUE)
        exposure_base[i] <- if(sum_pop > 0) sum_product / sum_pop else NA
      }
      geo_data_proj$exposure_base <- exposure_base
      
      # If alternate
      if (!is.null(rv$alternate_raster)) {
        exposure_alt <- numeric(nrow(geo_data_proj))
        for(i in 1:nrow(geo_data_proj)){
          poly <- geo_data_proj[i,]
          conc_masked <- mask(crop(rv$resampled_conc_raster_alt, extent(poly)), poly)
          pop_masked <- mask(crop(rv$pop_raster_resampled, extent(poly)), poly)
          sum_product <- cellStats(conc_masked * pop_masked, 'sum', na.rm = TRUE)
          sum_pop <- cellStats(pop_masked, 'sum', na.rm = TRUE)
          exposure_alt[i] <- if(sum_pop > 0) sum_product / sum_pop else NA
        }
        geo_data_proj$exposure_alt <- exposure_alt
      }
      
      rv$exposure_by_geo <- geo_data_proj
      
    }, "Calculating Population Exposure...")
  })
  
  output$exposure_results_display_ui <- renderUI({
    req(rv$exposure_conc_base)
    
    # Value box style
    value_box_style <- "border: 1px solid #ddd; padding: 15px; border-radius: 5px; text-align: center; background-color: #f5f5f5;"
    
    tagList(
      if (is.null(rv$exposure_conc_alt)) {
        fluidRow(
          column(12, style = value_box_style,
                 h4("Population-Weighted Exposure Concentration"),
                 h3(sprintf("%.4f µg/m³", rv$exposure_conc_base))
          )
        )
      } else {
        fluidRow(
          column(6, style = value_box_style,
                 h4(paste(input$scenario_name_base, "Exposure")),
                 h3(sprintf("%.4f µg/m³", rv$exposure_conc_base))
          ),
          column(6, style = value_box_style,
                 h4(paste(input$scenario_name_alt, "Exposure")),
                 h3(sprintf("%.4f µg/m³", rv$exposure_conc_alt))
          )
        )
      },
      hr(),
      if (is.null(rv$exposure_conc_alt)) {
        fluidRow(
          column(6, h4("Concentration (30m)"), leafletOutput("exposure_conc_map_base")),
          column(6, h4("Population (30m)"), leafletOutput("exposure_pop_map"))
        )
      } else {
        fluidRow(
          column(4, h4(paste(input$scenario_name_base, "Concentration (30m)")), leafletOutput("exposure_conc_map_base")),
          column(4, h4(paste(input$scenario_name_alt, "Concentration (30m)")), leafletOutput("exposure_conc_map_alt")),
          column(4, h4("Population (30m)"), leafletOutput("exposure_pop_map"))
        )
      },
      hr(),
      h4(paste("Population-Weighted Mean Exposure by", input$exposure_geo)),
      if (is.null(rv$exposure_conc_alt)) {
        leafletOutput("exposure_map_by_geo_base")
      } else {
        fluidRow(
          column(6, leafletOutput("exposure_map_by_geo_base")),
          column(6, leafletOutput("exposure_map_by_geo_alt"))
        )
      }
    )
  })
  
  # Leaflet map for baseline resampled concentration
  output$exposure_conc_map_base <- renderLeaflet({
    req(rv$resampled_conc_raster_base)
    raster_wgs84 <- projectRaster(rv$resampled_conc_raster_base, crs = "EPSG:4326", method = "bilinear")
    pal <- colorNumeric("viridis", domain = values(raster_wgs84), na.color = "transparent", reverse = TRUE)
    leaflet() %>% addTiles() %>%
      addRasterImage(raster_wgs84, colors = pal, opacity = 0.7) %>%
      addLegend(pal = pal, values = values(raster_wgs84), title = "Concentration")
  })
  
  # Leaflet map for alternate resampled concentration
  output$exposure_conc_map_alt <- renderLeaflet({
    req(rv$resampled_conc_raster_alt)
    raster_wgs84 <- projectRaster(rv$resampled_conc_raster_alt, crs = "EPSG:4326", method = "bilinear")
    pal <- colorNumeric("viridis", domain = values(raster_wgs84), na.color = "transparent", reverse = TRUE)
    leaflet() %>% addTiles() %>%
      addRasterImage(raster_wgs84, colors = pal, opacity = 0.7) %>%
      addLegend(pal = pal, values = values(raster_wgs84), title = "Concentration")
  })
  
  # Leaflet map for population
  output$exposure_pop_map <- renderLeaflet({
    req(rv$pop_raster_resampled)
    raster_wgs84 <- projectRaster(rv$pop_raster_resampled, crs = "EPSG:4326", method = "bilinear")
    pal <- colorNumeric("YlGnBu", domain = values(raster_wgs84), na.color = "transparent")
    leaflet() %>% addTiles() %>%
      addRasterImage(raster_wgs84, colors = pal, opacity = 0.7) %>%
      addLegend(pal = pal, values = values(raster_wgs84), title = "Population")
  })
  
  # Aggregated exposure map base
  output$exposure_map_by_geo_base <- renderLeaflet({
    req(rv$exposure_by_geo)
    map_data <- st_transform(rv$exposure_by_geo, 4326)
    pal <- colorNumeric("viridis", domain = map_data$exposure_base, na.color = "transparent")
    
    leaflet(map_data) %>%
      addTiles() %>%
      addPolygons(fillColor = ~pal(exposure_base),
                  fillOpacity = 0.7,
                  color = "#BDBDC3",
                  weight = 1,
                  label = ~paste0(GEOID, ": ", round(exposure_base, 4), " µg/m³")) %>%
      addLegend(pal = pal, values = ~exposure_base, title = "Mean Exposure (µg/m³)",
                position = "bottomright")
  })
  
  # Aggregated exposure map alt
  output$exposure_map_by_geo_alt <- renderLeaflet({
    req(rv$exposure_by_geo, rv$exposure_conc_alt)
    map_data <- st_transform(rv$exposure_by_geo, 4326)
    pal <- colorNumeric("viridis", domain = map_data$exposure_alt, na.color = "transparent")
    
    leaflet(map_data) %>%
      addTiles() %>%
      addPolygons(fillColor = ~pal(exposure_alt),
                  fillOpacity = 0.7,
                  color = "#BDBDC3",
                  weight = 1,
                  label = ~paste0(GEOID, ": ", round(exposure_alt, 4), " µg/m³")) %>%
      addLegend(pal = pal, values = ~exposure_alt, title = "Mean Exposure (µg/m³)",
                position = "bottomright")
  })
  
  # --- Health Impact Assessment (HIA) Logic ---
  
  # Dynamic UI for RR selection
  output$hia_rr_ui <- renderUI({
    req(input$hia_pollutant)
    choices <- rr_values[[input$hia_pollutant]]
    if (any(is.na(choices))) {
      p("HIA for long-term CO exposure is not recommended due to scientific uncertainty.")
    } else {
      selectInput("hia_rr_value", "Select Relative Risk (RR) per 10 µg/m³:",
                  choices = choices)
    }
  })
  
  # Main HIA calculation observer
  observeEvent(input$run_hia, {
    req(rv$baseline_raster, input$hia_pollutant, input$hia_rr_value, input$hia_mortality_geo)
    
    if(is.na(input$hia_rr_value)) {
      showNotification("Cannot run HIA with the selected pollutant.", type = "warning")
      return()
    }
    
    with_modal({
      # Ensure pop and resampled conc are available
      if (is.null(rv$pop_raster_resampled)) {
        rv$pop_raster <- get_population_raster(rv$baseline_raster)
        results_base <- calculate_weighted_exposure(rv$baseline_raster, rv$pop_raster)
        rv$resampled_conc_raster_base <- results_base$conc_resampled
        rv$pop_raster_resampled <- results_base$pop_resampled
        if (!is.null(rv$alternate_raster)) {
          results_alt <- calculate_weighted_exposure(rv$alternate_raster, rv$pop_raster)
          rv$resampled_conc_raster_alt <- results_alt$conc_resampled
        }
      }
      
      # Select geo data
      geo_data <- switch(input$hia_mortality_geo,
                         "County" = county_mr_sf,
                         "Tract" = tract_mr_sf,
                         "Block Group" = blkgp_mr_sf)
      
      # Crop to study area
      raster_bbox <- st_bbox(rv$baseline_raster)
      bbox_sf <- st_as_sfc(raster_bbox, crs = crs(rv$baseline_raster))
      bbox_sf <- st_transform(bbox_sf, st_crs(geo_data))
      geo_data_cropped <- geo_data[st_intersects(geo_data, bbox_sf, sparse = FALSE)[,1], ]
      geo_data_proj <- st_transform(geo_data_cropped, crs(rv$baseline_raster))
      
      # Calculate weighted conc for base
      base_conc <- numeric(nrow(geo_data_proj))
      for(i in 1:nrow(geo_data_proj)){
        poly <- geo_data_proj[i,]
        conc_masked <- mask(crop(rv$resampled_conc_raster_base, extent(poly)), poly)
        pop_masked <- mask(crop(rv$pop_raster_resampled, extent(poly)), poly)
        sum_product <- cellStats(conc_masked * pop_masked, 'sum', na.rm = TRUE)
        sum_pop <- cellStats(pop_masked, 'sum', na.rm = TRUE)
        base_conc[i] <- if(sum_pop > 0) sum_product / sum_pop else NA
      }
      geo_data_proj$base_conc <- base_conc
      
      # If alternate
      if (!is.null(rv$alternate_raster)) {
        alt_conc <- numeric(nrow(geo_data_proj))
        for(i in 1:nrow(geo_data_proj)){
          poly <- geo_data_proj[i,]
          conc_masked <- mask(crop(rv$resampled_conc_raster_alt, extent(poly)), poly)
          pop_masked <- mask(crop(rv$pop_raster_resampled, extent(poly)), poly)
          sum_product <- cellStats(conc_masked * pop_masked, 'sum', na.rm = TRUE)
          sum_pop <- cellStats(pop_masked, 'sum', na.rm = TRUE)
          alt_conc[i] <- if(sum_pop > 0) sum_product / sum_pop else NA
        }
        geo_data_proj$alt_conc <- alt_conc
        delta_x <- geo_data_proj$base_conc - geo_data_proj$alt_conc
      } else {
        delta_x <- geo_data_proj$base_conc
      }
      
      # Calculate Beta from RR
      beta <- log(as.numeric(input$hia_rr_value)) / 10
      
      # Health Impact Function
      geo_data_proj$health_impact <- geo_data_proj$Mortality_Rate * (1 - exp(-beta * delta_x)) * geo_data_proj$Total_Pop
      rv$total_hia_impact <- sum(geo_data_proj$health_impact, na.rm = TRUE)
      rv$hia_results <- geo_data_proj
      
    }, "Running Health Impact Assessment...")
  })
  
  # UI to display HIA results
  output$hia_results_display_ui <- renderUI({
    req(rv$hia_results, rv$total_hia_impact)
    
    impact_title <- if (!is.null(rv$alternate_raster)) "Annual Deaths Avoided" else "Annual Attributable Deaths"
    
    tagList(
      leafletOutput("hia_map"),
      div(class = "total-impact-box",
          h3(impact_title),
          h2(round(rv$total_hia_impact, 2))
      )
    )
  })
  
  # Leaflet map for HIA results
  output$hia_map <- renderLeaflet({
    req(rv$hia_results)
    map_data <- st_transform(rv$hia_results, 4326)
    
    # Handle potential Inf/-Inf values for color palette
    finite_impacts <- map_data$health_impact[is.finite(map_data$health_impact)]
    if (length(finite_impacts) == 0) return()
    
    pal_domain <- range(finite_impacts, na.rm = TRUE)
    
    # Use a diverging palette for avoided deaths, sequential for attributable
    pal_colors <- if (!is.null(rv$alternate_raster)) "RdYlGn" else "YlOrRd"
    pal <- colorNumeric(pal_colors, domain = pal_domain, na.color = "transparent")
    
    leaflet(map_data) %>%
      addTiles() %>%
      addPolygons(fillColor = ~pal(health_impact),
                  fillOpacity = 0.8,
                  color = "#BDBDC3",
                  weight = 1,
                  label = ~paste0(GEOID, ": ", round(health_impact, 3))) %>%
      addLegend(pal = pal, values = ~health_impact, title = "Health Impact",
                position = "bottomright")
  })
  
}

# Launch the app
shinyApp(ui, server)