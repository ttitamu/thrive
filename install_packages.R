# install_packages.R

# Install the 'remotes' package, which is needed to install from GitHub
install.packages("remotes", repos = "https://packagemanager.posit.co/cran/__linux__/jammy/latest")

# Install esri2sf from GitHub
remotes::install_github("yonghah/esri2sf")

# List of all other packages from CRAN
packages <- c(
  "shiny", "shinyWidgets", "leaflet", "sf", "tigris", "leaflet.extras",
  "geojsonio", "jsonlite", "ggplot2", "ggpubr", "dplyr", "tidyr", "zip",
  "stringr", "shinycssloaders", "raster", "gstat", "stars",
  "data.table", "lubridate", "httr"
)

# Install the rest of the packages from the Posit repository
install.packages(packages, repos = "https://packagemanager.posit.co/cran/__linux__/jammy/latest")