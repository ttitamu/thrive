# Use a specific, stable version of the Rocker image
FROM rocker/shiny:4.3.2

# Install system dependencies required for R packages like sf, raster, etc.
# This is a crucial step for spatial packages.
RUN apt-get update && apt-get install -y \
    git-core \
	libgdal-dev \
    libproj-dev \
	libjq1 \
    libgeos-dev \
    libudunits2-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy the package installation script into the container first
COPY install_packages.R /

# Run the package installation script. This is more efficient than multiple RUN commands.
RUN Rscript /install_packages.R

# Copy all other application files (app.R, data/, www/) into the server directory
COPY . /srv/shiny-server/


# Expose the default Shiny Server port
EXPOSE 3838

# Define the command to run the Shiny Server
CMD ["/usr/bin/shiny-server"]
