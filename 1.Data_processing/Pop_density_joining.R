library(raster)
library(terra)
library(dplyr)
library(stringr)
library(sf)
library(lubridate)
library(sf)
library(terra)
library(dplyr)
library(stringr)

popdensityextract <- function() {

  # Load input data
  spatiotemporaloutbreaks <- read.csv("joined_outbreak_vaccine_data.csv")
  cat("Loaded input data with columns:\n")
  print(names(spatiotemporaloutbreaks))
  
  # Create continuous numeric year_month variable
  spatiotemporaloutbreaks$year_month <- 
    (spatiotemporaloutbreaks$year.x - min(spatiotemporaloutbreaks$year.x, na.rm = TRUE)) * 12 + 
    spatiotemporaloutbreaks$month
  
  # Select relevant variables
  myvars <- c("year.x", "month", "outbreak_2", "name_2", "district_country.x", "country",
              "code.x", "rainfall", "eastward_wind", "north_wind", "aod", "humidity", 
              "vaccine", "windspeed", "temp", "cropland", "forest", "barren", "year_month")
  
  final_data <- spatiotemporaloutbreaks[, myvars]
  
  # Load shapefile
  shape2 <- st_read("Shapefile_improved.shp")
  cat("Loaded shapefile with", nrow(shape2), "features\n")
  
  # Drop unwanted countries
  countries_to_drop <- c("Cabo Verde", "Mauritius", "Seychelles")
  shape2 <- shape2[!shape2$COUNTRY %in% countries_to_drop, ]
  cat("After dropping countries, shapefile has", nrow(shape2), "features\n")
  
  # List all .tif files in pop_density folder
  pop_density_folder <- file.path(getwd(), "pop_density")
  files <- list.files(pop_density_folder, pattern = "\\.tif$", full.names = TRUE)
  cat("Population density raster files found:\n")
  print(files)
  
  # Initialize pop_density column in final_data
  final_data$pop_density <- NA_real_
  
  # Loop through population density rasters and extract values
  for (f in files) {
    cat("Processing file:", f, "\n")
    
    r <- rast(f)  # Load raster using terra
    file_name <- basename(f)
    
    # Extract year from filename (assumes 4-digit year in filename)
    year <- as.numeric(str_extract(file_name, "\\d{4}"))
    if (is.na(year)) {
      warning("Could not extract year from file name:", file_name)
      next
    }
    
    # Extract mean population density per polygon
    extracted <- terra::extract(r, vect(shape2), fun = max, na.rm = TRUE)
    
    # Add extracted pop_density to shapefile data
    shape2$pop_density <- extracted[, 2]  # second column is extracted value
    
    # Prepare temporary dataframe to merge with final_data
    temp_df <- shape2 %>%
      st_drop_geometry() %>%
      select(name_2 = NAME_2, country = COUNTRY, pop_density) %>%
      mutate(year.x = year)
    
    # Left join to final_data by name_2, country, and year.x
    final_data <- final_data %>%
      left_join(temp_df, by = c("name_2", "country", "year.x"), suffix = c("", "_new")) %>%
      mutate(pop_density = ifelse(is.na(pop_density), pop_density_new, pop_density)) %>%
      select(-pop_density_new)
  }
  
  # Export final dataset
  write.csv(final_data, "final_popdensity_vaccine.csv", row.names = FALSE)
  cat("Final dataset saved as final_popdensity_vaccine.csv\n")
}
