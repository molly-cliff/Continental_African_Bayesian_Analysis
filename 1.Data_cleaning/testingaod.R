library(raster)
library(terra)
library(dplyr)
library(stringr)
library(sf)
library(lubridate)

testing_extract_aod <- function() {
  # Load input data
  weekly_data <- read.csv("full-outbreak-matched.csv")
  shape2 <- st_read("Shapefile_improved.shp")
  
  # Drop unwanted countries
  countries_to_drop <- c("Cabo Verde", "Mauritius", "Seychelles")
  shape2 <- shape2[!shape2$COUNTRY %in% countries_to_drop, ]
  
  # Summarize outbreak data
  weekly_data <- weekly_data %>%
    group_by(year, month, district_country) %>%
    summarise(outbreak_occur = as.integer(any(outbreak > 0, na.rm = TRUE)), .groups = 'drop')
  print(head(weekly_data))
  
  # Add unique code for merging
  final_data <- weekly_data %>%
    mutate(
      year_month = paste(year, sprintf("%02d", month), sep = " "),
      code = paste(year_month, district_country, sep = " ")
    )
  
  # Set file path and read in NetCDF files
  path <- "aod"
  nc_files <- list.files(path, pattern = "\\.nc4$", full.names = TRUE)
  
print(nc_files)
# Initialize an empty list to hold the monthly rasters
monthly_rasters <- list()

# Loop through each file
for (file in nc_files) {
  message("Processing: ", basename(file))
  
  r <- rast(file)
  
  # Select only AODANA layers
  aod_layers <- r[[grep("^AODANA", names(r))]]
  
  # Average across the 8 time steps
  monthly_mean <- mean(aod_layers)
  
  # Optional: name it by date (e.g., 2003_01)
  date_str <- gsub(".*\\.(\\d{6})\\.nc4$", "\\1", basename(file))  # e.g., "200301"
  names(monthly_mean) <- paste0("AOD_", date_str)
  
  # Optional: Save individual GeoTIFFs
  writeRaster(monthly_mean, paste0("AOD_", date_str, ".tif"), overwrite = TRUE)
  
  # Store in list
  monthly_rasters[[date_str]] <- monthly_mean
}

# Stack all monthly rasters
aod_stack <- rast(monthly_rasters)

r_stack<-aod_stack

print(r_stack)
plot(r_stack[[1]])  # First month's AOD



r2<-r_stack
# Reproject if CRS mismatch
if (!crs(r2) == crs(shape2)) {
  message("CRS mismatch - reprojecting")
  r2 <- project(r2, crs(shape2))
} else {
  message("CRS aligned")
}

plot(r2[[1]])



# Crop and mask the raster stack
r2 <- mask(crop(r2, shape2), shape2)

# Plot first few layers for sanity check
plot(r2[[1]])
plot(r2[[2]])

# Prepare final_data for merging
final_data <- weekly_data %>%
  mutate(
    year_month = paste(year, sprintf("%02d", month), sep = " "),
    code = paste(year_month, district_country, sep = " "),
    aod = NA  # initialize
  )

# Loop through each raster layer
for (i in 1:nlyr(r2)) {
  lyr <- r2[[i]]
  lyr_name <- names(r2)[i]
  cat("Processing:", lyr_name, "\n")
  
  # Extract year and month from layer name like "AOD_200301"
  ym_parts <- str_match(lyr_name, "(\\d{4})(\\d{2})")
  year <- as.integer(ym_parts[, 2])
  month <- as.integer(ym_parts[, 3])
  
  # Extract mean AOD for each district in shape2
  extracted <- data.frame(
    shape2,
    aod_temp = terra::extract(lyr, shape2, fun = mean, na.rm = TRUE)[, 2]
  ) %>%
    mutate(
      year = year,
      month = month,
      year_month = paste(year, sprintf("%02d", month), sep = " "),
      district_country = paste(NAME_2, COUNTRY, sep = " "),
      code = paste(year_month, district_country, sep = " ")
    )
  
  # Merge with final_data and fill missing aod
  final_data <- merge(final_data, extracted[, c("code", "aod_temp")], by = "code", all.x = TRUE)
  final_data$aod[is.na(final_data$aod)] <- final_data$aod_temp[is.na(final_data$aod)]
  final_data$aod_temp <- NULL
}

print(head(final_data))


na_aod <- final_data[is.na(final_data$aod), ]
print(na_aod)

# Write output
write.csv(final_data, "final_aod_data.csv", row.names = FALSE)
}
