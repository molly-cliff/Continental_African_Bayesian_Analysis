
library(raster)
library(terra)
library(dplyr)
library(stringr)
library(sf)
library(lubridate)  # For ymd, year, month

testing_extract_hum <- function() {
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
path <- "humidity"
files <- list.files(path, pattern = "\\.nc$", full.names = TRUE)
r <- rast(files[1])

print(names(r))
r <- rotate(r)

r2<- r

# Reproject if CRS mismatch
if (!crs(r2) == crs(shape2)) {
  message("CRS mismatch - reprojecting")
  r2 <- project(r2, crs(shape2))
} else {
  message("CRS aligned")
}

# Crop and mask raster
r2 <- mask(crop(r2, shape2), shape2)
plot(r2[[1]])
# Create time-based names for layers
start_date <- ymd("2003-01-01")
date_seq <- seq(start_date, by = "1 month", length.out = nlyr(r2))
layer2_names <- paste0("humidity_", year(date_seq), "_", sprintf("%02d", month(date_seq)))
names(r2) <- layer2_names

print(names(r2))


final_data <- final_data %>%
  mutate(
    year_month = paste(year, sprintf("%02d", month), sep = " "),
    code = paste(year_month, district_country, sep = " "),
    humidity = NA  
  )

# Loop through each raster layer and extract humidity data
for (i in 1:nlyr(r2)) {
  lyr <- r2[[i]]
  lyr_name <- names(r2)[i]
  cat("Processing:", lyr_name, "\n")
  
  ym_parts <- str_match(lyr_name, "humidity_(\\d{4})_(\\d{2})")
  year <- as.integer(ym_parts[, 2])
  month <- as.integer(ym_parts[, 3])
  
  extracted <- data.frame(
    shape2,
    humidity_temp = terra::extract(lyr, shape2, fun = mean, na.rm = TRUE)[, 2]
  ) %>%
    mutate(
      year = year,
      month = month,
      year_month = paste(year, sprintf("%02d", month), sep = " "),
      district_country = paste(NAME_2, COUNTRY, sep = " "),
      code = paste(year_month, district_country, sep = " ")
    )
  
  # Merge with final_data and update humidity_test2 where missing
  final_data <- merge(final_data, extracted[, c("code", "humidity_temp")], by = "code", all.x = TRUE)
  final_data$humidity[is.na(final_data$humidity)] <- final_data$humidity_temp[is.na(final_data$humidity)]
  final_data$humidity_temp <- NULL
}

# Final preview
print(head(final_data))


na_humidity <- final_data[is.na(final_data$humidity), ]
print(na_humidity)

# Write output
write.csv(final_data, "final_humidity_data.csv", row.names = FALSE)
}