library(raster)
library(terra)
library(dplyr)
library(stringr)
library(sf)
library(lubridate)  # For ymd, year, month

testing_extract_datawind <- function() {
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
  weekly_data <- weekly_data %>%
    mutate(
      year_month = paste(year, sprintf("%02d", month), sep = " "),
      code = paste(year_month, district_country, sep = " ")
    )
  
  # Prepare file paths
  subfolder <- "windspeed"
  path <- file.path(getwd(), subfolder)
  files <- list.files(path, pattern = "\\.nc$", full.names = TRUE)
  
  r <- rast(files[1])
  print(names(r))
  r <- rotate(r)
  
  if (!crs(r) == crs(shape2)) {
    message("CRS mismatch - Reprojecting")
    r <- project(r, crs(shape2))
  } else {
    message("CRS aligned")
  }
  
  plot(r[[1]], main = "Rotated raster")
  
  r <- mask(crop(r, shape2), shape2)
  start_date <- ymd("2003-01-01") 
  date_seq <- seq(start_date, by = "1 month", length.out = nlyr(r))
  layer_names <- paste0("windspeed_", year(date_seq), "_", sprintf("%02d", month(date_seq)))
  names(r) <- layer_names
  print(names(r))
  
  # Filter by time
  start_filter <- ymd("2003-01-01")
  end_filter <- ymd("2022-12-31")
  keep_indices <- which(date_seq >= start_filter & date_seq <= end_filter)
  r <- r[[keep_indices]]
  
  # Initialize final data
  final_data <- weekly_data %>%
    mutate(windspeed = NA)
  
  # Loop through raster layers
  for (i in 1:nlyr(r)) {
    lyr <- r[[i]]
    lyr_name <- names(r)[i]
    cat("Processing:", lyr_name, "\n")
    ym_parts <- str_match(lyr_name, "windspeed_(\\d{4})_(\\d{2})")
    lyr_year <- as.integer(ym_parts[, 2])
    lyr_month <- as.integer(ym_parts[, 3])
    
    extracted <- data.frame(
      shape2,
      windspeed_temp = terra::extract(lyr, shape2, fun = mean, na.rm = TRUE)[, 2]
    )
    
    extracted <- extracted %>%
      mutate(
        year = lyr_year,
        month = lyr_month,
        year_month = paste(year, sprintf("%02d", month), sep = " "),
        district_country = paste(NAME_2, COUNTRY, sep = " "),
        code = paste(year_month, district_country, sep = " ")
      )
    
    final_data <- merge(final_data, extracted[, c("code", "windspeed_temp")], by = "code", all.x = TRUE)
    final_data$windspeed[is.na(final_data$windspeed)] <- final_data$windspeed_temp[is.na(final_data$windspeed)]
    final_data$windspeed_temp <- NULL
  }
  
  # Check missing windspeed
  na_windspeed <- final_data[is.na(final_data$windspeed), ]
  print(na_windspeed)
  
  # Write output
  write.csv(final_data, "final_windspeed_data.csv", row.names = FALSE)
}
