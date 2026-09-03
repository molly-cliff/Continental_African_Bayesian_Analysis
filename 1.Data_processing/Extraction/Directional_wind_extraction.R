library(raster)
library(terra)
library(dplyr)
library(stringr)
library(sf)
library(lubridate)

testing_extract_datameridonal_zonal <- function(raster_file) {
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
  
  # Add unique code for merging
  final_data <- weekly_data %>%
    mutate(
      year_month = paste(year, sprintf("%02d", month), sep = " "),
      code = paste(year_month, district_country, sep = " "),
      eastward_wind = NA_real_,
      north_wind = NA_real_
    )
  
  # Read raster files
  path <- "meriodonal_zonal"
  files <- list.files(path, pattern = "\\.nc$", full.names = TRUE)
  r <- rast(files[1])
  
  # Split raster into two sets
  r1 <- r[[1:240]]  # Eastward wind
  r2 <- r[[241:480]]  # Northward wind
  
  # Process first raster set (eastward)
  r1 <- rotate(r1)
  if (!crs(r1) == crs(shape2)) {
    message("CRS mismatch - reprojecting r1")
    r1 <- project(r1, crs(shape2))
  } else {
    message("CRS aligned (r1)")
  }
  r1 <- mask(crop(r1, shape2), shape2)
  
  start_date1 <- ymd("2003-01-01")
  date_seq1 <- seq(start_date1, by = "1 month", length.out = nlyr(r1))
  names(r1) <- paste0("windspeed_", year(date_seq1), "_", sprintf("%02d", month(date_seq1)))
  
  for (i in 1:nlyr(r1)) {
    lyr <- r1[[i]]
    lyr_name <- names(r1)[i]
    cat("Processing eastward wind:", lyr_name, "\n")
    
    ym_parts <- str_match(lyr_name, "windspeed_(\\d{4})_(\\d{2})")
    yr <- as.integer(ym_parts[, 2])
    mo <- as.integer(ym_parts[, 3])
    
    extracted_vals <- terra::extract(lyr, shape2, fun = mean, na.rm = TRUE)[, 2]
    
    extracted_df <- shape2 %>%
      st_drop_geometry() %>%
      mutate(
        windspeed_temp = extracted_vals,
        year = yr,
        month = mo,
        year_month = paste(yr, sprintf("%02d", mo), sep = " "),
        district_country = paste(NAME_2, COUNTRY, sep = " "),
        code = paste(year_month, district_country, sep = " ")
      ) %>%
      select(code, windspeed_temp)
    
    final_data <- left_join(final_data, extracted_df, by = "code")
    update_idx <- is.na(final_data$eastward_wind) & !is.na(final_data$windspeed_temp)
    final_data$eastward_wind[update_idx] <- final_data$windspeed_temp[update_idx]
    final_data$windspeed_temp <- NULL
  }
  
  # Process second raster set (northward)
  r2 <- rotate(r2)
  if (!crs(r2) == crs(shape2)) {
    message("CRS mismatch - reprojecting r2")
    r2 <- project(r2, crs(shape2))
  } else {
    message("CRS aligned (r2)")
  }
  r2 <- mask(crop(r2, shape2), shape2)
  
  start_date2 <- ymd("2003-01-01")
  date_seq2 <- seq(start_date2, by = "1 month", length.out = nlyr(r2))
  names(r2) <- paste0("windspeed_", year(date_seq2), "_", sprintf("%02d", month(date_seq2)))
  
  for (i in 1:nlyr(r2)) {
    lyr <- r2[[i]]
    lyr_name <- names(r2)[i]
    cat("Processing northward wind:", lyr_name, "\n")
    
    ym_parts <- str_match(lyr_name, "windspeed_(\\d{4})_(\\d{2})")
    yr <- as.integer(ym_parts[, 2])
    mo <- as.integer(ym_parts[, 3])
    
    extracted_vals <- terra::extract(lyr, shape2, fun = mean, na.rm = TRUE)[, 2]
    
    extracted_df <- shape2 %>%
      st_drop_geometry() %>%
      mutate(
        windspeed_temp = extracted_vals,
        year = yr,
        month = mo,
        year_month = paste(yr, sprintf("%02d", mo), sep = " "),
        district_country = paste(NAME_2, COUNTRY, sep = " "),
        code = paste(year_month, district_country, sep = " ")
      ) %>%
      select(code, windspeed_temp)
    
    final_data <- left_join(final_data, extracted_df, by = "code")
    update_idx <- is.na(final_data$north_wind) & !is.na(final_data$windspeed_temp)
    final_data$north_wind[update_idx] <- final_data$windspeed_temp[update_idx]
    final_data$windspeed_temp <- NULL
  }
  
  # Final preview
  print(head(final_data))
  
  # Write output to CSV
  write.csv(final_data, "final_meridonal_zonal_data.csv", row.names = FALSE)
  
  return(final_data)
}
