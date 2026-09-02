library(raster)
library(terra)
library(dplyr)
library(stringr)
library(sf)

testing_extract_data <- function() {
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
  weekly_data$year_month <- paste(weekly_data$year, weekly_data$month, sep = " ")
  weekly_data$code <- paste(weekly_data$year_month, weekly_data$district_country, sep = " ")
  
  # Prepare file paths
  subfolder <- "CHIRPS_data_monthly"
  path <- file.path(getwd(), subfolder)
  files <- list.files(path, pattern = "tif$", full.names = TRUE)
  print(files)
  
  # Initialize final data frame
  final_data <- weekly_data
  final_data$rainfall <- NA
  
  # Loop through rainfall files
  for (f in files) {
    cat("Processing file:", f, "\n")
    r <- raster(f)
    file_name <- basename(f)
    year <- str_extract(file_name, "\\d{4}")
    month <- str_extract(file_name, "(?<=\\d{4}\\.)\\d{2}")
    year_month_label <- paste0("rainfall_", year, "_", month)
    
    # Extract and format rainfall data
    extracted_data <- data.frame(
      shape2,
      rainfall = terra::extract(r, shape2, fun = mean, na.rm = TRUE)
    ) %>%
      mutate(
        year = as.integer(year),
        month = as.integer(month),
        year_month = paste(year, month, sep = " "),
        district_country = paste(NAME_2, COUNTRY, sep = " "),
        code = paste(year_month, district_country, sep = " ")
      )
    
    extracted_data <- extracted_data[, c("code", "rainfall")]
    names(extracted_data)[2] <- year_month_label
    
    # Merge with final data
    final_data <- merge(final_data, extracted_data, by = "code", all.x = TRUE)
    new_vals <- final_data[[year_month_label]]
    final_data$rainfall[is.na(final_data$rainfall) & !is.na(new_vals)] <- new_vals[is.na(final_data$rainfall) & !is.na(new_vals)]
    final_data[[year_month_label]] <- NULL
  }
  
  print(head(final_data))
  
  # Check for NAs
  na_in_X_based_on_Y <- final_data[is.na(final_data$rainfall), ]
  print(na_in_X_based_on_Y)
  
  # Write output to CSV
  write.csv(final_data, "final_rainfall_data.csv", row.names = FALSE)
}
