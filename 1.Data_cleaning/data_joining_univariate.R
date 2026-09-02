library(dplyr)
library(INLA)
library(sf)
library(spdep)
library(sp)
library(caret)
library(pROC)
library(tidyr)
library(broom)
library(corrplot)
library(car)
library(ggplot2)
library(reshape2)

data_join <-function() {
  library(dplyr)
  dataset_files <- list(
    rainfall = "final_rainfall_data.csv",
    meridonal_zonal = "final_meridonal_zonal_data.csv",
    aod = "final_aod_data.csv",
    humidity = "final_humidity_data.csv",
    windspeed = "final_windspeed_data.csv",
    temp = "final_temp_data.csv"
  )
  
  # Load outbreak data
  spatiotemporaloutbreaks <- read.csv("full-outbreak-matched.csv", stringsAsFactors = FALSE)
  
  # Check columns in outbreak data
  print("Column names in outbreak data:")
  print(names(spatiotemporaloutbreaks))
  
  # Select and summarise outbreak data
  myvars <- c("week", "month", "outbreak", "year", "name_2", "district_country", "country")
  newdata <- spatiotemporaloutbreaks[, myvars]
  
  bf1 <- newdata %>%
    group_by(year, month, name_2) %>%
    summarise(
      outbreak_2 = as.integer(any(outbreak > 0, na.rm = TRUE)),
      district_country = first(district_country),
      country = first(country),
      .groups = 'drop'
    )
  
  print(table(bf1$country, bf1$outbreak_2))
  
  # Create padded code in outbreak data
  bf1$year_month <- paste(bf1$year, sprintf("%02d", bf1$month), sep = " ")
  bf1$code <- paste(bf1$year_month, bf1$district_country, sep = " ")
  print(names(bf1))
  
  # Start with outbreak summary as the base
  merged_data <- bf1
  
  # Merge each environmental dataset
  for (name in names(dataset_files)) {
    env_data <- read.csv(dataset_files[[name]], stringsAsFactors = FALSE)
    
    # If the env dataset has year/month/district_country, rebuild the code correctly
    if (all(c("year", "month", "district_country") %in% names(env_data))) {
      env_data$year_month <- paste(env_data$year, sprintf("%02d", env_data$month), sep = " ")
      env_data$code <- paste(env_data$year_month, env_data$district_country, sep = " ")
    }
    
    # Deduplicate by code
    if ("code" %in% names(env_data)) {
      env_data <- env_data %>% distinct(code, .keep_all = TRUE)
    } else {
      warning(paste("Column 'code' not found in", name, "dataset. Deduplication skipped."))
    }
    
    # Left join onto merged_data
    merged_data <- merged_data %>%
      left_join(env_data, by = "code")
  }
  
  # Final column names
  print("Final column names after merging:")
  print(names(merged_data))
  print(table(merged_data$outbreak_2, merged_data$country))
  
  # Select & rename useful columns
  myvars <- c("year", "month", "name_2", "outbreak_2", "district_country", "country", "code",
              "rainfall", "eastward_wind", "north_wind", "aod", "humidity", "windspeed", "temp")
  
  # Keep only columns that exist (protects against missing columns)
  existing_vars <- intersect(myvars, names(merged_data))
  merged_data <- merged_data[, existing_vars]
  
  # Data structure and summary
  print("---- Structure of merged data ----")
  str(merged_data)
  
  print("---- Summary statistics ----")
  print(summary(merged_data))
  
  print("---- Missing values per column ----")
  na_counts <- sapply(merged_data, function(x) sum(is.na(x)))
  print(na_counts)
  
  print("---- Proportion of missing values per column ----")
  na_proportion <- round(na_counts / nrow(merged_data), 3)
  print(na_proportion)
  
  # Summary for numeric columns
  numeric_cols <- sapply(merged_data, is.numeric)
  print("---- Means and ranges for numeric columns ----")
  print(summary(merged_data[, numeric_cols]))
  
  # Unique counts for character columns
  categorical_cols <- sapply(merged_data, is.character)
  print("---- Unique values for character columns ----")
  unique_counts <- sapply(merged_data[, categorical_cols, drop = FALSE], function(x) length(unique(x)))
  print(unique_counts)
  
  # Example rows
  print("---- Example rows ----")
  print(head(merged_data, 10))
  
  
  landcover_data <- read.csv("final_landcover_data.csv", stringsAsFactors = FALSE)
  
  # Create merge key in landcover data
  landcover_data$year_district <- paste(landcover_data$year, landcover_data$district_country, sep = " ")
  
  # Create same merge key in merged_data
  merged_data$year_district <- paste(merged_data$year, merged_data$district_country, sep = " ")
  
  # Deduplicate landcover data (optional, in case there are duplicates)
  landcover_data <- landcover_data %>% distinct(year_district, .keep_all = TRUE)
  
  # Merge landcover onto merged_data
  merged_data <- merged_data %>%
    left_join(landcover_data, by = "year_district")
  
  # Optional: Remove the temporary merge key
  merged_data$year_district <- NULL
  merged_data<- na.omit(merged_data)
  # Check final structure
  print("---- Structure after merging landcover ----")
  str(merged_data)
  
  # Data structure and summary
  print("---- Structure of merged data ----")
  str(merged_data)
  
  print("---- Summary statistics ----")
  print(summary(merged_data))
  
  print("---- Missing values per column ----")
  na_counts <- sapply(merged_data, function(x) sum(is.na(x)))
  print(na_counts)
  
  print("---- Proportion of missing values per column ----")
  na_proportion <- round(na_counts / nrow(merged_data), 3)
  print(na_proportion)
  
  # Summary for numeric columns
  numeric_cols <- sapply(merged_data, is.numeric)
  print("---- Means and ranges for numeric columns ----")
  print(summary(merged_data[, numeric_cols]))
  
  # Unique counts for character columns
  categorical_cols <- sapply(merged_data, is.character)
  print("---- Unique values for character columns ----")
  unique_counts <- sapply(merged_data[, categorical_cols, drop = FALSE], function(x) length(unique(x)))
  print(unique_counts)
  
  # Example rows
  print("---- Example rows ----")
  print(head(merged_data, 10))
  
  write.csv(merged_data, "joined_outbreak_data.csv", row.names = FALSE)
}
