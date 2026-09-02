library(raster)
library(terra)
library(dplyr)
library(stringr)
library(sf)
library(lubridate)

testing_extract_landcover <- function() {
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
      code = paste(year_month, district_country, sep = " "))
  
  # Set file path and read in NetCDF files
  
  
  path <- "landcover"
  nc_files <- list.files(path,pattern = "\\.nc$", full.names = TRUE)
  
  # Initialize an empty list to hold the monthly rasters
  monthly_rasters <- list()
  
  # Loop through each file
  
  for (file in nc_files) {
    message("Processing: ", basename(file))
    
    r <- rast(file)
    var_names <- names(r)
    print(var_names)
  }
  
  plot(r[[1]])
  varnames(r)
  
  
  r2<-subset(r,1:1166)
  r2 <- r2[[ (nlyr(r2) - 12) : nlyr(r2) ]]
  r3<-subset(r,1167:2332)
  r3 <- r3[[ (nlyr(r3) - 12) : nlyr(r3) ]]
  
  r4<-subset(r,2333:3498)
  r4 <- r4[[ (nlyr(r4) - 12) : nlyr(r4) ]]
  r5<-subset(r,3499:4664)
  r5 <- r5[[ (nlyr(r5) - 12) : nlyr(r5) ]]
  r6<-subset(r,4665:5830)
  r6 <- r6[[ (nlyr(r6) - 12) : nlyr(r6) ]]
  r7<-subset(r,5831:6996)
  r7 <- r7[[ (nlyr(r7) - 12) : nlyr(r7) ]]
  r8<-subset(r,6997:8162)
  r8 <- r8[[ (nlyr(r8) - 12) : nlyr(r8) ]]
  r9<-subset(r,8163:9328)
  r9 <- r9[[ (nlyr(r9) - 12) : nlyr(r9) ]]
  r10<-subset(r,9329:10494)
  r10 <- r10[[ (nlyr(r10) - 12) : nlyr(r10) ]]
  r11<-subset(r,10495:11660)
  r11 <- r11[[ (nlyr(r11) - 12) : nlyr(r11) ]]
  r12<-subset(r,11661:13992)
  r12 <- r12[[ (nlyr(r12) - 12) : nlyr(r12) ]]
  r13<-subset(r,11662:15158)
  r13 <- r13[[ (nlyr(r13) - 12) : nlyr(r13) ]]
  r14<-subset(r,15158:16324)
  r14 <- r14[[ (nlyr(r14) - 12) : nlyr(r14) ]]
  
  
  
  align_crs <- function(raster_obj, vector_obj) {
    if (!crs(raster_obj) == crs(vector_obj)) {
      message("CRS mismatch - Reprojecting")
      raster_obj <- project(raster_obj, crs(vector_obj))
    } else {
      message("CRS aligned")
    }
    return(raster_obj)
  }
  
  align_crs(r2,shape2)
  align_crs(r3,shape2)
  align_crs(r4,shape2)
  align_crs(r5,shape2)
  align_crs(r6,shape2)
  align_crs(r7,shape2)
  align_crs(r8,shape2)
  align_crs(r9,shape2)
  align_crs(r10,shape2)
  align_crs(r11,shape2)
  align_crs(r12,shape2)
  align_crs(r13,shape2)
  align_crs(r14,shape2)
  
  
  
  r2 <- mask(crop(r2, shape2), shape2)
  r3 <- mask(crop(r3, shape2), shape2)
  r4 <- mask(crop(r4, shape2), shape2)
  r5 <- mask(crop(r5, shape2), shape2)
  r6 <- mask(crop(r6, shape2), shape2)
  r7 <- mask(crop(r7, shape2), shape2)
  r8 <- mask(crop(r8, shape2), shape2)
  r9 <- mask(crop(r9, shape2), shape2)
  r10 <- mask(crop(r10, shape2), shape2)
  r11 <- mask(crop(r11, shape2), shape2)
  r12 <- mask(crop(r12, shape2), shape2)
  r13 <- mask(crop(r13, shape2), shape2)
  r14 <- mask(crop(r14, shape2), shape2)
  
  
  final_data$code <- paste(final_data$year, final_data$district_country, sep = " ")
  
  
  # Prepare final_data for merging
  final_data <- final_data %>%
    mutate(
      code = paste(year, district_country, sep = " "),
      landcover = NA  # initialize
    )
  
  
  # Assign the correct names (2003 to 2015) to the layers
  years <- 2003:2015
  
  # Function to rename layers for each raster stack
  rename_raster_layers <- function(raster_stack, years) {
    # Make sure the number of layers matches the number of years
    if (nlyr(raster_stack) == length(years)) {
      names(raster_stack) <- as.character(years)
    } else {
      stop("The number of layers does not match the number of years.")
    }
    return(raster_stack)
  }
  
  # Now, apply this to all your raster stacks
  r2 <- rename_raster_layers(r2, years)
  r3 <- rename_raster_layers(r3, years)
  r4 <- rename_raster_layers(r4, years)
  r5 <- rename_raster_layers(r5, years)
  r6 <- rename_raster_layers(r6, years)
  r7 <- rename_raster_layers(r7, years)
  r8 <- rename_raster_layers(r8, years)
  r9 <- rename_raster_layers(r9, years)
  r10 <- rename_raster_layers(r10, years)
  r11 <- rename_raster_layers(r11, years)
  r12 <- rename_raster_layers(r12, years)
  r13 <- rename_raster_layers(r13, years)
  r14 <- rename_raster_layers(r14, years)
  
  
  
  
  #### LANDCOVER DATA####  
  
  path <- "landcover2"
  
  
  nc_files <- list.files(path,pattern = "\\.nc$", full.names = TRUE)
  
  for (file in nc_files) {
    message("Processing: ", basename(file))
    
    r <- rast(file)
    var_names <- names(r)
    print(var_names)
  }
  
  plot(r[[1]])
  varnames(r)
  
  
  
  
  r2a<-subset(r,1:86)
  r2a <- r2a[[2:8]]
  r3a<-subset(r,87:172)
  r3a <- r3a[[2:8]]
  
  r4a<-subset(r,173:258)
  r4a <- r4a[[2:8]]
  r5a<-subset(r,259:344)
  r5a <- r5a[[2:8]]
  r6a<-subset(r,345:430)
  r6a <- r6a[[2:8]]
  r7a<-subset(r,431:516)
  r7a <- r7a[[2:8]]
  r8a<-subset(r,517:602)
  r8a <- r8a[[2:8]]
  r9a<-subset(r,603:774)
  r9a <- r9a[[2:8]]
  r10a<-subset(r,775:860)
  r10a <- r10a[[2:8]]
  r11a<-subset(r,861:946)
  r11a <- r11a[[2:8]]
  r12a<-subset(r,947:1033)
  r12a <- r12a[[2:8]]
  r13a<-subset(r,1033:1120)
  r13a <- r13a[[2:8]]
  r14a<-subset(r,1119:1205)
  r14a <- r14a[[2:8]]
  
  
  
  align_crs <- function(raster_obj, vector_obj) {
    if (!crs(raster_obj) == crs(vector_obj)) {
      message("CRS mismatch - Reprojecting")
      raster_obj <- project(raster_obj, crs(vector_obj))
    } else {
      message("CRS aligned")
    }
    return(raster_obj)
  }
  
  align_crs(r2a,shape2)
  align_crs(r3a,shape2)
  align_crs(r4a,shape2)
  align_crs(r5a,shape2)
  align_crs(r6a,shape2)
  align_crs(r7a,shape2)
  align_crs(r8a,shape2)
  align_crs(r9a,shape2)
  align_crs(r10a,shape2)
  align_crs(r11a,shape2)
  align_crs(r12a,shape2)
  align_crs(r13a,shape2)
  align_crs(r14a,shape2)
  
  
  
  r2a <- mask(crop(r2a, shape2), shape2)
  r3a <- mask(crop(r3a, shape2), shape2)
  r4a <- mask(crop(r4a, shape2), shape2)
  r5a <- mask(crop(r5a, shape2), shape2)
  r6a <- mask(crop(r6a, shape2), shape2)
  r7a <- mask(crop(r7a, shape2), shape2)
  r8a <- mask(crop(r8a, shape2), shape2)
  r9a <- mask(crop(r9a, shape2), shape2)
  r10a <- mask(crop(r10a, shape2), shape2)
  r11a <- mask(crop(r11a, shape2), shape2)
  r12a <- mask(crop(r12a, shape2), shape2)
  r13a <- mask(crop(r13a, shape2), shape2)
  r14a <- mask(crop(r14a, shape2), shape2)
  
  
  years <- 2016:2022
  
  r2a <- rename_raster_layers(r2a, years)
  r3a <- rename_raster_layers(r3a, years)
  r4a <- rename_raster_layers(r4a, years)
  r5a <- rename_raster_layers(r5a, years)
  r6a <- rename_raster_layers(r6a, years)
  r7a <- rename_raster_layers(r7a, years)
  r8a <- rename_raster_layers(r8a, years)
  r9a <- rename_raster_layers(r9a, years)
  r10a <- rename_raster_layers(r10a, years)
  r11a <- rename_raster_layers(r11a, years)
  r12a <- rename_raster_layers(r12a, years)
  r13a <- rename_raster_layers(r13a, years)
  r14a <- rename_raster_layers(r14a, years)
  
  
  
  r2 <- c(r2, r2a)
  r3 <- c(r3, r3a)
  r4 <- c(r4, r4a)
  r5 <- c(r5, r5a)
  r6 <- c(r6, r6a)
  r7 <- c(r7, r7a)
  r8 <- c(r8, r8a)
  r9 <- c(r9, r9a)
  r10 <- c(r10, r10a)
  r11 <- c(r11, r11a)
  r12 <- c(r12, r12a)
  r13 <- c(r13, r13a)
  r14 <- c(r14, r14a)
  
  
  # Year-mapping aware processing function
  process_raster_stack <- function(raster_stack, shape2, final_data, value_column = "landcover") {
    
    for (i in 1:nlyr(raster_stack)) {
      lyr <- raster_stack[[i]]
      lyr_name <- names(raster_stack)[i]
      cat("Processing:", lyr_name, "\n")
      
      # Directly use the layer name as the year
      year <- as.integer(lyr_name)
      
      # Extract mean values per district
      extracted <- data.frame(
        shape2,
        value_temp = terra::extract(lyr, shape2, fun = mean, na.rm = TRUE)[, 2]
      ) %>%
        mutate(
          year = year,
          district_country = paste(NAME_2, COUNTRY, sep = " "),
          code = paste(year, district_country, sep = " ")
        ) %>%
        dplyr::select(code, value_temp)
      
      # Merge and fill NA values only
      final_data <- merge(final_data, extracted, by = "code", all.x = TRUE)
      final_data[[value_column]][is.na(final_data[[value_column]])] <- final_data$value_temp[is.na(final_data[[value_column]])]
      
      # Remove the temporary column
      final_data$value_temp <- NULL
    }
    
    return(final_data)
  }
  
  
  final_data$landcover <- NULL
  final_data$landcover <- NA
  final_data <- process_raster_stack(r2, shape2, final_data, value_column = "landcover")
  final_data$"primf"<-final_data$landcover
  final_data$landcover <- NULL
  final_data$landcover <- NA
  final_data <- process_raster_stack(r3, shape2, final_data, value_column = "landcover")
  final_data$"primn"<-final_data$"landcover"
  final_data$landcover <- NULL
  final_data$landcover <- NA
  final_data <- process_raster_stack(r4, shape2, final_data, value_column = "landcover")
  final_data$"secdf"<-final_data$"landcover"
  final_data$landcover <- NULL
  final_data$landcover <- NA
  final_data <- process_raster_stack(r5, shape2, final_data, value_column = "landcover")
  final_data$"secdn"<-final_data$"landcover"
  final_data$landcover <- NULL
  final_data$landcover <- NA
  final_data <- process_raster_stack(r6, shape2, final_data, value_column = "landcover")
  final_data$"urban"<-final_data$"landcover"
  final_data$landcover <- NULL
  final_data$landcover <- NA
  final_data <- process_raster_stack(r7, shape2, final_data, value_column = "landcover")
  final_data$"c3ann"<-final_data$"landcover"
  final_data$landcover <- NULL
  final_data$landcover <- NA
  final_data <- process_raster_stack(r8, shape2, final_data, value_column = "landcover")
  final_data$"c4ann"<-final_data$"landcover"
  final_data$landcover <- NULL
  final_data$landcover <- NA
  final_data <- process_raster_stack(r9, shape2, final_data, value_column = "landcover")
  final_data$"c3per"<-final_data$"landcover"
  final_data$landcover <- NULL
  final_data$landcover <- NA
  final_data <- process_raster_stack(r10, shape2, final_data, value_column = "landcover")
  final_data$"c4per"<-final_data$"landcover"
  final_data$landcover <- NULL
  final_data$landcover <- NA
  final_data <- process_raster_stack(r11, shape2, final_data, value_column = "landcover")
  final_data$"c3nfx"<-final_data$"landcover"
  final_data$landcover <- NULL
  final_data$landcover <- NA
  final_data <- process_raster_stack(r12, shape2, final_data, value_column = "landcover")
  final_data$"range"<-final_data$"landcover"
  final_data$landcover <- NULL
  final_data$landcover <- NA
  final_data <- process_raster_stack(r13, shape2, final_data, value_column = "landcover")
  final_data$"secmb"<-final_data$"landcover"
  final_data$landcover <- NULL
  final_data$landcover <- NA
  final_data <- process_raster_stack(r14, shape2, final_data, value_column = "landcover")
  final_data$"secma"<-final_data$"landcover"
  final_data$landcover <- NULL
  
  
  write_xlsx(final_data, "landcovertesting.xlsx")
  
}