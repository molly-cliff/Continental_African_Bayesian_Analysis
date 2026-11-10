
library(lubridate)
library(terra)
library(dplyr)
library(stringr)
library(sf)

data_dir <- getwd()
landcover_dir <- file.path(data_dir, "landcover")
outbreak_file <- file.path(data_dir, "full-outbreak-matched.csv")
shapefile <- file.path(data_dir, "Shapefile_improved.shp")
output_file <- file.path(data_dir, "final_landcover_data.csv")
years <- 2003:2022

weekly_data <- read.csv(outbreak_file, stringsAsFactors = FALSE)
shape2 <- st_read(shapefile, quiet = TRUE)

required_outbreak_columns <- c("year", "month", "district_country", "outbreak")
missing_outbreak_columns <- setdiff(required_outbreak_columns, names(weekly_data))
if (length(missing_outbreak_columns) > 0) {
  stop("Missing outbreak columns: ", paste(missing_outbreak_columns, collapse = ", "))
}
if (!all(c("NAME_2", "COUNTRY") %in% names(shape2))) {
  stop("The shapefile must contain NAME_2 and COUNTRY fields.")
}

shape2 <- shape2[!shape2$COUNTRY %in%
                   c("Cabo Verde", "Mauritius", "Seychelles"), ]
shape2$district_country <- paste(shape2$NAME_2, shape2$COUNTRY)

outbreak_data <- weekly_data %>%
  group_by(year, month, district_country) %>%
  summarise(
    outbreak_occur = as.integer(any(outbreak > 0, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(code = paste(year, month, district_country, sep = " "))

nc_files <- list.files(landcover_dir, pattern = "\\.nc$", full.names = TRUE)
if (length(nc_files) == 0) {
  stop("No NetCDF (.nc) files found in: ", landcover_dir)
}

landcover <- rast(nc_files)
message("Reading ", nlyr(landcover), " landcover layers from ", length(nc_files), " file(s)")

select_layers <- function(pattern) {
  layer_numbers <- which(grepl(pattern, names(landcover)))
  if (length(layer_numbers) != length(years)) {
    stop(
      "Expected ", length(years), " layers matching /", pattern,
      "/ but found ", length(layer_numbers), "."
    )
  }
  landcover[[layer_numbers]]
}

variables <- list(
  primf = select_layers("^primf"),
  secdf = select_layers("^secdf"),
  urban = select_layers("^urban"),
  c3ann = select_layers("^c3ann"),
  c4ann = select_layers("^c4ann"),
  c3per = select_layers("^c3per"),
  c4per = select_layers("^c4per"),
  c3nfx = select_layers("^c3nfx"),
  pastr = select_layers("^pastr"),
  range = select_layers("^range")
)

target_crs <- st_crs(shape2)$wkt
variables <- lapply(variables, function(raster_data) {
  if (!same.crs(raster_data, vect(shape2))) {
    project(raster_data, target_crs)
  } else {
    raster_data
  }
})

cropland <- variables$c3ann + variables$c4ann + variables$c3per +
  variables$c4per + variables$c3nfx
forest <- variables$primf + variables$secdf
barren <- ((variables$urban < 0.01) &
             (variables$pastr < 0.01) &
             (variables$range < 0.01) &
             (cropland < 0.01) &
             (forest < 0.01)) * 1

extract_annual_means <- function(raster_data, variable_name) {
  raster_data <- mask(crop(raster_data, vect(shape2)), vect(shape2))
  extracted <- terra::extract(raster_data, vect(shape2), fun = mean, na.rm = TRUE)
  extracted <- extracted[, -1, drop = FALSE]
  names(extracted) <- as.character(years)

  data.frame(
    district_country = shape2$district_country,
    extracted,
    check.names = FALSE
  ) %>%
    tidyr::pivot_longer(
      cols = all_of(as.character(years)),
      names_to = "year",
      values_to = variable_name
    ) %>%
    mutate(year = as.integer(year))
}

landcover_data <- list(
  extract_annual_means(cropland, "cropland"),
  extract_annual_means(forest, "forest"),
  extract_annual_means(barren, "barren")
) %>%
  Reduce(function(left, right) {
    left %>% left_join(right, by = c("district_country", "year"))
  }, x = .)

final_data <- outbreak_data %>%
  left_join(landcover_data, by = c("year", "district_country"))

write.csv(final_data, output_file, row.names = FALSE)

