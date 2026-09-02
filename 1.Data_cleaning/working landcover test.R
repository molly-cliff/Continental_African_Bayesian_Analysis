setwd("C:/Users/mvc32/OneDrive - University of Cambridge/Documents/Climate_meningitis_belt")
# Load shapefile
shape2 <- st_read("Shapefile_improved.shp")


country_list <- c(
  "Egypt", "Libya", "Algeria",
  "Tunisia", "Mauritania",  "Western Sahara", "Niger",
  "Chad", "Mali", "Sudan"
)


shape2 <- shape2  %>% 
  filter(COUNTRY %in% country_list)
plot(st_geometry(shape2))

# Plot geometry
windows(record = TRUE)

landcover_dir=("landcover")

setwd("C:/Users/mvc32/OneDrive - University of Cambridge/Documents/INLA project/landcovertest")

# ---- Read landcover file ----

nc_files <- list.files(pattern = "\\.nc$", full.names = TRUE)
if (length(nc_files) == 0) stop("No NetCDF (.nc) files found in the directory!")

message("Reading: ", basename(nc_files[1]))
r <- rast(nc_files[1])

# ---- Helper: Extract last 12 layers for a variable ----
extract_last_layers <- function(pattern) {
  x <- subset(r, grepl(pattern, names(r)))
  x[[ (nlyr(x) - 12 + 1):nlyr(x) ]]
}


# ---- Extract required variables ----
primf <- extract_last_layers("^primf")
primn <- extract_last_layers("^primn")
secdf <- extract_last_layers("^secdf")
secdn <- extract_last_layers("^secdn")
urban <- extract_last_layers("^urban")
c3ann <- extract_last_layers("^c3ann")
c4ann <- extract_last_layers("^c4ann")
c3per <- extract_last_layers("^c3per")
c4per <- extract_last_layers("^c4per")
c3nfx <- extract_last_layers("^c3nfx")
pastr <- extract_last_layers("^pastr")
range <- extract_last_layers("^range")

setwd("C:/Users/mvc32/OneDrive - University of Cambridge/Documents/INLA project/landcovertest3")


nc_files <- list.files(pattern = "\\.nc$", full.names = TRUE)
if (length(nc_files) == 0) stop("No NetCDF (.nc) files found in the directory!")

message("Reading: ", basename(nc_files[1]))
r <- rast(nc_files[1])

extract_first_layers <- function(pattern) {
  x <- subset(r, grepl(pattern, names(r)))
  x[[ 1:8 ]]
}
# ---- Extract required variables ----
primf2 <- extract_first_layers("^primf")
primn2 <- extract_first_layers("^primn")
secdf2 <- extract_first_layers("^secdf")
secdn2 <- extract_first_layers("^secdn")
urban2 <- extract_first_layers("^urban")
c3ann2 <- extract_first_layers("^c3ann")
c4ann2 <- extract_first_layers("^c4ann")
c3per2 <- extract_first_layers("^c3per")
c4per2 <- extract_first_layers("^c4per")
c3nfx2 <- extract_first_layers("^c3nfx")
pastr2 <- extract_first_layers("^pastr")
range2 <- extract_first_layers("^range")


primf<- c(primf,primf2)
primn<- c(primn,primn2)
secdf<- c(secdf,secdf2)
secdn<- c(secdn,secdn2)
urban<- c(urban,urban2)
c3ann<- c(c3ann,c3ann2)
c4ann<- c(c4ann,c4ann2)
c3per<- c(c3per,c3per2)
c4per<- c(c4per,c4per2)
c3nfx<- c(c3nfx,c3nfx2)
pastr<- c(pastr,pastr2)
range<- c(range,range2)

# ---- Align CRS for all rasters at once ----
target_crs <- crs(shape2)
all_rasters <- list(primf, primn, secdf, secdn, urban, c3ann, c4ann, c3per, c4per, c3nfx, pastr, range)
all_rasters <- lapply(all_rasters, function(ras) {
  if (!terra::crs(ras) == st_crs(shape2)$wkt) {
    message("CRS mismatch - Reprojecting raster")
    project(ras, st_crs(shape2)$wkt)
  } else {
    ras
  }
  
})
list2env(setNames(all_rasters, c("primf", "primn", "secdf", "secdn", "urban", 
                                 "c3ann", "c4ann", "c3per", "c4per", "c3nfx", 
                                 "pastr", "range")), envir = environment())

# ---- Compute cropland, forest ----
cropland_stack <- c3ann + c4ann + c3per + c4per + c3nfx
forest_stack   <- primf + secdf

# ---- Compute barren ----
barren_stack <- ((urban < 0.01) &
                   (pastr < 0.01) &
                   (range < 0.01) &
                   (cropland_stack < 0.01) &
                   (forest_stack < 0.01)) * 1
# ---- Crop & mask once for each final raster ----
mask_crop <- function(x) mask(crop(x, shape2), shape2)
cropland_stack <- mask_crop(cropland_stack)
forest_stack   <- mask_crop(forest_stack)
barren_stack   <- mask_crop(barren_stack)

# ---- Extract mean values ----
years <- 2003:2022

extract_to_df <- function(stack, varname) {
  vals <- terra::extract(stack, shape2, fun = mean, na.rm = TRUE)
  colnames(vals)[-1] <- years
  
  vals %>%
    mutate(district_country = paste(shape2$NAME_2, shape2$COUNTRY)) %>%
    pivot_longer(-district_country, names_to = "year", values_to = varname) %>%
    mutate(year = as.integer(year),
           code = paste(year, district_country)) %>%
    select(code, year, district_country, !!varname)
}

cropland_df <- extract_to_df(cropland_stack, "cropland") %>% 
  filter(!is.na(year))

forest_df <- extract_to_df(forest_stack, "forest") %>% 
  filter(!is.na(year))

barren_df <- extract_to_df(barren_stack, "barren") %>% 
  filter(!is.na(year))
# ---- Merge into one data.frame ----
final_df <- cropland_df %>%
  left_join(forest_df, by = c("code", "year", "district_country")) %>%
  left_join(barren_df, by = c("code", "year", "district_country")) %>%
  distinct(code, .keep_all = TRUE)

return(final_df)
