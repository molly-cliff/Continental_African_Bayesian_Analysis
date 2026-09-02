# Create output folder (if it doesn't exist yet)
library(dplyr)
library(INLA)
library(sf)
library(spdep)
library(sp)
library(caret)
library(pROC)
library(PRROC)
library(tidyr)
library(broom)
library(Hmisc)
library(magick)
library(corrplot)
library(car)
library(ggplot2)
library(reshape2)
library(cowplot)
testing_map<- function() {
spatiotemporaloutbreaks <- read.csv("final_popdensity_vaccine.csv")
print(names(spatiotemporaloutbreaks))


spatiotemporaloutbreaks$year_month <- (spatiotemporaloutbreaks$year.x - min(spatiotemporaloutbreaks$year.x)) * 12 + spatiotemporaloutbreaks$month

myvars <- c("year.x", "month", "outbreak_2", "name_2", "district_country.x", "country",
            "code.x","rainfall","eastward_wind","north_wind", "aod","humidity","vaccine"          
            ,"windspeed", "temp","cropland",          
            "forest", "barren","year_month", "pop_density" )
merged_data <- spatiotemporaloutbreaks[, myvars]
merged_data <- merged_data[merged_data$country %in% c("Benin", "Burkina Faso", "Burundi", "Cameroon", "Central African Republic",
                                                      "Chad", "Côte d'Ivoire", "Eritrea", "Ethiopia", 
                                                      "Gambia", "Ghana", "Guinea", "Guinea-Bissau", "Kenya", "Mali", "Mauritania",
                                                      "Niger", "Nigeria", "Rwanda", "Senegal", "South Sudan", "Sudan", "Tanzania",
                                                      "Togo", "Uganda"), ]



shape2 <- st_read("Shapefile_improved.shp")

shape2<- shape2[shape2$COUNTRY %in% c("Benin", "Burkina Faso", "Burundi", "Cameroon", "Central African Republic",
                                      "Chad", "Côte d'Ivoire", "Eritrea", "Ethiopia",
                                      "Gambia", "Ghana", "Guinea", "Guinea-Bissau", "Kenya", "Mali", "Mauritania",
                                      "Niger", "Nigeria", "Rwanda", "Senegal", "South Sudan", "Sudan", "Tanzania",
                                      "Togo", "Uganda"), ]

# 2. Check and fix invalid geometries
shape2 <- st_make_valid(shape2)  # Optionally fix geometries

# Optional: Remove invalid ones instead
# shape2 <- shape2[st_is_valid(shape2), ]

# 3. Convert to Spatial* object
shapefile_spatial <- as(shape2, "Spatial")

# 4. Create neighborhood structure
nb <- poly2nb(shapefile_spatial, queen = TRUE)

nb2INLA("full_africa_map.adj", nb)
ken.adj <- paste0(getwd(), "/full_africa_map.adj")

shape2$district_country.x <- paste(shape2$NAME_2, shape2$COUNTRY, sep = " ")

merged_data$area<-as.numeric(as.factor(merged_data$district_country.x))

missing_in_shape <- setdiff(merged_data$district_country.x, shape2$district_country.x)
if (length(missing_in_shape) > 0) {
  stop("❌ These districts are in merged_data but not in shapefile: ",
       paste(missing_in_shape, collapse = ", "))
}

# --- 2. In shape but not in data ---
missing_in_data <- setdiff(shape2$district_country.x, merged_data$district_country.x)

if (length(missing_in_data) > 0) {
  message("⚠ Dropping ", length(missing_in_data), 
          " districts from shapefile that are not in merged_data:\n",
          paste(missing_in_data, collapse = ", "))
  
  # Keep only those in both
  shape2 <- shape2[shape2$district_country.x %in% merged_data$district_country.x, ]
}


# --- 3. Check sorted lists match ---
if (!identical(sort(unique(merged_data$district_country.x)),
               sort(unique(shape2$district_country.x)))) {
  stop("❌ District-country lists differ even after sorting!")
} else {
  message("✅ All district-country matches verified.")
}

# Ensure consistent ordering
shape2 <- shape2[order(shape2$district_country.x), ]
merged_data <- merged_data[order(merged_data$district_country.x), ]

# Assign IDs
merged_data$area <- match(merged_data$district_country.x, shape2$district_country.x)


shape2 <- shape2[order(shape2$district_country.x), ]
merged_data <- merged_data[order(merged_data$district_country.x), ]
merged_data$area <- match(merged_data$district_country.x, shape2$district_country.x)



test<-table(merged_data$outbreak_2)
test<-as.data.frame(test)
print("OUTBREAK OCCURANCE:")
print(test)
print("NUMBER OF OUTBREAKS V NO OUTBREAKS:")
outbreak<-test[2,2]
no_outbreak<-test[1,2]

print(outbreak)
print(no_outbreak)

outbreak<-as.numeric(outbreak)
no_outbreak<-as.numeric(no_outbreak)
total<-no_outbreak+outbreak
non_outbreak_weight <- outbreak/total
outbreak_weight<-1-non_outbreak_weight

print("OUTBREAK WEIGHTS:")
print(outbreak_weight)

weights <- ifelse(merged_data$outbreak_2 == 1, outbreak_weight, non_outbreak_weight) 
# ---- Covariates / scaling ----
  merged_data$eastward_wind <- as.numeric(merged_data$eastward_wind)
  
  merged_data$eastward_wind_scale <- as.numeric(scale(merged_data$eastward_wind))
  merged_data$aod_scale <- scale(merged_data$aod)
  merged_data$humidity_scale <- scale(merged_data$humidity)
  merged_data$temp_scale <- scale(merged_data$temp)
  merged_data$windspeed_scale <- scale(merged_data$windspeed)
  merged_data$pop_density_scale <- scale(merged_data$pop_density)
  
  merged_data$temp_index <- as.numeric(cut(merged_data$temp_scale, breaks = 50))
  merged_data$pop_index <- as.numeric(cut(merged_data$pop_density_scale, breaks = 50))
  merged_data$aod_index <- as.numeric(cut(merged_data$aod_scale, breaks = 50))
  
  merged_data$wind_strength <- abs(merged_data$eastward_wind)
  merged_data$wind_strength_scale <- as.numeric(scale(merged_data$wind_strength))
  
  merged_data$wind_dir <- ifelse(merged_data$eastward_wind < 0, 1, 0)
  
  merged_data$transport_index <- merged_data$aod_scale * pmax(0, -merged_data$eastward_wind_scale)
  
  # ---- Model ----
  formula1 <- outbreak_2 ~ aod_scale +
    wind_strength_scale +
    humidity_scale +
    wind_dir +
    aod_scale:wind_dir +
    f(year_month, model = "rw2",
      hyper = list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))) +
    f(temp_index, model = "rw2",
      hyper = list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))) +
    f(area, model = "bym2", graph = ken.adj, scale.model = TRUE, constr = TRUE,
      hyper = list(
        prec = list(prior = "pc.prec", param = c(1, 0.01)),
        phi  = list(prior = "pc", param = c(0.5, 2/3))
      ))
  
  model <- inla(
    formula1,
    data = merged_data,
    family = "binomial",
    weights = weights,
    control.fixed = list(mean = 0, prec = 1),
    control.predictor = list(compute = TRUE, link = 1),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)
  )
  
  print(summary(model))
# ---- 2. Model Summary ----
cat("\n==============================\n")
cat("Full Model Summary\n")
cat("==============================\n")
print(summary(model))
cat("DIC:", model$dic$dic, "\n")
cat("WAIC:", model$waic$waic, "\n")
log_cpo <- -mean(log(model$cpo$cpo), na.rm = TRUE)
cat("Mean Log CPO:", log_cpo, "\n")

merged_data$M5 <- model$summary.fitted.values[, 1]

# Create year_month identifier
merged_data$year_month <- paste(
  merged_data$year.x,
  sprintf("%02d", as.integer(merged_data$month)),
  sep = "_"
)

# Merge with shapefile
merged_sf <- shape2 %>%
  left_join(merged_data, by = "district_country.x")

# Split by year_month
split_sf <- split(merged_sf, merged_sf$year_month)

# Create list of sf objects (optionally assign globally)
shapefiles_list <- lapply(names(split_sf), function(ym) {
  df <- split_sf[[ym]]
  df_sf <- st_as_sf(df)  # Ensure sf class
  assign(paste0("shapefile_", ym), df_sf, envir = .GlobalEnv)  # Optional
  return(df_sf)
})
names(shapefiles_list) <- names(split_sf)

# Output directory
outdir <- "monthly_predictions_plots"
if (!dir.exists(outdir)) dir.create(outdir)
# Loop over year_month shapefiles
for (year_month in names(shapefiles_list)) {
  current_shapefile <- shapefiles_list[[year_month]]
  current_shapefile$M5 <- as.numeric(current_shapefile$M5)
  
  # Convert year_month index to actual year and month
  # Assuming year_month starts from 1 and increments by 1 for each month
  year_month_index <- as.numeric(year_month)
  
  # Manually set starting point
  start_year <- 2003
  start_month <- 1
  
  # Calculate offset from start
  total_months_offset <- year_month_index - 1
  actual_year <- start_year + floor(total_months_offset / 12)
  actual_month <- ((total_months_offset) %% 12) + start_month
  
  # Handle month overflow
  if (actual_month > 12) {
    actual_year <- actual_year + 1
    actual_month <- actual_month - 12
  }
  
  # Format as "YYYY MM"
  date_label <- sprintf("%04d %02d", actual_year, actual_month)
  
  plot <- ggplot(data = current_shapefile) +
    geom_sf(aes(fill = M5), color = "black", size = 0.1) +
    scale_fill_gradient(
      low = "white", high = "red", na.value = "gray"
    ) +
    labs(
      title = paste("Outbreak Probability in", date_label),
      fill = "Probability"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5)
    )
  
  # Save PNG with date label
  ggsave(
    filename = file.path(outdir, paste0("outbreak_prob_", date_label, ".png")),
    plot = plot,
    width = 8, height = 6, dpi = 300
  )
}

# Create GIF from PNGs
png_files <- list.files(outdir, pattern = "outbreak_prob_.*\\.png$", full.names = TRUE)
# Sort files by date
png_files <- sort(png_files)
img_list <- image_read(png_files)
animation <- image_animate(img_list, fps = 5)

# Preview and save
print(animation)
image_write(animation, "animation_nodrc.gif")

}