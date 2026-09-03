library(dplyr)
library(INLA)
library(sf)
library(spdep)
library(sp)
library(caret)
library(pROC)
library(tidyr)
library(PRROC)
library(broom)
library(corrplot)
library(car)
library(ggplot2)
library(reshape2)
library(Hmisc)
library(writexl)

univariate_analysis <-function() {
  # Left join all environmental datasets onto spatiotemporaloutbreaks
  
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
  
  nb2INLA("map.adj", nb)
  ken.adj <- paste0(getwd(), "/map.adj")
  
  shape2$district_country.x <- paste(shape2$NAME_2, shape2$COUNTRY, sep = " ")
  
  merged_data$area<-as.numeric(as.factor(merged_data$district_country.x))
  
  missing_in_shape <- setdiff(merged_data$district_country.x, shape2$district_country.x)
  if (length(missing_in_shape) > 0) {
    stop("districts are in merged_data but not in shapefile: ",
         paste(missing_in_shape, collapse = ", "))
  }
  
  # --- 2. In shape but not in data ---
  missing_in_data <- setdiff(shape2$district_country.x, merged_data$district_country.x)
  
  if (length(missing_in_data) > 0) {
    message("dropping ", length(missing_in_data), 
            " districts from shapefile that are not in merged_data:\n",
            paste(missing_in_data, collapse = ", "))
    
    # Keep only those in both
    shape2 <- shape2[shape2$district_country.x %in% merged_data$district_country.x, ]
  }
  
  
  # --- 3. Check sorted lists match ---
  if (!identical(sort(unique(merged_data$district_country.x)),
                 sort(unique(shape2$district_country.x)))) {
    stop("dstrict-country lists differ even after sorting!")
  } else {
    message("matches verified.")
  }
  
  # Ensure consistent ordering
  shape2 <- shape2[order(shape2$district_country.x), ]
  merged_data <- merged_data[order(merged_data$district_country.x), ]
  
  # Assign IDs
  merged_data$area <- match(merged_data$district_country.x, shape2$district_country.x)
  
  
  shape2 <- shape2[order(shape2$district_country.x), ]
  merged_data <- merged_data[order(merged_data$district_country.x), ]
  merged_data$area <- match(merged_data$district_country.x, shape2$district_country.x)
  
  

  
  # Step 2: Calculate monthly climatology for multiple variables
  climatology <- merged_data %>%
    group_by(month) %>%
    summarise(across(c(temp, rainfall,windspeed,eastward_wind, north_wind, humidity, pop_density, aod, cropland, forest, barren), ~ mean(.x, na.rm = TRUE), .names = "clim_{.col}"))
  
  # Step 3: Join climatology to full dataset and compute anomalies
  merged_data <- merged_data %>%
    left_join(climatology, by = "month") %>%
    mutate(across(c(temp, rainfall,windspeed,eastward_wind, north_wind, humidity, pop_density, aod, cropland, forest, barren),
                  ~ .x - get(paste0("clim_", cur_column())),
                  .names = "{.col}_anomaly"))
  
  # Inspect results
  head(merged_data)
  
  
  
  
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
  
  
  # Select variables
  vars <- c("rainfall",
            "north_wind", "humidity",      
            "windspeed",
            "cropland","eastward_wind",
            "barren",  "vaccine","temp","forest","aod","pop_density")
  
  # Compute correlation matrix
  corr_matrix <- cor(merged_data[, vars], use = "complete.obs")
  
  
  # Save corrplot as PNG
  png("my_corrplot.png", width = 1200, height = 1000, res = 150)  # Adjust size/resolution if needed
  corrplot(corr_matrix, method = "color", type = "upper", 
           addCoef.col = "black", number.cex = 1.2, tl.cex = 1.2,
           tl.col = "black", col = colorRampPalette(c("#ff0000", "white", "#195696"))(100))
  dev.off()  # Close the PNG device
  
  
  # ==== VIF Analysis ====
  

  merged_data$rainfall_scale<- scale(merged_data$rainfall)
  merged_data$windspeed_scale<- scale(merged_data$windspeed)
  merged_data$eastward_wind_scale<- scale(merged_data$eastward_wind)
  merged_data$north_wind_scale<- scale(merged_data$north_wind)
  merged_data$humidity_scale<- scale(merged_data$humidity)
  merged_data$aod_scale<- scale(merged_data$aod)
  merged_data$pop_density_scale<- scale(merged_data$pop_density)
  merged_data$cropland_scale<- scale(merged_data$cropland)
  merged_data$forest_scale<- scale(merged_data$forest)
  merged_data$barren_scale<- scale(merged_data$barren)
  merged_data$temp_scale<- scale(merged_data$temp)
  
  test_model<-glm(outbreak_2 ~ aod + eastward_wind + north_wind + humidity + windspeed + cropland + barren +forest + vaccine + temp + pop_density_scale, data = merged_data)
  
  vif<-vif(test_model)
  print(vif)
 
  
  
  vars_inla <- c(  "windspeed_scale","north_wind_scale", "humidity_scale",
                 "vaccine", "pop_density_scale", "temp_scale", "temp_lag1",
                 "eastward_wind","eastward_wind_scale","aod_scale", "cropland",
                 "barren", "forest"
  )
  
  
  
  
  
  posterior_table_m2 <- data.frame()
  
  for (var in vars_inla) {
    formula <- as.formula(paste(
      "outbreak_2 ~", var, 
      "+ f(year_month, model = 'rw2') ",
      "+ f(area, model = 'bym2', graph = ken.adj, scale.model = TRUE, constr = TRUE)"
    ))
    
    
    M5 <- inla(
      formula,
      data = merged_data,
      family = "binomial",
      weights=weights, 
      control.predictor = list(compute = TRUE),
      control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)
    )
    
    effect_row <- M5$summary.fixed[grepl(var, rownames(M5$summary.fixed)), ]
    lcpocpo_value <- -mean(log(M5$cpo$cpo), na.rm = TRUE)
    
    posterior_table_m2 <- rbind(posterior_table_m2, data.frame(
      covariate = var,
      estimate = effect_row$mean,
      CI_1 = effect_row$`0.025quant`,
      CI_2 = effect_row$`0.975quant`,
      waic = M5$waic$waic,
      dic = M5$dic$dic,
      cpo = lcpocpo_value
    ))
    
    cat("\n--- INLA Model Summary for:", var, "---\n")
    cat("Mean:", effect_row$mean, "\n")
    cat("95% CI:", effect_row$`0.025quant`, "-", effect_row$`0.975quant`, "\n")
    cat("Significant?:", ifelse(effect_row$`0.025quant` > 0 | effect_row$`0.975quant` < 0, "Yes", "No"), "\n")
    cat("DIC:", M5$dic$dic, ", WAIC:", M5$waic$waic, "\n")
  }
  
  # Final cleanup and export
  rownames(posterior_table_m2) <- 1:nrow(posterior_table_m2)
  print(posterior_table_m2)
  
  
  # Compute the odds ratios (OR) by exponentiating the estimates and the CI bounds
  posterior_table_m2_or <- posterior_table_m2
  posterior_table_m2_or[, c("estimate", "CI_1", "CI_2")] <- exp(posterior_table_m2_or[, c("estimate", "CI_1", "CI_2")])
  print(posterior_table_m2_or)
  
  write_xlsx(posterior_table_m2_or, path = "univariate_analysis.xlsx")
  
}
  