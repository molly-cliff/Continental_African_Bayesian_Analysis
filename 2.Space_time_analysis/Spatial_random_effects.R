# Load required libraries
library(dplyr)
library(INLA)
library(sf)
library(spdep)
library(sp)
library(caret)
library(pROC)
library(Hmisc)
library(tidyr)
library(PRROC)
library(broom)
library(corrplot)
library(car)
library(ggplot2)
library(writexl)
library(reshape2)



monthly_models_test3 <- function() {
  # Read data
  
  spatiotemporaloutbreaks <- read.csv("final_popdensity_vaccine.csv")
  print(names(spatiotemporaloutbreaks))
  
  
  spatiotemporaloutbreaks$year_month <- (spatiotemporaloutbreaks$year.x - min(spatiotemporaloutbreaks$year.x)) * 12 + spatiotemporaloutbreaks$month
  
  myvars <- c("year.x", "month", "outbreak_2", "name_2", "district_country.x", "country",
              "code.x","year_month")
  merged_data <- spatiotemporaloutbreaks[, myvars]
  merged_data <- merged_data[merged_data$country %in% c("Benin", "Burkina Faso", "Burundi", "Cameroon", "Central African Republic",
                                                        "Chad", "Côte d'Ivoire", "Eritrea", "Ethiopia", "Democratic Republic of the Congo",
                                                        "Gambia", "Ghana", "Guinea", "Guinea-Bissau", "Kenya", "Mali", "Mauritania",
                                                        "Niger", "Nigeria", "Rwanda", "Senegal", "South Sudan", "Sudan", "Tanzania",
                                                        "Togo", "Uganda"), ]
  
  
  
  shape2 <- st_read("Shapefile_improved.shp")
  
  shape2<- shape2[shape2$COUNTRY %in% c("Benin", "Burkina Faso", "Burundi", "Cameroon", "Central African Republic",
                                        "Chad", "Côte d'Ivoire", "Eritrea", "Ethiopia", "Democratic Republic of the Congo",
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
  
  nb2INLA("full_africa_map2_test.adj", nb)
  ken.adj <- paste0(getwd(), "/full_africa_map2_test.adj")
  
  shape2$district_country.x <- paste(shape2$NAME_2, shape2$COUNTRY, sep = " ")
  
  merged_data$area<-as.numeric(as.factor(merged_data$district_country.x))
  
  missing_in_shape <- setdiff(merged_data$district_country.x, shape2$district_country.x)
  if (length(missing_in_shape) > 0) {
    stop("These districts are in merged_data but not in shapefile: ",
         paste(missing_in_shape, collapse = ", "))
  }
  
  # --- 2. In shape but not in data ---
  missing_in_data <- setdiff(shape2$district_country.x, merged_data$district_country.x)
  
  if (length(missing_in_data) > 0) {
    message("Dropping ", length(missing_in_data), 
            " districts from shapefile that are not in merged_data:\n",
            paste(missing_in_data, collapse = ", "))
    
    # Keep only those in both
    shape2 <- shape2[shape2$district_country.x %in% merged_data$district_country.x, ]
  }
  
  
  # --- 3. Check sorted lists match ---
  if (!identical(sort(unique(merged_data$district_country.x)),
                 sort(unique(shape2$district_country.x)))) {
    stop("District-country lists differ even after sorting")
  } else {
    message("All district-country matches verified.")
  }
  
  # Ensure consistent ordering
  shape2 <- shape2[order(shape2$district_country.x), ]
  merged_data <- merged_data[order(merged_data$district_country.x), ]
  
  # Assign IDs
  merged_data$area <- match(merged_data$district_country.x, shape2$district_country.x)
  
  
  shape2 <- shape2[order(shape2$district_country.x), ]
  merged_data <- merged_data[order(merged_data$district_country.x), ]
  merged_data$area <- match(merged_data$district_country.x, shape2$district_country.x)
  
  
  # -------------------------
  # 4. Hyperparameters for spatial models
  # -------------------------
  # Note: BYM2 typically uses a specific parameterization; here we set a PC prior on precision.
  shyper <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
  # adjust the param values to suit prior beliefs (scale, tail prob)
  
  # -------------------------
  # 5. Model formulas and fits
  # -------------------------
  # Use an index for the AR1 time effect (time_id)
  # Use family = "binomial" with outbreak_2 (0/1)
  # Random effect for area uses the graph ' ken.adj'
  
  # Helper function to compute -mean(log(cpo)) safely
  compute_logcpo <- function(model) {
    cpo <- model$cpo$cpo
    cpo <- cpo[!is.na(cpo) & cpo > 0]
    if (length(cpo) == 0) return(NA_real_)
    return(-mean(log(cpo)))
  }
  
  # Model 1: BYM2
  formula_bym2 <- outbreak_2 ~
    f(year_month, model = "rw2") +
    f(area, model = "bym2", graph =  ken.adj, hyper = shyper, scale.model = TRUE)
  
  fit_bym2 <- inla(
    formula_bym2,
    data = merged_data,
    family = "binomial",
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
  )
  
  # Model 2: BYM (classic)
  formula_bym <- outbreak_2 ~
    f(year_month, model = "rw2") +
    f(area, model = "bym", graph =  ken.adj, scale.model = TRUE)
  
  fit_bym <- inla(
    formula_bym,
    data = merged_data,
    family = "binomial",
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
  )
  
  # Model 3: Besag (ICAR)
  formula_besag <- outbreak_2 ~
    f(year_month, model = "rw2") +
    f(area, model = "besag", graph =  ken.adj, scale.model = TRUE)
  
  fit_besag <- inla(
    formula_besag,
    data = merged_data,
    family = "binomial",
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
  )
  
  # Model 4: Besag proper
  formula_besagproper <- outbreak_2 ~
    f(year_month, model = "rw2") +
    f(area, model = "besagproper", graph =  ken.adj)
  
  fit_besagproper <- inla(
    formula_besagproper,
    data = merged_data,
    family = "binomial",
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
  )
  
  # -------------------------
  # 6. Model diagnostics and comparison
  # -------------------------
  results <- data.frame(
    model = c("bym", "bym2", "besag", "besagproper"),
    DIC = c(fit_bym$dic$dic, fit_bym2$dic$dic, fit_besag$dic$dic, fit_besagproper$dic$dic),
    WAIC = c(fit_bym$waic$waic, fit_bym2$waic$waic, fit_besag$waic$waic, fit_besagproper$waic$waic),
    LogCPO = c(
      compute_logcpo(fit_bym),
      compute_logcpo(fit_bym2),
      compute_logcpo(fit_besag),
      compute_logcpo(fit_besagproper)
    ),
    stringsAsFactors = FALSE
  )
  
  print(results)
  write_xlsx(results, path = "model_comparison_spatial.xlsx")
  
  # Print summaries (fixed and random (area) summaries)
  print("BYM2 fixed effects:")
  print(fit_bym2$summary.fixed)
  print("BYM2 area random summary (first rows):")
  print(head(fit_bym2$summary.random$area))
  
  print("BYM fixed effects:")
  print(fit_bym$summary.fixed)
  print("BESAG fixed effects:")
  print(fit_besag$summary.fixed)
  print("BESAGPROPER fixed effects:")
  print(fit_besagproper$summary.fixed)
  
}
