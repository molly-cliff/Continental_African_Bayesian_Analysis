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
library(gridExtra)
library(reshape2)

non_linear_dust_wind_iter <-function() {
    
    # ============================================================
    # 1. LOAD AND PREPARE SPATIOTEMPORAL DATA
    # ============================================================
    
    spatiotemporaloutbreaks <- read.csv("final_popdensity_vaccine.csv")
    print(names(spatiotemporaloutbreaks))
    
    # Create continuous year-month index
    spatiotemporaloutbreaks$year_month <-
      (spatiotemporaloutbreaks$year.x -
         min(spatiotemporaloutbreaks$year.x)) * 12 +
      spatiotemporaloutbreaks$month
    
    # Select required variables
    myvars <- c(
      "year.x", "month", "outbreak_2", "name_2",
      "district_country.x", "country", "code.x",
      "rainfall", "eastward_wind", "north_wind",
      "aod", "humidity", "vaccine", "windspeed",
      "temp", "cropland", "forest", "barren",
      "year_month", "pop_density"
    )
    
    merged_data <- spatiotemporaloutbreaks[, myvars]
    
    
    # ============================================================
    # 2. RESTRICT DATA TO STUDY COUNTRIES
    # ============================================================
    
    study_countries <- c(
      "Benin", "Burkina Faso", "Burundi", "Cameroon",
      "Central African Republic", "Chad", "Côte d'Ivoire",
      "Eritrea", "Ethiopia", "Gambia", "Ghana", "Guinea",
      "Guinea-Bissau", "Kenya", "Mali", "Mauritania",
      "Niger", "Nigeria", "Rwanda", "Senegal",
      "South Sudan", "Sudan", "Tanzania", "Togo", "Uganda"
    )
    
    merged_data <- merged_data[
      merged_data$country %in% study_countries,
    ]
    
    
    # ============================================================
    # 3. PREPARE SPATIAL DATA AND ADJACENCY STRUCTURE
    # ============================================================
    
    shape2 <- st_read("Shapefile_improved.shp")
    
    shape2 <- shape2[
      shape2$COUNTRY %in% study_countries,
    ]
    
    # Fix invalid geometries
    shape2 <- st_make_valid(shape2)
    
    # Convert to Spatial* object
    shapefile_spatial <- as(shape2, "Spatial")
    
    # Create neighbourhood structure
    nb <- poly2nb(shapefile_spatial, queen = TRUE)
    
    nb2INLA("dust_wind.adj", nb)
    ken.adj <- paste0(getwd(), "/dust_wind.adj")
    
    
    # ============================================================
    # 4. MATCH DISTRICTS BETWEEN DATA AND SHAPEFILE
    # ============================================================
    
    shape2$district_country.x <-
      paste(shape2$NAME_2, shape2$COUNTRY, sep = " ")
    
    missing_in_shape <- setdiff(
      merged_data$district_country.x,
      shape2$district_country.x
    )
    
    if (length(missing_in_shape) > 0) {
      stop(
        "districts are in merged_data but not in shapefile: ",
        paste(missing_in_shape, collapse = ", ")
      )
    }
    
    # Check for districts present in shapefile but absent from data
    missing_in_data <- setdiff(
      shape2$district_country.x,
      merged_data$district_country.x
    )
    
    if (length(missing_in_data) > 0) {
      message(
        "Dropping ", length(missing_in_data),
        " districts from shapefile that are not in merged_data:\n",
        paste(missing_in_data, collapse = ", ")
      )
      
      shape2 <- shape2[
        shape2$district_country.x %in%
          merged_data$district_country.x,
      ]
    }
    
    # Confirm district-country lists match
    if (!identical(
      sort(unique(merged_data$district_country.x)),
      sort(unique(shape2$district_country.x))
    )) {
      stop("dstrict-country differ")
    } else {
      message("ll district-country matches verified.")
    }
    
    # Ensure consistent ordering
    shape2 <- shape2[order(shape2$district_country.x), ]
    merged_data <- merged_data[order(merged_data$district_country.x), ]
    
    # Assign spatial IDs
    merged_data$area <-
      match(
        merged_data$district_country.x,
        shape2$district_country.x
      )
    
    
    # ============================================================
    # 5. OUTBREAK FREQUENCY AND OBSERVATION WEIGHTS
    # ============================================================
    
    test <- table(merged_data$outbreak_2)
    test <- as.data.frame(test)
    
    print("OUTBREAK OCCURANCE:")
    print(test)
    
    print("NUMBER OF OUTBREAKS V NO OUTBREAKS:")
    
    outbreak <- test[2, 2]
    no_outbreak <- test[1, 2]
    
    print(outbreak)
    print(no_outbreak)
    
    outbreak <- as.numeric(outbreak)
    no_outbreak <- as.numeric(no_outbreak)
    
    total <- no_outbreak + outbreak
    
    non_outbreak_weight <- outbreak / total
    outbreak_weight <- 1 - non_outbreak_weight
    
    print("OUTBREAK WEIGHTS:")
    print(outbreak_weight)
    
    weights <- ifelse(
      merged_data$outbreak_2 == 1,
      outbreak_weight,
      non_outbreak_weight
    )
    
    
    # ============================================================
    # 6. SCALE AND DERIVE ENVIRONMENTAL VARIABLES
    # ============================================================
    
    merged_data$eastward_wind <-
      as.numeric(merged_data$eastward_wind)
    
    merged_data$eastward_wind_scale <-
      scale(merged_data$eastward_wind)
    
    merged_data$eastward_wind_scale <-
      as.numeric(merged_data$eastward_wind_scale)
    
    merged_data$aod_scale <-
      scale(merged_data$aod)
    
    merged_data$humidity_scale <-
      scale(merged_data$humidity)
    
    merged_data$temp_scale <-
      scale(merged_data$temp)
    
    merged_data$windspeed_scale <-
      scale(merged_data$windspeed)
    
    merged_data$pop_density_scale <-
      scale(merged_data$pop_density)
    
    
    # Create indices for non-linear effects
    merged_data$temp_index <-
      as.numeric(cut(merged_data$temp_scale, breaks = 50))
    
    
    merged_data$aod_index <-
      as.numeric(cut(merged_data$aod_scale, breaks = 50))
    
    
    # ============================================================
    # 7. WIND VARIABLES
    # ============================================================
    
    # Wind strength = absolute eastward wind magnitude
    merged_data$wind_strength <-
      abs(merged_data$eastward_wind)
    
    merged_data$wind_strength_scale <-
      as.numeric(scale(merged_data$wind_strength))
    
    # Wind direction indicator
    # 1 = east -> west
    merged_data$wind_dir <-
      ifelse(merged_data$eastward_wind < 0, 1, 0)
    

    formula1 <- outbreak_2 ~
      aod_scale +
      wind_strength_scale +
      humidity_scale +
      wind_dir +
      aod_scale:wind_dir +
      f(year_month, model = "rw2") +
      f(temp_index, model = "rw2") +
      f(
        area,
        model = "bym2",
        graph = ken.adj,
        scale.model = TRUE,
        constr = TRUE,
        hyper = list(
          prec = list(
            prior = "pc.prec",
            param = c(1, 0.01)
          ),
          phi = list(
            prior = "pc",
            param = c(0.5, 2/3)
          )
        )
      )
    
    model1 <- inla(
      formula1,
      data = merged_data,
      family = "binomial",
      weights = weights,
      control.fixed = list(mean = 0, prec = 1),
      control.predictor = list(
        compute = TRUE,
        link = 1
      ),
      control.compute = list(
        dic = TRUE,
        waic = TRUE,
        cpo = TRUE,
        config = TRUE
      )
    )
    
    
    print(summary(model1))
    
 
    formula2 <- outbreak_2 ~
      aod_scale +
      humidity_scale +
      wind_dir +
      aod_scale:wind_dir +
      f(year_month, model = "rw2") +
      f(temp_index, model = "rw2") +
      f(
        area,
        model = "bym2",
        graph = ken.adj,
        scale.model = TRUE,
        constr = TRUE,
        hyper = list(
          prec = list(
            prior = "pc.prec",
            param = c(1, 0.01)
          ),
          phi = list(
            prior = "pc",
            param = c(0.5, 2/3)
          )
        )
      )
    
    model2 <- inla(
      formula2,
      data = merged_data,
      family = "binomial",
      weights = weights,
      control.fixed = list(mean = 0, prec = 1),
      control.predictor = list(
        compute = TRUE,
        link = 1
      ),
      control.compute = list(
        dic = TRUE,
        waic = TRUE,
        cpo = TRUE,
        config = TRUE
      )
    )
    

    
    print(summary(model2))
    
    
    # ============================================================
    # MODEL COMPARISON
    # ============================================================
    

    
    comparison <- data.frame(
      Model = c(
        "Model 1: AOD + wind strength",
        "Model 2: AOD without wind strength"
      ),
      DIC = c(
        model1$dic$dic,
        model2$dic$dic
      ),
      WAIC = c(
        model1$waic$waic,
        model2$waic$waic
      ),
      CPO = c(
        -mean(log(model1$cpo$cpo), na.rm = TRUE),
        -mean(log(model2$cpo$cpo), na.rm = TRUE)
      )
    )
    
    print(comparison)
    

    
    best_dic <-
      comparison$Model[which.min(comparison$DIC)]
    
    best_waic <-
      comparison$Model[which.min(comparison$WAIC)]
    
    best_cpo <-
      comparison$Model[which.min(comparison$CPO)]
    
    cat("\nBest model by DIC: ", best_dic, "\n")
    cat("Best model by WAIC:", best_waic, "\n")
    cat("Best model by CPO: ", best_cpo, "\n")
    

  
  
}
  