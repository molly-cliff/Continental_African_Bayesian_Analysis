library(dplyr)
library(INLA)
library(sf)
library(spdep)
library(sp)
library(pROC)
library(PRROC)


  # ============================================================
  # 1. LOAD AND PREPARE OUTBREAK DATA
  # ============================================================

  # Read the joined outbreak and environmental dataset.
  spatiotemporaloutbreaks <- read.csv("final_popdensity_vaccine.csv")
  print(names(spatiotemporaloutbreaks))
  
  # Create a continuous month index for the temporal RW2 effect.
  spatiotemporaloutbreaks$year_month <- (spatiotemporaloutbreaks$year.x - min(spatiotemporaloutbreaks$year.x)) * 12 + spatiotemporaloutbreaks$month
  
  # Retain variables used by the spatial and statistical models.
  myvars <- c("year.x", "month", "outbreak_2", "name_2", "district_country.x", "country",
              "code.x","rainfall","eastward_wind","north_wind", "aod","humidity","vaccine"          
              ,"windspeed", "temp","cropland",          
              "forest", "barren","year_month", "pop_density" )
  merged_data <- spatiotemporaloutbreaks[, myvars]

  # Restrict observations to the study countries.
  merged_data <- merged_data[merged_data$country %in% c("Benin", "Burkina Faso", "Burundi", "Cameroon", "Central African Republic",
                                                        "Chad", "Côte d'Ivoire", "Eritrea", "Ethiopia",
                                                        "Gambia", "Ghana", "Guinea", "Guinea-Bissau", "Kenya", "Mali", "Mauritania",
                                                        "Niger", "Nigeria", "Rwanda", "Senegal", "South Sudan", "Sudan", "Tanzania",
                                                        "Togo", "Uganda"), ]
  
  
  
  # ============================================================
  # 2. PREPARE SPATIAL DATA AND ADJACENCY STRUCTURE
  # ============================================================

  # Read and restrict the district shapefile to the study countries.
  shape2 <- st_read("Shapefile_improved.shp")
  
  shape2<- shape2[shape2$COUNTRY %in% c("Benin", "Burkina Faso", "Burundi", "Cameroon", "Central African Republic",
                                        "Chad", "Côte d'Ivoire", "Eritrea", "Ethiopia",
                                        "Gambia", "Ghana", "Guinea", "Guinea-Bissau", "Kenya", "Mali", "Mauritania",
                                        "Niger", "Nigeria", "Rwanda", "Senegal", "South Sudan", "Sudan", "Tanzania",
                                        "Togo", "Uganda"), ]
  
  # Repair invalid geometries before building neighbourhoods.
  shape2 <- st_make_valid(shape2)  # Optionally fix geometries
  
  # Optional: Remove invalid ones instead
  # shape2 <- shape2[st_is_valid(shape2), ]
  
  # Convert to a Spatial object for neighbourhood construction.
  shapefile_spatial <- as(shape2, "Spatial")
  
  # Create queen-contiguity neighbourhoods and write the INLA graph.
  nb <- poly2nb(shapefile_spatial, queen = TRUE)
  
  nb2INLA("platt-scale.adj", nb)
  ken.adj <- paste0(getwd(), "/platt-scale.adj")
  
  # ============================================================
  # 3. VALIDATE AND ALIGN DISTRICTS
  # ============================================================

  # Build the district-country key used to align both data sources.
  shape2$district_country.x <- paste(shape2$NAME_2, shape2$COUNTRY, sep = " ")
  
  # Check for observations that have no matching spatial polygon.
  missing_in_shape <- setdiff(merged_data$district_country.x, shape2$district_country.x)
  if (length(missing_in_shape) > 0) {
    stop("districts in merged_data but not in shapefile: ",
         paste(missing_in_shape, collapse = ", "))
  }
  
  # Drop districts present in the shapefile but absent from the data.
  missing_in_data <- setdiff(shape2$district_country.x, merged_data$district_country.x)
  
  if (length(missing_in_data) > 0) {
    message("dropping ", length(missing_in_data), 
            " districts from shapefile that are not in merged_data:\n",
            paste(missing_in_data, collapse = ", "))
    
    # Keep only districts represented in both sources.
    shape2 <- shape2[shape2$district_country.x %in% merged_data$district_country.x, ]
  }
  
  
  # Confirm that district-country identifiers match.
  if (!identical(sort(unique(merged_data$district_country.x)),
                 sort(unique(shape2$district_country.x)))) {
    stop("district-country lists differ even after sorting!")
  } else {
    message("all district-country matches verified.")
  }
  
  # Sort both sources identically before assigning spatial IDs.
  shape2 <- shape2[order(shape2$district_country.x), ]
  merged_data <- merged_data[order(merged_data$district_country.x), ]
  
  # Assign the integer area index used by the BYM2 effect.
  merged_data$area <- match(merged_data$district_country.x, shape2$district_country.x)
  
  
  
  # ============================================================
  # 4. CALCULATE OUTCOME WEIGHTS
  # ============================================================

  # Count outbreak and non-outbreak observations.
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
  # Give the minority class a larger contribution to the likelihood.
  total<-no_outbreak+outbreak
  non_outbreak_weight <- outbreak/total
  outbreak_weight<-1-non_outbreak_weight
  
  print("OUTBREAK WEIGHTS:")
  print(outbreak_weight)
  weights <- ifelse(merged_data$outbreak_2 == 1, outbreak_weight, non_outbreak_weight) 
  
  
  # ============================================================
  # 5. SCALE COVARIATES AND DERIVE WIND VARIABLES
  # ============================================================

  # Standardise environmental predictors used as fixed effects.
  merged_data$eastward_wind_scale <- scale(merged_data$eastward_wind)
  merged_data$aod_scale <- scale(merged_data$aod)
  merged_data$humidity_scale <- scale(merged_data$humidity)
  merged_data$temp_scale <- scale(merged_data$temp)

  # Convert temperature into bins for its nonlinear RW2 effect.
  merged_data$temp_index <- as.numeric(cut(merged_data$temp_scale, breaks = 50))
  
  # Derive wind magnitude and direction predictors.
  merged_data$wind_strength <- abs(merged_data$eastward_wind)
  merged_data$wind_strength_scale <- as.numeric(scale(merged_data$wind_strength))
  # Wind direction indicator (1 = east->west)
  merged_data$wind_dir <- ifelse(merged_data$eastward_wind < 0, 1, 0)
  
 
  
  # ============================================================
  # 6. FIT FULL MODEL WITH AOD-WIND INTERACTION
  # ============================================================

  # Specify the temporal, nonlinear temperature, and spatial effects.
  formula1 <- outbreak_2 ~ wind_strength_scale + humidity_scale +
    aod_scale * wind_dir +
    f(year_month, model = "rw2",
      hyper = list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))) +
    f(temp_index, model = "rw2",
      hyper = list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))) +
    f(area, model = "bym2", graph = ken.adj, scale.model = TRUE, constr = TRUE,
      hyper = list(
        prec = list(prior = "pc.prec", param = c(1, 0.01)),
        phi  = list(prior = "pc", param = c(0.5, 2/3))
      ))
  
  # Fit the full weighted binomial INLA model.
  model <- INLA::inla(
    formula1, data = merged_data, family = "binomial",
    weights = weights,
    control.fixed = list(mean = 0, prec = 1),
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)
  )
  

  # Inspect model estimates and fit criteria.
  print(summary(model))
  

  
  # ============================================================
  # 7. CREATE OUTBREAK EPISODES AND ASSIGN FOLDS
  # ============================================================
  k <- 5  # number of folds
  set.seed(123)
  
  # Identify consecutive outbreak episodes within each district.
  merged_data <- merged_data  %>%
    arrange(district_country.x, year.x, month) %>%
    group_by(district_country.x) %>%
    mutate(
      outbreak_prev = lag(outbreak_2, default = 0),                     # previous month outbreak
      episode_start = (outbreak_2 == 1 & outbreak_prev == 0),          # TRUE at start of new episode
      episode_counter = cumsum(as.integer(episode_start)),                  # running count of episodes per district
      episode_id = if_else(outbreak_2 == 1,                             
                           paste0(district_country.x, "_ep", episode_counter),
                           NA_character_)                                   # assign ID only to outbreak months
    ) %>%
    ungroup()
  
  # Assign the same fold to all months in an outbreak episode.
  episode_folds <- merged_data %>%
    filter(outbreak_2 == 1) %>%
    distinct(episode_id) %>% 
    mutate(fold_space = sample(rep(1:k, length.out = n())))  # random k-fold assignment
  
  # Assign non-outbreak months to folds independently.
  non_outbreak_folds <- merged_data %>%
    filter(outbreak_2 == 0) %>%
    mutate(row_id = row_number(),
           fold_space_non = sample(rep(1:k, length.out = n())))
  
  # Merge fold assignments back into the analysis data.
  merged_data <- merged_data %>%
    left_join(episode_folds, by = "episode_id") %>%                # outbreak folds
    mutate(row_id = row_number()) %>%                               # temporary row identifier
    left_join(non_outbreak_folds %>% dplyr::select(row_id, fold_space_non), by = "row_id") %>%
    mutate(
      fold_space = ifelse(!is.na(fold_space), fold_space, fold_space_non)   # combine folds
    ) %>%
    dplyr::select(-fold_space_non, -row_id, -outbreak_prev, -episode_start, -episode_counter)
  
  
  # Validate that outbreak episodes retained their assigned folds.
  check_episode_folds <- merged_data %>%
    filter(!is.na(episode_id)) %>%                       # only outbreak months
    left_join(episode_folds, by = "episode_id") %>%      # bring in original assigned fold
    mutate(fold_match = (fold_space.x == fold_space.y))  # compare merged_data fold vs episode_folds
  
  # Record any mismatches for inspection.
  mismatches <- check_episode_folds %>%
    filter(fold_match == FALSE) %>%
    dplyr::select(district_country.x, year.x, month, episode_id, fold_space.x, fold_space.y)
  
  
  correct_match <- check_episode_folds %>%
    filter(fold_match ==TRUE) %>%
    dplyr::select(district_country.x, year.x, month, episode_id, fold_space.x, fold_space.y)
  
  
  # ============================================================
  # 8. INITIALISE CROSS-VALIDATION OUTPUTS
  # ============================================================
  all_roc_objs <- list()
  metrics_list <- data.frame()
  
  eps <- 1e-15
  
  all_metrics <- list()
  all_test_preds <- list()
  # ============================================================
  # 9. FIT AND EVALUATE EACH TEST FOLD
  # ============================================================
  for (k in 1:5) {
    
    message("==== Fold ", k, " / 5 ====")
    
    # Split the data into training and test indices.
    data_copy <- merged_data
    test_indices  <- which(data_copy$fold_space == k)
    train_indices <- which(data_copy$fold_space != k)
    
    # Mask test outcomes so INLA fits using training observations only.
    data_copy$outbreak_2[test_indices] <- NA
    
    # Fit INLA on the training data by masking test outcomes.
    message("Fitting INLA with fold_space == ", k, " masked...")
    M_fold <- inla(
      outbreak_2 ~ aod_scale +
        wind_strength_scale +
        humidity_scale +
        wind_dir +
        aod_scale:wind_dir +
        f(year_month, model = "rw2",
          hyper = list(
            prec = list(prior = "pc.prec", param = c(0.5, 0.01))
          )
        ) +
        f(temp_index, model = "rw2",
          hyper = list(
            prec = list(prior = "pc.prec", param = c(0.5, 0.01))
          )
        ) +
        f(area, model = "bym2", graph = ken.adj, scale.model = TRUE, constr = TRUE,
          hyper = list(
            prec = list(prior = "pc.prec", param = c(1, 0.01)),
            phi  = list(prior = "pc", param = c(0.5, 2/3))
          )
        ),
      data = data_copy,
      family = "binomial",
      weights = data_copy$weights, 
      control.fixed = list(mean = 0, prec = 1),
      control.predictor = list(compute = TRUE, link = 1),
      control.compute = list(dic = FALSE, waic = FALSE)
    )
    
    # ------------------------------------------------------------
    # 9.1 EXTRACT HELD-OUT FOLD PREDICTIONS
    # ------------------------------------------------------------
    # INLA returns fitted probabilities for every row. Keep only the
    # observations belonging to the current test fold for evaluation.
    p_all <- M_fold$summary.fitted.values$mean
    
    p_test  <- p_all[test_indices]
    y_test  <- merged_data$outbreak_2[test_indices]
    w_test  <- merged_data$weights[test_indices]
    
    # the idea is here for log loss, that predicted values that are either extremely close to 0 or 1 will be clipped to avoid log(0) issues.
    # important here due to the rareity of outbreaks, which can lead to extreme predicted probabilities.
    p0_test  <- pmin(pmax(as.numeric(p_test),  eps), 1 - eps)
    
    # ------------------------------------------------------------
    # 9.2 CALCULATE WEIGHTED TEST-FOLD PERFORMANCE
    # ------------------------------------------------------------

    # class weights should be comparable across folds. This is important for the Brier score and log loss.
    wN <- as.numeric(w_test) / sum(as.numeric(w_test))
    
    # Brier score: weighted squared error of the predicted probabilities.
    brier <- sum(wN * (p0_test - y_test)^2)

    # Log loss: weighted negative log likelihood for binary outcomes.
    logloss <- -sum(wN * (y_test * log(p0_test) +
                            (1 - y_test) * log(1 - p0_test)))
    
    # ------------------------------------------------------------
    # 9.3 CALCULATE AND STORE ROC PERFORMANCE
    # ------------------------------------------------------------
    # Build a weighted ROC curve and extract its area under the curve.
    roc_obj <- pROC::roc(y_test, p0_test, weights = wN, quiet = TRUE)
    auc_val <- as.numeric(pROC::auc(roc_obj))

    # Retain the complete ROC object for the cross-fold plot later.
    all_roc_objs[[paste0("Fold", k)]] <- roc_obj

    # ------------------------------------------------------------
    # 9.4 CALCULATE HIGH-RISK ALERT PERFORMANCE
    # ------------------------------------------------------------
    # Evaluate how well the model identifies outbreaks among the
    # highest-risk 1%, 5%, and 10% of test observations.
    alert_metrics <- function(risk, y, ks = c(0.01, 0.05, 0.10)) {
      # Remove invalid predictions or missing outcomes before scoring.
      valid <- is.finite(risk) & !is.na(y)
      risk <- as.numeric(risk[valid])
      y <- as.numeric(y[valid])

      # Establish the observed outbreak prevalence and number of positives.
      prevalence <- mean(y == 1)
      n_pos <- sum(y == 1)

      dplyr::bind_rows(lapply(ks, function(alert_fraction) {
        # Select the requested top-risk fraction using its risk threshold.
        threshold <- as.numeric(
          stats::quantile(risk, 1 - alert_fraction, na.rm = TRUE)
        )
        alert <- risk >= threshold

        # Compare alerts with observed outbreaks.
        true_positive <- sum(alert & y == 1)
        n_alert <- sum(alert)
        precision <- if (n_alert > 0) true_positive / n_alert else NA_real_
        recall <- if (n_pos > 0) true_positive / n_pos else NA_real_

        # Return alert threshold, classification metrics, and lift.
        data.frame(
          alert_fraction = alert_fraction,
          threshold = threshold,
          prevalence = prevalence,
          precision = precision,
          recall = recall,
          lift = precision / prevalence,
          stringsAsFactors = FALSE
        )
      }))
    }

    # Apply the alert analysis to the current test fold, risk refers to the predicted probabilities and y to actual outcomes.
    alert_df <- alert_metrics(
      risk = p0_test,
      y = y_test,
      ks = c(0.01, 0.05, 0.10)
    )

    # Add fold-level model metrics to every alert-fraction row, this is later seen in output.
    metrics_df <- dplyr::mutate(
      alert_df,
      fold = k,
      roc_auc = auc_val,
      brier = brier,
      logloss = logloss,
      .before = 1
    )
    
    # Store all metrics for later cross-fold aggregation.
    all_metrics[[k]] <- metrics_df
    
    # ------------------------------------------------------------
    # 9.5 SAVE TEST-FOLD PREDICTIONS AND METRICS
    # ------------------------------------------------------------
    # Keep the observed outcome, predicted probability, and weight
    # for every observation in the held-out fold.
    preds_df <- data.frame(
      fold = k,
      truth = as.numeric(y_test),
      probability = p0_test,
      weight = as.numeric(w_test),
      stringsAsFactors = FALSE
    )
    all_test_preds[[k]] <- preds_df

    # Write fold-specific prediction and metric files.
    out_prefix <- paste0("fold", k, "_cv_ipw")
    write.csv(preds_df,    paste0(out_prefix, "_preds.csv"),       row.names = FALSE)
    write.csv(metrics_df,  paste0(out_prefix, "_metrics.csv"),     row.names = FALSE)
  }
  
  # ============================================================
  # 10. SUMMARISE ROC CURVES ACROSS FOLDS
  # ============================================================

  # Extract the AUC for each fold and calculate the mean AUC.
  fold_auc <- vapply(all_roc_objs, function(roc_obj) {
    as.numeric(pROC::auc(roc_obj))
  }, numeric(1))
  mean_auc <- mean(fold_auc, na.rm = TRUE)

  # Plot all fold-specific ROC curves in one image.
  png("roc_curves_all_folds.png", width = 8, height = 6, units = "in", res = 300)
  plot(all_roc_objs[[1]], col = 1, lwd = 1.5,
       main = "ROC Curves - All Folds")
  abline(a = 0, b = 1, lty = 2, col = "gray")

  # Add each remaining fold curve and label it with its AUC.
  fold_colors <- seq_along(all_roc_objs)
  for (i in 2:length(all_roc_objs)) {
    lines(all_roc_objs[[i]], col = fold_colors[i], lwd = 1.5)
  }
  legend("bottomright",
         legend = c(paste0(names(all_roc_objs), " (AUC = ",
                            round(fold_auc, 3), ")"),
                    paste0("Mean AUC = ", round(mean_auc, 3))),
         col = c(fold_colors, NA),
         lwd = c(rep(1.5, length(all_roc_objs)), NA),
         cex = 0.8)
  dev.off()
  
  # ============================================================
  # 11. CALCULATE OVERALL CROSS-VALIDATION SUMMARIES
  # ============================================================

  # Combine alert metrics from all folds into one table.
  all_metrics_df <- do.call(rbind, all_metrics)
  rownames(all_metrics_df) <- NULL

  # Save the combined cross-validation results.
  write.csv(all_metrics_df, "cv_all_folds_metricstest.csv", row.names = FALSE)
  
  