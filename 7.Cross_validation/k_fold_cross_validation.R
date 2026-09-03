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
library(cowplot)
library(ggplot2)
library(reshape2)
library(stringr)
kfold_cross_val <-function() {
  
  spatiotemporaloutbreaks <- read.csv("final_popdensity_vaccine.csv")
  spatiotemporaloutbreaks$year_month <-
    (spatiotemporaloutbreaks$year.x - min(spatiotemporaloutbreaks$year.x, na.rm = TRUE)) * 12 +
    spatiotemporaloutbreaks$month
  
  myvars <- c(
    "year.x", "month", "outbreak_2", "name_2", "district_country.x", "country",
    "code.x", "rainfall", "eastward_wind", "north_wind", "aod", "humidity", "vaccine",
    "windspeed", "temp", "cropland", "forest", "barren", "year_month", "pop_density"
  )
  
  merged_data <- spatiotemporaloutbreaks[, myvars]
  
  keep_countries <- c(
    "Benin", "Burkina Faso", "Burundi", "Cameroon", "Central African Republic",
    "Chad", "Côte d'Ivoire", "Eritrea", "Ethiopia",
    "Gambia", "Ghana", "Guinea", "Guinea-Bissau", "Kenya", "Mali", "Mauritania",
    "Niger", "Nigeria", "Rwanda", "Senegal", "South Sudan", "Sudan", "Tanzania",
    "Togo", "Uganda"
  )
  merged_data <- merged_data[merged_data$country %in% keep_countries, ]
  
  # ---- Read shapefile & build adjacency ----
  shape2 <- st_read("Shapefile_improved.shp", quiet = TRUE)
  shape2 <- shape2[shape2$COUNTRY %in% keep_countries, ]
  shape2 <- st_make_valid(shape2)
  
  shapefile_spatial <- as(shape2, "Spatial")
  nb <- poly2nb(shapefile_spatial, queen = TRUE)
  
  ken.adj <- file.path(tempdir(), "dust_wind2.adj")
  spdep::nb2INLA(ken.adj, nb)
  stopifnot(file.exists(ken.adj))
  
  shape2$district_country.x <- paste(shape2$NAME_2, shape2$COUNTRY, sep = " ")
  
  # ---- Align districts between data and shape ----
  missing_in_shape <- setdiff(merged_data$district_country.x, shape2$district_country.x)
  if (length(missing_in_shape) > 0) {
    stop("These districts are in merged_data but not in shapefile: ",
         paste(missing_in_shape, collapse = ", "))
  }
  
  missing_in_data <- setdiff(shape2$district_country.x, merged_data$district_country.x)
  if (length(missing_in_data) > 0) {
    message("Dropping ", length(missing_in_data), " districts from shapefile not in merged_data.")
    shape2 <- shape2[shape2$district_country.x %in% merged_data$district_country.x, ]
  }
  
  if (!identical(sort(unique(merged_data$district_country.x)),
                 sort(unique(shape2$district_country.x)))) {
    stop("District-country lists differ even after filtering/sorting.")
  }
  
  shape2 <- shape2[order(shape2$district_country.x), ]
  merged_data <- merged_data[order(merged_data$district_country.x), ]
  
  merged_data$area <- match(merged_data$district_country.x, shape2$district_country.x)
  
  # ---- Class weights ----
  tab <- as.data.frame(table(merged_data$outbreak_2))
  outbreak <- as.numeric(tab[tab$Var1 == 1, "Freq"])
  no_outbreak <- as.numeric(tab[tab$Var1 == 0, "Freq"])
  total <- outbreak + no_outbreak
  
  non_outbreak_weight <- outbreak / total
  outbreak_weight <- 1 - non_outbreak_weight
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

  k <- 5  # number of folds
  set.seed(123)
  
  #  Identify outbreak episodes per district
  # Each episode is a run of consecutive months with outbreak_2 == 1
  merged_data <- merged_data  %>%
    arrange(district_country.x, year.x, month) %>%
    group_by(district_country.x) %>%
    mutate(
      outbreak_prev = lag(outbreak_2, default = 0),                    
      episode_start = (outbreak_2 == 1 & outbreak_prev == 0),        
      episode_counter = cumsum(as.integer(episode_start)),                
      episode_id = if_else(outbreak_2 == 1,                             
                           paste0(district_country.x, "_ep", episode_counter),
                           NA_character_)                            
    ) %>%
    ungroup()
  
  # ================================
  # Assign k-folds to outbreak episodes
  # All months in the same outbreak episode get the same fold
  # ================================
  
  
  episode_folds <- merged_data %>%
    filter(outbreak_2 == 1) %>%
    distinct(episode_id) %>% 
    mutate(fold_space = sample(rep(1:k, length.out = n())))  # random k-fold assignment
  
  # ================================
  #Assign k-folds to non-outbreak months independently
  # ================================
  non_outbreak_folds <- merged_data %>%
    filter(outbreak_2 == 0) %>%
    mutate(row_id = row_number(),
           fold_space_non = sample(rep(1:k, length.out = n())))
  
  # ================================
  #  Merge fold assignments back into merged_data
  # ================================
  merged_data <- merged_data %>%
    left_join(episode_folds, by = "episode_id") %>%             
    mutate(row_id = row_number()) %>%                             
    left_join(non_outbreak_folds %>% dplyr::select(row_id, fold_space_non), by = "row_id") %>%
    mutate(
      fold_space = ifelse(!is.na(fold_space), fold_space, fold_space_non)   # combine folds
    ) %>%
    dplyr::select(-fold_space_non, -row_id, -outbreak_prev, -episode_start, -episode_counter)
  
  
  # ================================
  # Check fold assignments against original episode_folds
  # ================================
  check_episode_folds <- merged_data %>%
    filter(!is.na(episode_id)) %>%                      
    left_join(episode_folds, by = "episode_id") %>%     
    mutate(fold_match = (fold_space.x == fold_space.y))  
  
  # Inspect any mismatches
  mismatches <- check_episode_folds %>%
    filter(fold_match == FALSE) %>%
    dplyr::select(district_country.x, year.x, month, episode_id, fold_space.x, fold_space.y)
  
  
  correct_match <- check_episode_folds %>%
    filter(fold_match ==TRUE) %>%
    dplyr::select(district_country.x, year.x, month, episode_id, fold_space.x, fold_space.y)
  
  
  all_roc_objs <- list()
  metrics_list <- data.frame()
  
  #################FOLD 1 ##################################### 
  
  
  k <- 1
  
  data_copy <- merged_data
  test_indices <- which(data_copy$fold_space == k)
  data_copy$outbreak_2[test_indices] <- NA  

  M_fold2 <- inla(
    
    
    outbreak_2 ~ aod_scale +
      wind_strength_scale +
      humidity_scale +
      wind_dir +                 
      aod_scale:wind_dir +         
      f(year_month, model="rw2",
        hyper = list(
          prec = list(
            prior = "pc.prec",
            param = c(0.5, 0.01)
          )
        )
      ) +
      f(temp_index, model="rw2",
        hyper = list(
          prec = list(
            prior = "pc.prec",
            param = c(0.5, 0.01)
          )
        )
      ) +
      f(area, model="bym2", graph = ken.adj, scale.model = TRUE, constr = TRUE,
        hyper = list(
          prec = list(prior = "pc.prec", param = c(1, 0.01)),
          phi  = list(prior = "pc", param = c(0.5, 2/3))
        )),
    data = data_copy,
    family = "binomial",
    weights = weights,
    control.fixed = list(mean = 0, prec = 1),
    control.predictor = list(compute = TRUE, link = 1),
    control.compute = list(dic = FALSE, waic = FALSE)
  )
  
  # Evaluate model on test fold (fold_time == 2)
  predicted_probs <- M_fold2$summary.fitted.values$mean[test_indices]
  actual_vals <- merged_data$outbreak_2[test_indices]
  predicted_classes <- ifelse(predicted_probs > 0.2, 1, 0)
  
  # Store predictions
  predictions_fold2 <- data.frame(
    row_id = test_indices,
    actual = actual_vals,
    predicted_prob = predicted_probs,
    predicted_class = predicted_classes,
    district = merged_data$district_country.x[test_indices]
  )
  
  # Create confusion matrix
  confusion_matrix_fold2 <- table(
    Actual = actual_vals,
    Predicted = predicted_classes
  )
  
  # Output

  print(confusion_matrix_fold2)
  
  # Precision, Recall, F1
  tp <- sum(predicted_classes == 1 & actual_vals == 1, na.rm = TRUE)
  fp <- sum(predicted_classes == 1 & actual_vals == 0, na.rm = TRUE)
  fn <- sum(predicted_classes == 0 & actual_vals == 1, na.rm = TRUE)
  tn <- sum(predicted_classes == 0 & actual_vals == 0, na.rm = TRUE)
  specificity <- tn / (tn + fp)
  recall <- tp / (tp + fn)
  correlation_coefficient <- cor(actual_vals, predicted_probs)


  # ROC AUC
  roc_obj <- roc(actual_vals, predicted_probs)
  auc_val <- auc(roc_obj)
  
  # PR AUC
  pr_obj <- pr.curve(scores.class0 = predicted_probs[actual_vals == 1],
                     scores.class1 = predicted_probs[actual_vals == 0],
                     curve = FALSE)
  pr_auc <- pr_obj$auc.integral
  
  # -----------------------
  # Brier Score
  # -----------------------
  brier_score <- mean((predicted_probs - actual_vals)^2)
  brier_score
 
  
  # -----------------------
  # LogLoss
  # -----------------------
  # Prevent log(0)
  eps <- 1e-15
  y_pred_clipped <- pmin(pmax(predicted_probs, eps), 1 - eps)  
  
  logloss <- -mean(actual_vals * log(y_pred_clipped) + (1 - actual_vals) * log(1 - y_pred_clipped))
  logloss
  
  # Output
  cat("\nConfusion Matrix for Fold 2:\n")
  print(confusion_matrix_fold2)
  cat("\nRecall (Sensitivity):", round(recall, 3))
  cat("\nROC AUC:", round(auc_val, 3))
  cat("\brier:", round (brier_score, 3), "\n")
  cat("logloss:", round (logloss, 3), "\n")
  
  
  
  all_roc_objs[[paste0("Fold", k)]] <- roc_obj
  metrics_list <- rbind(metrics_list,
                        data.frame(Fold = k,
                                   AUC = auc_val,
                                 specificity = specificity,
                                   Recall = recall,
                                   Brier = brier_score,
                                   LogLoss = logloss))
  
  
  alert_metrics <- function(  M_fold2, merged_data, ks = c(0.01, 0.05, 0.10)) {
    risk <-   M_fold2$summary.linear.predictor$mean
    y <- data_copy$outbreak_2
    
    prev <- mean(y == 1, na.rm = TRUE)
    n_pos <- sum(y == 1, na.rm = TRUE)
    
    dplyr::bind_rows(lapply(ks, function(k) {
      thr <- as.numeric(stats::quantile(risk, 1 - k, na.rm = TRUE))
      alert <- risk >= thr
      
      tp <- sum(alert & y == 1, na.rm = TRUE)
      n_alert <- sum(alert, na.rm = TRUE)
      
      precision <- tp / n_alert
      recall <- tp / n_pos
      
      data.frame(
        top_frac = k,
        threshold = thr,
        prevalence = prev,
        precision = precision,
        recall = recall,
        lift = precision / prev
      )
    }))
  }
  
  met <- alert_metrics(M_fold2, data_copy)
  utils::write.csv(met, "alert_tradeoff_metricstestfold1.csv", row.names = FALSE)
  

  
  #################FOLD 2 ##################################### 
  k <- 2
  

  data_copy <- merged_data
  test_indices <- which(data_copy$fold_space == k)
  data_copy$outbreak_2[test_indices]

  M_fold2 <- inla(
    outbreak_2 ~ aod_scale +
      wind_strength_scale +
      humidity_scale +
      wind_dir +                    # Main effect
      aod_scale:wind_dir +          # Interaction: AOD effect varies by wind direction
      f(year_month, model="rw2",
        hyper = list(
          prec = list(
            prior = "pc.prec",
            param = c(0.5, 0.01)
          )
        )
      ) +
      f(temp_index, model="rw2",
        hyper = list(
          prec = list(
            prior = "pc.prec",
            param = c(0.5, 0.01)
          )
        )
      ) +
      f(area, model="bym2", graph = ken.adj, scale.model = TRUE, constr = TRUE,
        hyper = list(
          prec = list(prior = "pc.prec", param = c(1, 0.01)),
          phi  = list(prior = "pc", param = c(0.5, 2/3))
        )),
    data = data_copy,
    family = "binomial",
    weights = weights,
    control.fixed = list(mean = 0, prec = 1),
    control.predictor = list(compute = TRUE, link = 1),
    control.compute = list(dic = FALSE, waic = FALSE)
  )
  
  # Evaluate model on test fold (fold_time == 2)
  predicted_probs <- M_fold2$summary.fitted.values$mean[test_indices]
  actual_vals <- merged_data$outbreak_2[test_indices]
  predicted_classes <- ifelse(predicted_probs > 0.2, 1, 0)

  # Store predictions
  predictions_fold2 <- data.frame(
    row_id = test_indices,
    actual = actual_vals,
    predicted_prob = predicted_probs,
    predicted_class = predicted_classes,
    district = merged_data$district_country.x[test_indices]
  )
  
  # Create confusion matrix
  confusion_matrix_fold2 <- table(
    Actual = actual_vals,
    Predicted = predicted_classes
  )
  
  # Output

  print(confusion_matrix_fold2)
  
  # Precision, Recall, F1
  tp <- sum(predicted_classes == 1 & actual_vals == 1, na.rm = TRUE)
  fp <- sum(predicted_classes == 1 & actual_vals == 0, na.rm = TRUE)
  fn <- sum(predicted_classes == 0 & actual_vals == 1, na.rm = TRUE)
  
  tn <- sum(predicted_classes == 0 & actual_vals == 0, na.rm = TRUE)
  specificity <- tn / (tn + fp)
  recall <- tp / (tp + fn)
  correlation_coefficient <- cor(actual_vals, predicted_probs)

  
  # ROC AUC
  roc_obj <- roc(actual_vals, predicted_probs)
  auc_val <- auc(roc_obj)
  optimal_index <- which.max(roc_obj$sensitivities + roc_full$specificities - 1)
  optimal_threshold <- roc_full$obj[optimal_index]
  
  # PR AUC
  pr_obj <- pr.curve(scores.class0 = predicted_probs[actual_vals == 1],
                     scores.class1 = predicted_probs[actual_vals == 0],
                     curve = FALSE)
  pr_auc <- pr_obj$auc.integral
  
  # -----------------------
  # Brier Score
  # -----------------------
  brier_score <- mean((predicted_probs - actual_vals)^2)
  brier_score
  # Output: a single numeric value
  
  # -----------------------
  # LogLoss
  # -----------------------
  # Prevent log(0)
  eps <- 1e-15
  y_pred_clipped <- pmin(pmax(predicted_probs, eps), 1 - eps)  
  
  logloss <- -mean(actual_vals * log(y_pred_clipped) + (1 - actual_vals) * log(1 - y_pred_clipped))
  logloss
  
  # Output
  cat("\nConfusion Matrix for Fold 2:\n")
  print(confusion_matrix_fold2)
  cat("\nRecall (Sensitivity):", round(recall, 3))
  cat("\nROC AUC:", round(auc_val, 3))
  cat("\brier:", round (brier_score, 3), "\n")
  cat("logloss:", round (logloss, 3), "\n")
  
  
  
  all_roc_objs[[paste0("Fold", k)]] <- roc_obj
  metrics_list <- rbind(metrics_list,
                        data.frame(Fold = k,
                                   AUC = auc_val,
                                   specificity = specificity,
                                   Recall = recall,
                                   Brier = brier_score,
                                   LogLoss = logloss))
  
  
  alert_metrics <- function(  M_fold2, merged_data, ks = c(0.01, 0.05, 0.10)) {
    risk <-   M_fold2$summary.linear.predictor$mean
    y <- data_copy$outbreak_2
    
    prev <- mean(y == 1, na.rm = TRUE)
    n_pos <- sum(y == 1, na.rm = TRUE)
    
    dplyr::bind_rows(lapply(ks, function(k) {
      thr <- as.numeric(stats::quantile(risk, 1 - k, na.rm = TRUE))
      alert <- risk >= thr
      
      tp <- sum(alert & y == 1, na.rm = TRUE)
      n_alert <- sum(alert, na.rm = TRUE)
      
      precision <- tp / n_alert
      recall <- tp / n_pos
      
      data.frame(
        top_frac = k,
        threshold = thr,
        prevalence = prev,
        precision = precision,
        recall = recall,
        lift = precision / prev
      )
    }))
  }
  
  met <- alert_metrics(M_fold2, data_copy)
  utils::write.csv(met, "alert_tradeoff_metricstestfold2.csv", row.names = FALSE)
  
  
  #################FOLD 3 ##################################### 
  
  k <- 3
  
  data_copy <- merged_data
  test_indices <- which(data_copy$fold_space == k)
  data_copy$outbreak_2[test_indices] <- NA 
  
  M_fold2 <- inla(
    outbreak_2 ~ aod_scale +
      wind_strength_scale +
      humidity_scale +
      wind_dir +                 
      aod_scale:wind_dir +        
      f(year_month, model="rw2",
        hyper = list(
          prec = list(
            prior = "pc.prec",
            param = c(0.5, 0.01)
          )
        )
      ) +
      f(temp_index, model="rw2",
        hyper = list(
          prec = list(
            prior = "pc.prec",
            param = c(0.5, 0.01)
          )
        )
      ) +
      f(area, model="bym2", graph = ken.adj, scale.model = TRUE, constr = TRUE,
        hyper = list(
          prec = list(prior = "pc.prec", param = c(1, 0.01)),
          phi  = list(prior = "pc", param = c(0.5, 2/3))
        )),
    data = data_copy,
    family = "binomial",
    weights = weights,
    control.fixed = list(mean = 0, prec = 1),
    control.predictor = list(compute = TRUE, link = 1),
    control.compute = list(dic = FALSE, waic = FALSE)
  )
  
  # Evaluate model on test fold (fold_time == 2)
  predicted_probs <- M_fold2$summary.fitted.values$mean[test_indices]
  actual_vals <- merged_data$outbreak_2[test_indices]
  predicted_classes <- ifelse(predicted_probs > 0.2, 1, 0)

  # Store predictions
  predictions_fold2 <- data.frame(
    row_id = test_indices,
    actual = actual_vals,
    predicted_prob = predicted_probs,
    predicted_class = predicted_classes,
    district = merged_data$district_country.x[test_indices]
  )
  
  # Create confusion matrix
  confusion_matrix_fold2 <- table(
    Actual = actual_vals,
    Predicted = predicted_classes
  )

  
  # Precision, Recall, F1
  tp <- sum(predicted_classes == 1 & actual_vals == 1, na.rm = TRUE)
  fp <- sum(predicted_classes == 1 & actual_vals == 0, na.rm = TRUE)
  fn <- sum(predicted_classes == 0 & actual_vals == 1, na.rm = TRUE)
  
  tn <- sum(predicted_classes == 0 & actual_vals == 0, na.rm = TRUE)
  specificity <- tn / (tn + fp)
  recall <- tp / (tp + fn)
  correlation_coefficient <- cor(actual_vals, predicted_probs)

  
  # ROC AUC
  roc_obj <- roc(actual_vals, predicted_probs)
  auc_val <- auc(roc_obj)
  
  # PR AUC
  pr_obj <- pr.curve(scores.class0 = predicted_probs[actual_vals == 1],
                     scores.class1 = predicted_probs[actual_vals == 0],
                     curve = FALSE)
  pr_auc <- pr_obj$auc.integral
  
  
  # -----------------------
  # Brier Score
  # -----------------------
  brier_score <- mean((predicted_probs - actual_vals)^2)
  brier_score

  # -----------------------
  # LogLoss 
  # -----------------------

  eps <- 1e-15
  y_pred_clipped <- pmin(pmax(predicted_probs, eps), 1 - eps)  
  
  logloss <- -mean(actual_vals * log(y_pred_clipped) + (1 - actual_vals) * log(1 - y_pred_clipped))
  logloss
  
  # Output
  cat("\nConfusion Matrix for Fold 2:\n")
  print(confusion_matrix_fold2)
  cat("\nRecall (Sensitivity):", round(recall, 3))
  cat("\nROC AUC:", round(auc_val, 3))
  cat("\brier:", round (brier_score, 3), "\n")
  cat("logloss:", round (logloss, 3), "\n")
  

  
  all_roc_objs[[paste0("Fold", k)]] <- roc_obj
  metrics_list <- rbind(metrics_list,
                        data.frame(Fold = k,
                                   AUC = auc_val,
                                   specificity = specificity,
                                   Recall = recall,
                                   Brier = brier_score,
                                   LogLoss = logloss))
  
  
  alert_metrics <- function(  M_fold2, merged_data, ks = c(0.01, 0.05, 0.10)) {
    risk <-   M_fold2$summary.linear.predictor$mean
    y <- data_copy$outbreak_2
    
    prev <- mean(y == 1, na.rm = TRUE)
    n_pos <- sum(y == 1, na.rm = TRUE)
    
    dplyr::bind_rows(lapply(ks, function(k) {
      thr <- as.numeric(stats::quantile(risk, 1 - k, na.rm = TRUE))
      alert <- risk >= thr
      
      tp <- sum(alert & y == 1, na.rm = TRUE)
      n_alert <- sum(alert, na.rm = TRUE)
      
      precision <- tp / n_alert
      recall <- tp / n_pos
      
      data.frame(
        top_frac = k,
        threshold = thr,
        prevalence = prev,
        precision = precision,
        recall = recall,
        lift = precision / prev
      )
    }))
  }
  
  met <- alert_metrics(M_fold2, data_copy)
  utils::write.csv(met, "alert_tradeoff_metricstestfold3.csv", row.names = FALSE)
  
  
  #################FOLD 4#####################################
  k <- 4
  
  data_copy <- merged_data
  test_indices <- which(data_copy$fold_space == k)
  data_copy$outbreak_2[test_indices]

  M_fold2 <- inla(
    outbreak_2 ~ aod_scale +
      wind_strength_scale +
      humidity_scale +
      wind_dir +                  
      aod_scale:wind_dir +        
      f(year_month, model="rw2",
        hyper = list(
          prec = list(
            prior = "pc.prec",
            param = c(0.5, 0.01)
          )
        )
      ) +
      f(temp_index, model="rw2",
        hyper = list(
          prec = list(
            prior = "pc.prec",
            param = c(0.5, 0.01)
          )
        )
      ) +
      f(area, model="bym2", graph = ken.adj, scale.model = TRUE, constr = TRUE,
        hyper = list(
          prec = list(prior = "pc.prec", param = c(1, 0.01)),
          phi  = list(prior = "pc", param = c(0.5, 2/3))
        )),
    data = data_copy,
    family = "binomial",
    weights = weights,
    control.fixed = list(mean = 0, prec = 1),
    control.predictor = list(compute = TRUE, link = 1),
    control.compute = list(dic = FALSE, waic = FALSE)
  )
  
  # Evaluate model on test fold (fold_time == 2)
  predicted_probs <- M_fold2$summary.fitted.values$mean[test_indices]
  actual_vals <- merged_data$outbreak_2[test_indices]
  predicted_classes <- ifelse(predicted_probs > 0.2, 1, 0)

  # Store predictions
  predictions_fold2 <- data.frame(
    row_id = test_indices,
    actual = actual_vals,
    predicted_prob = predicted_probs,
    predicted_class = predicted_classes,
    district = merged_data$district_country.x[test_indices]
  )
  
  # Create confusion matrix
  confusion_matrix_fold2 <- table(
    Actual = actual_vals,
    Predicted = predicted_classes
  )
  
  # Output
  print(head(predictions_fold2))
  print(confusion_matrix_fold2)
  
  # Precision, Recall, F1
  tp <- sum(predicted_classes == 1 & actual_vals == 1, na.rm = TRUE)
  fp <- sum(predicted_classes == 1 & actual_vals == 0, na.rm = TRUE)
  fn <- sum(predicted_classes == 0 & actual_vals == 1, na.rm = TRUE)
  
  tn <- sum(predicted_classes == 0 & actual_vals == 0, na.rm = TRUE)
  specificity <- tn / (tn + fp)
  recall <- tp / (tp + fn)
  correlation_coefficient <- cor(actual_vals, predicted_probs)
  
  # ROC AUC
  roc_obj <- roc(actual_vals, predicted_probs)
  auc_val <- auc(roc_obj)
  
  # PR AUC
  pr_obj <- pr.curve(scores.class0 = predicted_probs[actual_vals == 1],
                     scores.class1 = predicted_probs[actual_vals == 0],
                     curve = FALSE)
  pr_auc <- pr_obj$auc.integral
  
  # -----------------------
  # Brier Score
  # -----------------------
  brier_score <- mean((predicted_probs - actual_vals)^2)
  brier_score
  # Output: a single numeric value
  
  # -----------------------
  # LogLoss 
  # -----------------------
  # Prevent log(0)
  eps <- 1e-15
  y_pred_clipped <- pmin(pmax(predicted_probs, eps), 1 - eps)  
  
  logloss <- -mean(actual_vals * log(y_pred_clipped) + (1 - actual_vals) * log(1 - y_pred_clipped))
  logloss
  
  # Output
 
  print(confusion_matrix_fold2)
  cat("\nRecall (Sensitivity):", round(recall, 3))
  cat("\nROC AUC:", round(auc_val, 3))
  cat("\brier:", round (brier_score, 3), "\n")
  cat("logloss:", round (logloss, 3), "\n")
  

  
  all_roc_objs[[paste0("Fold", k)]] <- roc_obj
  metrics_list <- rbind(metrics_list,
                        data.frame(Fold = k,
                                   AUC = auc_val,
                                   specificity = specificity,
                                   Recall = recall,
                                   Brier = brier_score,
                                   LogLoss = logloss))
  
  
  alert_metrics <- function(  M_fold2, merged_data, ks = c(0.01, 0.05, 0.10)) {
    risk <-   M_fold2$summary.linear.predictor$mean
    y <- data_copy$outbreak_2
    
    prev <- mean(y == 1, na.rm = TRUE)
    n_pos <- sum(y == 1, na.rm = TRUE)
    
    dplyr::bind_rows(lapply(ks, function(k) {
      thr <- as.numeric(stats::quantile(risk, 1 - k, na.rm = TRUE))
      alert <- risk >= thr
      
      tp <- sum(alert & y == 1, na.rm = TRUE)
      n_alert <- sum(alert, na.rm = TRUE)
      
      precision <- tp / n_alert
      recall <- tp / n_pos
      
      data.frame(
        top_frac = k,
        threshold = thr,
        prevalence = prev,
        precision = precision,
        recall = recall,
        lift = precision / prev
      )
    }))
  }
  
  met <- alert_metrics(M_fold2, data_copy)
  utils::write.csv(met, "alert_tradeoff_metricstestfold4.csv", row.names = FALSE)
  
  
  #################FOLD 5#####################################

  
  k <- 5
  
  # Copy data and mask outcomes only for fold_time == 2
  data_copy <- merged_data
  test_indices <- which(data_copy$fold_space == k)
  data_copy$outbreak_2[test_indices] <- NA  
  
 
  M_fold2 <- inla(
    outbreak_2 ~ aod_scale +
      wind_strength_scale +
      humidity_scale +
      wind_dir +                    # Main effect
      aod_scale:wind_dir +          # Interaction: AOD effect varies by wind direction
      f(year_month, model="rw2",
        hyper = list(
          prec = list(
            prior = "pc.prec",
            param = c(0.5, 0.01)
          )
        )
      ) +
      f(temp_index, model="rw2",
        hyper = list(
          prec = list(
            prior = "pc.prec",
            param = c(0.5, 0.01)
          )
        )
      ) +
      f(area, model="bym2", graph = ken.adj, scale.model = TRUE, constr = TRUE,
        hyper = list(
          prec = list(prior = "pc.prec", param = c(1, 0.01)),
          phi  = list(prior = "pc", param = c(0.5, 2/3))
        )),
    data = data_copy,
    family = "binomial",
    weights = weights,
    control.fixed = list(mean = 0, prec = 1),
    control.predictor = list(compute = TRUE, link = 1),
    control.compute = list(dic = FALSE, waic = FALSE)
  )
  
  
  predicted_probs <- M_fold2$summary.fitted.values$mean[test_indices]
  actual_vals <- merged_data$outbreak_2[test_indices]
  predicted_classes <- ifelse(predicted_probs > 0.2, 1, 0)

  # Store predictions
  predictions_fold2 <- data.frame(
    row_id = test_indices,
    actual = actual_vals,
    predicted_prob = predicted_probs,
    predicted_class = predicted_classes,
    district = merged_data$district_country.x[test_indices]
  )
  
  # Create confusion matrix
  confusion_matrix_fold2 <- table(
    Actual = actual_vals,
    Predicted = predicted_classes
  )
  
  # Output
  print(confusion_matrix_fold2)
  
  # Precision, Recall, F1
  tp <- sum(predicted_classes == 1 & actual_vals == 1, na.rm = TRUE)
  fp <- sum(predicted_classes == 1 & actual_vals == 0, na.rm = TRUE)
  fn <- sum(predicted_classes == 0 & actual_vals == 1, na.rm = TRUE)
  
  tn <- sum(predicted_classes == 0 & actual_vals == 0, na.rm = TRUE)
  specificity <- tn / (tn + fp)
  recall <- tp / (tp + fn)
  correlation_coefficient <- cor(actual_vals, predicted_probs)
  
  # ROC AUC
  roc_obj <- roc(actual_vals, predicted_probs)
  auc_val <- auc(roc_obj)
  
  # PR AUC
  pr_obj <- pr.curve(scores.class0 = predicted_probs[actual_vals == 1],
                     scores.class1 = predicted_probs[actual_vals == 0],
                     curve = FALSE)
  pr_auc <- pr_obj$auc.integral
  
  # -----------------------
  # Brier Score
  # -----------------------
  brier_score <- mean((predicted_probs - actual_vals)^2)
  brier_score
 
  # -----------------------
  # LogLoss 
  # -----------------------
  # Prevent log(0)
  eps <- 1e-15
  y_pred_clipped <- pmin(pmax(predicted_probs, eps), 1 - eps)  
  
  logloss <- -mean(actual_vals * log(y_pred_clipped) + (1 - actual_vals) * log(1 - y_pred_clipped))
  logloss
  
  # Output
  print(confusion_matrix_fold2)
  cat("\nRecall (Sensitivity):", round(recall, 3))
  cat("\nROC AUC:", round(auc_val, 3))
  cat("\brier:", round (brier_score, 3), "\n")
  cat("logloss:", round (logloss, 3), "\n")
  
  alert_metrics <- function(  M_fold2, merged_data, ks = c(0.01, 0.05, 0.10)) {
    risk <-   M_fold2$summary.linear.predictor$mean
    y <- data_copy$outbreak_2
    
    prev <- mean(y == 1, na.rm = TRUE)
    n_pos <- sum(y == 1, na.rm = TRUE)
    
    dplyr::bind_rows(lapply(ks, function(k) {
      thr <- as.numeric(stats::quantile(risk, 1 - k, na.rm = TRUE))
      alert <- risk >= thr
      
      tp <- sum(alert & y == 1, na.rm = TRUE)
      n_alert <- sum(alert, na.rm = TRUE)
      
      precision <- tp / n_alert
      recall <- tp / n_pos
      
      data.frame(
        top_frac = k,
        threshold = thr,
        prevalence = prev,
        precision = precision,
        recall = recall,
        lift = precision / prev
      )
    }))
  }
  
  met <- alert_metrics(M_fold2, data_copy)
  utils::write.csv(met, "alert_tradeoff_metricstestfold5.csv", row.names = FALSE)
  
  
  all_roc_objs[[paste0("Fold", k)]] <- roc_obj
  metrics_list <- rbind(metrics_list,
                        data.frame(Fold = k,
                                   AUC = auc_val,
                                   specificity = specificity,
                                   Recall = recall,
                                   Brier = brier_score,
                                   LogLoss = logloss))
  
  
  # Combined ROC plot
  test<-plot(all_roc_objs[[1]], col = 1, lwd = 1.5, main = "ROC Curves - All Folds")
  for (i in 2:length(all_roc_objs)) plot(all_roc_objs[[i]], col = i, add = TRUE, lwd = 1.5)
  mean_auc <- mean(metrics_list$AUC)
  legend("bottomright", legend = paste0("Fold ", 1:5, " (AUC=", round(metrics_list$AUC, 3), ")"),
         col = 1:5, lwd = 2, cex = 0.8)
  abline(a = 0, b = 1, lty = 2, col = "gray")
  title(sub = paste0("Mean AUC = ", round(mean_auc, 3)))
  
  save_plot( test, " test_ROC.png", 
            base_width = 14, height = 8)
  print( test)
  

  # Summary statistics
  summary_stats <- data.frame(
    Metric = c("AUC", "specificity", "Recall", "F1", "Brier", "LogLoss"),
    Mean = sapply(metrics_list[, -1], mean, na.rm = TRUE),
    SD = sapply(metrics_list[, -1], sd, na.rm = TRUE)
  )
  
  print(summary_stats)
  write.csv(summary_stats, "cv_summary_metrics.csv", row.names = FALSE)
  
}
  
  
  