
library(INLA)
library(ggplot2)
library(mgcv)
library(dplyr)
library(DAAG)
library(reshape2)
library(writexl)


# =====================================================================
# Function: monthly_models_test3
# Description: 
#   - Aggregates outbreak data by month and district
#   - Fits several temporal INLA models (RW1, RW2, AR1)
#   - Compares model performance (DIC, WAIC, Log CPO)
#   - Exports model comparison results to Excel
# =====================================================================

monthly_models_test3 <- function() {
  
  # ----------------------------------------------------------
  # STEP 1: Read and inspect data
  # ----------------------------------------------------------
  spatiotemporaloutbreaks <- read.csv("full-outbreak-matched.csv")
  print(head(spatiotemporaloutbreaks))
  
  # ----------------------------------------------------------
  # STEP 2: Select relevant variables
  # ----------------------------------------------------------
  myvars <- c("district_country", "month", "outbreak", "year", "name_2")
  newdata <- spatiotemporaloutbreaks[, myvars]
  
  # ----------------------------------------------------------
  # STEP 3: Aggregate data by year, month, and district
  #          -> outbreak_occur = 1 if any outbreak > 0
  # ----------------------------------------------------------
  monthly_data <- newdata %>%
    group_by(year, month, district_country) %>%
    summarise(outbreak_occur = as.integer(any(outbreak > 0, na.rm = TRUE)),
              .groups = 'drop')
  
  # ----------------------------------------------------------
  # STEP 4: Create continuous month variable (time index)
  # ----------------------------------------------------------
  monthly_data <- monthly_data %>%
    arrange(year, month) %>%
    mutate(continuous_month = (year - min(year)) * 12 + month)
  
  monthly_data$continuous_month <- as.numeric(monthly_data$continuous_month)
  
  cat("Monthly aggregated data:\n")
  print(head(monthly_data))
  
  
  # =====================================================================
  # SECTION A: TEMPORAL EFFECTS FOR MONTH – CYCLIC MODELS
  # =====================================================================
  # These models capture repeating seasonal effects (months loop cyclically)
  # =====================================================================
  
  # ----------------------------------------------------------
  # MODEL 1A: RW2 (Second-order Random Walk, Cyclic)
  # ----------------------------------------------------------
  formula <- outbreak_occur ~ f(month, model = "rw2", cyclic = TRUE)
  
  rw2_month_cyclic <- inla(
    formula,
    data = monthly_data,
    family = "binomial",
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
  )
  
  dic_rw2_month_cyclic <- rw2_month_cyclic$dic$dic
  waic_rw2_month_cyclic <- rw2_month_cyclic$waic$waic
  cpo_valuesrw2_month_cyclic <- -mean(log(rw2_month_cyclic$cpo$cpo))
  
  print(rw2_month_cyclic$summary.fixed)
  print(rw2_month_cyclic$summary.random$week)
  summary(rw2_month_cyclic)
  
  
  # ----------------------------------------------------------
  # MODEL 1B: RW1 (First-order Random Walk, Cyclic)
  # ----------------------------------------------------------
  formula <- outbreak_occur ~ f(month, model = "rw1", cyclic = TRUE)
  
  rw1_month_cyclic <- inla(
    formula,
    data = monthly_data,
    family = "binomial",
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
  )
  
  dic_rw1_month_cyclic <- rw1_month_cyclic$dic$dic
  waic_rw1_month_cyclic <- rw1_month_cyclic$waic$waic
  cpo_valuesrw1_month_cyclic <- -mean(log(rw1_month_cyclic$cpo$cpo))
  
  print(rw1_month_cyclic$summary.fixed)
  print(rw1_month_cyclic$summary.random$week)
  summary(rw1_month_cyclic)
  
  
  # ----------------------------------------------------------
  # MODEL 1C: AR1 (Autoregressive, Cyclic)
  # ----------------------------------------------------------
  formula <- outbreak_occur ~ f(month, model = "ar1", cyclic = TRUE)
  
  ar1_month_cyclic <- inla(
    formula,
    data = monthly_data,
    family = "binomial",
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
  )
  
  dic_ar1_month_cyclic <- ar1_month_cyclic$dic$dic
  waic_ar1_month_cyclic <- ar1_month_cyclic$waic$waic
  cpo_valuesar1_month_cyclic <- -mean(log(ar1_month_cyclic$cpo$cpo))
  
  print(ar1_month_cyclic$summary.fixed)
  print(ar1_month_cyclic$summary.random$week)
  summary(ar1_month_cyclic)
  
  
  # =====================================================================
  # SECTION B: TEMPORAL EFFECTS FOR CONTINUOUS MONTH
  # =====================================================================
  # These models capture long-term temporal trends (non-cyclic)
  # =====================================================================
  
  # ----------------------------------------------------------
  # MODEL 2A: RW2 on continuous month
  # ----------------------------------------------------------
  formula <- outbreak_occur ~ f(continuous_month, model = "rw2")
  
  rw2_month_cont <- inla(
    formula,
    data = monthly_data,
    family = "binomial",
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
  )
  
  dic_rw2_month_cont <- rw2_month_cont$dic$dic
  waic_rw2_month_cont <- rw2_month_cont$waic$waic
  cpo_valuesrw2_month_cont <- -mean(log(rw2_month_cont$cpo$cpo))
  
  print(rw2_month_cont$summary.fixed)
  print(rw2_month_cont$summary.random$week)
  summary(rw2_month_cont)
  
  
  # ----------------------------------------------------------
  # MODEL 2B: RW1 on continuous month
  # ----------------------------------------------------------
  formula <- outbreak_occur ~ f(continuous_month, model = "rw1")
  
  rw1_month_cont <- inla(
    formula,
    data = monthly_data,
    family = "binomial",
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
  )
  
  dic_rw1_month_cont <- rw1_month_cont$dic$dic
  waic_rw1_month_cont <- rw1_month_cont$waic$waic
  cpo_valuesrw1_month_cont <- -mean(log(rw1_month_cont$cpo$cpo))
  
  
  # ----------------------------------------------------------
  # MODEL 2C: AR1 on continuous month
  # ----------------------------------------------------------
  formula <- outbreak_occur ~ f(continuous_month, model = "ar1")
  
  ar1_month_cont <- inla(
    formula,
    data = monthly_data,
    family = "binomial",
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
  )
  
  dic_ar1_month_cont <- ar1_month_cont$dic$dic
  waic_ar1_month_cont <- ar1_month_cont$waic$waic
  cpo_valuesar1_month_cont <- -mean(log(ar1_month_cont$cpo$cpo))
  
  print(ar1_month_cont$summary.fixed)
  print(ar1_month_cont$summary.random$week)
  summary(ar1_month_cont)
  
  
  # =====================================================================
  # SECTION C: MODEL PERFORMANCE COMPARISON
  # =====================================================================
  # Collect DIC, WAIC, and LogCPO from all models
  # =====================================================================
  dic <- data.frame(
    criteria = c("DIC", "WAIC", "Logcpo"),
    rw2_month_cyclic = c(dic_rw2_month_cyclic, waic_rw2_month_cyclic, cpo_valuesrw2_month_cyclic),
    rw1_month_cyclic = c(dic_rw1_month_cyclic, waic_rw1_month_cyclic, cpo_valuesrw1_month_cyclic),
    ar1_month_cyclic = c(dic_ar1_month_cyclic, waic_ar1_month_cyclic, cpo_valuesar1_month_cyclic),
    ar1_month_cont = c(dic_ar1_month_cont, waic_ar1_month_cont, cpo_valuesar1_month_cont),
    rw2_month_cont = c(dic_rw2_month_cont, waic_rw2_month_cont, cpo_valuesrw2_month_cont),
    rw1_month_cont = c(dic_rw1_month_cont, waic_rw1_month_cont, cpo_valuesrw1_month_cont)
  )
  
  # Save results
  write_xlsx(dic, path = "model_comparison_temporal.xlsx")
}
