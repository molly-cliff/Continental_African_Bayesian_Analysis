
code_cleaning <- function() {
  
  # ---- Libraries ----
  library(dplyr)
  library(INLA)
  library(sf)
  library(spdep)
  library(sp)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(stringr)
  library(lubridate)
  
  # ---- Output folder ----
  out_dir <- "results"
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # ---- Helpers ----
  invlogit <- function(x) 1 / (1 + exp(-x))
  
  get_intercept <- function(model) {
    if (!is.null(model$summary.fixed) && "(Intercept)" %in% rownames(model$summary.fixed)) {
      return(model$summary.fixed["(Intercept)", "mean"])
    }
    0
  }
  
  parse_time_id <- function(id) {
    id_chr <- as.character(id)
    
    if (grepl("^\\d{4}[-_\\s]\\d{1,2}$", id_chr, perl = TRUE)) {
      parts <- strsplit(id_chr, "[-_\\s]")[[1]]
      return(lubridate::ymd(paste0(parts[1], "-", sprintf("%02d", as.integer(parts[2])), "-01")))
    }
    if (grepl("^\\d{6}$", id_chr)) {
      return(lubridate::ymd(paste0(substr(id_chr, 1, 4), "-", substr(id_chr, 5, 6), "-01")))
    }
    if (grepl("^\\d+$", id_chr)) return(as.Date(NA))
    
    try_date <- suppressWarnings(lubridate::ymd(id_chr))
    if (!is.na(try_date)) return(try_date)
    
    as.Date(NA)
  }
  
  # ---- Load data ----
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
  
  
  # -----------------------------
  # Temperature non-linear plot 
  # -----------------------------
  
  # 1) Get scaling parameters (so we can convert scaled values back to original units)
  temp_center <- attr(merged_data$temp_scale, "scaled:center")
  temp_scale_val <- attr(merged_data$temp_scale, "scaled:scale")
  
  # 2) Recreate the same binning used for temp_index
  temp_breaks <- pretty(range(merged_data$temp_scale, na.rm = TRUE), n = 50)
  merged_data$temp_index <- as.integer(
    cut(merged_data$temp_scale, temp_breaks, include.lowest = TRUE)
  )
  
  # 3) Compute bin midpoints (scaled) then convert midpoints to original temperature units
  temp_mid_original <- ((temp_breaks[-1] + temp_breaks[-length(temp_breaks)]) / 2) *
    temp_scale_val + temp_center
  
  # 4) Only plot if INLA produced a random-effect summary for temp_index
  if (!is.null(model$summary.random$temp_index)) {
    
    # 5) Pull out the random-effect summary table for temp_index
    temp_eff <- as.data.frame(model$summary.random$temp_index)
    
    # 6) Convert the random-effect ID labels into integer bin indices (with a fallback)
    temp_eff$idx <- suppressWarnings(as.integer(as.character(temp_eff$ID)))
    if (anyNA(temp_eff$idx)) {
      temp_eff$idx <- suppressWarnings(as.integer(gsub("\\D+", "", temp_eff$ID)))
    }
    
    # 7) Safety checks to prevent misaligned indexing
    if (anyNA(temp_eff$idx)) stop("Could not parse temp_index IDs into integers.")
    stopifnot(all(temp_eff$idx >= 1 & temp_eff$idx <= length(temp_mid_original)))
    
    # 8) Attach the original-unit temperature value for each bin (x-axis for plotting)
   temp_eff$temp_original <- temp_mid_original[temp_eff$idx] - 273.15
    
    p_temp <- ggplot(temp_eff, aes(temp_original, mean)) +
      geom_line(linewidth = 1) +
      geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), alpha = 0.25) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(
        x = "Temperature (°C)",  # change to °C if appropriate
        y = "Effect on log-odds of outbreak"
      ) +
      theme_minimal()
    

    ggsave(
      file.path(out_dir, "temp_nonlin_effect_test.png"),
      p_temp,
      width = 8, height = 5, units = "in", dpi = 300, bg = "white"
    )
  }
  
  
  
  
  
  # -----------------------------
  # Temporal random effect: year_month (simplified + labeled)
  # -----------------------------
  
  # 0) Default output (in case we skip plotting)
  year_eff <- NULL
  
  # 1) Only run if the model actually has a random-effect summary for year_month
  if (!is.null(model$summary.random$year_month)) {
    
    # 2) Extract year_month random-effect summary into a data.frame
    time_eff <- as.data.frame(model$summary.random$year_month)
    
    # 3) Add 95% interval columns (for plotting ribbon)
    time_eff$lower <- time_eff$`0.025quant`
    time_eff$upper <- time_eff$`0.975quant`
    
    # 4) Convert the INLA ID values to actual dates using your helper
    time_eff$time_date <- vapply(time_eff$ID, parse_time_id, FUN.VALUE = as.Date(NA))
    
   
    
    # 6) Plot mean temporal effect over time + 95% interval ribbon
    p_time <- ggplot(time_eff, aes(x = time_date, y = mean, group = 1)) +
      geom_line(linewidth = 1) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(
        x = "Time",
        y = "Temporal effect (log-odds)"
      ) +
      theme_minimal()
    
    # 7) Save the plot to disk
    ggsave(
      filename = file.path(out_dir, "temporal_effect_time_series_test.png"),
      plot = p_time,
      width = 10, height = 5, units = "in", dpi = 300, bg = "white"
    )
    
    # 8) Create a clean table (year/month columns) for later seasonality analysis
    year_eff <- time_eff |>
      dplyr::mutate(
        year = lubridate::year(time_date),
        month = lubridate::month(time_date)
      ) |>
      dplyr::filter(!is.na(year), !is.na(month)) |>
      dplyr::select(ID, mean, lower, upper, time_date, year, month)
    
  }
  
  merged_data$pred_prob <- model$summary.fitted.values[, 1]
  
  if ("month" %in% names(merged_data)) {
    seasonal_df <- merged_data %>%
      group_by(month) %>%
      summarise(mean_prob = mean(pred_prob, na.rm = TRUE), .groups = "drop")
    
    p_season <- ggplot(seasonal_df, aes(x = month, y = mean_prob, group = 1)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      scale_x_continuous(breaks = 1:12) +
      labs(x = "Month", y = "Mean outbreak probability") +
      theme_minimal()
    
    ggsave(
      filename = file.path(out_dir, "seasonal_cycle_mean_prob_test.png"),
      plot = p_season,
      width = 10, height = 5, units = "in", dpi = 300, bg = "white"
      )}
  
  
  # Exceedance threshold on the probability scale.
  # We will compute Pr(p > tau) for each observation, then average by year-month.
  tau <- 0.2
  
  # Compute exceedance probability per observation.
  # Prefer INLA marginals if available (most accurate), otherwise fall back to a
  # normal approximation using the linear predictor summary.
  ex_probs <- if (!is.null(model$marginals.fitted.values)) {
    
    # model$marginals.fitted.values is a list of marginal distributions (one per obs).
    # inla.pmarginal(tau, mg) returns Pr(X <= tau); so 1 - that is Pr(X > tau).
    vapply(
      model$marginals.fitted.values,
      \(mg) 1 - INLA::inla.pmarginal(tau, mg),
      numeric(1)
    )
    
  } else {
    
    # summary.linear.predictor contains mean/sd for the linear predictor (eta) per obs.
    # Convert tau from probability to logit scale, then compute Pr(eta > logit(tau))
    # assuming eta ~ Normal(mean, sd).
    eta <- model$summary.linear.predictor
    1 - pnorm(qlogis(tau), mean = eta$mean, sd = eta$sd)
  }
  
  # Attach exceedance probabilities to the data, then aggregate to monthly means.
  ex_df <- merged_data %>%
    mutate(exceed = ex_probs, year = year.x) %>%    # year.x comes from a prior join; month already exists
    group_by(year, month) %>%
    summarise(pr_exceed = mean(exceed, na.rm = TRUE), .groups = "drop")
  
  # Plot a year (rows) by month (columns) heatmap of the aggregated exceedance probabilities.
  p_exceed_heat <- ggplot(ex_df, aes(x = factor(month), y = factor(year), fill = pr_exceed)) +
    geom_tile() +
    scale_fill_gradient(limits = c(0, 1), low = "white", high = "navy") +
    labs(
      x = "Month",
      y = "Year",
      fill = paste0("Pr(p > ", tau * 100, "%)"),
      title = "Exceedance probability of monthly outbreak risk"
    ) +
    theme_minimal()
  
  # Save the exceedance heatmap to disk.
  ggsave(
    filename = file.path(out_dir, "exceedance_probability_heatmap_test.png"),
    plot = p_exceed_heat,
    width = 10, height = 5, units = "in", dpi = 300, bg = "white"
  )
  
  
  # ------------------------------------------------------------
  # Risk-score summaries + exceedance heatmap (top decile)
  # ------------------------------------------------------------
  
  # Risk score = posterior mean of linear predictor (log-odds).
  merged_data <- merged_data %>%
    dplyr::mutate(
      year = year.x,
      risk_score = model$summary.linear.predictor$mean
    )
  
  # Helper to save plots consistently
  safe_ggsave <- function(filename, plot, width = 10, height = 5) {
    ggplot2::ggsave(
      filename = file.path(out_dir, filename),
      plot = plot,
      width = width, height = height, units = "in", dpi = 300, bg = "white"
    )
  }
  
  # 1) Seasonal cycle: average risk by month
  seasonal_risk <- merged_data %>%
    dplyr::group_by(month) %>%
    dplyr::summarise(mean_risk = mean(risk_score, na.rm = TRUE), .groups = "drop")
  
  p_season_risk <- ggplot2::ggplot(seasonal_risk, ggplot2::aes(month, mean_risk)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_continuous(breaks = 1:12) +
    ggplot2::labs(
      x = "Month", y = "Mean risk score (log-odds)",
      title = "Seasonal cycle (relative risk)"
    ) +
    ggplot2::theme_minimal()
  
  safe_ggsave("seasonal_cycle_mean_risk.png", p_season_risk)
  
  # 2) Monthly time series: average risk by year-month
  monthly_risk <- merged_data %>%
    dplyr::group_by(year, month) %>%
    dplyr::summarise(mean_risk = mean(risk_score, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(year, month) %>%
    dplyr::mutate(ym = as.Date(sprintf("%d-%02d-01", year, month)))
  
  p_monthly_risk <- ggplot2::ggplot(monthly_risk, ggplot2::aes(ym, mean_risk)) +
    ggplot2::geom_line() +
    ggplot2::labs(
      x = "Year-Month", y = "Mean risk score (log-odds)",
      title = "Monthly relative risk (ranking)"
    ) +
    ggplot2::theme_minimal()
  
  safe_ggsave("monthly_mean_risk.png", p_monthly_risk)
  
  # 3) Exceedance heatmap: fraction of districts in the top decile of risk_score
  thr <- as.numeric(stats::quantile(merged_data$risk_score, 0.90, na.rm = TRUE))
  
  ex_df <- merged_data %>%
    dplyr::mutate(exceed = as.numeric(risk_score > thr)) %>%
    dplyr::group_by(year, month) %>%
    dplyr::summarise(pr_exceed = mean(exceed, na.rm = TRUE), .groups = "drop")
  
  p_exceed_heat <- ggplot2::ggplot(ex_df, ggplot2::aes(factor(month), factor(year), fill = pr_exceed)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(limits = c(0, 1), low = "white", high = "navy") +
    ggplot2::labs(
      x = "Month", y = "Year",
      fill = "Frac. in top 10%",
      title = "Exceedance of outbreak risk (risk-score, top decile)"
    ) +
    ggplot2::theme_minimal()
  
  safe_ggsave("exceedance_riskscore_top_decile_heatmap.png", p_exceed_heat, width = 10, height = 5)
  
  # ------------------------------------------------------------
  # AOD × wind_dir interaction plot (approx 95% band; no cov.fixed)
  # ------------------------------------------------------------
  plot_aod_winddir_interaction_no_covfixed <- function(model, aod_range = c(-2.5, 2.5), n = 200) {
    
    fe <- model$summary.fixed
    need <- c("(Intercept)", "aod_scale", "wind_dir", "aod_scale:wind_dir")
    if (!all(need %in% rownames(fe))) {
      stop("Missing required fixed effects: ", paste(setdiff(need, rownames(fe)), collapse = ", "))
    }
    
    # Convenience extractor for mean + 95% CrI endpoints
    get3 <- function(term) {
      c(
        mean = fe[term, "mean"],
        lo   = fe[term, "0.025quant"],
        hi   = fe[term, "0.975quant"]
      )
    }
    
    b0  <- get3("(Intercept)")
    ba  <- get3("aod_scale")
    bw  <- get3("wind_dir")
    bi  <- get3("aod_scale:wind_dir")
    
    aod_seq <- seq(aod_range[1], aod_range[2], length.out = n)
    
    make_df <- function(wd) {
      # lp = (b0 + wd*bw) + (ba + wd*bi) * aod
      int  <- b0 + wd * bw
      slope <- ba + wd * bi
      
      data.frame(
        wind_dir = factor(wd),
        aod_scale = aod_seq,
        lp_mean  = int[1] + slope[1] * aod_seq,
        lp_lower = int[2] + slope[2] * aod_seq,
        lp_upper = int[3] + slope[3] * aod_seq
      )
    }
    
    df <- rbind(make_df(0), make_df(1))
    
    wind_labels <- c('0' = "West to East", '1' = "East to West")
    
    ggplot2::ggplot(df, ggplot2::aes(aod_scale, lp_mean, color = factor(wind_dir), fill = factor(wind_dir))) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lp_lower, ymax = lp_upper), alpha = 0.15, color = NA) +
      ggplot2::labs(
        x = "AOD (scaled, SD units)",
        y = "Linear predictor (log-odds)") +
      ggplot2::scale_color_discrete(labels = wind_labels, name = "Wind direction") +
      ggplot2::scale_fill_discrete(labels = wind_labels, name = "Wind direction") +
      ggplot2::theme_minimal()
  }
  
  p_int <- plot_aod_winddir_interaction_no_covfixed(model)
  safe_ggsave("aod_by_winddir_logodds.png", p_int, width = 8, height = 5)
  }