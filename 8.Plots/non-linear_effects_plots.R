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
  out_dir <- "result_test"
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
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
  
  merged_data$temp_index <- as.numeric(cut(merged_data$temp_scale, breaks = 50))
  merged_data$aod_index <- as.numeric(cut(merged_data$aod_scale, breaks = 50))
  
  merged_data$wind_strength <- abs(merged_data$eastward_wind)
  merged_data$wind_strength_scale <- as.numeric(scale(merged_data$wind_strength))
  
  merged_data$wind_dir <- ifelse(merged_data$eastward_wind < 0, 1, 0)
  
  temp_k <- merged_data$temp
  
  # Use one consistent binning scheme on the original temperature scale
  breaks_k <- pretty(range(temp_k, na.rm = TRUE), n = 50)
  
  # Create factor bins with explicit levels
  merged_data$temp_bin <- cut(
    temp_k,
    breaks = breaks_k,
    include.lowest = TRUE,
    right = TRUE
  )
  
  # Refit model using temperature bins
  formula1 <- outbreak_2 ~ aod_scale +
    wind_strength_scale +
    humidity_scale +
    wind_dir +
    aod_scale:wind_dir +
    f(year_month, model = "rw2",
      hyper = list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))) +
    f(temp_bin, model = "rw2",
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
  
  library(dplyr)
  library(ggplot2)
  
  # ----------------------------------
  # Extract INLA random effect
  # ----------------------------------
  temp_eff <- as.data.frame(model$summary.random$temp_bin)
  
  # Factor levels used by the model
  bin_levels <- levels(merged_data$temp_bin)
  
  # Check extracted random-effect dimensions
  cat("Number of temp bins:", length(bin_levels), "\n")
  cat("Rows in random effect:", nrow(temp_eff), "\n")
  
  # ----------------------------------
  # Parse temperature interval labels
  # ----------------------------------
  get_midpoint <- function(interval_label) {
    
    # remove brackets/parentheses
    clean <- gsub("\\[|\\]|\\(|\\)", "", interval_label)
    
    # split lower and upper bound
    bounds <- strsplit(clean, ",")[[1]]
    
    mean(as.numeric(bounds))
  }
  
  # Convert bin labels to midpoints in Kelvin
  bin_mid_k <- sapply(bin_levels, get_midpoint)
  
  # ----------------------------------
  # Match random-effect rows directly
  # No regex extraction of IDs
  # ----------------------------------
  
  # Keep only valid temperature bins
  n_keep <- min(nrow(temp_eff), length(bin_mid_k))
  
  temp_eff <- temp_eff %>%
    mutate(
      temp_k = bin_mid_k[seq_len(n_keep)],
      temp_c = temp_k - 273.15
    ) %>%
    filter(
      !is.na(temp_c),
      !is.na(mean),
      !is.na(`0.025quant`),
      !is.na(`0.975quant`)
    ) %>%
    arrange(temp_c)
  
  # ----------------------------------
  # Diagnostics
  # ----------------------------------
  cat(
    "\nTemperature range plotted:",
    round(min(temp_eff$temp_c), 1),
    "to",
    round(max(temp_eff$temp_c), 1),
    "°C\n"
  )
  
  print(
    tail(
      temp_eff[, c("temp_c", "mean")],
      10
    )
  )
  
  # ----------------------------------
  # Plot temperature effect
  # ----------------------------------
  p_temp <- ggplot(
    temp_eff,
    aes(x = temp_c, y = mean)
  ) +
    geom_ribbon(
      aes(
        ymin = `0.025quant`,
        ymax = `0.975quant`
      ),
      alpha = 0.25
    ) +
    geom_line(linewidth = 1) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed"
    ) +
    scale_x_continuous(
      limits = range(temp_eff$temp_c),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      x = "Temperature (°C)",
      y = "Effect on log-odds"
    ) +
    theme_minimal()
  
  print(p_temp)
  
  ggsave(
    file.path(
      out_dir,
      "temp_nonlin_effect_fixed.png"
    ),
    p_temp,
    width = 8,
    height = 5,
    dpi = 300,
    bg = "white"
  )
  ```
  
  library(dplyr)
  library(ggplot2)
  
  # Extract temperature random effect after model fitting
  
  temp_eff <- as.data.frame(model$summary.random$temp_bin)
  
  # 1. Parse the bin index from the INLA ID
  temp_eff$bin_index <- as.integer(gsub("\\D+", "", as.character(temp_eff$ID)))
  
  # 2. Get the ordered factor levels used by the model
  bin_levels <- levels(merged_data$temp_bin)
  
  # 3. Parse bin midpoints safely
  get_mid <- function(lbl) {
    nums <- as.numeric(unlist(regmatches(lbl, gregexpr("-?[0-9]*\\.?[0-9]+", lbl))))
    mean(nums)
  }
  bin_mid_k <- vapply(bin_levels, get_mid, numeric(1))
  
  # 4. Map bin index to midpoint and remove invalid rows
  temp_eff <- temp_eff %>%
    mutate(
      temp_k = bin_mid_k[bin_index],
      temp_c = temp_k - 273.15
    ) %>%
    filter(!is.na(temp_c), !is.na(mean), !is.na(`0.025quant`), !is.na(`0.975quant`)) %>%
    arrange(temp_c)
  
  # 5. Quick diagnostic checks
  cat("Rows in temp_eff:", nrow(temp_eff), "\n")
  cat("Non-NA temp_c:", sum(!is.na(temp_eff$temp_c)), "\n")
  print(head(temp_eff[, c("ID", "bin_index", "temp_c", "mean")]))
  
  # 6. Plot temperature effect
  p_temp <- ggplot(temp_eff, aes(x = temp_c, y = mean)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), alpha = 0.25) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = "Temperature (°C)", y = "Effect on log-odds") +
    theme_minimal()
  
  print(p_temp)
  ggsave(
    file.path(
      out_dir,
      "temp_nonlin_effect_testingtest.png"
    ),
    p_temp,
    width=8,
    height=5,
    dpi=300,
    bg="white"
  )
  
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
  # -----------------------------
  # 1. Recover scaling information
  # ------------------------------------------------------------
  
  center <- attr(merged_data$temp_scale, "scaled:center")
  
  # This is the ORIGINAL mean temperature used in scale()
  
  scale_val <- attr(merged_data$temp_scale, "scaled:scale")
  
  # ------------------------------------------------------------
  # 2. Create temperature bins, this is for ease of plotting
  # ------------------------------------------------------------
  
  breaks <- pretty(
    range(merged_data$temp_scale, na.rm = TRUE),
    n = 50
  )
  
  # ------------------------------------------------------------
  # 3. Assign each temperature (scaled) to a bin
  # ------------------------------------------------------------
  
  merged_data$temp_index <- as.integer(
    cut(
      merged_data$temp_scale,
      breaks,
      include.lowest = TRUE
    )
  )
  
  # temp index essentially becomes interger of the different bins as opposed to raw scaled value
  
  # ------------------------------------------------------------
  # 4. Calculate midpoint of each bin
  # ------------------------------------------------------------
  
  mid_scaled <-
    (
      head(breaks,-1) +
        tail(breaks,-1)
    )/2
  
  # midpoint of each bin
  # Still in SCALED units, again ease of plotting due to large number of data points
  
  # ------------------------------------------------------------
  # 5. Convert midpoints back
  #    to original temperatures
  # ------------------------------------------------------------
  
  mid_temp <-
    mid_scaled *
    scale_val +
    center
  
  # Undo z-score scaling:
  #
  # original = scaled × sd + mean
  #
  
  # ------------------------------------------------------------
  # 6. Continue if INLA estimated
  # temp effect
  # ------------------------------------------------------------
  
  if (!is.null(model$summary.random$temp_index)) {
    
    temp_eff <-
      as.data.frame(
        model$summary.random$temp_index
      )
    
    # Extract estimated effect
    # for each temperature bin
    
    # ------------------------------------------------------------
    # 7. Extract numeric IDs
    # ------------------------------------------------------------
    
    temp_eff$idx <-
      as.integer(
        gsub("\\D+","",temp_eff$ID)
      )
    
    # ------------------------------------------------------------
    # 9. Match bins to temperature
    # ------------------------------------------------------------
    
    temp_eff$temp_original <-
      mid_temp[
        temp_eff$idx
      ] -
      273.15
    
    # ------------------------------------------------------------
    # 10. Plot effect curve
    # ------------------------------------------------------------
    
    p_temp <-
      ggplot(
        temp_eff,
        aes(
          temp_original,
          mean
        )
      )+
      
      geom_line(
        linewidth=1
      )+
      
      geom_ribbon(
        aes(
          ymin=`0.025quant`,
          ymax=`0.975quant`
        ),
        alpha=.25
      )+
      
      geom_hline(
        yintercept=0,
        linetype="dashed"
      )+
      
      labs(
        x="Temperature (°C)",
        y="Effect on log-odds"
      )+
      
      theme_minimal()
    
    # ------------------------------------------------------------
    # 11. Save figure
    # ------------------------------------------------------------
    
    ggsave(
      file.path(
        out_dir,
        "temp_nonlin_effect.png"
      ),
      p_temp,
      width=8,
      height=5,
      dpi=300,
      bg="white"
    )
    
  }
  
  # -----------------------------
  # Temporal random effect: year_month (robust mapping to real dates)
  # -----------------------------
  
  # Build lookup: year_month index -> first day of month
  ym_lookup <- merged_data %>%
    dplyr::distinct(year_month, year.x, month) %>%
    dplyr::mutate(
      time_date = as.Date(sprintf("%04d-%02d-01", as.integer(year.x), as.integer(month)))
    ) %>%
    dplyr::arrange(year_month)
  
  year_eff <- NULL
  
  if (!is.null(model$summary.random$year_month)) {
    
    time_eff <- as.data.frame(model$summary.random$year_month)
    
    # INLA rw2 IDs are typically numeric indices
    time_eff$year_month <- suppressWarnings(as.integer(as.character(time_eff$ID)))
    
    time_eff <- time_eff %>%
      dplyr::left_join(
        ym_lookup %>% dplyr::select(year_month, time_date),
        by = "year_month"
      ) %>%
      dplyr::mutate(
        lower = `0.025quant`,
        upper = `0.975quant`
      ) %>%
      dplyr::arrange(time_date)
    
    # Keep a clean table if you want to inspect/export later
    year_eff <- time_eff %>%
      dplyr::select(year_month, time_date, mean, lower, upper, sd, kld)
    
    # Safety checks
    if (all(is.na(time_eff$time_date))) {
      warning("All temporal dates are NA after mapping. Check year_month construction and INLA IDs.")
    } else {
      p_time <- ggplot(time_eff, aes(x = time_date, y = mean, group = 1)) +
        geom_line(linewidth = 1) +
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        labs(
          x = "Time (monthly over 20 years)",
          y = "Effect on log odds of outbreak"
        ) +
        theme_minimal()
      
      ggsave(
        filename = file.path(out_dir, "temporal_effect_time_series.png"),
        plot = p_time,
        width = 10, height = 5, units = "in", dpi = 300, bg = "white"
      )
    }
    
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
          x = "AOD (scaled-zscore)",
          y = "Effect on log odds of outbreak") +
        ggplot2::scale_color_discrete(labels = wind_labels, name = "Wind direction") +
        ggplot2::scale_fill_discrete(labels = wind_labels, name = "Wind direction") +
        ggplot2::theme_minimal()
    }
    
    p_int <- plot_aod_winddir_interaction_no_covfixed(model)
    ggsave("aod_by_winddir_logodds.png", p_int, width = 8, height = 5)
  }
}
}
