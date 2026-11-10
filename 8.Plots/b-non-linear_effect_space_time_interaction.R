
non_linear_plots <- function() {
  
  # ---- Libraries ----
  library(dplyr)
  library(INLA)
  library(sf)
  library(spdep)
  library(sp)
  library(ggplot2)
  library(lubridate)
  
  # ============================================================
  # 1. SET UP OUTPUTS AND TIME-PARSING HELPERS
  # ============================================================

  # Create the directory used for generated plots.
  out_dir <- "results"
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # Convert supported year-month identifiers into Date objects.
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
  
  # ============================================================
  # 2. LOAD AND FILTER OUTBREAK DATA
  # ============================================================

  # Load the joined outbreak and environmental dataset.
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
  
  # ============================================================
  # 3. PREPARE SPATIAL DATA AND ADJACENCY STRUCTURE
  # ============================================================

  # Read the district shapefile and retain the study countries.
  shape2 <- st_read("Shapefile_improved.shp", quiet = TRUE)
  shape2 <- shape2[shape2$COUNTRY %in% keep_countries, ]
  shape2 <- st_make_valid(shape2)
  
  shapefile_spatial <- as(shape2, "Spatial")
  nb <- poly2nb(shapefile_spatial, queen = TRUE)
  
  ken.adj <- file.path(tempdir(), "dust_wind2.adj")
  spdep::nb2INLA(ken.adj, nb)
  stopifnot(file.exists(ken.adj))
  
  shape2$district_country.x <- paste(shape2$NAME_2, shape2$COUNTRY, sep = " ")
  
  # ============================================================
  # 4. VALIDATE AND ALIGN DISTRICTS
  # ============================================================

  # Check that every district in the data has a matching polygon.
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
  
  # ============================================================
  # 5. CALCULATE OUTCOME WEIGHTS
  # ============================================================

  # Calculate class weights to account for outbreak imbalance.
  tab <- as.data.frame(table(merged_data$outbreak_2))
  outbreak <- as.numeric(tab[tab$Var1 == 1, "Freq"])
  no_outbreak <- as.numeric(tab[tab$Var1 == 0, "Freq"])
  total <- outbreak + no_outbreak
  
  non_outbreak_weight <- outbreak / total
  outbreak_weight <- 1 - non_outbreak_weight
  weights <- ifelse(merged_data$outbreak_2 == 1, outbreak_weight, non_outbreak_weight)
  
  # ============================================================
  # 6. SCALE COVARIATES AND DERIVE MODEL PREDICTORS
  # ============================================================

  # Standardise environmental predictors and derive wind variables.

  merged_data$eastward_wind_scale <- (scale(merged_data$eastward_wind))
  merged_data$aod_scale <- scale(merged_data$aod)
  merged_data$humidity_scale <- scale(merged_data$humidity)
  merged_data$temp_scale <- scale(merged_data$temp)
  merged_data$temp_index <- as.numeric(cut(merged_data$temp_scale, breaks = 50))

  
  merged_data$wind_strength <- abs(merged_data$eastward_wind)
  merged_data$wind_dir <- ifelse(merged_data$eastward_wind < 0, 1, 0)
  
  # ============================================================
  # 7. FIT SPATIOTEMPORAL INLA MODEL
  # ============================================================

  # Specify fixed effects, temporal RW2 effect, nonlinear temperature
  # effect, and spatial BYM2 effect.
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
  
  # Fit the weighted binomial model.
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
  
 
  

  # ============================================================
  # 8. EXTRACT AND PLOT THE TEMPORAL RANDOM EFFECT
  # ============================================================
  
  # Map each year-month index to the first day of its calendar month.
  ym_lookup <- merged_data %>%
    dplyr::distinct(year_month, year.x, month) %>%
    dplyr::mutate(
      time_date = as.Date(sprintf("%04d-%02d-01", as.integer(year.x), as.integer(month)))
    ) %>%
    dplyr::arrange(year_month)
  
  year_eff <- NULL
  
  if (!is.null(model$summary.random$year_month)) {
    
    time_eff <- as.data.frame(model$summary.random$year_month)
    
    # Convert INLA random-effect IDs into the year-month index.
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
    
    # Keep a clean table for inspection or later export.
    year_eff <- time_eff %>%
      dplyr::select(year_month, time_date, mean, lower, upper, sd, kld)
    
    # Plot the temporal effect if the date mapping succeeded.
   
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
   
  

  # ============================================================
  # 9. PLOT THE AOD AND WIND-DIRECTION INTERACTION
  # ============================================================

  # Construct an approximate 95% credible band from fixed-effect
  # posterior summaries, without conditioning on other covariates.
  plot_aod_winddir_interaction_no_covfixed <- function(model, aod_range = c(-2.5, 2.5), n = 200) {
    
    fe <- model$summary.fixed
    
    # Extract posterior means and 95% credible interval endpoints.
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
      # Calculate the linear predictor for one wind-direction category.
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

    # ============================================================
    # 10. PREPARE AND PLOT THE SPATIAL BYM2 RANDOM EFFECT
    # ============================================================

    # Extract posterior spatial effects from the fitted BYM2 model.
    area_eff <- as.data.frame(model$summary.random$area)

    # Convert INLA area IDs into the integer indices used by the shapefile.
    area_eff$area_id <- as.integer(gsub("[^0-9]", "", as.character(area_eff$ID)))

    # Recreate the area-index-to-district mapping used by the model.
    shape2_map <- shape2 %>%
      st_drop_geometry() %>%
      select(district_country.x) %>%
      mutate(area_id = row_number())

    # Join posterior means and credible intervals to district identifiers.
    area_eff2 <- area_eff %>%
      select(area_id, mean, `0.025quant`, `0.975quant`) %>%
      rename(lower = `0.025quant`, upper = `0.975quant`) %>%
      left_join(shape2_map, by = "area_id")

    # Attach the spatial effects to the district geometries.
    map_df <- shape2 %>%
      left_join(area_eff2, by = "district_country.x")

    # Plot and save posterior mean spatial effects by district.
    p_spatial_mean <- ggplot(map_df) +
      geom_sf(aes(fill = mean), color = NA) +
      scale_fill_gradient(low = "white", high = "navy", na.value = "grey90") +
      labs(
        title = "Spatial effect (BYM2 area random effect, posterior mean)",
        fill = "log-odds"
      ) +
      theme_minimal()

    ggsave(
      "spatial_bym2_mean_thresholds.png",
      p_spatial_mean,
      width = 9,
      height = 7
    )

    print(p_spatial_mean)
}
