library(dplyr)
library(INLA)
library(sf)
library(spdep)
library(sp)

non_linear_testing <- function(
    data_path  = "final_popdensity_vaccine.csv",
    shape_path = "Shapefile_improved.shp",
    output_dir = "nonlinearity_plots"
) {
  countries <- c(
    "Benin", "Burkina Faso", "Burundi", "Cameroon", "Central African Republic",
    "Chad", "Côte d'Ivoire", "Eritrea", "Ethiopia", "Gambia", "Ghana", "Guinea",
    "Guinea-Bissau", "Kenya", "Mali", "Mauritania", "Niger", "Nigeria", "Rwanda",
    "Senegal", "South Sudan", "Sudan", "Tanzania", "Togo", "Uganda"
  )
  
  keep_vars <- c(
    "year.x", "month", "outbreak_2", "name_2", "district_country.x", "country",
    "code.x", "rainfall", "eastward_wind", "north_wind", "aod", "humidity",
    "vaccine", "windspeed", "temp", "cropland", "forest", "barren", "pop_density"
  )
  
  # -----------------------------
  # Load + prep data
  # -----------------------------
  dat <- read.csv(data_path) %>%
    mutate(
      year_month = (year.x - min(year.x, na.rm = TRUE)) * 12 + month
    ) %>%
    select(all_of(c(keep_vars, "year_month"))) %>%
    filter(country %in% countries)
  
  shp <- st_read(shape_path, quiet = TRUE) %>%
    filter(COUNTRY %in% countries) %>%
    st_make_valid() %>%
    mutate(district_country.x = paste(NAME_2, COUNTRY, sep = " "))
  
  # Ensure district alignment
  missing_in_shape <- setdiff(unique(dat$district_country.x), unique(shp$district_country.x))
  if (length(missing_in_shape) > 0) {
    stop("Districts in data but not shapefile: ", paste(missing_in_shape, collapse = ", "))
  }
  
  shp <- shp %>% filter(district_country.x %in% dat$district_country.x)
  
  if (!identical(sort(unique(dat$district_country.x)), sort(unique(shp$district_country.x)))) {
    stop("District-country mismatch after filtering.")
  }
  
  shp <- shp %>% arrange(district_country.x)
  dat <- dat %>%
    arrange(district_country.x) %>%
    mutate(area = match(district_country.x, shp$district_country.x))
  
  # -----------------------------
  # Spatial adjacency for BYM2
  # -----------------------------
  nb <- poly2nb(as(shp, "Spatial"), queen = TRUE)
  adj_file <- file.path(getwd(), "map_test_non_lin.adj")
  nb2INLA(adj_file, nb)
  
  # -----------------------------
  # Class weights
  # -----------------------------
  dat$outbreak_2 <- as.numeric(dat$outbreak_2)
  p_outbreak <- mean(dat$outbreak_2 == 1, na.rm = TRUE)
  w <- ifelse(dat$outbreak_2 == 1, 1 - p_outbreak, p_outbreak)
  
  # -----------------------------
  # Scaled vars + fixed temp RW2
  # -----------------------------
  dat <- dat %>%
    mutate(
      aod_scale = as.numeric(scale(aod)),
      eastward_wind_scale = as.numeric(scale(eastward_wind)),
      humidity_scale = as.numeric(scale(humidity)),
      temp_scale = as.numeric(scale(temp)),
      temp_scale2 = as.numeric(cut(temp_scale, breaks = 50))
    )
  
  test_vars <- c("aod_scale", "eastward_wind_scale", "humidity_scale")
  
  # -----------------------------
  # Base model (linear test vars)
  # -----------------------------
  formula_base <- outbreak_2 ~
    aod_scale + eastward_wind_scale + humidity_scale +
    f(temp_scale2, model = "rw2") +
    f(year_month, model = "rw2") +
    f(area, model = "bym2", graph = adj_file, scale.model = TRUE, constr = TRUE)
  
  M_base <- inla(
    formula_base,
    data = dat,
    family = "binomial",
    weights = w,
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)
  )
  
  base_dic  <- M_base$dic$dic
  base_waic <- M_base$waic$waic
  base_cpo  <- -mean(log(M_base$cpo$cpo), na.rm = TRUE)
  
  # -----------------------------
  # Test non-linearity per variable
  # -----------------------------
  results <- lapply(test_vars, function(var) {
    var_disc <- paste0(var, "_disc")
    dat[[var_disc]] <- as.numeric(cut(dat[[var]], breaks = 50))
    other_vars <- setdiff(test_vars, var)
    
    formula_nl <- as.formula(paste(
      "outbreak_2 ~",
      paste(other_vars, collapse = " + "),
      "+ f(", var_disc, ", model = 'rw2')",
      "+ f(temp_scale2, model = 'rw2')",
      "+ f(year_month, model = 'rw2')",
      "+ f(area, model = 'bym2', graph = adj_file, scale.model = TRUE, constr = TRUE)"
    ))
    
    M_nl <- inla(
      formula_nl,
      data = dat,
      family = "binomial",
      weights = w,
      control.predictor = list(compute = TRUE),
      control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)
    )
    
    nl_dic  <- M_nl$dic$dic
    nl_waic <- M_nl$waic$waic
    nl_cpo  <- -mean(log(M_nl$cpo$cpo), na.rm = TRUE)
    
    d_dic  <- base_dic - nl_dic
    d_waic <- base_waic - nl_waic
    d_cpo  <- base_cpo - nl_cpo
    
    rec <- if (d_dic > 5 && d_waic > 5) {
      "YES - SUBSTANTIAL"
    } else if (d_dic > 2 && d_waic > 2) {
      "Maybe - Moderate"
    } else {
      "No - Minimal"
    }
    
    data.frame(
      variable = var,
      dic_linear = base_dic,
      dic_nonlinear = nl_dic,
      dic_improvement = d_dic,
      waic_linear = base_waic,
      waic_nonlinear = nl_waic,
      waic_improvement = d_waic,
      cpo_linear = base_cpo,
      cpo_nonlinear = nl_cpo,
      cpo_improvement = d_cpo,
      recommendation = rec
    )
  })
  
  comparison_results <- bind_rows(results)
  
  # -----------------------------
  # Save + return
  # -----------------------------
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  write.csv(
    comparison_results,
    file.path(output_dir, "nonlinearity_comparison.csv"),
    row.names = FALSE
  )
  
  print(comparison_results)
  invisible(comparison_results)
}