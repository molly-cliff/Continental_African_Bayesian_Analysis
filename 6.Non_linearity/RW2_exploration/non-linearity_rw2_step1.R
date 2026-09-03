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
  # -----------------------------
  # 1) Load and filter data
  # -----------------------------
  countries <- c(
    "Benin", "Burkina Faso", "Burundi", "Cameroon", "Central African Republic",
    "Chad", "Côte d'Ivoire", "Eritrea", "Ethiopia", "Gambia", "Ghana", "Guinea",
    "Guinea-Bissau", "Kenya", "Mali", "Mauritania", "Niger", "Nigeria", "Rwanda",
    "Senegal", "South Sudan", "Sudan", "Tanzania", "Togo", "Uganda"
  )
  
  vars_keep <- c(
    "year.x", "month", "outbreak_2", "name_2", "district_country.x", "country",
    "code.x", "rainfall", "eastward_wind", "north_wind", "aod", "humidity",
    "vaccine", "windspeed", "temp", "cropland", "forest", "barren", "pop_density"
  )
  
  dat <- read.csv(data_path) %>%
    mutate(
      year_month = (year.x - min(year.x, na.rm = TRUE)) * 12 + month
    ) %>%
    select(all_of(c(vars_keep, "year_month"))) %>%
    filter(country %in% countries)
  
  shp <- st_read(shape_path, quiet = TRUE) %>%
    filter(COUNTRY %in% countries) %>%
    st_make_valid() %>%
    mutate(district_country.x = paste(NAME_2, COUNTRY, sep = " "))
  
  # -----------------------------
  # 2) Align shapefile and data
  # -----------------------------
  missing_in_shape <- setdiff(unique(dat$district_country.x), unique(shp$district_country.x))
  if (length(missing_in_shape) > 0) {
    stop(
      "These districts are in data but not shapefile: ",
      paste(missing_in_shape, collapse = ", ")
    )
  }
  
  shp <- shp %>% filter(district_country.x %in% dat$district_country.x)
  
  if (!identical(sort(unique(dat$district_country.x)), sort(unique(shp$district_country.x)))) {
    stop("District-country lists differ after filtering.")
  }
  
  shp <- shp %>% arrange(district_country.x)
  dat <- dat %>%
    arrange(district_country.x) %>%
    mutate(area = match(district_country.x, shp$district_country.x))
  
  # -----------------------------
  # 3) Build adjacency for BYM2
  # -----------------------------
  nb <- poly2nb(as(shp, "Spatial"), queen = TRUE)
  adj_file <- file.path(getwd(), "map_test_non_lin.adj")
  nb2INLA(adj_file, nb)
  
  # -----------------------------
  # 4) Weights for class imbalance
  # -----------------------------
  dat$outbreak_2 <- as.numeric(dat$outbreak_2)
  n_outbreak <- sum(dat$outbreak_2 == 1, na.rm = TRUE)
  n_total <- nrow(dat)
  
  w_non_outbreak <- n_outbreak / n_total
  w_outbreak <- 1 - w_non_outbreak
  w <- ifelse(dat$outbreak_2 == 1, w_outbreak, w_non_outbreak)
  
  # -----------------------------
  # 5) Scale candidate variables
  # -----------------------------
  test_vars_raw <- c("temp", "aod", "eastward_wind", "humidity")
  test_vars <- paste0(test_vars_raw, "_scale")
  
  for (v in test_vars_raw) {
    dat[[paste0(v, "_scale")]] <- as.numeric(scale(dat[[v]]))
  }
  
  # -----------------------------
  # 6) Fit base linear model
  # -----------------------------
  base_formula <- as.formula(paste(
    "outbreak_2 ~", paste(test_vars, collapse = " + "),
    "+ f(year_month, model='rw2')",
    "+ f(area, model='bym2', graph=adj_file, scale.model=TRUE, constr=TRUE)"
  ))
  
  M_base <- inla(
    base_formula,
    data = dat,
    family = "binomial",
    weights = w,
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)
  )
  
  base_cpo <- -mean(log(M_base$cpo$cpo), na.rm = TRUE)
  
  # -----------------------------
  # 7) Test each variable as RW2
  # -----------------------------
  results <- lapply(test_vars, function(v) {
    v_disc <- paste0(v, "_disc")
    dat[[v_disc]] <- as.numeric(cut(dat[[v]], breaks = 50))
    
    other_vars <- setdiff(test_vars, v)
    
    f_nl <- as.formula(paste(
      "outbreak_2 ~", paste(other_vars, collapse = " + "),
      "+ f(", v_disc, ", model='rw2')",
      "+ f(year_month, model='rw2')",
      "+ f(area, model='bym2', graph=adj_file, scale.model=TRUE, constr=TRUE)"
    ))
    
    M_nl <- inla(
      f_nl,
      data = dat,
      family = "binomial",
      weights = w,
      control.predictor = list(compute = TRUE),
      control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)
    )
    
    dic_diff  <- M_base$dic$dic - M_nl$dic$dic
    waic_diff <- M_base$waic$waic - M_nl$waic$waic
    cpo_diff  <- base_cpo - (-mean(log(M_nl$cpo$cpo), na.rm = TRUE))
    
    rec <- if (dic_diff > 5 && waic_diff > 5) {
      "YES - SUBSTANTIAL"
    } else if (dic_diff > 2 && waic_diff > 2) {
      "Maybe - Moderate"
    } else {
      "No - Minimal"
    }
    
    data.frame(
      variable         = v,
      dic_linear       = M_base$dic$dic,
      dic_nonlinear    = M_nl$dic$dic,
      dic_improvement  = dic_diff,
      waic_linear      = M_base$waic$waic,
      waic_nonlinear   = M_nl$waic$waic,
      waic_improvement = waic_diff,
      cpo_linear       = base_cpo,
      cpo_nonlinear    = -mean(log(M_nl$cpo$cpo), na.rm = TRUE),
      cpo_improvement  = cpo_diff,
      recommendation   = rec
    )
  })
  
  comparison_results <- bind_rows(results)
  
  # -----------------------------
  # 8) Save + return
  # -----------------------------
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  out_csv <- file.path(output_dir, "nonlinearity_comparison.csv")
  write.csv(comparison_results, out_csv, row.names = FALSE)
  
  print(comparison_results)
  
}