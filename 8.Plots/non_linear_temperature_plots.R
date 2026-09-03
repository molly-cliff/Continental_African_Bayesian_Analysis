temp_debug <- function() {
  
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


  # Convert temperature from Kelvin to Celsius
merged_data$temp_c <- merged_data$temp - 273.15

# Keep rows with non-missing values
plot_data <- merged_data[!is.na(merged_data$temp_c) & !is.na(merged_data$outbreak_2), ]

# Histogram-style bar plot:
# counts of months, split by outbreak status (0 = no outbreak, 1 = outbreak)
library(ggplot2)

outbreak_data <- merged_data[
  merged_data$outbreak_2 == 1 & !is.na(merged_data$temp_c),
]


plot<-ggplot(outbreak_data, aes(x = temp_c)) +
  geom_histogram(
    binwidth = 1,      # 1°C bins
    fill = "#2085f0",
  ) +
  labs(
    x = "Temperature (°C)",
    y = "Total district-months experiencing bacterial meningitis outbreaks"
  ) +
  theme_minimal()

  
  ggsave("histo.png", plot, width = 8, height = 5)
  
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
  
  # Use ONE consistent binning scheme on the ORIGINAL scale
  breaks_k <- pretty(range(temp_k, na.rm = TRUE), n = 50)
  
  # Make factor bins with explicit levels
  merged_data$temp_bin <- cut(
    temp_k,
    breaks = breaks_k,
    include.lowest = TRUE,
    right = TRUE
  )
  
  # Refit model using temp_bin instead of temp_index
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
 
  
  # ----------------------------------
  # Extract INLA random effect
  # ----------------------------------
  temp_eff <- as.data.frame(model$summary.random$temp_bin)
  
  # Factor levels actually used in model
  bin_levels <- levels(merged_data$temp_bin)
  
  # Check dimensions
  cat("Number of temp bins:", length(bin_levels), "\n")
  cat("Rows in random effect:", nrow(temp_eff), "\n")
  
  # ----------------------------------
  # Parse interval labels safely
  # ----------------------------------
  get_midpoint <- function(interval_label) {
    
    # remove brackets/parentheses
    clean <- gsub("\\[|\\]|\\(|\\)", "", interval_label)
    
    # split lower and upper bound
    bounds <- strsplit(clean, ",")[[1]]
    
    mean(as.numeric(bounds))
  }
  
  # Convert all bin labels → midpoint in Kelvin
  bin_mid_k <- sapply(bin_levels, get_midpoint)
  
  # ----------------------------------
  # Match random effect rows directly
  # (NO regex extraction of IDs)
  # ----------------------------------
  
  # Keep only available bins
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
  # Plot
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
      "temp_nonlin_effect_fixed.png"
    ),
    p_temp,
    width = 8,
    height = 5,
    dpi = 300,
    bg = "white"
  )

    }
    
