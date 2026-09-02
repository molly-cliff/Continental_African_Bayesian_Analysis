burnett_func <- function() {
  # =========================================================
  # 0) LOAD LIBRARIES
  # =========================================================
  library(dplyr)
  library(openxlsx)
  library(sf)
  library(sp)
  library(spdep)
  library(INLA)
  library(igraph)
  library(ggplot2)
  library(scales)
  library(viridis)
  library(rnaturalearth)
  
  # =========================================================
  # 1) LOAD AND PREPARE DATA
  # =========================================================
  spatiotemporaloutbreaks <- read.csv("final_popdensity_vaccine.csv", stringsAsFactors = FALSE)
  
  required_cols <- c("year.x", "month", "outbreak_2", "name_2", "district_country.x", "country", "code.x")
  if (!all(required_cols %in% names(spatiotemporaloutbreaks))) {
    stop("CSV is missing some columns!")
  }
  
  spatiotemporaloutbreaks$year_month <-
    (spatiotemporaloutbreaks$year.x - min(spatiotemporaloutbreaks$year.x)) * 12 +
    spatiotemporaloutbreaks$month
  
  myvars <- c("year.x", "month", "outbreak_2", "name_2", "district_country.x", "country", "code.x", "year_month")
  merged_data <- spatiotemporaloutbreaks[, myvars]
  
  countries_of_interest <- c(
    "Benin", "Burkina Faso", "Burundi", "Cameroon", "Central African Republic", "Chad",
    "Côte d'Ivoire", "Eritrea", "Ethiopia", "Gambia", "Ghana", "Guinea", "Guinea-Bissau",
    "Kenya", "Mali", "Mauritania", "Niger", "Nigeria", "Rwanda", "Senegal", "South Sudan",
    "Sudan", "Tanzania", "Togo", "Uganda"
  )
  merged_data <- merged_data[merged_data$country %in% countries_of_interest, ]
  
  # =========================================================
  # 2) LOAD SHAPEFILE + GEOMETRY CLEANING
  # =========================================================
  shape2 <- st_read("Shapefile_improved.shp", quiet = TRUE)
  shape2$district_country.x <- paste(shape2$NAME_2, shape2$COUNTRY, sep = " ")
  shape2 <- shape2[shape2$COUNTRY %in% countries_of_interest, ]
  shape2 <- st_make_valid(shape2)
  
  # =========================================================
  # 3) OUTBREAKS OVER TIME (LINE PLOT)
  # =========================================================
  merged_data <- merged_data %>%
    mutate(year_month_label = sprintf("%s %d", month.abb[month], year.x))
  
  outbreaks_over_time <- merged_data %>%
    group_by(year.x, month, year_month_label) %>%
    summarise(total_outbreaks = sum(outbreak_2, na.rm = TRUE), .groups = "drop") %>%
    arrange(year.x, month) %>%
    mutate(year_month_label = factor(year_month_label, levels = unique(year_month_label)))
  
  breaks_x <- outbreaks_over_time$year_month_label[seq(1, nrow(outbreaks_over_time), 6)]
  
  p_outbreaks <- ggplot(outbreaks_over_time, aes(x = year_month_label, y = total_outbreaks, group = 1)) +
    geom_line(color = "#3D0C5F", linewidth = 1.2) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) +
    scale_x_discrete(breaks = breaks_x, guide = guide_axis(angle = 45)) +
    labs(
      x = NULL,
      y = "Total district-months with a bacterial meningitis outbreak"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  ggsave("outbreaks_over_time.png", p_outbreaks, width = 12, height = 6, dpi = 300)
  
  # =========================================================
  # 4) CREATE NEIGHBOR MATRIX (SPATIAL ADJACENCY)
  # =========================================================
  shape2 <- st_make_valid(shape2)
  shapefile_spatial <- as(shape2, "Spatial")
  nb <- poly2nb(shapefile_spatial, queen = TRUE)
  W_matrix <- as.matrix(nb2mat(nb, style = "B", zero.policy = TRUE))
  rownames(W_matrix) <- shape2$district_country.x
  colnames(W_matrix) <- shape2$district_country.x
  
  # Check for merge consistency
  missing_in_shape <- setdiff(merged_data$district_country.x, shape2$district_country.x)
  if (length(missing_in_shape) > 0) {
    stop("districts are in merged_data but not in shapefile: ", paste(missing_in_shape, collapse = ", "))
  }
  
  missing_in_data <- setdiff(shape2$district_country.x, merged_data$district_country.x)
  if (length(missing_in_data) > 0) {
    message("Dropping ", length(missing_in_data), " districts from shapefile:\n",
            paste(missing_in_data, collapse = ", "))  # <- typo fix from "ropping"
    shape2 <- shape2[shape2$district_country.x %in% merged_data$district_country.x, ]
    W_matrix <- W_matrix[shape2$district_country.x, shape2$district_country.x, drop = FALSE]
  }
  
  if (!identical(sort(unique(merged_data$district_country.x)),
                 sort(unique(shape2$district_country.x)))) {
    stop("District-country lists differ")
  } else {
    message("Matches verified.")
  }
  
  # Ensure consistent ordering
  shape2 <- shape2[order(shape2$district_country.x), ]
  merged_data <- merged_data[order(merged_data$district_country.x), ]
  merged_data$area <- match(merged_data$district_country.x, shape2$district_country.x)
  
  # =========================================================
  # 5) IDENTIFY FIRST OUTBREAKS (GAP > 2 MONTHS RULE)
  # =========================================================
  first_outbreaks <- merged_data %>%
    filter(outbreak_2 == 1) %>%
    arrange(district_country.x, year_month) %>%
    group_by(district_country.x) %>%
    mutate(
      prev_year_month = lag(year_month),
      prev_year = lag(year.x),
      months_since_last = year_month - prev_year_month,
      is_first = case_when(
        is.na(prev_year_month) ~ TRUE,
        months_since_last > 2 ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>%
    filter(is_first == TRUE) %>%
    select(-prev_year_month, -prev_year, -months_since_last, -is_first) %>%
    ungroup()
  
  print(head(first_outbreaks))
  print(names(first_outbreaks))
  
  # =========================================================
  # 6) CLUSTER OUTBREAKS BY SPACE + TIME
  # =========================================================
  identify_clusters <- function(outbreaks_df, W_matrix) {
    all_clusters <- list()
    cluster_id <- 1
    processed_rows <- integer(0)
    
    outbreaks_df <- outbreaks_df %>%
      arrange(year_month) %>%
      mutate(row_id = row_number())
    
    for (i in 1:nrow(outbreaks_df)) {
      if (i %in% processed_rows) next
      
      current_district <- outbreaks_df$district_country.x[i]
      current_year <- outbreaks_df$year.x[i]
      current_year_month <- outbreaks_df$year_month[i]
      
      # same year, within ±2 months
      concurrent_outbreaks <- outbreaks_df %>%
        filter(year.x == current_year, abs(year_month - current_year_month) <= 2)
      
      affected_districts <- concurrent_outbreaks$district_country.x
      valid_districts <- affected_districts[affected_districts %in% rownames(W_matrix)]
      if (length(valid_districts) == 0) next
      
      W_subset <- W_matrix[valid_districts, valid_districts, drop = FALSE]
      g <- igraph::graph_from_adjacency_matrix(W_subset, mode = "undirected")
      comp <- igraph::components(g)
      
      if (current_district %in% names(comp$membership)) {
        component_id <- comp$membership[current_district]
        cluster_districts <- names(comp$membership[comp$membership == component_id])
        
        cluster_data <- concurrent_outbreaks %>%
          filter(district_country.x %in% cluster_districts) %>%
          mutate(outbreak_cluster_id = cluster_id)
        
        all_clusters[[cluster_id]] <- cluster_data
        processed_rows <- c(processed_rows, cluster_data$row_id)
        cluster_id <- cluster_id + 1
      }
    }
    
    dplyr::bind_rows(all_clusters)
  }
  
  outbreak_clusters <- identify_clusters(first_outbreaks, W_matrix)
  
  # =========================================================
  # 7) COUNT UNIQUE OUTBREAK CLUSTERS PER DISTRICT
  # =========================================================
  district_unique_outbreaks <- outbreak_clusters %>%
    distinct(district_country.x, outbreak_cluster_id)
  
  district_totals <- shape2 %>%
    st_drop_geometry() %>%
    select(district_country.x) %>%
    left_join(
      district_unique_outbreaks %>% count(district_country.x, name = "n_outbreaks"),
      by = "district_country.x"
    ) %>%
    mutate(
      n_outbreaks = ifelse(is.na(n_outbreaks), 0, n_outbreaks),
      fill_var = ifelse(n_outbreaks == 0, NA, n_outbreaks)
    )
  
  map_data <- shape2 %>%
    left_join(district_totals, by = "district_country.x")
  
  # =========================================================
  # 8) MAP OUTBREAK INTENSITY + AFRICA BACKGROUND
  # =========================================================
  africa_bg <- ne_countries(continent = "Africa", returnclass = "sf")
  
  p_map <- ggplot() +
    geom_sf(data = africa_bg, fill = "grey96", color = "grey85", size = 0.15) +
    geom_sf(data = map_data, aes(fill = fill_var), color = "grey30", size = 0.1) +
    scale_fill_viridis_c(
      option = "B",
      name = "Total outbreaks",
      na.value = "grey90",
      direction = -1
    ) +
    coord_sf(
      xlim = st_bbox(africa_bg)[c("xmin", "xmax")],
      ylim = st_bbox(africa_bg)[c("ymin", "ymax")]
    ) +
    theme_minimal(base_family = "sans", base_size = 14) +
    theme(
      panel.grid = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      plot.title = element_text(size = 17, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
  
  ggsave("01_outbreak_intensity_overview_nodrc.pdf", p_map, width = 13, height = 9, dpi = 300)


}

