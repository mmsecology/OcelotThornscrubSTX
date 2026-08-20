library(landscapemetrics)
library(terra)
library(purrr)
library(tidyverse)

# =================================================================================
# Characterize and model thornscrub metrics within ocelot UDs
# =================================================================================

ud_polygons <- readRDS("output/ctmm_occurrence_ud_polygons_sf.rds")

# combines separate ud polygons at each level into single spatial object
ud_analysis_areas <- purrr::imap(
  ud_polygons, ~.x %>% group_by(Deployment_ID, level) %>%
    summarise(n_components = n(), total_area_km2 = sum(area_km2), 
              geometry = sf::st_union(geometry), .groups = "drop")
            )

ud_analysis_areas$EO19M_1

ud_components <- ud_polygons$EO19M_1 %>% arrange(level, area_km2)
ggplot(ud_components,
       aes(x = area_km2)) +
  geom_histogram() +
  facet_wrap(~ level, scales = "free_y") +
  scale_x_log10() + theme_bw(base_size = 20)


library(dplyr)
ud_polygons_df <- purrr::map_dfr(ud_polygons, identity)

ud_component_rank <- ud_polygons_df %>%
  group_by(Deployment_ID, level) %>%
  arrange(desc(area_km2), .by_group = TRUE) %>%
  mutate(
    n_components = n(),
    component_rank = row_number(),
    total_area_km2 = sum(area_km2),
    cumulative_area_km2 = cumsum(area_km2),
    cumulative_prop = cumulative_area_km2 / total_area_km2
  ) %>%
  ungroup()

area_thresholds <- c(0.90, 0.95, 0.99)

ud_component_summary <- purrr::map_dfr(
  area_thresholds,
  function(threshold) {
    
    ud_component_rank %>%
      group_by(Deployment_ID, level) %>%
      filter(cumulative_prop >= threshold) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      mutate(area_threshold = threshold) %>%
      select(
        Deployment_ID,
        level,
        area_threshold,
        n_components = component_rank,
        total_components = n_components,
        total_area_km2,
        retained_area_km2 = cumulative_area_km2,
        retained_prop = cumulative_prop,
        smallest_component_km2 = area_km2
      )
  }
)

ud_component_summary <- ud_component_summary %>%
  rename(
    components_to_threshold = n_components,
    original_n_components = total_components
  )

ud_component_summary

ud_component_summary %>%
  group_by(level, area_threshold) %>%
  summarise(
    median_components = median(components_to_threshold),
    mean_components = mean(components_to_threshold),
    min_components = min(components_to_threshold),
    max_components = max(components_to_threshold),
    median_smallest_component = median(smallest_component_km2),
    .groups = "drop"
  )

ggplot(
  ud_component_summary,
  aes(
    x = level,
    y = components_to_threshold,
    group = Deployment_ID
  )
) +
  geom_line(alpha = 0.2) +
  geom_point(alpha = 0.4) +
  facet_wrap(~area_threshold) +
  scale_y_log10() +
  labs(
    x = "UD level",
    y = "Components required to retain cumulative UD area",
    title = "Number of occurrence components required to retain UD area"
  ) +
  theme_minimal()

ud_component_summary %>%
  mutate(
    prop_components_retained =
      components_to_threshold / original_n_components
  ) %>%
  select(
    Deployment_ID,
    level,
    area_threshold,
    components_to_threshold,
    original_n_components,
    prop_components_retained,
    retained_prop
  )

ud_component_summary %>%
  mutate(
    prop_components_retained =
      components_to_threshold / original_n_components
  ) %>%
  group_by(level, area_threshold) %>%
  summarise(
    median_components = median(components_to_threshold),
    median_prop_components = median(prop_components_retained),
    median_prop_area = median(retained_prop),
    .groups = "drop"
  )


retain_ud_components <- function(x, thresholds = c(0.90, 0.95)) {
  
  x %>%
    group_by(Deployment_ID, level) %>%
    arrange(desc(area_km2), .by_group = TRUE) %>%
    mutate(
      total_ud_area_km2 = sum(area_km2),
      cumulative_area_km2 = cumsum(area_km2),
      cumulative_prop = cumulative_area_km2 / total_ud_area_km2
    ) %>%
    group_modify(~ {
      
      map_dfr(thresholds, function(threshold) {
        
        # First component that reaches the desired cumulative area
        cutoff_rank <- which(.x$cumulative_prop >= threshold)[1]
        
        .x %>%
          slice_head(n = cutoff_rank) %>%
          mutate(
            area_threshold = threshold
          )
      })
      
    }) %>%
    ungroup()
}

ud_retained <- retain_ud_components(ud_polygons_df)


ud_retained %>%
  count(Deployment_ID, level, area_threshold)

ud_component_rank %>%
  filter(
    Deployment_ID == "EO19M_1",
    level == 0.5
  ) %>%
  ggplot(aes(component_rank, cumulative_prop)) +
  geom_line() +
  geom_hline(
    yintercept = c(0.90, 0.95),
    linetype = "dashed"
  ) +
  labs(
    x = "Component rank",
    y = "Cumulative proportion of UD area"
  ) +
  theme_minimal()

retain_ud_components <- function(x, thresholds = c(0.90, 0.95)) {
  
  x %>%
    group_by(Deployment_ID, level) %>%
    arrange(desc(area_km2), .by_group = TRUE) %>%
    mutate(
      total_ud_area_km2 = sum(area_km2),
      original_n_components = n(),
      component_rank = row_number(),
      cumulative_area_km2 = cumsum(area_km2),
      cumulative_prop = cumulative_area_km2 / total_ud_area_km2
    ) %>%
    group_modify(~ {
      
      map_dfr(thresholds, function(threshold) {
        
        cutoff_rank <- which(.x$cumulative_prop >= threshold)[1]
        
        .x %>%
          slice_head(n = cutoff_rank) %>%
          mutate(
            area_threshold = threshold,
            retained_n_components = cutoff_rank,
            retained_area_km2 = sum(area_km2),
            retained_prop = sum(area_km2) / first(total_ud_area_km2)
          )
      })
      
    }) %>%
    ungroup()
}

ud_retained <- retain_ud_components(ud_polygons_df)

ud_retained_multipolygon <- ud_retained %>%
  group_by(Deployment_ID, level, area_threshold) %>%
  summarise(
    n_retained_components = n(),
    retained_area_km2 = sum(area_km2),
    original_n_components = first(original_n_components),
    total_ud_area_km2 = first(total_ud_area_km2),
    .groups = "drop"
  )

ud_retained <- st_as_sf(
  ud_retained,
  sf_column_name = "geometry",
  crs = st_crs(ud_polygons_df)
)
class(ud_retained)
st_geometry_type(ud_retained)

ud_retained_multipolygon <- ud_retained %>%
  group_by(
    Deployment_ID,
    level,
    area_threshold
  ) %>%
  summarise(
    n_retained_components = n(),
    retained_area_km2 = sum(area_km2),
    original_n_components = first(original_n_components),
    total_ud_area_km2 = first(total_ud_area_km2),
    .groups = "drop"
  )

plot(ud_retained_multipolygon)
ud_retained_multipolygon %>%
  filter(
    Deployment_ID == "EO19M_1",
    level == 0.5
  ) %>%
  select(
    Deployment_ID,
    level,
    area_threshold,
    n_retained_components,
    retained_area_km2
  )


# Function to crop and mask thornscrub raster by ud
get_landscape <- function(r, ud_area) {  
  r_crop <- terra::crop(r, terra::vect(ud_area))
  r_mask <- terra::mask(r_crop, terra::vect(ud_area))
  r_mask
}


thornscrub_binary <- rast("output/south_texas_thornscrub_binary.tif")
strict_thornscrub_binary <- rast("output/strict_south_texas_thornscrub_binary.tif")
terra::same.crs(thornscrub_binary, strict_thornscrub_binary)

ud_analysis_areas <- lapply(ud_analysis_areas, function(x) st_transform(x, crs = crs(thornscrub_binary)))
ud_retained_multipolygon <- st_transform(ud_retained_multipolygon, crs = crs(thornscrub_binary))

########## Test #################################################################
r_19_50 <- get_landscape(
  thornscrub_binary,
  ud_analysis_areas$EO19M_1 %>% filter(level == 0.50)
)

r_19_50
terra::plot(r_19_50)
freq(r_19_50)

r_19_50 <- get_landscape(
  strict_thornscrub_binary,
  ud_analysis_areas$EO19M_1 %>% filter(level == 0.50)
)

r_19_50
terra::plot(r_19_50)
freq(r_19_50)

test <- ud_analysis_areas$EO19M_1 %>% filter(level == 0.50)
class_metrics <- calculate_lsm(r_19_50, what = c("lsm_c_ai", "lsm_c_contig_mn", "lsm_c_cohesion", "lsm_c_gyrate_mn", "lsm_c_pland"))
print(class_metrics)

class.metrics <- lapply(lc.use.areas, function(x) calculate_lsm(x, what = c("lsm_c_ai", "lsm_c_contig_mn", "lsm_c_cohesion", "lsm_c_gyrate_mn", 
                                                                            "lsm_c_pland", "lsm_c_pd", "lsm_c_area_mn", "lsm_c_shape_mn", "lsm_c_lpi")))
metrics.list <- lapply(class.metrics, function(x) filter(x, class == 1))
metrics <- lapply(metrics.list, function(x) pivot_wider(x, names_from = metric, values_from = value))
class.metrics.df <-  bind_rows(metrics, .id = "column_label")

lsm_metrics <- c(
  "lsm_c_ai",
  "lsm_c_contig_mn",
  "lsm_c_cohesion",
  "lsm_c_gyrate_mn",
  "lsm_c_pland",
  "lsm_c_pd",
  "lsm_c_area_mn",
  "lsm_c_shape_mn",
  "lsm_c_lpi"
)

metrics <- purrr::imap_dfr(
  ud_retained_multipolygon ,
  function(ud_df, deployment_id) {
    purrr::map_dfr(
      seq_len(nrow(ud_df)),
      function(i) {
        r <- get_landscape(
          thornscrub_binary,
          ud_df[i, ]
        )
        calculate_lsm(
          r,
          what = lsm_metrics
        ) %>%
          filter(class == 1) %>%
          mutate(
            Deployment_ID = deployment_id,
            UD_level = ud_df$level[i]
          )
      }
    )
  }
)

metrics %>% count(Deployment_ID, UD_level)

metrics_wide <- metrics %>%
  select(
    Deployment_ID,
    UD_level,
    metric,
    value
  ) %>%
  pivot_wider(
    names_from = metric,
    values_from = value
  )

ggplot(data = metrics) + 
  geom_boxplot(aes(x = as.factor(UD_level), y = value)) + 
  facet_wrap(~metric, scales = "free_y") +
  theme_bw(base_size = 20)

ggplot(metrics_wide,
       aes(x = UD_level, y = pland,
           group = Deployment_ID)) +
  geom_line(alpha = 0.3) +
  geom_point() +
  scale_x_continuous(
    breaks = c(.1, .25, .5, .75, .95)
  )

ud_areas <- purrr::imap_dfr(
  ud_analysis_areas,
  function(ud_df, deployment_id) {
    ud_df %>%
      sf::st_drop_geometry() %>%
      dplyr::select(level, total_area_km2) %>%
      dplyr::mutate(
        Deployment_ID = deployment_id,
        UD_level = level
      ) %>%
      dplyr::select(Deployment_ID, UD_level, total_area_km2)
  }
)

metrics_long <- metrics %>%
  left_join(
    ud_areas,
    by = c("Deployment_ID", "UD_level")
  )

ggplot(metrics_wide,
       aes(x = total_area_km2, y = cohesion,
           group = Deployment_ID)) +
  geom_line(alpha = 0.3) +
  geom_point() +
  scale_x_log10()

ggplot(data = metrics_long, aes(x = total_area_km2, y = value,
           group = Deployment_ID)) + 
  geom_line(alpha = 0.3) + 
  geom_point() + 
  facet_wrap(~metric, scales = "free_y") +
  theme_bw(base_size = 20) +
  scale_x_log10()


library(brms)

m_pland <- brm(
  pland ~ UD_level + (1 | Deployment_ID),
  data = metrics_wide,
  family = gaussian(),
  chains = 4,
  cores = 4,
  iter = 4000
)

summary(m_pland)
conditional_effects(m_pland)

m_pland_area <- brm(
  pland ~ log(total_area_km2) + (1 | Deployment_ID),
  data = metrics_wide,
  family = gaussian(),
  chains = 4,
  iter = 4000
)

summary(m_pland_area)
conditional_effects(m_pland_area)

cor(
  metrics_wide$UD_level,
  log(metrics_wide$total_area_km2),
  use = "complete.obs"
)




thorn_lsm_areas <- purrr::pmap(
  ud_retained_multipolygon %>%
    st_drop_geometry() %>%
    select(Deployment_ID, level, area_threshold),
  
  function(Deployment_ID, level, area_threshold) {
    
    area <- ud_retained_multipolygon %>%
      filter(
        .data$Deployment_ID == Deployment_ID,
        .data$level == level,
        .data$area_threshold == area_threshold
      )
    
    r_crop <- terra::crop(thornscrub_binary, terra::vect(area))
    r_mask <- terra::mask(r_crop, terra::vect(area))
    
    r_mask
  }
)




metrics <- purrr::map_dfr(
  seq_len(nrow(ud_retained_multipolygon)),
  function(i) {
    
    ud_area <- ud_retained_multipolygon[i, ]
    
    r <- get_landscape(
      thornscrub_binary,
      ud_area
    )
    
    calculate_lsm(
      r,
      what = lsm_metrics
    ) %>%
      filter(class == 1) %>%
      mutate(
        Deployment_ID = ud_area$Deployment_ID,
        UD_level = ud_area$level,
        area_threshold = ud_area$area_threshold
      )
  }
)

ggplot(data = metrics[metrics$area_threshold == 0.9,]) + 
  geom_boxplot(aes(x = as.factor(UD_level), y = value)) + 
  facet_wrap(~metric, scales = "free_y") +
  theme_bw(base_size = 20)

ggplot(data = metrics) + 
  geom_boxplot(aes(x = as.factor(UD_level), y = value, color = as.factor(area_threshold))) + 
  facet_wrap(~metric, scales = "free_y") +
  theme_bw(base_size = 20)



patches_r <- terra::patches(thornscrub_binary, directions = 8, zeroAsNA = TRUE)

patch_metrics <- landscapemetrics::calculate_lsm(
  thornscrub_binary,
  level = "patch",
  metric = c("area", "shape", "core", "contig")
)