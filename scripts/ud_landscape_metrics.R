library(sf); library(terra); library(ctmm)
library(landscapemetrics); library(dplyr)
library(purrr); library(ggplot2)
library(tidyterra); library(patchwork)
library(tictoc); library(Hmisc)


# ---------------------------------------------------------------
# Read in objects
# ---------------------------------------------------------------

ud_list <- readRDS("output/ctmm_occurrence_uds.rds")
ecomap_stx <- rast("output/south_texas_thornscrub_binary.tif")
nlcd_stx <- rast("output/nlcd_stx.tif")
thornscrub_binary <- rast("output/south_texas_thornscrub_binary.tif")
envelope_available <- st_read("output/available_area.shp")
ud_retained_sf <- readRDS("output/ud_retained_sf.rds")

# ---------------------------------------------------------------
# Get available landsacpe from locations
# ---------------------------------------------------------------

get_ud_polygons <- function(ud, levels = c(0.99)) {
  purrr::map_dfr(levels, function(lvl) {
    x <- as.sf(ud, level.UD = lvl)
    x %>% sf::st_cast("POLYGON") %>%
      dplyr::mutate(level = lvl, area_km2 = as.numeric(sf::st_area(geometry)) / 1e6)
  })
}

target_crs <- crs(thornscrub_binary)  # your working Albers CRS, already used for the patch mosaic

ud_polygons <- purrr::imap(ud_list, ~get_ud_polygons(.x) %>% mutate(Deployment_ID = .y))
ud_polygons_sf <- purrr::map_dfr(ud_polygons, identity) %>% st_transform(target_crs)

envelope_raw <- ud_polygons_sf %>% st_union() %>% st_sf(geometry = .)

mapview::mapview(envelope_raw)

# Buffer to close small gaps
buffer_dist <- 2000  # meters, placeholder - justify from your data
envelope <- st_buffer(envelope_raw, dist = buffer_dist) %>% st_union() %>% st_sf(geometry = .)
mapview::mapview(envelope)

envelope_hull <- st_union(envelope_raw) %>% st_convex_hull() %>% st_sf(geometry = .)

envelope <- st_buffer(envelope_hull, dist = 1000) %>% st_union() %>% st_sf(geometry = .)
mapview::mapview(envelope)

terra::freq(ecomap_stx)  # look for a 0 category, and check if NA cells exist separately
water_mask <- terra::ifel(nlcd_stx == 11, 1, NA)
terra::plot(water_mask)

exclude_poly <- terra::as.polygons(water_mask, dissolve = TRUE) %>%
  st_as_sf() %>% st_make_valid() %>% st_transform(st_crs(envelope))

envelope_available <- st_difference(envelope, st_union(exclude_poly))

mapview::mapview(envelope_available)
st_write(envelope_available, "output/available_area.shp", append = FALSE)

crs(thornscrub_binary)
crs(envelope_raw)
crs(envelope)
crs(envelope_available)

# --------------------------------------------------------------------
# Patch attribute full landscape (binary thornscrub raster)
# --------------------------------------------------------------------

# thornscrub_binary: 1 = thornscrub, 0/NA = non-thornscrub
patch_metrics <- landscapemetrics::calculate_lsm(
  thornscrub_binary,
  directions = 8,
  level = "patch"
)

class_metrics <- landscapemetrics::calculate_lsm(
  thornscrub_binary,
  directions = 8,
  level = "class"
) 

# reshape wide - one row per patch, one column per metric
patch_table <- patch_metrics |>
  tidyr::pivot_wider(id_cols = c(class, id), names_from = metric, values_from = value)

saveRDS(patch_table, "output/thornscub_patch_table.rds")

# get patch polygons for spatial joins (id matches patch_table$id)
patch_polys <- landscapemetrics::get_patches(thornscrub_binary, directions = 8, class = 1)[[1]][[1]] |>
  terra::as.polygons(dissolve = TRUE) |>
  st_as_sf() |>
  rename(id = 1) |>
  left_join(patch_table, by = "id")

# re-verify clean join
nrow(patch_polys)                          # should still be ~7978
sum(is.na(patch_polys$area))               # should be 0
range(sf::st_drop_geometry(patch_polys)$id, na.rm = TRUE)  # should now be 2657-10634, matching table exactly

patch_centroids <- sf::st_centroid(patch_polys) %>% st_transform(crs = st_crs(envelope_available))
in_envelope <- sf::st_within(patch_centroids, envelope_available, sparse = FALSE)[, 1]
patch_polys_available <- patch_polys[in_envelope, ]

nrow(patch_polys_available)
summary(patch_polys_available$area)
patch_polys_available |> sf::st_drop_geometry() |> dplyr::arrange(desc(area)) |> head(10)

ggplot() +
  geom_sf(data = patch_polys_available, aes(fill = area), color = NA) +
  geom_sf(data = envelope_available, fill = NA, color = "black", linewidth = 0.6) +
  scale_fill_viridis_c(name = "Patch area (ha)", trans = "log10") +
  theme_minimal(base_size = 25) +
  labs(title = NULL)

# MEGA_PATCH_ID <- 2782

# get_captured_summary <- function(ud_poly, patches_sf, patch_id_col = "id") {
#   patch_centroids <- sf::st_centroid(patches_sf)
#   captured <- patches_sf[sf::st_within(patch_centroids, ud_poly, sparse = FALSE)[,1], ]
#   captured_df <- sf::st_drop_geometry(captured)
  
#   if (nrow(captured_df) == 0) {
#     return(tibble::tibble(n_patches = 0, total_captured_area = 0, contains_mega_patch = FALSE))
#   }
  
#   metric_cols <- c("area", "cai", "circle", "contig", "core", "enn", 
#                     "frac", "gyrate", "ncore", "para", "perim", "shape")
#   metric_cols <- intersect(metric_cols, names(captured_df))  # only summarize what exists
  
#   summary_stats <- captured_df |>
#     summarise(across(all_of(metric_cols),
#                       list(mean = ~mean(.x, na.rm = TRUE),
#                            median = ~median(.x, na.rm = TRUE)),
#                       .names = "{.col}_{.fn}"))
  
#   summary_stats |>
#     mutate(
#       n_patches = nrow(captured_df),
#       total_captured_area = sum(captured_df$area, na.rm = TRUE),
#       contains_mega_patch = MEGA_PATCH_ID %in% captured_df[[patch_id_col]]
#     )
# }

# --------------------------------------------------------------------
# Get 95% and 50% UD and retain polygons up to 90% of total area
# --------------------------------------------------------------------

ud_polygons <- purrr::imap(ud_list, ~get_ud_polygons(.x, levels = c(0.50, 0.95)) %>% 
                              mutate(Deployment_ID = .y))

ud_polygons_sf <- purrr::map_dfr(ud_polygons, identity)

retain_ud_components <- function(x, component_retention = 0.90) {
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
      cutoff <- which(.x$cumulative_prop >= component_retention)[1]
      .x %>%
        slice_head(n = cutoff) %>%
        mutate(
          component_retention = component_retention,
          retained_n_components = cutoff,
          retained_area_km2 = sum(area_km2),
          retained_prop = retained_area_km2 / first(total_ud_area_km2)
        )
    }) %>%
    ungroup()
}

ud_retained <- retain_ud_components(ud_polygons_sf, component_retention = 0.90)

ud_retained_sf <- ud_retained %>%
  group_by(Deployment_ID, level) %>%
  summarise(
    geometry = st_union(geometry),
    retained_n_components = n(),
    retained_area_km2 = sum(area_km2),
    total_ud_area_km2 = first(total_ud_area_km2),
    original_n_components = first(original_n_components),
    retained_prop = sum(area_km2) / first(total_ud_area_km2),
    .groups = "drop"
  ) %>%
  st_as_sf(sf_column_name = "geometry", crs = st_crs(ud_polygons_sf)) %>%
  st_transform(crs(thornscrub_binary))

saveRDS(ud_retained_sf, "output/ud_retained_sf.rds")

# --------------------------------------------------------------------
# Create null distribution of UDs
# --------------------------------------------------------------------

# random_shift() unchanged from before - works on any single-row sf polygon
random_shift <- function(poly, envelope, rotate = TRUE, max_attempts = 500) {
  poly_geom <- sf::st_geometry(poly)
  centroid <- sf::st_centroid(poly_geom)
  env_bbox <- sf::st_bbox(envelope)
  
  for (i in seq_len(max_attempts)) {
    ang <- if (rotate) runif(1, 0, 360) else 0
    
    if (ang != 0) {
      ang_rad <- ang * pi / 180
      rot_mat <- matrix(c(cos(ang_rad), sin(ang_rad),
                          -sin(ang_rad), cos(ang_rad)), 2, 2)
      centroid_coords <- sf::st_coordinates(centroid)[1, ]
      rotated <- (poly_geom - centroid_coords) * rot_mat + centroid_coords
      rotated <- sf::st_set_crs(rotated, sf::st_crs(poly))
    } else {
      rotated <- poly_geom
    }
    
    target_x <- runif(1, env_bbox["xmin"], env_bbox["xmax"])
    target_y <- runif(1, env_bbox["ymin"], env_bbox["ymax"])
    target_pt <- sf::st_sfc(sf::st_point(c(target_x, target_y)), crs = sf::st_crs(envelope))
    
    if (!sf::st_within(target_pt, envelope, sparse = FALSE)[1, 1]) next
    
    rotated_centroid_coords <- sf::st_coordinates(sf::st_centroid(rotated))[1, ]
    shift_vec <- c(target_x, target_y) - rotated_centroid_coords
    shifted <- rotated + shift_vec
    shifted_sf <- sf::st_sfc(shifted, crs = sf::st_crs(poly))
    
    if (sf::st_within(shifted_sf, envelope, sparse = FALSE)[1, 1]) {
      return(sf::st_sf(geometry = shifted_sf))
    }
  }
  warning(paste("max_attempts reached -", poly$Deployment_ID[1], "level", poly$level[1]))
  return(NULL)
}

# ---------------------------------------------------------------
# Generate null placements directly from ud_retained_sf rows
# ---------------------------------------------------------------

ud_retained_sf <- st_transform(ud_retained_sf, crs = crs(envelope_available))

n_null <- 999

tictoc::tic()
null_uds <- pmap_dfr(
  list(ud_retained_sf$Deployment_ID, seq_len(nrow(ud_retained_sf))),
  function(dep_id, row_i) {
    poly <- ud_retained_sf[row_i, ]
    
    map_dfr(1:n_null, function(rep_i) {
      shifted <- random_shift(poly, envelope_available, rotate = TRUE)
      if (is.null(shifted)) return(NULL)
      shifted %>%
        mutate(
          Deployment_ID = dep_id,
          level = poly$level,
          rep = rep_i,
          source = "null"
        )
    })
  }
)
tictoc::toc()

saveRDS(null_uds, "output/null.rds")

crs(null_uds)

# check placement success rate per individual
null_uds %>% st_drop_geometry() %>% count(Deployment_ID, level) %>% filter(n < n_null * 0.8)

# Visual test - one individual, one level, a handful of replicates
test_poly <- ud_retained_sf %>% filter(Deployment_ID == "EO34M_2")
test_nulls <- map(1:10, ~ random_shift(test_poly, envelope_available, rotate = TRUE)) |> compact()

envelope_available <- st_transform(envelope_available, st_crs(ud_retained_sf))

test_nulls_sf <- test_nulls %>% bind_rows() %>% mutate(rep = row_number())

ggplot() +
  geom_sf(data = envelope_available, fill = NA) +
  geom_sf(data = test_poly, fill = "red", alpha = 0.4) +
  geom_sf(data = test_nulls_sf, fill = "blue", alpha = 0.4) +
  coord_sf()

# ---------------------------------------------------------------
# Simplified landscape metrics run
# ---------------------------------------------------------------

class_metrics <- c("lsm_c_pland", "lsm_c_area_mn", "lsm_c_shape_mn", "lsm_c_cai_mn", "lsm_c_enn_mn", "lsm_c_clumpy")

ud_retained_sf <- st_transform(ud_retained_sf, crs = crs(thornscrub_binary))
crs(ud_retained_sf)

get_landscape <- function(r, ud_area) {  
  r_crop <- terra::crop(r, terra::vect(ud_area))
  r_mask <- terra::mask(r_crop, terra::vect(ud_area))
  r_mask
}

calculate_ud_metrics <- function(ud, landscape, metrics) {
  r <- get_landscape(landscape, ud)
  calculate_lsm(r, what = metrics) %>%
    filter(class == 1) %>%
    mutate(Deployment_ID = ud$Deployment_ID, UD_level = ud$level)
}

observed_metrics <- purrr::map_dfr(
  seq_len(nrow(ud_retained_sf)),
  ~ calculate_ud_metrics(
      ud = ud_retained_sf[.x, ],
      landscape = thornscrub_binary,
      metrics = class_metrics
    )
)

null_uds <- st_transform(null_uds, st_crs(ud_retained_sf))
same.crs(null_uds, ud_retained_sf)

null_metrics <- purrr::map_dfr(
  seq_len(nrow(null_uds)),
  function(i) {
    calculate_ud_metrics(
      ud = null_uds[i, ],
      landscape = thornscrub_binary,
      metrics = class_metrics
    ) %>%
      mutate(rep = null_uds$rep[i])
  }
)

saveRDS(null_metrics, "output/null_ud_metrics.rds")

# ---------------------------------------------------------------
# Compare observed to null
# ---------------------------------------------------------------

comparison <- null_metrics %>%
  left_join(
    observed_metrics %>% select(Deployment_ID, UD_level, metric, observed_value = value),
    by = c("Deployment_ID", "UD_level", "metric")
  ) %>%
  group_by(Deployment_ID, UD_level, metric) %>%
  summarize(
    observed_value = first(observed_value),
    null_mean = mean(value, na.rm = TRUE),
    null_sd = sd(value, na.rm = TRUE),
    percentile_rank = mean(value <= observed_value, na.rm = TRUE),
    ses = (first(observed_value) - mean(value, na.rm = TRUE)) / sd(value, na.rm = TRUE),
    .groups = "drop"
  )
saveRDS(comparison, "output/obs_null_comparison_uds.rds")

percentile_ranks <- null_metrics %>%
  left_join(
    observed_metrics %>% select(Deployment_ID, UD_level, metric, observed_value = value),
    by = c("Deployment_ID", "UD_level", "metric")
  ) %>%
  group_by(Deployment_ID, UD_level, metric) %>%
  summarize(
    observed_value = first(observed_value),
    percentile_rank = mean(value <= observed_value, na.rm = TRUE),
    .groups = "drop"
  )

percentile_ranks

ggplot(comparison, aes(x = ses)) +
  geom_histogram(bins = 30) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  facet_grid(UD_level ~ metric, scales = "free_x") +
  theme_bw(base_size = 25) +
  labs(x = "Standardized effect size", y = "Number of individuals")



ggplot(comparison, aes(x = metric, y = ses)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_jitter(aes(color = ses > 0), width = 0.15, size = 3, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", size = 5, color = "black") +
  facet_wrap(~ UD_level) +
  coord_flip() +
  scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "firebrick"), guide = "none") +
  theme_bw(base_size = 18) +
  labs(x = NULL, y = "Standardized effect size")

ggplot(comparison, aes(x = ses, y = reorder(Deployment_ID, ses))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3) +
  facet_grid(UD_level ~ metric, scales = "free_x") +
  theme_bw(base_size = 14) +
  labs(x = "Standardized effect size", y = NULL)

ggplot(comparison, aes(x = ses)) +
  geom_density(fill = "grey80", alpha = 0.6) +
  geom_rug(sides = "b") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  facet_grid(UD_level ~ metric, scales = "free_x") +
  theme_bw(base_size = 18)

ggplot(comparison, aes(x = ses)) +
  geom_density(aes(y = after_stat(scaled)), fill = "grey80", color = "grey30", alpha = 0.6) +
  geom_rug(sides = "b", alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  facet_grid(UD_level ~ metric, scales = "free_x") +
  theme_bw(base_size = 18) +
  labs(x = "Standardized effect size", y = "Scaled density")

ggplot(comparison, aes(x = metric, y = ses)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_jitter(aes(color = ses > 0), width = 0.15, size = 3, alpha = 0.8) +
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange", size = 0.8, color = "black") +
  facet_wrap(~ UD_level) +
  coord_flip() +
  scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "firebrick"), guide = "none") +
  theme_bw(base_size = 18) +
  labs(x = NULL, y = "Standardized effect size")

# Pre-compute bootstrapped mean + CI per metric/level
summary_df <- comparison %>%
  group_by(metric, UD_level) %>%
  dplyr::summarize(
    boot = list(Hmisc::smean.cl.boot(ses)),
    .groups = "drop"
  ) %>%
  mutate(
    mean_ses = map_dbl(boot, ~ .x[["Mean"]]),
    lower    = map_dbl(boot, ~ .x[["Lower"]]),
    upper    = map_dbl(boot, ~ .x[["Upper"]])
  ) %>%
  select(-boot) %>%
  mutate(
    sig = case_when(
      lower > 0 ~ "positive",
      upper < 0 ~ "negative",
      TRUE ~ "ns"
    )
  )

  ## GO WITH THIS ONE; still need to change a few things
ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_jitter(
    data = comparison,
    aes(x = metric, y = ses, color = ses > 0),
    width = 0.15, size = 3, alpha = 0.5
  ) +
  geom_pointrange(
    data = summary_df,
    aes(x = metric, y = mean_ses, ymin = lower, ymax = upper, fill = sig),
    shape = 23, size = 1, linewidth = 1, color = "black"
  ) +
  facet_wrap(~ UD_level) +
  coord_flip() +
  scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "firebrick"), guide = "none") +
  scale_fill_manual(
    values = c("positive" = "firebrick", "negative" = "steelblue", "ns" = "grey40"),
    guide = "none"
  ) +
  theme_bw(base_size = 18) +
  labs(x = NULL, y = "Standardized effect size")
#e.g. "points are per-individual SES (observed relative to n = 999 individual-specific null UDs); black points/ranges are bootstrapped mean ± 95% CI across individuals (n = 28)."