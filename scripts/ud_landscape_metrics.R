library(sf)
library(terra)
library(ctmm)
library(landscapemetrics)
library(dplyr)
library(purrr)
library(ggplot2)
library(tidyterra)
library(patchwork)

# ---------------------------------------------------------------
# Get available landsacpe from locations
# ---------------------------------------------------------------

ud_list <- readRDS("output/ctmm_occurrence_uds.rds")

get_ud_polygons <- function(ud, levels = c(0.99)) {
  purrr::map_dfr(levels, function(lvl) {
    x <- as.sf(ud, level.UD = lvl)
    x %>% sf::st_cast("POLYGON") %>%
      dplyr::mutate(level = lvl, area_km2 = as.numeric(sf::st_area(geometry)) / 1e6)
  })
}

ud_polygons <- purrr::imap(ud_list, ~get_ud_polygons(.x) %>% mutate(Deployment_ID = .y))

# Union across individuals into a single multipolygon
envelope_raw <- do.call(rbind, ud_polygons) |>
  st_union() |>
  st_sf(geometry = _)

mapview::mapview(envelope_raw)

# Buffer to close small gaps and give a bit of "reachable but unsampled" margin
# choose buffer_dist based on something defensible - e.g. median individual home range diameter
buffer_dist <- 2000  # meters, placeholder - justify from your data
envelope <- st_buffer(envelope_raw, dist = buffer_dist) |>
  st_union() |>
  st_sf(geometry = _)
mapview::mapview(envelope)

envelope_hull <- st_union(envelope_raw) |>
  st_convex_hull() |>
  st_sf(geometry = _)

envelope <- st_buffer(envelope_hull, dist = 1000) |>  # small margin only, hull already closes gaps
  st_union() |>
  st_sf(geometry = _)
mapview::mapview(envelope)

terra::freq(ecomap_stx)  # look for a 0 category, and check if NA cells exist separately
water_mask <- terra::ifel(nlcd_stx == 11, 1, NA)
terra::plot(water_mask)

exclude_poly <- terra::as.polygons(water_mask, dissolve = TRUE) |>
  st_as_sf() |>
  st_make_valid() |>
  st_transform(st_crs(envelope))

envelope_available <- st_difference(envelope, st_union(exclude_poly))

mapview::mapview(envelope_available)

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
# Get 95% UD and retain polygons up to 90% of total area
# --------------------------------------------------------------------

ud_list <- readRDS("output/ctmm_occurrence_uds.rds")

get_ud_polygons <- function(ud, levels) {
  purrr::map_dfr(levels, function(lvl) {
    x <- as.sf(ud, level.UD = lvl)
    x %>% sf::st_cast("POLYGON") %>%
      dplyr::mutate(level = lvl, area_km2 = as.numeric(sf::st_area(geometry)) / 1e6)
  })
}

ud95_polygons <- purrr::imap(ud_list, ~get_ud_polygons(.x, levels = 0.95) %>% mutate(Deployment_ID = .y))
#str(ud90_polygons)

ud95_polygons_sf <- purrr::map_dfr(ud95_polygons, identity)

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

ud_retained <- retain_ud_components(
  ud95_polygons_sf,
  component_retention = 0.90
)

# Check 1 individual
ud_retained %>% filter(Deployment_ID == "EO19M_1", level == 0.95) %>%
  select(
    area_km2,
    component_rank,
    cumulative_prop,
    component_retention,
    retained_n_components,
    retained_prop
  )
ud_retained

ud_retained_multipolygon <- ud_retained %>%
  group_by(Deployment_ID, level) %>%
  summarise(
    geometry = st_union(geometry),
    retained_n_components = n(),
    retained_area_km2 = sum(area_km2),
    total_ud_area_km2 = first(total_ud_area_km2),
    original_n_components = first(original_n_components),
    retained_prop = sum(area_km2) / first(total_ud_area_km2),
    .groups = "drop"
  )

summary(ud_retained_multipolygon)

top_fragmented <- ud_retained_multipolygon %>%
  arrange(desc(retained_n_components)) %>%
  slice_head(n = 4)

ud_retained_sf <- st_as_sf(
  ud_retained_multipolygon,
  sf_column_name = "geometry",
  crs = st_crs(ud95_polygons_sf)
)

class(ud_retained_sf)
st_geometry_type(ud_retained_sf)

mapview::mapview(ud_retained_sf[ud_retained_sf$Deployment_ID == "EO34M_2",])
mapview::mapview(ud95_polygons_sf[ud95_polygons_sf$Deployment_ID == "EO34M_2",])

# --------------------------------------------------------------------
# Landscape metrics within UD + buffers
# --------------------------------------------------------------------
class_metrics <- c("lsm_c_pland", "lsm_c_area_mn", "lsm_c_shape_mn", "lsm_c_cai_mn", "lsm_c_enn_mn", "lsm_c_clumpy")

buffer_sizes <- c(0, 250, 500, 1000, 2000, 3000, 5000)

thornscrub_binary <- rast("output/south_texas_thornscrub_binary.tif")
terra::plot(thornscrub_binary)

ud_retained_sf <- st_transform(ud_retained_sf, crs = crs(thornscrub_binary))

# Function to crop and mask thornscrub raster by ud
get_landscape <- function(r, ud_area) {  
  r_crop <- terra::crop(r, terra::vect(ud_area))
  r_mask <- terra::mask(r_crop, terra::vect(ud_area))
  r_mask
}

calculate_ud_metrics <- function(ud, landscape, buffer_sizes, metrics) {
  purrr::map_dfr(
    buffer_sizes,
    function(buffer_dist) {
      ud_buffer <- st_buffer(
        ud,
        dist = buffer_dist
      )
      
      r <- get_landscape(
        landscape,
        ud_buffer
      )
      
      calculate_lsm(
        r,
        what = metrics
      ) %>%
        filter(class == 1) %>%
        mutate(
          Deployment_ID = ud$Deployment_ID,
          UD_level = ud$level,
          buffer_m = buffer_dist
        )
    }
  )
}

metrics <- purrr::map_dfr(
  seq_len(nrow(ud_retained_sf)),
  ~ calculate_ud_metrics(
      ud = ud_retained_sf[.x, ],
      landscape = thornscrub_binary,
      buffer_sizes = buffer_sizes,
      metrics = class_metrics
    )
)

ggplot() + geom_sf(data = ud_retained_sf, aes(color = Deployment_ID), fill = NA, color = "black", linewidth = 0.8)
mapview::mapview(ud_retained_sf, zcol = "Deployment_ID")

plots <- purrr::map(
  unique(ud_retained_sf$Deployment_ID),
  function(id) {

    ud_i <- ud_retained_sf |>
      filter(Deployment_ID == id)

    bb <- st_bbox(
      st_buffer(ud_i, dist = 1000)
    )

    ggplot() +
      geom_spatraster(
        data = thornscrub_binary
      ) +
      geom_sf(
        data = ud_i,
        fill = "black",
        alpha = 0.25,
        color = "black",
        linewidth = 0.4
      ) +
      coord_sf(
        xlim = c(bb["xmin"], bb["xmax"]),
        ylim = c(bb["ymin"], bb["ymax"])
      ) +
      labs(title = id) +
      theme_classic()
  }
)

plots[[4]]

wrap_plots(plots[1:15])
wrap_plots(plots[16:31])

# --------------------------------------------------------------------
# Landscape patch metrics within UD + buffers
# --------------------------------------------------------------------

patch_metrics <- c("lsm_p_area")

calculate_patch_ud_metrics <- function(ud, landscape, buffer_sizes, metrics) {
  purrr::map_dfr(
    buffer_sizes,
    function(buffer_dist) {
      ud_buffer <- st_buffer(
        ud,
        dist = buffer_dist
      )
     
      r <- get_landscape(
        landscape,
        ud_buffer
      )
      
      calculate_lsm(
        r,
        what = metrics
      ) %>%
        mutate(
          Deployment_ID = ud$Deployment_ID,
          UD_level = ud$level,
          buffer_m = buffer_dist
        )
    }
  )
}

metrics_patch <- purrr::map_dfr(
  seq_len(nrow(ud_retained_sf)),
  ~ calculate_patch_ud_metrics(
      ud = ud_retained_sf[.x, ],
      landscape = thornscrub_binary,
      buffer_sizes = buffer_sizes,
      metrics = patch_metrics
    )
)

