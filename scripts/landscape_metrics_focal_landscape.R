# =====================================================================
# FULL LANDSCAPE METRICS PIPELINE -- 164m radius
# =====================================================================
library(terra); library(dplyr); library(sf); library(tidyr); library(ggplot2)
library(tidyterra); library(patchwork)

# -----------------------------------------------------------------
# Setup
# -----------------------------------------------------------------
radius <- 164
cell_size <- res(thornscrub_cropped_final)[1]

w <- focalMat(thornscrub_cropped_final, d = radius, type = "circle")
w_binary <- ifelse(w > 0, 1, 0)
window_area_ha <- (sum(w_binary) * cell_size^2) / 10000
window_area_100ha <- window_area_ha / 100

base_rast <- thornscrub_cropped_final

# -----------------------------------------------------------------
# Patch IDs, core cells, edge layer
# -----------------------------------------------------------------
thornscrub_patches <- base_rast
thornscrub_patches[thornscrub_patches != 1] <- NA
patch_id_rast <- patches(thornscrub_patches, directions = 8)
writeRaster(patch_id_rast, "output/objects/patch_id_rast_164m.tif", overwrite = TRUE)

rook_kernel <- matrix(c(0,1,0, 1,0,1, 0,1,0), nrow = 3)
neighbor_sum <- focal(base_rast, w = rook_kernel, fun = sum, na.rm = TRUE)
core_cells <- ifel(base_rast == 1 & neighbor_sum == 4, 1, 0)
names(core_cells) <- "is_core"

m <- as.matrix(base_rast, wide = TRUE)
nr <- nrow(m); nc <- ncol(m)
horiz_diff <- matrix(0, nr, nc); horiz_diff[, 1:(nc-1)] <- (m[, 1:(nc-1)] != m[, 2:nc]) * 1
vert_diff  <- matrix(0, nr, nc); vert_diff[1:(nr-1), ]  <- (m[1:(nr-1), ] != m[2:nr, ]) * 1
horiz_diff[is.na(horiz_diff)] <- 0; vert_diff[is.na(vert_diff)] <- 0
total_edges <- rast(horiz_diff + vert_diff, extent = ext(base_rast), crs = crs(base_rast))
edges_masked <- total_edges * (base_rast == 1)

patch_df <- as.data.frame(patch_id_rast, xy = TRUE, na.rm = TRUE)
names(patch_df)[3] <- "patch_id"

# -----------------------------------------------------------------
# PLAND, ED
# -----------------------------------------------------------------
tictoc::tic("pland + ed")
pland_surface <- focal(base_rast, w = w_binary, fun = mean, na.rm = TRUE) * 100

edge_count_sum <- focal(total_edges, w = w_binary, fun = sum, na.rm = TRUE)
ed_surface <- (edge_count_sum * cell_size) / window_area_ha
tictoc::toc()
writeRaster(c(pland_surface, ed_surface), "output/objects/surfaces_pland_ed_164m.tif", overwrite = TRUE)

# -----------------------------------------------------------------
# NP, PD
# -----------------------------------------------------------------
tictoc::tic("np + pd")
np_surface <- focal(patch_id_rast, w = w_binary, fun = function(x, ...) {
  length(unique(na.omit(x)))
}, na.rm = TRUE)
pd_surface <- np_surface / window_area_100ha
tictoc::toc()
writeRaster(c(np_surface, pd_surface), "output/objects/surfaces_np_pd_164m.tif", overwrite = TRUE)

# -----------------------------------------------------------------
# AREA_MN, LPI (patch-weighted, clipped to in-window cell counts)
# -----------------------------------------------------------------
tictoc::tic("area_mn + lpi")
area_mn_surface <- focal(patch_id_rast, w = w_binary, fun = function(x, ...) {
  ids <- na.omit(x)
  if (length(ids) == 0) return(0)
  areas_ha <- table(ids) * cell_size^2 / 10000
  mean(areas_ha)
}, na.rm = TRUE)

lpi_surface <- focal(patch_id_rast, w = w_binary, fun = function(x, ...) {
  ids <- na.omit(x)
  if (length(ids) == 0) return(0)
  areas_ha <- table(ids) * cell_size^2 / 10000
  win_area_ha <- length(x) * cell_size^2 / 10000
  max(areas_ha) / win_area_ha * 100
}, na.rm = TRUE)
tictoc::toc()
writeRaster(c(area_mn_surface, lpi_surface), "output/objects/surfaces_area_lpi_164m.tif", overwrite = TRUE)

# -----------------------------------------------------------------
# SHAPE_MN (windowed perimeter/area -- validated fix)
# -----------------------------------------------------------------
tictoc::tic("shape_mn")
thorn_cell_count <- focal(base_rast, w = w_binary, fun = sum, na.rm = TRUE)
edge_sum_window <- focal(edges_masked, w = w_binary, fun = sum, na.rm = TRUE)
area_m2 <- thorn_cell_count * cell_size^2
perimeter_m <- edge_sum_window * cell_size
shape_mn_surface <- ifel(area_m2 > 0, 0.25 * perimeter_m / sqrt(area_m2), 0)
tictoc::toc()
writeRaster(shape_mn_surface, "output/objects/surface_shape_mn_164m.tif", overwrite = TRUE)

# -----------------------------------------------------------------
# CAI_MN (patch-weighted, via encoded patch_id/core layer -- validated v3)
# -----------------------------------------------------------------
tictoc::tic("cai_mn")
encoded <- patch_id_rast * 2 + core_cells

cai_mn_surface <- focal(encoded, w = w_binary, fun = function(x, ...) {
  x <- na.omit(x)
  if (length(x) == 0) return(0)
  pid <- floor(x / 2); is_core <- x %% 2
  mean(tapply(is_core, pid, mean)) * 100
}, na.rm = TRUE)
tictoc::toc()
writeRaster(cai_mn_surface, "output/objects/surface_cai_mn_164m.tif", overwrite = TRUE)

# -----------------------------------------------------------------
# INTERIOR DEPTH
# -----------------------------------------------------------------
tictoc::tic("interior_depth")
edge_targets <- base_rast
edge_targets[edge_targets == 1] <- NA
edge_targets[edge_targets == 0] <- 1
dist_to_edge <- distance(edge_targets)
dist_to_edge[base_rast == 0] <- NA

interior_depth_surface <- focal(dist_to_edge, w = w_binary, fun = mean, na.rm = TRUE)
tictoc::toc()
writeRaster(interior_depth_surface, "output/objects/surface_interior_depth_164m.tif", overwrite = TRUE)

# -----------------------------------------------------------------
# PLADJ
# -----------------------------------------------------------------
tictoc::tic("pladj")
h_left <- m[, 1:(nc-1)]; h_right <- m[, 2:nc]
like_h <- matrix(0, nr, nc); total_h <- matrix(0, nr, nc)
like_h[, 1:(nc-1)] <- ifelse(h_left==1 & h_right==1, 1, 0)
like_h[, 2:nc]     <- like_h[, 2:nc] + ifelse(h_left==1 & h_right==1, 1, 0)
total_h[, 1:(nc-1)] <- ifelse(h_left==1 | h_right==1, 1, 0)
total_h[, 2:nc]     <- total_h[, 2:nc] + ifelse(h_left==1 | h_right==1, 1, 0)

v_top <- m[1:(nr-1), ]; v_bottom <- m[2:nr, ]
like_v <- matrix(0, nr, nc); total_v <- matrix(0, nr, nc)
like_v[1:(nr-1), ] <- ifelse(v_top==1 & v_bottom==1, 1, 0)
like_v[2:nr, ]     <- like_v[2:nr, ] + ifelse(v_top==1 & v_bottom==1, 1, 0)
total_v[1:(nr-1), ] <- ifelse(v_top==1 | v_bottom==1, 1, 0)
total_v[2:nr, ]     <- total_v[2:nr, ] + ifelse(v_top==1 | v_bottom==1, 1, 0)

like_rast <- rast(like_h + like_v, extent = ext(base_rast), crs = crs(base_rast))
adj_rast  <- rast(total_h + total_v, extent = ext(base_rast), crs = crs(base_rast))

like_sum <- focal(like_rast, w = w_binary, fun = sum, na.rm = TRUE)
adj_sum  <- focal(adj_rast,  w = w_binary, fun = sum, na.rm = TRUE)

pladj_surface <- ifel(adj_sum == 0, 0, (like_sum / adj_sum) * 100)
tictoc::toc()
writeRaster(pladj_surface, "output/objects/surface_pladj_164m.tif", overwrite = TRUE)

# -----------------------------------------------------------------
# ED_PER_PLAND (derived; NA where pland <= 5%)
# -----------------------------------------------------------------
ed_per_pland_surface <- ifel(pland_surface > 5, ed_surface / (pland_surface / 100), NA)
writeRaster(ed_per_pland_surface, "output/objects/surface_ed_per_pland_164m.tif", overwrite = TRUE)

# =====================================================================
# Stack, visualize, extract, and screen
# =====================================================================
all_surfaces <- c(pland_surface, ed_surface, np_surface, pd_surface,
                   area_mn_surface, lpi_surface, shape_mn_surface, cai_mn_surface,
                   interior_depth_surface, pladj_surface, ed_per_pland_surface)

names(all_surfaces) <- c("pland", "ed", "np", "pd", "area_mn", "lpi", "shape_mn",
                          "cai_mn", "interior_depth", "pladj", "ed_per_pland")

writeRaster(all_surfaces, "output/objects/all_surfaces_164m.tif", overwrite = TRUE)

# visualize all at once
plot_list <- lapply(names(all_surfaces), function(nm) {
  ggplot() +
    geom_spatraster(data = all_surfaces[[nm]]) +
    scale_fill_viridis_c(na.value = "white", name = nm) +
    labs(title = nm) +
    theme_minimal(base_size = 9) +
    theme(axis.title = element_blank(), axis.text = element_blank(),
          legend.key.size = unit(0.3, "cm"))
})
wrap_plots(plot_list, ncol = 3)

# extract at all real points
tictoc::tic("full extraction")
all_metric_values <- terra::extract(all_surfaces, vect(pts_sf))
tictoc::toc()
all_metric_values$point_id <- pts_sf$point_id

# apply fill-zero convention for metrics where 0 = "no thornscrub present"
fill_zero_metrics <- c("cai_mn", "interior_depth")
all_metric_values <- all_metric_values %>%
  mutate(across(all_of(fill_zero_metrics), ~ coalesce(., 0)))

saveRDS(all_metric_values, "output/objects/all_metric_values_164m.rds")

# summary, NA/CV screen, correlation matrix
all_metric_values %>% select(-ID, -point_id) %>% summary()

all_metric_values %>% select(-ID, -point_id) %>%
  summarise(across(everything(), list(
    pct_na = ~mean(is.na(.)),
    cv = ~sd(., na.rm = TRUE) / mean(., na.rm = TRUE)
  ))) %>%
  pivot_longer(everything(), names_sep = "_(?=[^_]+$)", names_to = c("metric", "stat")) %>%
  pivot_wider(names_from = stat, values_from = value)

cor_matrix_full <- all_metric_values %>% select(-ID, -point_id) %>%
  cor(use = "pairwise.complete.obs")
print(round(cor_matrix_full, 2))

# =====================================================================
# Run lsm on full dataset to compare
# =====================================================================
class_metrics <- c("lsm_c_pland", "lsm_c_area_mn", "lsm_c_shape_mn", "lsm_c_cai_mn",
                    "lsm_c_ed", "lsm_c_lpi", "lsm_c_pd", "lsm_c_np")
batch_size <- 1500

run_batched_checkpointed <- function(radius, pts, raster, metrics, batch_size, checkpoint_every = 50,
                                      checkpoint_path = "output/objects/run_batched_164_checkpoint.rds") {
  n_batches <- ceiling(nrow(pts) / batch_size)
  out <- vector("list", n_batches)

  for (b in 1:n_batches) {
    idx <- ((b - 1) * batch_size + 1):min(b * batch_size, nrow(pts))
    batch_result <- sample_lsm(
      landscape = raster, y = pts[idx, ],
      size = radius, what = metrics, shape = "circle", return_raster = FALSE
    ) %>% filter(class == 1)

    batch_result$point_id <- idx[batch_result$plot_id]
    out[[b]] <- batch_result

    gc()
    if (b %% checkpoint_every == 0 || b == n_batches) {
      saveRDS(list(completed_through = b, data = bind_rows(out[1:b])), checkpoint_path)
      cat("  Checkpoint saved at batch", b, "of", n_batches, "\n")
    }
    cat("  Radius", radius, "- batch", b, "of", n_batches, "done at", format(Sys.time()), "\n")
  }
  bind_rows(out) %>% mutate(buffer_radius = radius)
}

full_run <- run_batched_checkpointed(164, pts_sf, thornscrub_cropped_final, class_metrics, batch_size)
