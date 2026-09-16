library(terra); library(sf); library(ggplot2)
library(tidyterra); library(patchwork)
library(amt); library(dplyr); library(lubridate); 
library(momentuHMM); library(landscapemetrics)
library(tidyverse)

# # -------------------------------------------------------------------
# # Used vs available points from AMT for local only
# #--------------------------------------------------------------------

# # Build the amt track 
# trk <- local_only %>%
#   arrange(Deployment_ID.x, timestamp) %>%
#   nest(data = -Deployment_ID.x) %>%
#   mutate(trk = lapply(data, function(d) {
#     make_track(d, .x = x_proj, .y = y_proj, .t = timestamp, crs = target_crs)
#   }))

# issf_local <- trk %>%
#   mutate(stp = lapply(trk, function(x) {
#     x %>%
#       track_resample(rate = hours(1), tolerance = minutes(10)) %>%
#       filter_min_n_burst(min_n = 3) %>%
#       steps_by_burst() %>%
#       random_steps(n_control = 10) %>%
#       remove_incomplete_strata() %>%
#       mutate(log_sl_ = log(sl_ + 0.1), cos_ta_ = cos(ta_))
#   })) %>%
#   select(-data, -trk) %>%
#   unnest(cols = stp)

# nrow(issf_local)

# # how many bursts survive per individual, and how long are they?
# burst_summary <- issf_local %>% 
#   distinct(Deployment_ID.x, burst_) %>%   # or whatever the burst column is named post-steps_by_burst()
#   count(Deployment_ID.x, name = "n_bursts")

# step_counts <- issf_local %>% filter(case_ == TRUE) %>% count(Deployment_ID.x)
# print(step_counts, n = 35)

# # how many individuals/bursts triggered this?
# step_counts_raw <- local_only %>% count(Deployment_ID.x, name = "n_raw_local")
# step_counts_final <- issf_local %>% filter(case_ == TRUE) %>% count(Deployment_ID.x, name = "n_final_steps")

# retained <- step_counts_raw %>% 
#   left_join(step_counts_final, by = c("Deployment_ID.x")) %>%
#   mutate(n_final_steps = coalesce(n_final_steps, 0),
#          pct_retained = n_final_steps / n_raw_local)
# retained

# # -------------------------------------------------------------------
# # Landscape metrics by points
# #--------------------------------------------------------------------

# library(landscapemetrics); library(terra); library(sf); library(dplyr); library(purrr)

# saveRDS(issf_local, "output/objects/issf_points_local.rds")
# #saveRDS(subsample_sf, "output/objects/issf_subsample_local.rds")

# thornscrub_raster_full <- rast("output/thornscrub_raster_full.tif")

# #issf_local <- issf_local %>% rename(id = Deployment_ID.x)

# # Points to extract covariates at -- typically the step ENDPOINT is where selection is evaluated
# pts_sf <- issf_local %>% 
#   st_as_sf(coords = c("x2_", "y2_"), crs = target_crs, remove = FALSE) %>%
#   mutate(point_id = row_number()) 

# thornscrub_raster_full <- rast("output/thornscrub_raster_full.tif")

# # ---------------------------------------------------------------
# # Crop tightly to points + max radius, trim, force into memory
# # ---------------------------------------------------------------
# buffer_ext <- st_buffer(st_as_sfc(st_bbox(pts_sf)), 4000)
# thornscrub_cropped <- crop(thornscrub_raster_full, vect(buffer_ext))
# thornscrub_cropped <- trim(thornscrub_cropped)

# # force fully into memory (avoids repeated disk reads per point)
# values(thornscrub_cropped) <- values(thornscrub_cropped)
# inMemory(thornscrub_cropped)  

# # ---------------------------------------------------------------
# # Batched extraction 
# # ---------------------------------------------------------------
# class_metrics <- c("lsm_c_pland", "lsm_c_area_mn", "lsm_c_shape_mn", "lsm_c_cai_mn",
#                     "lsm_c_enn_mn", "lsm_c_clumpy", "lsm_c_ed", "lsm_c_cohesion", "lsm_c_lpi", "lsm_c_gyrate_mn", "lsm_c_pd", "lsm_c_np")
# batch_size <- 500

# run_batched <- function(radius, pts, raster, metrics, batch_size) {
#   n_batches <- ceiling(nrow(pts) / batch_size)
#   out <- vector("list", n_batches)

#   for (b in 1:n_batches) {
#     idx <- ((b - 1) * batch_size + 1):min(b * batch_size, nrow(pts))
#     batch_result <- sample_lsm(
#       landscape = raster, y = pts[idx, ],
#       size = radius, what = metrics, shape = "circle", return_raster = FALSE
#     ) %>% filter(class == 1)

#     batch_result$point_id <- idx[batch_result$plot_id]
#     out[[b]] <- batch_result

#     gc()
#     cat("  Radius", radius, "- batch", b, "of", n_batches, "done at", format(Sys.time()), "\n")
#   }
#   bind_rows(out) %>% mutate(buffer_radius = radius)
# }

# # Get the precise value to base the radius on
# issf_local %>% filter(case_ == TRUE) %>% select(sl_) %>% 
#   summarise(quantile_95 = quantile(sl_, probs = 0.95, na.rm = TRUE), median_sl = median(sl_, na.rm = TRUE), mean_sl = mean(sl_, na.rm = TRUE))
# #  summarise(median_sl = median(sl_, na.rm = TRUE), mean_sl = mean(sl_, na.rm = TRUE))

# # radius = half of whichever you choose to anchor on
# final_radius <- round(207/ 2)  # or mean_sl / 2, your call
# landscape_metrics_all_extraction <- run_batched(final_radius, pts_sf, thornscrub_cropped, class_metrics, batch_size)
# saveRDS(landscape_metrics_all_extraction, "output/objects/landscape_metrics_local_95.rds") # added lpi, gyrate, pd, np

# #
# final_radius <- round(67/ 2)  # median
# landscape_metrics_all_extraction <- run_batched(final_radius, pts_sf, thornscrub_cropped, class_metrics, batch_size)
# saveRDS(landscape_metrics_all_extraction, "output/objects/landscape_metrics_local_median.rds") # added lpi, gyrate, pd, np


# traveling_median <- readRDS(landscape_metrics_all_extraction, "output/objects/final_landscape_metrics3.rds") # all metrics but at 95% of step length
# traveling_p95
# local_median

# local_p95 <- readRDS("output/objects/landscape_metrics_local_95.rds")
# local_median <- readRDS("output/objects/landscape_metrics_local_median.rds")

# local_p95_wide <- local_p95 %>%
#   select(point_id, metric, value) %>%
#   pivot_wider(names_from = metric, values_from = value)

# local_median_wide <- local_median %>%
#   select(point_id, metric, value) %>%
#   pivot_wider(names_from = metric, values_from = value)

# saveRDS(local_p95_wide, "output/objects/local_p95_wide.rds") # added lpi, gyrate, pd, np
# saveRDS(local_median_wide, "output/objects/local_median_wide.rds") # added lpi, gyrate, pd, np





# traveling_p95 <- readRDS("output/objects/final_landscape_metrics3.rds") # all metrics but at 95% of step length

# traveling_p95_wide <- traveling_p95 %>%
#   select(point_id, metric, value) %>%
#   pivot_wider(names_from = metric, values_from = value)
# saveRDS(traveling_p95_wide, "output/objects/traveling_p95_wide.rds")

# traveling_median1 <- readRDS("output/objects/final_landscape_metrics.rds") # added lpi, gyrate, pd, np
# traveling_median2 <- readRDS("output/objects/final_landscape_metrics2.rds") # added lpi, gyrate, pd, np

# traveling_median <- bind_rows(traveling_median1, traveling_median2)
# traveling_median_wide <- traveling_median %>%
#   select(point_id, metric, value) %>%
#   pivot_wider(names_from = metric, values_from = value)
# saveRDS(traveling_median_wide, "output/objects/traveling_median_wide.rds")


# colnames(traveling_median_wide)

# landscape_metrics_all_extraction <- readRDS("output/objects/final_landscape_metrics2.rds")
# # Confirm class filter held
# table(landscape_metrics_all_extraction$class)  # should show only 1

# summary(traveling_median_wide)
# summary(traveling_p95_wide)
# summary(local_median_wide)
# summary(local_p95_wide)

















# # -------------------------------------------------------------------
# # Pivot to wide: one row per point, one column per metric
# # -------------------------------------------------------------------
# metrics_wide <- landscape_metrics_all_extraction %>%
#   select(point_id, metric, value) %>%
#   pivot_wider(names_from = metric, values_from = value)

# # pland: missing = genuinely 0% thornscrub in buffer (no row extracted at all)
# # other metrics: missing = no patch present, genuinely undefined -- leave as NA
# all_points <- tibble(point_id = unique(pts_sf$point_id))




## GOOD CODE TO START HERE #################################################################################

# reattach state labels
data_adj$state <- viterbi(m3_cosinor)
data_adj$state_label <- recode(as.character(data_adj$state),
                                "1" = "local_movement", "2" = "encamped", "3" = "traveling")

# reproject lon/lat into meters-scale target CRS
data_adj <- data_adj %>%
  left_join(hmm_input %>% st_drop_geometry() %>% select(BurstID, timestamp, longitude, latitude),
            by = c("ID" = "BurstID", "timestamp" = "timestamp"))

data_adj_sf <- data_adj %>% st_as_sf(coords = c("longitude","latitude"), crs = 4326, remove=FALSE) %>%
  st_transform(target_crs)
coords <- st_coordinates(data_adj_sf)
data_adj$x_proj <- coords[,1]
data_adj$y_proj <- coords[,2]

trk_all <- data_adj %>%
  arrange(Deployment_ID, timestamp) %>%
  nest(data = -Deployment_ID) %>%
  mutate(trk = lapply(data, function(d) {
    make_track(d, .x = x_proj, .y = y_proj, .t = timestamp, crs = target_crs)
  }))

issf_all_states <- trk_all %>%
  mutate(stp = lapply(trk, function(x) {
    x %>%
      track_resample(rate = hours(1), tolerance = minutes(10)) %>%
      filter_min_n_burst(min_n = 3) %>%
      steps_by_burst() %>%
      random_steps(n_control = 10) %>%
      remove_incomplete_strata() %>%
      mutate(log_sl_ = log(sl_ + 0.1), cos_ta_ = cos(ta_))
  })) %>%
  select(-data, -trk) %>%
  unnest(cols = stp)

issf_all_states <- issf_all_states %>% mutate(step_id_global = paste0(Deployment_ID, "_", step_id_))
n_distinct(issf_all_states$step_id_global)
issf_all_states %>% count(step_id_global) %>% filter(n > 11)  # sanity check no runaway duplication (11 = 1 used + 10 available)

used_points_with_state <- issf_all_states %>%
  filter(case_ == TRUE) %>%
  left_join(
    data_adj %>% select(Deployment_ID, timestamp, state_label),
    by = c("Deployment_ID" = "Deployment_ID", "t2_" = "timestamp")
  )
# check match quality
nrow(used_points_with_state)
sum(is.na(used_points_with_state$state_label))

state_by_stratum <- used_points_with_state %>% select(step_id_global, state = state_label)

issf_all_states <- issf_all_states %>%
  select(-any_of("state")) %>%
  left_join(state_by_stratum, by = "step_id_global")

sum(is.na(issf_all_states$state))
table(issf_all_states$state[issf_all_states$case_])  # should now sum to ~130,488

issf_all_states <- issf_all_states %>%
  mutate(step_id_ = step_id_global) %>%
  select(-step_id_global)

issf_all_states %>% filter(case_ == TRUE) %>% 
  summarise(median_sl = median(sl_, na.rm=TRUE), mean_sl = mean(sl_, na.rm=TRUE), quantile_95 = quantile(sl_, probs = 0.95, na.rm = TRUE))

saveRDS(issf_all_states, "output/objects/issf_all_states.rds")

radii_test <- c(mean_half = round(164/2), mean_full = round(164), 
                 p95_half = round(610/2), p95_full = round(610))
radii_test

# Points to extract covariates at -- typically the step ENDPOINT is where selection is evaluated
pts_sf <- issf_all_states %>%
  st_as_sf(coords = c("x2_", "y2_"), crs = target_crs, remove = FALSE) %>%
  mutate(point_id = row_number()) 

# Need to include water for the radii
thornscrub_raster_full <- rast("output/thornscrub_raster_full.tif")

# ===============================================
# Crop landscape around points with buffer
# ===============================================
# makes downstream calculations faster

# margin around points -- covers max radius plus some buffer 
buffer_ext <- st_buffer(st_as_sfc(st_bbox(pts_sf)), 610 * 5)  # 5x radius margin
thornscrub_cropped_final <- crop(thornscrub_raster_full, vect(buffer_ext))
thornscrub_cropped_final <- trim(thornscrub_cropped_final)  # drop any leftover empty NA border

# # force into memory, same discipline as earlier in the session
# values(thornscrub_cropped_final) <- values(thornscrub_cropped_final)

# sanity check the new extent covers all your points comfortably
ext(thornscrub_cropped_final)
st_bbox(pts_sf)
dim(thornscrub_cropped_final)  

# quick base plot: raster + points
plot(thornscrub_cropped_final, main = "Thornscrub with points overlaid")
plot(st_geometry(pts_sf), add = TRUE, pch = 16, cex = 0.1, col = rgb(1,0,0,0.1))

# ===============================================
# Create pland raster using focal (max radius; 610 m; p95 of step length)
# ===============================================
w <- focalMat(thornscrub_cropped_final, d = 610, type = "circle")
w_binary <- ifelse(w > 0, 1, 0)
pland_surface <- focal(thornscrub_cropped_final, w = w_binary, fun = mean, na.rm = TRUE) * 100

pts_vect <- vect(pts_sf)
pland_values <- terra::extract(pland_surface, pts_vect)
sum(is.na(pland_values$focal_mean))
summary(pland_values$focal_mean)
plot(pland_surface)





# -------------------------------------------------------------------
# 2. Validate against sample_lsm() at the same test points
# -------------------------------------------------------------------
direct_check_pland <- sample_lsm(thornscrub_cropped_final, y = test_pts_sf, size = 610,
                                  what = "lsm_c_pland", shape = "circle", return_raster = FALSE) %>%
  filter(class == 1) %>%
  select(plot_id, value) %>%
  rename(point_id = plot_id, pland = value)

focal_check_pland <- terra::extract(pland_surface, vect(test_pts_sf))
names(focal_check_pland) <- c("ID", "pland_focal")

comparison_pland <- direct_check_pland %>%
  left_join(focal_check_pland, by = c("point_id" = "ID"))

print(comparison_pland)

# -------------------------------------------------------------------
# 3. Panel plots -- same pattern as the others
# -------------------------------------------------------------------
plot_pts <- test_pts_sf[1:6, ]
buffer_polys <- st_buffer(plot_pts, dist = 610)
crop_df <- as.data.frame(thornscrub_cropped_final, xy = TRUE)
names(crop_df)[3] <- "thornscrub"

plot_list_pland <- lapply(1:6, function(i) {
  pt <- plot_pts[i, ]
  buf <- buffer_polys[i, ]
  val <- terra::extract(pland_surface, vect(pt))[,2]
  direct_val <- comparison_pland$pland[comparison_pland$point_id == i]
  
  ggplot() +
    geom_raster(data = crop_df, aes(x = x, y = y, fill = factor(thornscrub))) +
    scale_fill_manual(values = c("0" = "grey85", "1" = "#4daf4a"), guide = "none") +
    geom_sf(data = buf, fill = NA, color = "red", linewidth = 0.8, inherit.aes = FALSE) +
    geom_sf(data = pt, color = "black", size = 2, inherit.aes = FALSE) +
    coord_sf(datum = NA) +
    labs(title = paste0("Pt ", i, ": pland=", round(val,1), " (lsm=", round(direct_val,1), ")")) +
    theme_minimal(base_size = 9) +
    theme(axis.title = element_blank(), axis.text = element_blank())
})

wrap_plots(plot_list_pland, ncol = 3)

# ===============================================
# Test and validation of raster wide metrics using focal approach
# ===============================================

# -------------------------------------------------------------------
# Small test crop -- pick an area with a mix of patch sizes
# -------------------------------------------------------------------
# crop to a manageable window, e.g. 3km x 3km, centered on a point you know has varied thornscrub
test_center <- st_centroid(st_as_sfc(st_bbox(pts_sf)))  
test_ext <- st_buffer(test_center, 1500) %>% st_bbox() %>% ext()

test_crop <- crop(thornscrub_raster_full, test_ext)
plot(test_crop, main = "Test crop: thornscrub binary")

# -------------------------------------------------------------------
# Identify patches ONCE across the test landscape
# -------------------------------------------------------------------
thornscrub_patches <- test_crop
thornscrub_patches[thornscrub_patches != 1] <- NA
patch_id_rast <- patches(thornscrub_patches, directions = 8)

plot(test_crop, main = "Thornscrub (binary)")
plot(patch_id_rast, main = "Patch IDs")  # visually confirm distinct patches got distinct IDs

# -------------------------------------------------------------------
# np and pd via focal window over patch IDs
# -------------------------------------------------------------------
w <- focalMat(test_crop, d = 610, type = "circle")
w_binary <- ifelse(w > 0, 1, 0)

tictoc::tic()
np_test <- focal(patch_id_rast, w = w_binary, fun = function(x, ...) {
  length(unique(na.omit(x)))
}, na.rm = TRUE)
tictoc::toc()

window_area_ha <- (sum(w_binary) * res(test_crop)[1]^2) / 10000
pd_test <- np_test / window_area_ha

par(mfrow = c(1,3))
plot(test_crop, main = "Thornscrub")
plot(np_test, main = "Number of patches (np)")
plot(pd_test, main = "Patch density (pd)")
par(mfrow = c(1,1))

# -------------------------------------------------------------------
# Validate against direct sample_lsm() at test points
# -------------------------------------------------------------------
set.seed(1)
test_pts <- spatSample(test_crop, size = 8, as.points = TRUE, na.rm = FALSE)
test_pts_sf <- st_as_sf(test_pts)
test_pts_sf$point_id <- 1:nrow(test_pts_sf)

direct_check <- sample_lsm(test_crop, y = test_pts_sf, size = 610,
                            what = c("lsm_c_np", "lsm_c_pd"),
                            shape = "circle", return_raster = FALSE) %>%
  filter(class == 1)

window_area_100ha <- (sum(w_binary) * res(test_crop)[1]^2) / 10000 / 100  # per 100 hectares
pd_test_corrected <- np_test / window_area_100ha

# recheck
focal_check <- terra::extract(c(np_test, pd_test_corrected), vect(test_pts_sf))
names(focal_check) <- c("ID", "np_focal", "pd_focal_corrected")

direct_check_wide <- direct_check %>%
  select(plot_id, metric, value) %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  rename(point_id = plot_id)

comparison <- direct_check_wide %>% left_join(focal_check, by = c("point_id"="ID"))
print(comparison)

# -------------------------------------------------------------------
# Visualize 6 test points
# -------------------------------------------------------------------

# use your existing 8 test points (or pick 6 specifically)
plot_pts <- test_pts_sf[1:6, ]

# build actual circular buffer polygons at the same radius used for np/pd
buffer_polys <- st_buffer(plot_pts, dist = 610)

# convert raster to df once for consistent plotting across panels
crop_df <- as.data.frame(test_crop, xy = TRUE)
names(crop_df)[3] <- "thornscrub"

plot_list <- lapply(1:6, function(i) {
  pt <- plot_pts[i, ]
  buf <- buffer_polys[i, ]
  
  # pull the np/pd focal values at this point for the panel title
  np_val <- terra::extract(np_test, vect(pt))[,2]
  pd_val <- terra::extract(pd_test_corrected, vect(pt))[,2]
  
  ggplot() +
    geom_raster(data = crop_df, aes(x = x, y = y, fill = factor(thornscrub))) +
    scale_fill_manual(values = c("0" = "grey85", "1" = "#4daf4a"), guide = "none") +
    geom_sf(data = buf, fill = NA, color = "red", linewidth = 0.8, inherit.aes = FALSE) +
    geom_sf(data = pt, color = "black", size = 2, inherit.aes = FALSE) +
    coord_sf(datum = NA) +
    labs(title = paste0("Point ", i, ": np=", round(np_val,1), ", pd=", round(pd_val,1))) +
    theme_minimal(base_size = 9) +
    theme(axis.title = element_blank(), axis.text = element_blank()) #+
    #labs(title = paste0("Pt ", i, ": np=", round(np_val,1), " (sample_lsm=", 
     #                 direct_check_wide$np[direct_check_wide$point_id==i], ")"))
})

wrap_plots(plot_list, ncol = 3)

# -------------------------------------------------------------------
# PLAND check 
# -------------------------------------------------------------------
direct_check_pland <- sample_lsm(thornscrub_cropped_final, y = test_pts_sf, size = 610,
                                  what = "lsm_c_pland", shape = "circle", return_raster = FALSE) %>%
  filter(class == 1) %>%
  select(plot_id, value) %>%
  rename(point_id = plot_id, pland = value)

focal_check_pland <- terra::extract(pland_surface, vect(test_pts_sf))
names(focal_check_pland) <- c("ID", "pland_focal")

comparison_pland <- direct_check_pland %>%
  left_join(focal_check_pland, by = c("point_id" = "ID"))

print(comparison_pland)

# -------------------------------------------------------------------
# Pland panel
# -------------------------------------------------------------------
plot_list_pland <- lapply(1:6, function(i) {
  pt <- plot_pts[i, ]
  buf <- buffer_polys[i, ]
  val <- terra::extract(pland_surface, vect(pt))[,2]
  direct_val <- comparison_pland$pland[comparison_pland$point_id == i]
  
  ggplot() +
    geom_raster(data = crop_df, aes(x = x, y = y, fill = factor(thornscrub))) +
    scale_fill_manual(values = c("0" = "grey85", "1" = "#4daf4a"), guide = "none") +
    geom_sf(data = buf, fill = NA, color = "red", linewidth = 0.8, inherit.aes = FALSE) +
    geom_sf(data = pt, color = "black", size = 2, inherit.aes = FALSE) +
    coord_sf(datum = NA) +
    labs(title = paste0("Pt ", i, ": pland=", round(val,1), " (lsm=", round(direct_val,1), ")")) +
    theme_minimal(base_size = 9) +
    theme(axis.title = element_blank(), axis.text = element_blank())
})

wrap_plots(plot_list_pland, ncol = 3)


# -------------------------------------------------------------------
# Compute patch areas across the whole test landscape (once)
# -------------------------------------------------------------------
patch_areas <- as.data.frame(freq(patch_id_rast)) %>%
  mutate(area_ha = count * res(test_crop)[1]^2 / 10000) %>%
  select(patch_id = value, area_ha)

# -------------------------------------------------------------------
# 2. area_mn and lpi via focal window, looking up whole-patch areas
# -------------------------------------------------------------------
window_area_ha <- (sum(w_binary) * res(test_crop)[1]^2) / 10000

area_mn_weighted <- focal(patch_id_rast, w = w_binary, fun = function(x, ...) {
  ids_present <- na.omit(x)
  if (length(ids_present) == 0) return(0)
  
  # count in-window cells per patch ID -- this IS the clipped area, no lookup needed
  cell_counts <- table(ids_present)
  areas_ha <- cell_counts * res(patch_id_rast)[1]^2 / 10000
  
  mean(areas_ha)
}, na.rm = TRUE)

lpi_weighted <- focal(patch_id_rast, w = w_binary, fun = function(x, ...) {
  ids_present <- na.omit(x)
  if (length(ids_present) == 0) return(0)
  
  cell_counts <- table(ids_present)
  areas_ha <- cell_counts * res(patch_id_rast)[1]^2 / 10000
  window_area_ha <- length(x) * res(patch_id_rast)[1]^2 / 10000
  
  max(areas_ha) / window_area_ha * 100
}, na.rm = TRUE)

par(mfrow = c(1,3))
plot(test_crop, main = "Thornscrub")
plot(area_mn_weighted, main = "Mean patch area (area_mn)")
plot(lpi_weighted, main = "Largest patch index (lpi)")
par(mfrow = c(1,1))

# -------------------------------------------------------------------
# patch area and lpi weighted against direct sample_lsm() at the same test points
# -------------------------------------------------------------------
direct_check_area <- sample_lsm(test_crop, y = test_pts_sf, size = 610,
                                 what = c("lsm_c_area_mn", "lsm_c_lpi"),
                                 shape = "circle", return_raster = FALSE) %>%
  filter(class == 1)

direct_check_area_wide <- direct_check_area %>%
  select(plot_id, metric, value) %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  rename(point_id = plot_id)

focal_check_area <- terra::extract(c(area_mn_weighted, lpi_weighted), vect(test_pts_sf))
names(focal_check_area) <- c("ID", "area_mn_focal", "lpi_focal")

comparison_area <- direct_check_area_wide %>%
  left_join(focal_check_area, by = c("point_id" = "ID"))

print(comparison_area)

# -------------------------------------------------------------------
# patch area and lpi panel
# -------------------------------------------------------------------
plot_list_area <- lapply(1:6, function(i) {
  pt <- plot_pts[i, ]
  buf <- buffer_polys[i, ]
  
  area_val <- terra::extract(area_mn_weighted, vect(pt))[,2]
  lpi_val <- terra::extract(lpi_weighted, vect(pt))[,2]
  
  area_direct <- comparison_area$area_mn[comparison_area$point_id == i]
  lpi_direct <- comparison_area$lpi[comparison_area$point_id == i]
  
  ggplot() +
    geom_raster(data = crop_df, aes(x = x, y = y, fill = factor(thornscrub))) +
    scale_fill_manual(values = c("0" = "grey85", "1" = "#4daf4a"), guide = "none") +
    geom_sf(data = buf, fill = NA, color = "red", linewidth = 0.8, inherit.aes = FALSE) +
    geom_sf(data = pt, color = "black", size = 2, inherit.aes = FALSE) +
    coord_sf(datum = NA) +
    labs(title = paste0("Pt ", i, ": area_mn=", round(area_val,2), " (lsm=", round(area_direct,2), ")\n",
                         "lpi=", round(lpi_val,1), " (lsm=", round(lpi_direct,1), ")")) +
    theme_minimal(base_size = 9) +
    theme(axis.title = element_blank(), axis.text = element_blank())
})

wrap_plots(plot_list_area, ncol = 3)

# -------------------------------------------------------------------
# Count true edges: cell differs from its right neighbor, or its neighbor below
# -------------------------------------------------------------------
r <- test_crop
cell_size <- res(r)[1]
m <- as.matrix(r, wide = TRUE)  # raw values as a matrix, rows/cols aligned to the grid

nr <- nrow(m); nc <- ncol(m)

# horizontal edges: compare each cell to the one to its right (shift columns)
horiz_diff <- matrix(0, nrow = nr, ncol = nc)
horiz_diff[, 1:(nc-1)] <- (m[, 1:(nc-1)] != m[, 2:nc]) * 1

# vertical edges: compare each cell to the one below (shift rows)
vert_diff <- matrix(0, nrow = nr, ncol = nc)
vert_diff[1:(nr-1), ] <- (m[1:(nr-1), ] != m[2:nr, ]) * 1

# handle NAs from the comparison (NA != NA gives NA, not TRUE/FALSE)
horiz_diff[is.na(horiz_diff)] <- 0
vert_diff[is.na(vert_diff)] <- 0

total_edges_m <- horiz_diff + vert_diff

total_edges <- rast(total_edges_m, extent = ext(r), crs = crs(r))
# don't touch res() separately -- it's already implied by extent/dimensions

res(total_edges)  # confirm this matches res(r) automatically
plot(total_edges)

# -------------------------------------------------------------------
# 2. Sum edge length within each focal window, divide by window area (ha)
# -------------------------------------------------------------------
w <- focalMat(test_crop, d = 610, type = "circle")
w_binary <- ifelse(w > 0, 1, 0)

edge_count_sum <- focal(total_edges, w = w_binary, fun = function(x, ...) sum(x, na.rm=TRUE), na.rm = TRUE)
total_edge_length_m <- edge_count_sum * cell_size

window_area_ha <- (sum(w_binary) * cell_size^2) / 10000
ed_test_v2 <- total_edge_length_m / window_area_ha

plot(ed_test_v2, main = "Edge density (ed), corrected")

par(mfrow = c(1,2))
plot(test_crop, main = "Thornscrub")
plot(ed_test_v2, main = "Edge density (ed), corrected")
par(mfrow = c(1,1))

direct_check_ed <- sample_lsm(test_crop, y = test_pts_sf, size = 610,
                               what = "lsm_c_ed", shape = "circle", return_raster = FALSE) %>%
  filter(class == 1)

direct_check_ed_wide <- direct_check_ed %>%
  select(plot_id, value) %>%
  rename(point_id = plot_id, ed = value)

focal_check_ed2 <- terra::extract(ed_test_v2, vect(test_pts_sf))
names(focal_check_ed2) <- c("ID", "ed_focal2")

comparison_ed <- direct_check_ed_wide %>%
  left_join(focal_check_ed2, by = c("point_id" = "ID"))

print(comparison_ed)

# -------------------------------------------------------------------
# 4. Panel plots
# -------------------------------------------------------------------
plot_list_ed <- lapply(1:6, function(i) {
  pt <- plot_pts[i, ]
  buf <- buffer_polys[i, ]
  ed_val <- terra::extract(ed_test_v2, vect(pt))[,2]
  ed_direct <- comparison_ed$ed[comparison_ed$point_id == i]
  
  ggplot() +
    geom_raster(data = crop_df, aes(x = x, y = y, fill = factor(thornscrub))) +
    scale_fill_manual(values = c("0" = "grey85", "1" = "#4daf4a"), guide = "none") +
    geom_sf(data = buf, fill = NA, color = "red", linewidth = 0.8, inherit.aes = FALSE) +
    geom_sf(data = pt, color = "black", size = 2, inherit.aes = FALSE) +
    coord_sf(datum = NA) +
    labs(title = paste0("Pt ", i, ": ed=", round(ed_val,1), " (lsm=", round(ed_direct,1), ")")) +
    theme_minimal(base_size = 9) +
    theme(axis.title = element_blank(), axis.text = element_blank())
})

wrap_plots(plot_list_ed, ncol = 3)













# # ---------------------------------------------------------------
# # Crop tightly to points + max radius, trim, force into memory
# # ---------------------------------------------------------------
# buffer_ext <- st_buffer(st_as_sfc(st_bbox(pts_sf)), 4000)
# thornscrub_cropped <- crop(thornscrub_raster_full, vect(buffer_ext))
# thornscrub_cropped <- trim(thornscrub_cropped)

# # force fully into memory (avoids repeated disk reads per point)
# values(thornscrub_cropped) <- values(thornscrub_cropped)
# inMemory(thornscrub_cropped)  

# # ---------------------------------------------------------------
# # Batched extraction 
# # ---------------------------------------------------------------
# class_metrics <- c("lsm_c_pland", "lsm_c_area_mn", "lsm_c_shape_mn", "lsm_c_cai_mn",
#                     "lsm_c_ed", "lsm_c_cohesion", "lsm_c_lpi", "lsm_c_gyrate_mn", "lsm_c_pd", "lsm_c_np")
# batch_size <- 500

# run_batched <- function(radius, pts, raster, metrics, batch_size) {
#   n_batches <- ceiling(nrow(pts) / batch_size)
#   out <- vector("list", n_batches)

#   for (b in 1:n_batches) {
#     idx <- ((b - 1) * batch_size + 1):min(b * batch_size, nrow(pts))
#     batch_result <- sample_lsm(
#       landscape = raster, y = pts[idx, ],
#       size = radius, what = metrics, shape = "circle", return_raster = FALSE
#     ) %>% filter(class == 1)

#     batch_result$point_id <- idx[batch_result$plot_id]
#     out[[b]] <- batch_result

#     gc()
#     cat("  Radius", radius, "- batch", b, "of", n_batches, "done at", format(Sys.time()), "\n")
#   }
#   bind_rows(out) %>% mutate(buffer_radius = radius)
# }

# ## --- Create subsample --- ##
# set.seed(1)
# subsample_pooled <- pts_sf_all %>%
#   st_drop_geometry() %>%
#   mutate(row_id = row_number()) %>%
#   group_by(state, case_) %>%
#   slice_sample(prop = 0.10) %>%  # 10% stratified by state and used/available
#   ungroup() %>%
#   pull(row_id)

# pts_sf_subsample <- pts_sf_all[subsample_pooled, ]
# nrow(pts_sf_subsample)  # sanity check on resulting size
# table(pts_sf_subsample$state, pts_sf_subsample$case_)  # confirm proportions preserved

# extraction_list <- list()
# for (nm in names(radii_test)) {
#   cat("Starting:", nm, "radius =", radii_test[[nm]], "at", format(Sys.time()), "\n")
#   extraction_list[[nm]] <- run_batched(radii_test[[nm]], pts_sf_subsample, thornscrub_cropped, class_metrics, batch_size)
#   saveRDS(extraction_list, "output/objects/radii_variation_test_subsample.rds")
#   cat("Completed:", nm, "at", format(Sys.time()), "\n\n")
# }

# extraction_wide <- lapply(extraction_list, function(x) {
# x %>% select(point_id, metric, value) %>% pivot_wider(names_from = metric, values_from = value)
# })


# pland_full_check <- run_batched(610, pts_sf_all, thornscrub_full, "lsm_c_pland", batch_size)

# nrow(pts_sf_all)
# n_distinct(pland_full_check$point_id)
# 1 - (n_distinct(pland_full_check$point_id) / nrow(pts_sf_all))  # % missing pland at 610m, full dataset


# library(terra)

# thornscrub_full <- rast("output/thornscrub_raster_full.tif")  # your original full, uncropped raster

# w2 <- focalMat(thornscrub_full, d = 82, type = "circle")
# w_binary2 <- ifelse(w2 > 0, 1, 0)

# pland_surface2 <- focal(thornscrub_full, w = w_binary2, fun = mean, na.rm = TRUE) * 100

# pts_vect <- vect(pts_sf_all)
# pland_values2 <- terra::extract(pland_surface2, pts_vect)

# sum(is.na(pland_values[,2]))
# mean(pland_values[,2] == 0, na.rm = TRUE)


# sum(is.na(pland_values$focal_mean))
# summary(pland_values$focal_mean)

# sum(is.na(pland_values2$focal_mean))
# summary(pland_values2$focal_mean)





# ext_check <- extraction_list[["mean_half"]]

# # does area_mn ever appear as a row with NA value, vs simply never appearing?
# ext_check %>% filter(metric == "area_mn") %>% summarise(n = n(), n_na_value = sum(is.na(value)))

# # compare point count for area_mn vs pland directly, same object
# ext_check %>% filter(metric == "area_mn") %>% pull(point_id) %>% n_distinct()
# ext_check %>% filter(metric == "pland") %>% pull(point_id) %>% n_distinct()

# wide_check <- ext_check %>% select(point_id, metric, value) %>% pivot_wider(names_from = metric, values_from = value)

# nrow(wide_check)  # how many rows after pivot?
# sum(is.na(wide_check$area_mn))
# sum(is.na(wide_check$pland))
# sum(!is.na(wide_check$area_mn) & is.na(wide_check$pland))  # any mismatched patterns?

# ext_check %>% filter(metric == "area_mn") %>% count(point_id) %>% filter(n > 1)  # any point with >1 row for area_mn?
# ext_check %>% filter(metric == "pland") %>% count(point_id) %>% filter(n > 1)
# ext_check %>% count(metric)  # total row count per metric -- should all be equal if coverage is truly identical

# out_check <- all_points %>% left_join(wide_check, by = "point_id")

# nrow(out_check)
# sum(is.na(out_check$area_mn))
# sum(is.na(out_check$pland))

# fill_zero_metrics
# intersect(fill_zero_metrics, names(wide_check))









# summary(extraction_wide[[1]])
# summary(extraction_wide[[2]])
# summary(extraction_wide[[3]])
# summary(extraction_wide[[4]])
# temp <- extraction_wide[[1]]

# nrow(pts_sf_subsample)                 # your known denominator, e.g. ~143,535
# names(extraction_list)                 # should be your 4 radius labels
# map_int(extraction_list, ~ n_distinct(.x$point_id))  # how many unique points got ANY row, per radius
# all_points <- tibble(point_id = unique(pts_sf_subsample$point_id))
# nrow(all_points)  # must equal nrow(pts_sf_subsample)
# fill_zero_metrics <- c("pland", "np", "pd", "lpi", "ed")

# wide_list <- map(extraction_list, function(ext) {
#   wide <- ext %>%
#     select(point_id, metric, value) %>%
#     pivot_wider(names_from = metric, values_from = value)
  
#   out <- all_points %>%
#     left_join(wide, by = "point_id") %>%
#     mutate(across(all_of(intersect(fill_zero_metrics, names(wide))), ~ coalesce(., 0)))
  
#   out
# })

# map_int(wide_list, nrow)  # should all equal nrow(all_points)

# check_na <- map_dfr(names(wide_list), function(nm) {
#   d <- wide_list[[nm]]
#   tibble(radius_set = nm, metric = setdiff(names(d), "point_id"),
#          pct_na = sapply(setdiff(names(d), "point_id"), function(m) mean(is.na(d[[m]]))))
# })
# print(check_na, n = 50)  # eyeball this before proceeding -- fill-zero metrics should show 0


# summary_across_radii <- map_dfr(names(wide_list), function(nm) {
#   d <- wide_list[[nm]]
#   map_dfr(setdiff(names(d), "point_id"), function(m) {
#     x_all <- d[[m]]
#     x_nonzero <- x_all[!is.na(x_all) & x_all != 0]
#     tibble(
#       radius_set = nm, metric = m,
#       pct_na = mean(is.na(x_all)),
#       pct_zero = mean(x_all == 0, na.rm = TRUE),
#       n_nonzero = length(x_nonzero),
#       cv_all = sd(x_all, na.rm = TRUE) / mean(x_all, na.rm = TRUE),
#       cv_nonzero = sd(x_nonzero) / mean(x_nonzero)
#     )
#   })
# })

# print(summary_across_radii %>% arrange(metric, radius_set), n = 100)


# na_points <- pts_sf_subsample %>% 
#   left_join(wide_check %>% select(point_id, pland_raw = pland), by = "point_id") %>%
#   filter(is.na(pland_raw))

# nrow(na_points)
# st_bbox(na_points)   # compare against your raster's extent
# ext(thornscrub_cropped)

# # quick visual check
# mapview::mapview(na_points %>% slice_sample(n = 2000))  # sample to keep it light

# issf_all_states_labeled <- pts_sf_subsample %>% st_drop_geometry() %>%
#   left_join(wide_check %>% select(point_id, pland_raw = pland), by = "point_id")

# issf_all_states_labeled %>% group_by(state) %>% 
#   summarise(pct_no_thornscrub = mean(is.na(pland_raw)), n = n())

# issf_all_states_labeled %>% group_by(state, case_) %>% 
#   summarise(pct_no_thornscrub = mean(is.na(pland_raw)), n = n())

# joined_temp <- left_join(temp, pts_sf_subsample, by = "point_id")
# nrow(temp)/nrow(pts_sf_subsample)

# fill_zero_metrics <- c("pland", "np", "pd", "lpi", "ed")
# all_points <- tibble(point_id = unique(pts_sf_subsample$point_id))

# wide_list <- map(extraction_list, function(ext) {
#   wide <- ext %>% select(point_id, metric, value) %>% pivot_wider(names_from=metric, values_from=value)
#   all_points %>% left_join(wide, by="point_id") %>%
#     mutate(across(all_of(intersect(fill_zero_metrics, names(wide))), ~coalesce(.,0)))
# })

# issf_all_states_labeled %>% filter(state == "traveling") %>% summarise(pct_na_82m = mean(is.na(pland_raw)))

# summary_across_radii <- map_dfr(names(wide_list), function(nm) {
#   d <- wide_list[[nm]]
#   map_dfr(setdiff(names(d), "point_id"), function(m) {
#     x <- d[[m]]
#     tibble(radius_set = nm, metric = m, n = sum(!is.na(x)), pct_na = mean(is.na(x)),
#            median = median(x, na.rm=TRUE), cv = sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE))
#   })
# })

# print(summary_across_radii %>% arrange(metric, radius_set), n = 100)

# ggplot(summary_across_radii, aes(x = radius_set, y = cv, fill = radius_set)) +
#   geom_col() + facet_wrap(~metric, scales = "free_y") +
#   theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   labs(title = "Metric variability across radius choices (10% subsample)")

# ggplot(summary_across_radii, aes(x = radius_set, y = pct_na, fill = radius_set)) +
#   geom_col() + facet_wrap(~metric, scales = "free_y") +
#   theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   labs(title = "Missing data rate across radius choices (10% subsample)")
# # step_counts_raw <- data_adj %>% count(Deployment_ID, name = "n_raw")
# # step_counts_final <- issf_all_states %>% filter(case_ == TRUE) %>% count(Deployment_ID, name = "n_final")

# # retention_pooled <- step_counts_raw %>%
# #   left_join(step_counts_final, by = "Deployment_ID") %>%
# #   mutate(n_final = coalesce(n_final, 0), pct_retained = n_final / n_raw)
# # print(retention_pooled)#, n = 35)


