library(amt); library(dplyr); library(lubridate); library(sf)
library(momentuHMM); library(terra); library(landscapemetrics)

# -------------------------------------------------------------------
# Read in needed data and outputs
#--------------------------------------------------------------------

hmm_input <- readRDS("output/objects/hmm_input.rds") %>% rename(BurstID = burst)
data_adj <- readRDS("output/objects/hmm_data_adj.rds")
hmm_model_list <- readRDS("output/objects/hmm_model_list.rds")
m3_cosinor <- hmm_model_list[[4]]
thornscrub_binary <- rast("output/south_texas_thornscrub_binary.tif")
target_crs <- crs(thornscrub_binary)

# -------------------------------------------------------------------
# Identify states for locations and subset traveling state
#--------------------------------------------------------------------

# Confirm model/data alignment before assigning states
print(m3_cosinor)  # spot check: loglik should be -57113.75
stopifnot(nrow(data_adj) == length(viterbi(m3_cosinor)))

# Attach decoded 3-state Viterbi output
data_adj$state <- viterbi(m3_cosinor)
data_adj$state_label <- recode(as.character(data_adj$state),
                                "1" = "local_movement", "2" = "encamped", "3" = "traveling")
table(data_adj$state_label)  # sanity check against known ~23%/40%/37% split

# Attach animal ID
animal_map <- hmm_input %>% st_drop_geometry() %>% distinct(BurstID, Deployment_ID)
data_adj <- data_adj %>% left_join(animal_map, by = c("ID" = "BurstID"))

# Reattach meters-scale x/y (undo km rescale) + timestamp -- check for row-count drift
n_before <- nrow(data_adj)
data_adj <- data_adj %>%
  left_join(hmm_input %>% st_drop_geometry() %>% select(BurstID, timestamp, x, y),
            by = c("ID" = "BurstID", "timestamp" = "timestamp"))
stopifnot(nrow(data_adj) == n_before)  # catches join-induced duplication

traveling_only <- data_adj %>% filter(state_label == "traveling")
nrow(traveling_only)

# Pull lon/lat from hmm_input and join onto traveling_only
traveling_only <- traveling_only %>%
  left_join(
    hmm_input %>% st_drop_geometry() %>% select(BurstID, timestamp, longitude, latitude),
    by = c("ID" = "BurstID", "timestamp" = "timestamp")
  )

# Confirm join worked and didn't duplicate rows
sum(is.na(traveling_only$longitude))  # should be 0
nrow(traveling_only)  # compare to pre-join row count -- should be unchanged


# Reproject lon/lat (WGS84) into the raster's CRS
traveling_pts_sf <- traveling_only %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
  st_transform(crs = target_crs)
mapview::mapview(traveling_pts_sf)

coords <- st_coordinates(traveling_pts_sf)
traveling_only$x_proj <- coords[, 1]
traveling_only$y_proj <- coords[, 2]

# -------------------------------------------------------------------
# Used vs available points from AMT for traveling only
#--------------------------------------------------------------------

# Build the amt track 
trk <- traveling_only %>%
  arrange(Deployment_ID.x, timestamp) %>%
  nest(data = -Deployment_ID.x) %>%
  mutate(trk = lapply(data, function(d) {
    make_track(d, .x = x_proj, .y = y_proj, .t = timestamp, crs = target_crs)
  }))

issf_traveling <- trk %>%
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

nrow(issf_traveling)

# how many bursts survive per individual, and how long are they?
burst_summary <- issf_traveling %>% 
  distinct(Deployment_ID.x, burst_) %>%   # or whatever the burst column is named post-steps_by_burst()
  count(Deployment_ID.x, name = "n_bursts")

step_counts <- issf_traveling %>% filter(case_ == TRUE) %>% count(Deployment_ID.x)
print(step_counts, n = 35)

# how many individuals/bursts triggered this?
step_counts_raw <- traveling_only %>% count(Deployment_ID.x, name = "n_raw_traveling")
step_counts_final <- issf_traveling %>% filter(case_ == TRUE) %>% count(Deployment_ID.x, name = "n_final_steps")

retained <- step_counts_raw %>% 
  left_join(step_counts_final, by = c("Deployment_ID.x")) %>%
  mutate(n_final_steps = coalesce(n_final_steps, 0),
         pct_retained = n_final_steps / n_raw_traveling)
retained

# -------------------------------------------------------------------
# Landscape metrics by points
#--------------------------------------------------------------------

library(landscapemetrics); library(terra); library(sf); library(dplyr); library(purrr)

saveRDS(issf_traveling, "output/objects/issf_points_traveling.rds")
saveRDS(subsample_sf, "output/objects/issf_subsample_traveling.rds")

issf_traveling <- readRDS("output/objects/issf_points_traveling.rds")
subsample_sf <- readRDS("output/objects/issf_subsample_traveling.rds")
thornscrub_raster_full <- rast("output/thornscrub_raster_full.tif")

issf_traveling <- issf_traveling %>% rename(id = Deployment_ID.x)

# Points to extract covariates at -- typically the step ENDPOINT is where selection is evaluated
pts_sf <- issf_traveling %>% 
  st_as_sf(coords = c("x2_", "y2_"), crs = target_crs, remove = FALSE) %>%
  mutate(point_id = row_number()) 

# Need to include water for the radii
thornscrub_raster_full <- thornscrub_binary
thornscrub_raster_full[is.na(thornscrub_raster_full)] <- 0

writeRaster(thornscrub_raster_full, "output/thornscrub_raster_full.tif")
thornscrub_raster_full <- rast("output/thornscrub_raster_full.tif")
freq(thornscrub_raster_full)  # confirm 0 and 1 only, no NA


# ---------------------------------------------------------------
# Crop tightly to points + max radius, trim, force into memory
# ---------------------------------------------------------------
buffer_ext <- st_buffer(st_as_sfc(st_bbox(pts_sf)), 4000)
thornscrub_cropped <- crop(thornscrub_raster_full, vect(buffer_ext))
thornscrub_cropped <- trim(thornscrub_cropped)

# force fully into memory (avoids repeated disk reads per point)
values(thornscrub_cropped) <- values(thornscrub_cropped)
inMemory(thornscrub_cropped)  

# ---------------------------------------------------------------
# Batched extraction 
# ---------------------------------------------------------------
class_metrics <- c("lsm_c_pland", "lsm_c_area_mn", "lsm_c_shape_mn", "lsm_c_cai_mn",
                    "lsm_c_enn_mn", "lsm_c_clumpy", "lsm_c_ed", "lsm_c_cohesion")
batch_size <- 500

run_batched <- function(radius, pts, raster, metrics, batch_size) {
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
    cat("  Radius", radius, "- batch", b, "of", n_batches, "done at", format(Sys.time()), "\n")
  }
  bind_rows(out) %>% mutate(buffer_radius = radius)
}


# Get the precise value to base the radius on
issf_traveling %>% filter(case_ == TRUE) %>% 
  summarise(median_sl = median(sl_, na.rm = TRUE), mean_sl = mean(sl_, na.rm = TRUE))

# radius = half of whichever you choose to anchor on
final_radius <- round(338 / 2)  # or mean_sl / 2, your call

landscape_metrics_all_extraction <- run_batched(final_radius, pts_sf, thornscrub_cropped, class_metrics, batch_size)
saveRDS(landscape_metrics_all_extraction, "output/objects/final_landscape_metrics.rds")