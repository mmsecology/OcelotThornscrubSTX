library(amt); library(dplyr); library(lubridate); library(sf)
library(momentuHMM); library(terra); library(landscapemetrics)
library(tidyverse); library(glmmTMB); library(mgcv)

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
table(data_adj$state_label)

# Attach animal ID
animal_map <- hmm_input %>% st_drop_geometry() %>% distinct(BurstID, Deployment_ID)
data_adj <- data_adj %>% left_join(animal_map, by = c("ID" = "BurstID"))

# Reattach meters-scale x/y (undo km rescale) + timestamp -- check for row-count drift
n_before <- nrow(data_adj)
data_adj <- data_adj %>%
  left_join(hmm_input %>% st_drop_geometry() %>% select(BurstID, timestamp, x, y),
            by = c("ID" = "BurstID", "timestamp" = "timestamp"))
stopifnot(nrow(data_adj) == n_before)  # catches join-induced duplication

## --- Create traveling only dataset --- ##
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

## --- Create local only dataset --- ##
local_only <- data_adj %>% filter(state_label == "local_movement")
nrow(local_only)

# Pull lon/lat from hmm_input and join onto local_only
local_only <- local_only %>%
  left_join(
    hmm_input %>% st_drop_geometry() %>% select(BurstID, timestamp, longitude, latitude),
    by = c("ID" = "BurstID", "timestamp" = "timestamp")
  )

# Confirm join worked and didn't duplicate rows
sum(is.na(local_only$longitude))  # should be 0
nrow(local_only)  # compare to pre-join row count -- should be unchanged

# Reproject lon/lat (WGS84) into the raster's CRS
local_pts_sf <- local_only %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
  st_transform(crs = target_crs)
mapview::mapview(local_pts_sf)

coords <- st_coordinates(local_pts_sf)
local_only$x_proj <- coords[, 1]
local_only$y_proj <- coords[, 2]

## --- Create encamped only dataset --- ##
encamped_only <- data_adj %>% filter(state_label == "encamped")
nrow(encamped_only)

# Pull lon/lat from hmm_input and join onto encamped_only
encamped_only <- encamped_only %>%
  left_join(
    hmm_input %>% st_drop_geometry() %>% select(BurstID, timestamp, longitude, latitude),
    by = c("ID" = "BurstID", "timestamp" = "timestamp")
  )

# Confirm join worked and didn't duplicate rows
sum(is.na(encamped_only$longitude))  # should be 0
nrow(encamped_only)  # compare to pre-join row count -- should be unchanged

# Reproject lon/lat (WGS84) into the raster's CRS
encamped_pts_sf <- encamped_only %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
  st_transform(crs = target_crs)
mapview::mapview(encamped_pts_sf)

coords <- st_coordinates(encamped_pts_sf)
encamped_only$x_proj <- coords[, 1]
encamped_only$y_proj <- coords[, 2]

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

#issf_traveling <- issf_traveling %>% rename(id = Deployment_ID.x)

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
                    "lsm_c_enn_mn", "lsm_c_clumpy", "lsm_c_ed", "lsm_c_cohesion", "lsm_c_lpi", "lsm_c_gyrate_mn", "lsm_c_pd", "lsm_c_np")
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
issf_traveling %>% filter(case_ == TRUE) %>% select(sl_) %>% summarise(quantile_95 = quantile(sl_, probs = 0.95, na.rm = TRUE))
#  summarise(median_sl = median(sl_, na.rm = TRUE), mean_sl = mean(sl_, na.rm = TRUE))

# radius = half of whichever you choose to anchor on
final_radius <- round(814/ 2)  # or mean_sl / 2, your call

landscape_metrics_all_extraction <- run_batched(final_radius, pts_sf, thornscrub_cropped, class_metrics, batch_size)
saveRDS(landscape_metrics_all_extraction, "output/objects/final_landscape_metrics2.rds") # added lpi, gyrate, pd, np
saveRDS(landscape_metrics_all_extraction, "output/objects/final_landscape_metrics3.rds") # all metrics but at 95% of step length

landscape_metrics_all_extraction <- readRDS("output/objects/final_landscape_metrics2.rds")
# Confirm class filter held
table(landscape_metrics_all_extraction$class)  # should show only 1

# -------------------------------------------------------------------
# Pivot to wide: one row per point, one column per metric
# -------------------------------------------------------------------
metrics_wide <- landscape_metrics_all_extraction %>%
  select(point_id, metric, value) %>%
  pivot_wider(names_from = metric, values_from = value)

# pland: missing = genuinely 0% thornscrub in buffer (no row extracted at all)
# other metrics: missing = no patch present, genuinely undefined -- leave as NA
all_points <- tibble(point_id = unique(pts_sf$point_id))

metrics_complete <- all_points %>%
  left_join(metrics_wide, by = "point_id") %>%
  mutate(pland = coalesce(pland, 0))

# quick NA audit per metric at this final radius
metrics_complete %>%
  summarise(across(-point_id, ~ mean(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "metric", values_to = "pct_na") %>%
  arrange(desc(pct_na))


nrow(pts_sf)  # total points
14906 / nrow(pts_sf)  # what % is this?

# is it concentrated among certain individuals?
missing_points <- pts_sf %>% st_drop_geometry() %>% 
  filter(!point_id %in% metrics_wide$point_id)

missing_points %>% count(id) %>% arrange(desc(n))  # or whatever the individual ID column is called here

# is it more common for used vs. available points?
missing_points %>% count(case_)
pts_sf %>% st_drop_geometry() %>% count(case_)  # compare proportions

all_points <- tibble(point_id = unique(pts_sf$point_id))

metrics_complete <- all_points %>%
  left_join(metrics_wide, by = "point_id") %>%
  mutate(pland = coalesce(pland, 0))

# confirm the fill worked and matches the expected count
sum(metrics_complete$pland == 0, na.rm = TRUE)  # should be >= 14906 (some non-missing points may also have pland genuinely near 0)

# -------------------------------------------------------------------
# Check correlations and run glmmTMB
# -------------------------------------------------------------------
issf_model_data <- pts_sf %>% st_drop_geometry() %>% left_join(metrics_complete, by = "point_id")
predictor_cols <- c("pland", "shape_mn", "cai_mn", "ed", "cohesion", "lpi", "gyrate_mn", "pd", "np")

cor_matrix <- metrics_complete %>% select(all_of(predictor_cols)) %>% cor(use = "pairwise.complete.obs")
print(round(cor_matrix, 2))

library(car)
vif_check <- lm(as.numeric(case_) ~ pland + shape_mn + cai_mn + ed + cohesion,
                 data = issf_model_data)
vif(vif_check)
nobs(vif_check)  # check how much listwise-deletion shrinks the sample given enn_mn/clumpy missingness

issf_model_data$step_id_ <- paste0(issf_model_data$id, "_", issf_model_data$step_id_)
n_distinct(issf_model_data$step_id_)

issf_model_data <- issf_model_data %>%
  mutate(across(c(pland, shape_mn, cai_mn, ed, cohesion), ~ as.numeric(scale(.)), .names = "{.col}_z"))

m_issf_final <- glmmTMB(
  case_ ~ -1 + 
    pland_z + (0 + pland_z | id) +
    shape_mn_z + (0 + shape_mn_z | id) +
    cai_mn_z + (0 + cai_mn_z | id) +
    ed_z + (0 + ed_z | id) +
    cohesion_z + (0 + cohesion_z | id) +
    log_sl_ + (0 + log_sl_ | id) +
    cos_ta_ +
    (1 | step_id_),
  family = poisson(), doFit = TRUE,
  data = issf_model_data,
  map = list(theta = factor(c(1:6, NA))),
  start = list(theta = c(rep(0, times = 6),log(10000)))
)

VarCorr(m_issf_final)  # confirm step_id_ variance is fixed large, others estimated
summary(m_issf_final)

# Pull fixed-effect coefficients directly
fixef_est <- fixef(m_issf_final)$cond
fixef_est

# Manual log-RSS function: log-RSS = sum(beta_i * (x1_i - x2_i)) for all covariates
compute_log_rss <- function(coefs, x1, x2, covs) {
  rowSums(sapply(covs, function(v) coefs[[v]] * (x1[[v]] - x2[[v]])))
}

make_rss_plot_manual <- function(model, data, covariate, xlab) {
  covs <- c("pland_z", "shape_mn_z", "cai_mn_z", "ed_z", "cohesion_z", "log_sl_")
  coefs <- fixef(model)$cond
  
  x1 <- as.data.frame(matrix(0, nrow = 100, ncol = length(covs)))
  names(x1) <- covs
  x1[[covariate]] <- seq(quantile(data[[covariate]], 0.02, na.rm = TRUE),
                          quantile(data[[covariate]], 0.98, na.rm = TRUE), length.out = 100)
  x1$log_sl_ <- mean(data$log_sl_, na.rm = TRUE)
  x1$cos_ta_ <- mean(data$cos_ta_, na.rm = TRUE)
  
  x2 <- x1
  x2[[covariate]] <- 0  # reference value
  
  log_rss_vals <- compute_log_rss(coefs, x1, x2, covs)
  
  ggplot(data.frame(x = x1[[covariate]], log_rss = log_rss_vals), aes(x = x, y = log_rss)) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    labs(x = xlab, y = "log-RSS", title = paste("Selection:", xlab)) +
    theme_minimal(base_size = 13)
}

make_rss_plot_manual(m_issf_final, issf_model_data, "pland_z", "Thornscrub % (standardized)")
make_rss_plot_manual(m_issf_final, issf_model_data, "shape_mn_z", "Patch shape (standardized)")
make_rss_plot_manual(m_issf_final, issf_model_data, "cai_mn_z", "Core area index (standardized)")
make_rss_plot_manual(m_issf_final, issf_model_data, "ed_z", "Edge density (standardized)")
make_rss_plot_manual(m_issf_final, issf_model_data, "cohesion_z", "Cohesion (standardized)")

# -------------------------------------------------------------------
# Run non-linear model
# -------------------------------------------------------------------

gam_data <- issf_model_data %>%
  group_by(step_id_) %>%
  filter(!any(case_ & if_any(all_of(predictor_cols), is.na))) %>%
  filter(!(!case_ & if_any(all_of(predictor_cols), is.na))) %>%
  ungroup()

length(unique(gam_data$step_id_))

gam_data$times <- 1
gam_data$case_binary <- as.numeric(gam_data$case_)
gam_data <- gam_data %>% 
  mutate(stratum = factor(step_id_), ID = factor(id))

sapply(gam_data %>% select(pland_z, shape_mn_z, cai_mn_z, ed_z, cohesion_z, log_sl_), 
       function(x) c(n_na = sum(is.na(x)), n_inf = sum(is.infinite(x)), n_nan = sum(is.nan(x))))

test_gam <- gam(
  cbind(times, stratum) ~ s(pland_z, k = 7),
  weights = case_binary, method = "REML", select = TRUE, family = cox.ph(),
  data = gam_data
)

gam_data %>% count(stratum) %>% count(n)  # distribution of rows-per-stratum -- should be one dominant value (e.g., 11)
gam_data %>% count(stratum) %>% filter(n == 1)  # any strata with only 1 row?

# Check how many strata would survive at different minimum-row thresholds
gam_data %>% count(stratum) %>% count(n >= 8)  # e.g., how many strata keep at least 8 of 11 rows?

# Rebuild the filter: keep only strata with the FULL original row count (cleanest, most conservative)
strata_row_counts <- gam_data %>% count(stratum)
full_strata <- strata_row_counts %>% filter(n >= 8) %>% pull(stratum)

gam_data_clean <- gam_data %>% filter(stratum %in% full_strata)

# confirm cleanup worked
gam_data_clean %>% count(stratum) %>% count(n)  # should show ONLY n=11
gam_data_clean %>% group_by(stratum) %>% summarise(n_used = sum(case_)) %>% count(n_used)

gam_test <- gam(
  cbind(times, stratum) ~ 
    s(pland_z, k = 7) + s(pland_z, ID, k = 5, bs = "fs") +
    log_sl_ + cos_ta_,
    data = gam_data_clean,
    method = "REML",
    family = cox.ph,
    weights = case_binary
)

summary(gam_test)

pop_smooth_pland <- smooth_estimates(gam_test, select = "s(pland_z)", n = 500) %>%
  add_confint()

head(pop_smooth_pland)  # check column names before plotting -- confirm .estimate, .se, and the covariate column name

smooths(gam_test)  # lists the exact smooth term labels available for `select`

ggplot(pop_smooth_pland, aes(x = pland_z, y = exp(.estimate))) +
  geom_ribbon(aes(ymin = exp(.lower_ci), ymax = exp(.upper_ci)), alpha = 0.3, fill = "steelblue4") +
  geom_line(linewidth = 1, color = "steelblue4") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  labs(x = "Thornscrub % (standardized)", y = "Relative selection strength",
       title = "Population-level selection for thornscrub cover") +
  theme_minimal(base_size = 13)

ggplot(pop_smooth_pland, aes(x = pland_z, y = .estimate)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.3, fill = "steelblue4") +
  geom_line(linewidth = 1, color = "steelblue4") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  labs(x = "Thornscrub % (standardized)", y = "Relative selection strength",
       title = "Population-level selection for thornscrub cover") +
  theme_minimal(base_size = 13)


# derivative to find where the curve flattens
deriv_pland <- derivatives(gam_test, select = "s(pland_z)", n = 500)

deriv_pland %>% filter(.lower_ci <= 0) %>% slice(1)  # first point where increase is no longer significant

# Pull ALL smooth estimates from the model (population + individual)
all_smooths <- smooth_estimates(gam_test, n = 500)

# check exact labels for the population smooth vs individual smooth
unique(all_smooths$.smooth)

# isolate population-level smooth
pop_smooth <- all_smooths %>% filter(.smooth == "s(pland_z)") %>% add_confint()

# isolate individual-level deviations
ind_smooth <- all_smooths %>% filter(.smooth == "s(pland_z,ID)")

# for each individual, add their deviation to the population estimate
# (requires matching on the same pland_z grid values between pop_smooth and ind_smooth)
ind_smooth_combined <- ind_smooth %>%
  left_join(pop_smooth %>% select(pland_z, pop_est = .estimate), by = "pland_z") %>%
  mutate(est_combined = .estimate + pop_est)

# Population plot
p_pop <- ggplot(pop_smooth, aes(x = pland_z, y = exp(.estimate))) +
  geom_ribbon(aes(ymin = exp(.lower_ci), ymax = exp(.upper_ci)), alpha = 0.3, fill = "grey60") +
  geom_line(linewidth = 0.75, color = "black") +
  ylim(0, max(exp(ind_smooth_combined$est_combined), na.rm = TRUE) * 1.1) +
  labs(x = "Thornscrub % (standardized)", y = "exp(f(pland))", title = "Population") +
  theme_bw()

# Individual-level plot
p_ind <- ggplot(pop_smooth, aes(x = pland_z, y = exp(.estimate))) +
  geom_line(data = ind_smooth_combined, aes(x = pland_z, y = exp(est_combined), group = ID),
            color = "steelblue", alpha = 0.5) +
  ylim(0, max(exp(ind_smooth_combined$est_combined), na.rm = TRUE) * 1.1) +
  labs(x = "Thornscrub % (standardized)", y = "exp(f(pland))", title = "Individual") +
  theme_bw()

library(cowplot)
plot_grid(p_pop, p_ind)

# Set your reference point -- e.g., pland_z = 0 (mean thornscrub cover, since standardized)
ref_pland <- 0

all_smooths <- smooth_estimates(gam_test, n = 500)
unique(all_smooths$.smooth)  # confirm exact labels

pop_smooth <- all_smooths %>% filter(.smooth == "s(pland_z)")

# find the fitted value AT the reference point (interpolate if ref isn't exactly on the grid)
ref_est <- approx(pop_smooth$pland_z, pop_smooth$.estimate, xout = ref_pland)$y

# log-RSS relative to the reference point
pop_smooth_rss <- pop_smooth %>%
  add_confint() %>%
  mutate(
    log_rss = .estimate - ref_est,
    log_rss_lower = .lower_ci - ref_est,
    log_rss_upper = .upper_ci - ref_est
  )

# Population RSS plot
p_pop_rss <- ggplot(pop_smooth_rss, aes(x = pland_z, y = exp(log_rss))) +
  geom_ribbon(aes(ymin = exp(log_rss_lower), ymax = exp(log_rss_upper)), alpha = 0.3, fill = "grey60") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  geom_line(linewidth = 0.75, color = "black") +
  labs(x = "Thornscrub % (standardized)", y = "Relative selection strength",
       title = paste("Population RSS (relative to pland_z =", ref_pland, ")")) +
  theme_bw()

# Individual RSS: same logic, but reference each individual's OWN combined curve at ref_pland
ind_smooth <- all_smooths %>% filter(.smooth == "s(pland_z,ID)")

ind_smooth_combined <- ind_smooth %>%
  left_join(pop_smooth %>% select(pland_z, pop_est = .estimate), by = "pland_z") %>%
  mutate(est_combined = .estimate + pop_est)

ref_est_by_id <- ind_smooth_combined %>%
  group_by(ID) %>%
  summarise(ref_val = approx(pland_z, est_combined, xout = ref_pland)$y)

ind_smooth_rss <- ind_smooth_combined %>%
  left_join(ref_est_by_id, by = "ID") %>%
  mutate(log_rss_ind = est_combined - ref_val)

p_ind_rss <- ggplot(pop_smooth_rss, aes(x = pland_z, y = exp(log_rss))) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  geom_line(data = ind_smooth_rss, aes(x = pland_z, y = exp(log_rss_ind), group = ID),
            color = "steelblue", alpha = 0.5) +
  labs(x = "Thornscrub % (standardized)", y = "Relative selection strength",
       title = "Individual RSS") +
  theme_bw()

library(cowplot)
plot_grid(p_pop_rss, p_ind_rss)




saveRDS(gam_data_clean, "output/objects/gam_data_clean.rds")

gam_data_clean <- readRDS("output/objects/gam_data_clean.rds")


tictoc::tic()
ocelot_gam_fs <- gam(
  cbind(times, stratum) ~ 
    s(pland_z, k = 7) + s(pland_z, ID, k = 10, bs = "re") +
    s(shape_mn_z, k = 7) + s(shape_mn_z, ID, k = 10, bs = "re") +
    s(cai_mn_z, k = 7) + s(cai_mn_z, ID, k = 10, bs = "re") +
    s(ed_z, k = 7) + s(ed_z, ID, k = 10, bs = "re") +
    log_sl_ + cos_ta_,
  data = gam_data_clean,
  method = "REML",
  family = cox.ph,
  weights = case_binary
)
tictoc::toc()

summary(ocelot_gam_fs)

library(gratia); library(dplyr); library(ggplot2); library(purrr); library(patchwork)

covariates <- c("pland_z", "shape_mn_z", "cai_mn_z", "ed_z", "cohesion_z")
xlabs <- c("Thornscrub cover (standardized)", "Patch shape (standardized)", 
           "Core area index (standardized)", "Edge density (standardized)", 
           "Cohesion (standardized)")

all_smooths <- smooth_estimates(ocelot_gam_fs, n = 500)
unique(all_smooths$.smooth)  # confirm exact population-term labels before filtering

# Build one plot per covariate
plot_list <- map2(covariates, xlabs, function(cov, xlab) {
  smooth_label <- paste0("s(", cov, ")")  # confirm this matches unique(all_smooths$.smooth) exactly
  
  pop_smooth <- all_smooths %>% 
    filter(.smooth == smooth_label) %>% 
    add_confint()
  
  ggplot(pop_smooth, aes(x = .data[[cov]], y = .estimate)) +
    geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.3, fill = "steelblue4") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_line(linewidth = 0.9, color = "steelblue4") +
    labs(x = xlab, y = "Partial effect") +
    theme_minimal(base_size = 11)
})

# Combine into one figure
wrap_plots(plot_list, ncol = 2) +
  plot_annotation(title = "Population-level selection curves: five landscape metrics")

hist(gam_data_clean$cohesion, breaks = 100)
quantile(gam_data_clean$cohesion, probs = c(0.01, 0.05, 0.5, 0.95, 0.99), na.rm = TRUE)


metric_spread <- sapply(c("pland", "shape_mn", "cai_mn", "ed", "cohesion"), function(x) 
  quantile(gam_data_clean[[x]], probs = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99), na.rm = TRUE))
print(round(metric_spread, 2))

# also worth a simple coefficient of variation check -- low CV = little relative spread
sapply(c("pland", "shape_mn", "cai_mn", "ed", "cohesion"), function(x) 
  sd(gam_data_clean[[x]], na.rm = TRUE) / mean(gam_data_clean[[x]], na.rm = TRUE))

quantile(gam_data_clean$shape_mn_z, probs = c(0.75, 0.90, 0.95, 0.99), na.rm = TRUE)
quantile(gam_data_clean$ed_z, probs = c(0.75, 0.90, 0.95, 0.99), na.rm = TRUE)

pop_smooth <- smooth_estimates(ocelot_gam_fs, select = "s(pland_z)", n = 500) %>%
  add_confint()

# Log scale (partial effect, 0 = landscape average)
ggplot(pop_smooth, aes(x = pland_z, y = .estimate)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.3, fill = "steelblue4") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_line(linewidth = 0.9, color = "steelblue4") +
  labs(x = "Thornscrub cover (standardized)", y = "Partial effect (log scale)") +
  theme_minimal(base_size = 13)

# Exp scale (partial effect, 1 = landscape average)
ggplot(pop_smooth, aes(x = pland_z, y = exp(.estimate))) +
  geom_ribbon(aes(ymin = exp(.lower_ci), ymax = exp(.upper_ci)), alpha = 0.3, fill = "steelblue4") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  geom_line(linewidth = 0.9, color = "steelblue4") +
  labs(x = "Thornscrub cover (standardized)", y = "exp(partial effect)") +
  theme_minimal(base_size = 13)
