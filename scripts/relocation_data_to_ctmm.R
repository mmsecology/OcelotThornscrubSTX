library(sf)
library(tidyverse)
library(mapview)
library(terra)
library(ctmm)
library(tictoc)

## ------ Read in location data for each set of collars ------ ##
## Collar locations have been processed using data_filtering_and_prep to clean up dataset
spec <- readRDS("E:/GitProjects/GitprojectDroughtSelection/data/all_gps_locations_sf.rds") %>% mutate(project = "spec") %>% distinct(Deployment_ID, DateUTC, .keep_all = TRUE)
min(spec$DateLocal)
max(spec$DateLocal)

unique(spec$Property)
spec <- spec %>% filter(Property %in% c("LANWR", "El Sauz"))

## ------ Summarise by deployment id for remove individuals with limited data ------ ##
(summary_by_species <- spec %>% group_by(Species) %>% summarise(n_individuals = n_distinct(Animal.ID), n_days_total = n_distinct(dayLocal), .groups = "drop"))

# Deployment length per individual
deployment_summary <- spec %>%
  group_by(Species, Deployment_ID) %>%
  summarise(
    n_fixes     = n(),
    date_start  = min(dayLocal),
    date_end    = max(dayLocal),
    n_days      = as.numeric(difftime(max(dayLocal), min(dayLocal), units = "days")),
    .groups     = "drop"
  )
print(deployment_summary)

spec <- spec %>% filter(!Animal.ID == "727051A") %>% mutate(dayUTC = date(DateUTC))
length(unique(spec$Animal.ID))

# MCP for study areas?
mcps <- st_read("E:/GitProjects/GitprojectDroughtSelection/data/buffered_mcps.shp")
mapview::mapview(mcps)

spec_ocelot <- spec %>% filter(Species == "Ocelot")
length(unique(spec_ocelot$Animal.ID))

# Deployment length per individual
deployment_summary <- spec_ocelot %>%
  group_by(Deployment_ID) %>%
  summarise(
    n_fixes     = n(),
    date_start  = min(dayLocal),
    date_end    = max(dayLocal),
    n_days      = as.numeric(difftime(max(dayLocal), min(dayLocal), units = "days")),
    short_deployment = n_days <= 10,
    .groups     = "drop"
  )
print(deployment_summary, n = 100)

relocation_summary <- deployment_summary %>%
  summarise(avg_days = mean(n_days), min_days = min(n_days), max_days = max(n_days), sd_days = sd(n_days))
print(relocation_summary)

summary(spec_ocelot$elapsed_time)

print(spec_ocelot %>% filter(elapsed_time < 30) %>% select(Deployment_ID, DateUTC, elapsed_time), n = 200)
print(spec_ocelot %>% filter(elapsed_time > 360) %>% select(Deployment_ID, DateUTC, elapsed_time), n = 200)

spec_ocelot %>%
  st_drop_geometry() %>%
  group_by(Deployment_ID) %>%
  summarize(median_interval = median(elapsed_time, na.rm = TRUE)) %>%
  filter(median_interval != 60)

# =======================================================================
# Fit ctmm model
# =======================================================================
data.ctmm <- spec_ocelot
  
data.ctmm.input <- data.ctmm %>% select(Deployment_ID, DateUTC, lon, lat, HDOP, x, y) %>% # Downstream UTC often works better with CTMM functions
  rename(timestamp = DateUTC, ID = Deployment_ID, longitude = lon, latitude = lat, GPS.HDOP = HDOP) %>% 
  arrange(ID, timestamp) %>% st_drop_geometry()

# create ctmm telemetry object
data.tele <- as.telemetry(data.ctmm.input, timezone="UTC")
ID <- names(data.tele)

###--- Find and fit best ctmm model for each individual ---###
data.guess <- lapply(data.tele, function(x) ctmm.guess(x, interactive = FALSE)) 
data.fit <- list()

tictoc::tic()
data.fit <- lapply(1:length(data.tele), function(i) {ctmm.select(data.tele[[i]], data.guess[[i]], method = 'pHREML')}) # Fleming et al. 2019
names(data.fit) <- names(data.tele)
tictoc::toc()

saveRDS(data.fit, "output/ctmm_fit_object.rds")

data.fit <- readRDS("output/ctmm_fit_object.rds")

model_summary <- purrr::imap_dfr(data.fit, function(fit, id) {
  s <- summary(fit)
  ci <- s$CI
  # Find position tau row
  tau_pos_row <- grep("^τ\\[position\\]", rownames(ci))
  tau_position = if (length(tau_pos_row) == 1) {
    ci[tau_pos_row, "est"]
  } else {
    NA_real_
  }
  tau_position_unit = if (length(tau_pos_row) == 1) {
    rownames(ci)[tau_pos_row]
  } else {
    NA_character_
  }
  # Convert tau[position] to days
  tau_position_days = case_when(
    grepl("\\(hours\\)", tau_position_unit) ~ tau_position / 24,
    grepl("\\(days\\)", tau_position_unit) ~ tau_position,
    grepl("\\(minutes\\)", tau_position_unit) ~ tau_position / (60 * 24),
    TRUE ~ NA_real_
  )
  # Velocity tau
  tau_vel_row <- grep("^τ\\[velocity\\]", rownames(ci))
  tau_velocity_min = if (length(tau_vel_row) == 1) {
    ci[tau_vel_row, "est"]
  } else {
    NA_real_
  }
  tibble(
    Deployment_ID = id,
    model = s$name,
    area_km2 = ci["area (square kilometers)", "est"],
    tau_position_days = tau_position_days,
    tau_position_unit = tau_position_unit,
    tau_velocity_min = tau_velocity_min,
    speed_km_day = if ("speed (kilometers/day)" %in% rownames(ci))
      ci["speed (kilometers/day)", "est"]
    else NA_real_,
    diffusion_km2_day = if ("diffusion (square kilometers/day)" %in% rownames(ci))
      ci["diffusion (square kilometers/day)", "est"]
    else NA_real_
  )
})

model_summary <- model_summary %>%
  left_join(
    deployment_summary %>%
      st_drop_geometry() %>%
      select(Deployment_ID, n_days),
    by = "Deployment_ID"
  ) %>%
  mutate(
    tau_fraction = tau_position_days / n_days
  )

print(model_summary, n = 100)

# =================================================================================
# Check EO34M_2: extreme outlier 
# =================================================================================
summary(data.fit$EO34M_2)

eo34m2 <- spec_ocelot %>%
  filter(Deployment_ID == "EO34M_2")

mapview::mapview(eo34m2)
ggplot(eo34m2) +
  geom_sf(size = 0.5) +
  theme_minimal() +
  labs(title = "EO34M_2 movement track")

eo34m2 %>% ggplot() + geom_sf(aes(color = DateLocal), size = 1) +
  scale_fill_viridis_b() +
  theme_minimal(base_size = 20) + labs(title = "EO34M_2 — Movement through time", color = "Date")

# =================================================================================
# Plot variograms
# =================================================================================
par(mfrow = c(3,3))
for(i in names(data.tele)){
  vg <- variogram(data.tele[[i]])
  ctmm::plot(vg, data.guess[[i]], main = i)
}
par(mfrow = c(1,1))


# =================================================================================
# Section to create UD occurrences for all individuals and get size of polygons
# =================================================================================

data_ud <- list()

tictoc::tic()
data_ud <- lapply(seq_along(data.tele), function(i) {
  occurrence(data.tele[[i]], data.fit[[i]])
})
names(data_ud) <- names(data.tele)
tictoc::toc()

saveRDS(data_ud, "output/ctmm_occurrence_uds.rds")

get_ud_polygons <- function(ud, levels = c(0.95, 0.75, 0.50, 0.25, 0.10)) {
  purrr::map_dfr(levels, function(lvl) {
    x <- as.sf(ud, level.UD = lvl)
    x %>% sf::st_cast("POLYGON") %>%
      dplyr::mutate(level = lvl, area_km2 = as.numeric(sf::st_area(geometry)) / 1e6)
  })
}

ud_polygons <- purrr::imap(data_ud, ~get_ud_polygons(.x) %>% mutate(Deployment_ID = .y))

saveRDS(ud_polygons, "output/ctmm_occurrence_ud_polygons_sf.rds")

# =================================================================================
# Interpolate missing steps
# =================================================================================

####--- Predict locations ---###
join.comp <- c() 
predicted.df <- c()

for(i in 1:length(data.fit)){
  ## Get predictions
  pred.temp <- predict(data.fit[[i]], data = data.tele[[i]], complete = TRUE)
  ## Convert to data frame from telemetry object ##
  temp1 <- do.call(cbind.data.frame, pred.temp)
  temp1$ID <- pred.temp@info$identity
  
  ## Convert to data frame from telemetry object ##
  data.temp <- do.call(cbind.data.frame, data.tele[[i]])
  data.temp$ID <- data.tele[[i]]@info$identity
  
  ## Join observed data to predicted data ##
  join.temp <- inner_join(data.temp, temp1, by = "timestamp")
  anti.temp <- anti_join(temp1, data.temp, by = "timestamp")
  
  ## Get difference from observed and predicted from same timestamp ##
  temp.join.comp <- join.temp %>% mutate(xdiff = x.x - x.y, ydiff = y.x - y.y) %>% select(ID.x, xdiff, ydiff)
  join.comp <- rbind(join.comp, temp.join.comp)
  
  ## Format data to output into single data frame for all IDs ##
  join.temp <- join.temp %>% select(ID.x, timestamp, x.x, y.x, longitude.x, latitude.x) %>% rename(ID = ID.x, x = x.x, y = y.x, longitude = longitude.x, latitude = latitude.x) %>% mutate(data = "Observed")
  anti.temp <- anti.temp %>% select(ID, timestamp, x, y, longitude, latitude) %>% mutate(data = "Predicted")
  temp.combined <- rbind(join.temp, anti.temp)
  
  predicted.df <- rbind(predicted.df, temp.combined)
  
}

gap_threshold_min <- 360  # 6 hours, in minutes (matches your elapsed_time units)

predicted_df_filtered <- predicted.df %>%
  arrange(ID, timestamp) %>%
  group_by(ID) %>%
  mutate(
    # time gap to the nearest OBSERVED point on either side determines whether 
    # a given point (real or interpolated) sits within a "trustworthy" span
    time_to_prev = as.numeric(difftime(timestamp, lag(timestamp), units = "mins")),
    time_to_next = as.numeric(difftime(lead(timestamp), timestamp, units = "mins"))
  ) %>%
  ungroup()

# get observed-only gap structure first
observed_gaps <- predicted_df_filtered %>%
  filter(data == "Observed") %>%
  group_by(ID) %>%
  arrange(timestamp) %>%
  mutate(
    gap_start = timestamp,
    gap_end = lead(timestamp),
    gap_length = as.numeric(difftime(gap_end, gap_start, units = "mins"))
  ) %>%
  filter(!is.na(gap_length), gap_length > gap_threshold_min) %>%
  select(ID, gap_start, gap_end) %>%
  ungroup()


# flag and drop any Predicted point that falls inside one of those oversized gaps
predicted_df_clean <- predicted_df_filtered %>%
  left_join(observed_gaps, by = "ID", relationship = "many-to-many") %>%
  mutate(
    in_oversized_gap = data == "Predicted" & !is.na(gap_start) & 
      timestamp > gap_start & timestamp < gap_end
  ) %>%
  group_by(ID, timestamp) %>%
  summarize(
    across(-c(gap_start, gap_end, in_oversized_gap), first),
    drop_point = any(in_oversized_gap),
    .groups = "drop"
  ) %>%
  filter(!drop_point) %>%
  select(-drop_point)
length(unique(predicted_df_clean$ID))

predicted_df_clean <- predicted_df_clean %>%
  arrange(ID, timestamp) %>%
  group_by(ID) %>%
  mutate(
    time_gap = as.numeric(difftime(timestamp, lag(timestamp), units = "mins")),
    new_burst = is.na(time_gap) | time_gap > gap_threshold_min,
    burst = paste0(ID, "_", cumsum(new_burst))
  ) %>%
  ungroup()
length(unique(predicted_df_clean$burst))

predicted_df_clean %>%
  distinct(ID, burst) %>%
  count(ID) %>%
  arrange(desc(n))

burst_lengths <- predicted_df_clean %>%
  count(burst)
summary(burst_lengths$n)
hist(burst_lengths$n, breaks = 50)

# check how much data a given threshold would drop
sapply(c(5, 10, 15, 20), function(thresh) {
  keep <- burst_lengths$burst[burst_lengths$n >= thresh]
  sum(predicted_df_clean$burst %in% keep) / nrow(predicted_df_clean)
})

valid_bursts <- burst_lengths %>%
  filter(n >= 20) %>%
  pull(burst)

hmm_ready_df <- predicted_df_clean %>%
  filter(burst %in% valid_bursts)

hmm_input <- resample_to_interval(hmm_ready_df, time_col = "timestamp", hours = 1, tolerance = minutes(9))
nrow(hmm_ready_df) - nrow(hmm_input)  # should land back near the small expected number

# =================================================================================
# Clean up intervals
# =================================================================================
resample_to_interval <- function(df, time_col, hours = 1, tolerance = minutes(9)) {
  unit_label <- paste(hours, "hours")
  
  df %>%
    mutate(
      # Round timestamp to nearest specified interval
      nearest_time = round_date(.data[[time_col]], unit = unit_label),
      # Time difference from nearest interval
      diff = abs(difftime(.data[[time_col]], nearest_time, units = "mins"))
    ) %>%
    # Keep observations close enough to target time
    filter(diff <= as.numeric(tolerance, units = "mins")) %>%
    select(-nearest_time, -diff)
}

hmm_input <- resample_to_interval(hmm_ready_df, time_col = "timestamp", hours = 1, tolerance = minutes(11))
nrow(hmm_ready_df) - nrow(hmm_input)  # should land back near the small expected number

hmm_ready_df <- hmm_ready_df %>% rename(Deployment_ID = ID)
hmm_input <- hmm_input %>% rename(Deployment_ID = ID)

removed <- anti_join(hmm_ready_df %>% st_drop_geometry(), 
                      hmm_input %>% st_drop_geometry(), 
                      by = c("Deployment_ID", "timestamp"))  # adjust join keys to your actual column names

removed %>% count(minute(timestamp))

saveRDS(hmm_input, "output/hmm_input.rds")

# =================================================================================
# Some checks on the prediction outputs
# =================================================================================

# for every gap in this table, confirm a burst change occurred across it
gap_check_full <- spec_ocelot %>%
  st_drop_geometry() %>%
  filter(elapsed_time > 360) %>%
  select(Deployment_ID, DateUTC, elapsed_time) %>%
  rename(ID = Deployment_ID, gap_end_time = DateUTC)

verify_cuts <- gap_check_full %>%
  rename(check_id = ID) %>%
  rowwise() %>%
  mutate(
    burst_before = predicted_df_clean %>%
      filter(ID == check_id, timestamp < gap_end_time) %>%
      slice_max(timestamp, n = 1, with_ties = FALSE) %>%
      pull(burst) %>% {if(length(.) == 0) NA_character_ else .},
    burst_after = predicted_df_clean %>%
      filter(ID == check_id, timestamp >= gap_end_time) %>%
      slice_min(timestamp, n = 1, with_ties = FALSE) %>%
      pull(burst) %>% {if(length(.) == 0) NA_character_ else .}
  ) %>%
  ungroup() %>%
  mutate(cut_confirmed = burst_before != burst_after)

# should be 0 rows if every long gap correctly produced a burst break
verify_cuts %>% filter(!cut_confirmed | is.na(cut_confirmed))

# also worth a positive check - confirm verify_cuts itself has sensible values, not all NA
verify_cuts %>% count(is.na(burst_before), is.na(burst_after))
summary(verify_cuts$cut_confirmed)

check_id <- "OF274_1"
check_gap_start <- as.POSIXct("2025-09-11 14:00:00", tz = "UTC")  # from your earlier gap table

predicted_df_clean %>%
  filter(ID == check_id, timestamp >= check_gap_start - 60*60*2, timestamp <= check_gap_start + 60*60*8) %>%
  select(timestamp, data, burst)

predicted_df_clean %>%
  arrange(ID, timestamp) %>%
  group_by(ID) %>%
  mutate(gap_check = as.numeric(difftime(timestamp, lag(timestamp), units = "mins"))) %>%
  filter(data == "Predicted") %>%
  pull(gap_check) %>%
  summary()

# confirm
length(unique(hmm_ready_df$ID))       # should be 30
length(unique(hmm_ready_df$burst))    # should be a bit under 246, minus OF274_1's bursts and any short fragments
nrow(hmm_ready_df) / nrow(predicted_df_clean)  # should be ~0.995ish, minus whatever OF274_1's rows accounted for



of274_data <- predicted_df_clean %>% filter(ID == "OF274_1")

ggplot(of274_data, aes(x = timestamp, y = 1, color = data)) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("Observed" = "#0083b3", "Predicted" = "#f78d2d")) +
  theme_minimal(base_size = 14) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank()) +
  labs(title = "OF274_1 - locations over time", x = "Date", color = NULL)

table(of274_data$data)

of274_data %>%
  filter(data == "Observed") %>%
  arrange(timestamp) %>%
  mutate(gap_to_next = as.numeric(difftime(lead(timestamp), timestamp, units = "hours"))) %>%
  filter(!is.na(gap_to_next)) %>%
  ggplot(aes(x = timestamp, y = gap_to_next)) +
  geom_point() +
  geom_hline(yintercept = 6, linetype = "dashed", color = "red") +
  theme_minimal(base_size = 14) +
  labs(title = "OF274_1 - gap length between fixes over time", x = "Date", y = "Gap (hours)")

of274_sf <- of274_data %>%
  filter(!is.na(x), !is.na(y)) %>%  
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)  # adjust CRS to match whatever x/y actually are

ggplot() +
  geom_sf(data = of274_sf, aes(color = data), size = 1) +
  scale_color_manual(values = c("Observed" = "#0083b3", "Predicted" = "#f78d2d")) +
  theme_minimal()

mapview::mapview(of274_sf, zcol = "data")
