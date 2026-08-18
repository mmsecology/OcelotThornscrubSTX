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
tictoc::toc() # 3-4 minutes


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

get_ud_polygons <- function(ud, levels = c(0.95, 0.75, 0.50, 0.25, 0.10)) {
  purrr::map_dfr(levels, function(lvl) {
    x <- as.sf(ud, level.UD = lvl)
    x %>% sf::st_cast("POLYGON") %>%
      dplyr::mutate(level = lvl, area_km2 = as.numeric(sf::st_area(geometry)) / 1e6)
  })
}

ud_polygons <- purrr::imap(
  data_ud,
  ~get_ud_polygons(.x) %>%
    mutate(Deployment_ID = .y)
)

ud_polygons_df <- ud_polygons |> list_rbind()

rank_area <- ud_polygons_df %>%
  group_by(Deployment_ID, level) %>%
  arrange(desc(area_km2), .by_group = TRUE) %>%
  mutate(
    rank = row_number()
  ) %>%
  ungroup()

ggplot(rank_area, aes(rank, area_km2, group = interaction(Deployment_ID, level))) +
  geom_line(alpha = 0.2) +
  scale_y_log10() +
  facet_wrap(~level, scales = "free") +
  labs(
    x = "Component rank",
    y = "Component area (km²)"
  )

rank_area_prop <- ud_polygons_df %>%
  group_by(Deployment_ID, level) %>%
  arrange(desc(area_km2), .by_group = TRUE) %>%
  mutate(
    rank = row_number(),
    area_prop = area_km2 / sum(area_km2)
  ) %>%
  ungroup()

ggplot(rank_area_prop, aes(rank, area_prop,
                           group = interaction(Deployment_ID, level))) +
  geom_line(alpha = 0.2) +
  scale_y_log10() +
  facet_wrap(~level, scales = "free_x") +
  labs(
    x = "Component rank",
    y = "Component area / total contour area"
  )


rank_area <- ud_polygons_df %>%
  group_by(Deployment_ID, level) %>%
  arrange(desc(area_km2), .by_group = TRUE) %>%
  mutate(
    rank = row_number(),
    n_components = n(),
    relative_rank = rank / n_components,
    area_prop = area_km2 / sum(area_km2)
  ) %>%
  ungroup()

ggplot(rank_area, aes(relative_rank, area_prop, group = interaction(Deployment_ID, level))) +
  geom_line(alpha = 0.2) +
  scale_y_log10() +
  facet_wrap(~level, scales = "free") +
  labs(
    x = "Component rank",
    y = "Component area (km²)"
  ) + theme_bw(base_size = 20)

ud_polygons_df %>%
  group_by(Deployment_ID, level) %>%
  summarise(
    n_components = n(),
    total_area_km2 = sum(area_km2),
    .groups = "drop"
  )

ud_polygons_df <- ud_polygons_df %>%
  group_by(Deployment_ID, level) %>%
  arrange(desc(area_km2), .by_group = TRUE) %>%
  mutate(
    rank = row_number(),
    log_area = log10(area_km2),
    cumulative_area_km2 = cumsum(area_km2),
    cumulative_area_prop = cumulative_area_km2 / sum(area_km2)
  ) %>%
  ungroup()

ud_component_summary <- ud_polygons_df %>%
  group_by(Deployment_ID, level) %>%
  arrange(desc(area_km2), .by_group = TRUE) %>%
  mutate(
    rank = row_number(),
    cumulative_area_prop = cumsum(area_km2) / sum(area_km2)
  ) %>%
  summarise(
    n_components = n(),
    total_area_km2 = sum(area_km2),
    n_for_50_area = which(cumulative_area_prop >= 0.50)[1],
    n_for_75_area = which(cumulative_area_prop >= 0.75)[1],
    n_for_90_area = which(cumulative_area_prop >= 0.90)[1],
    n_for_95_area = which(cumulative_area_prop >= 0.95)[1],
    .groups = "drop"
  )
print(ud_component_summary)

get_breakpoint <- function(dat) {
  
  dat <- dat %>%
    arrange(desc(area_km2)) %>%
    mutate(
      rank = row_number(),
      log_area = log10(area_km2)
    )
  
  if (nrow(dat) < 10) {
    return(tibble(
      breakpoint = NA_real_,
      breakpoint_relative = NA_real_,
      slope_1 = NA_real_,
      slope_2 = NA_real_
    ))
  }
  
  m <- lm(log_area ~ rank, data = dat)
  
  fit <- try(
    segmented(
      m,
      seg.Z = ~rank,
      psi = list(rank = max(3, round(nrow(dat) * 0.05)))
    ),
    silent = TRUE
  )
  
  if (inherits(fit, "try-error")) {
    return(tibble(
      breakpoint = NA_real_,
      breakpoint_relative = NA_real_,
      slope_1 = NA_real_,
      slope_2 = NA_real_
    ))
  }
  
  psi <- summary(fit)$psi
  
  if (is.null(psi) || nrow(psi) == 0) {
    return(tibble(
      breakpoint = NA_real_,
      breakpoint_relative = NA_real_,
      slope_1 = NA_real_,
      slope_2 = NA_real_
    ))
  }
  
  bp <- psi[1, "Est."]
  
  sl <- try(slope(fit)$rank, silent = TRUE)
  
  if (inherits(sl, "try-error") || nrow(sl) < 2) {
    slope_1 <- NA_real_
    slope_2 <- NA_real_
  } else {
    slope_1 <- sl[1, "Est."]
    slope_2 <- sl[2, "Est."]
  }
  
  tibble(
    breakpoint = bp,
    breakpoint_relative = bp / nrow(dat),
    slope_1 = slope_1,
    slope_2 = slope_2
  )
}

breakpoints <- ud_polygons_df %>%
  group_by(Deployment_ID, level) %>%
  group_modify(~get_breakpoint(.x)) %>%
  ungroup()
print(breakpoints, n = 200)

ggplot(
  breakpoints,
  aes(level, breakpoint_relative)
) +
#  geom_boxplot() +
  geom_jitter(width = 0.05, alpha = 0.5)

# =================================================================================
# Test out ctmm occurrence
# =================================================================================


ud_19_sf <- as.sf(ud_19, level.UD = 0.50)



ud_19 <- occurrence(
  data.tele[["EO19M_1"]],
  data.fit[["EO19M_1"]]
)

ud_34 <- occurrence(
  data.tele[["EO34M_2"]],
  data.fit[["EO34M_2"]]
)
par(mfrow = c(1,2))

plot(ud_19, main = "EO19M_1", level.UD = 0.90)
plot(ud_34, main = "EO34M_2")

ud_19_sf <- as.sf(ud_19, level.UD = 0.50)
patches <- terra::patches(ud_19_sf, directions = 8)

ud_19_rows <- ud_19_sf %>% summarise(geometry = st_union(geometry)) %>% st_cast("POLYGON") %>% mutate(area_km2 = as.numeric(st_area(geometry)) / 1e6)
ud_19_rows %>%
  st_drop_geometry() %>%
  arrange(desc(area_km2))
ggplot(
  ud_19_rows,
  aes(x = area_km2)
) +
  geom_histogram() +
  scale_x_log10()
ud_19_rows %>%
  st_drop_geometry() %>%
  arrange(desc(area_km2)) %>%
  mutate(rank = row_number()) %>%
  ggplot(aes(rank, area_km2)) +
  geom_line() +
  scale_y_log10()

print(ud_19_rows)
plot(ud_19_sf)
str(ud_19_sf)
mapview::mapview(ud_19_rows)


summary(ud_19)
summary(ud_34)

ud <- data_ud[["EO19M_1"]]

str(ud_19)
r <- raster(ud_19, DF = "PDF")
res(r)
ext(r)
ncell(r)
terra::plot(r)
summary(r)
