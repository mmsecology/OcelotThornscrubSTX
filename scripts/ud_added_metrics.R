library(sf); library(terra); library(ctmm)
library(landscapemetrics); library(dplyr)
library(purrr); library(ggplot2)
library(tidyterra); library(patchwork)
library(tictoc); library(Hmisc)

# ---------------------------------------------------------------
# Simplified landscape metrics run
# ---------------------------------------------------------------

ud_retained_sf <- readRDS("output/ud_retained_sf.rds")
thornscrub_binary <- rast("output/south_texas_thornscrub_binary.tif")
null_uds <- readRDS("output/null.rds")

class_metrics <- c("lsm_c_pland", "lsm_c_area_mn", "lsm_c_shape_mn", "lsm_c_cai_mn",
                    "lsm_c_enn_mn", "lsm_c_clumpy", "lsm_c_ed", "lsm_c_cohesion")

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
    calculate_ud_metrics(ud = null_uds[i, ], landscape = thornscrub_binary, metrics = class_metrics) %>%
      mutate(rep = null_uds$rep[i])
  }
)
saveRDS(null_metrics, "output/null_ud_metrics.rds")

#null_metrics <- readRDS("output/null_ud_metrics.rds")

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
comparison <- readRDS("output/obs_null_comparison_uds.rds")

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


#e.g. "points are per-individual SES (observed relative to n = 999 individual-specific null UDs); black points/ranges are bootstrapped mean ± 95% CI across individuals (n = 28)."
ud_labels <- c("0.5" = "Core-use area (50% UD)", "0.95" = "Overall-use area (95% UD)")

# Human-readable metric labels, in the order you want them to appear
metric_labels <- c(
  "area_mn"  = "Mean patch area",
  "cai_mn"   = "Core area index",
  "clumpy"   = "Aggregation (CLUMPY)",
  "cohesion" = "Patch cohesion",
  "ed"       = "Edge density",
  "enn_mn"   = "Mean nearest-neighbor\ndistance",
  "pland"    = "Percent thornscrub",
  "shape_mn" = "Mean shape complexity"
)

metric_labels_units <- c(
  "area_mn"  = "Mean patch area (ha)",
  "cai_mn"   = "Core area index (%)",
  "clumpy"   = "Aggregation (CLUMPY)",
  "cohesion" = "Patch cohesion (%)",
  "ed"       = "Edge density (m/ha)",
  "enn_mn"   = "Mean nearest-neighbor distance (m)",
  "pland"    = "Percent thornscrub (%)",
  "shape_mn" = "Mean shape complexity"
)

# Apply as an ordered factor so panel order is controlled explicitly
comparison <- comparison %>%
  mutate(metric = factor(metric, levels = names(metric_labels)))

summary_df <- summary_df %>%
  mutate(metric = factor(metric, levels = names(metric_labels)))

ud_comp_fig <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_jitter(
    data = comparison,
    aes(x = metric, y = ses, color = ses > 0),
    width = 0.15, size = 2.5, alpha = 0.45
  ) +
  geom_pointrange(
    data = summary_df,
    aes(x = metric, y = mean_ses, ymin = lower, ymax = upper, fill = sig),
    shape = 23, size = 1, linewidth = 1, color = "black"
  ) +
  facet_wrap(~ UD_level, labeller = as_labeller(ud_labels)) +
  scale_x_discrete(labels = metric_labels) +
  coord_flip() +
  scale_color_manual(values = c("FALSE" = "#002c4c", "TRUE" = "#4c3b00"), guide = "none") +
  scale_fill_manual(
    values = c("positive" = "#4c3b00", "negative" = "#002c4c", "ns" = "grey40"),
    guide = "none") +
  theme_bw(base_size = 25) +
  labs(x = NULL, y = "Standardized effect size")

ud_comp_fig

ggsave("output/fig_ud_comp.pdf", plot = ud_comp_fig)


## -----------------------------------------------------------------------
## 1. POPULATION-LEVEL SUMMARY: raw mean value by metric x UD level
## -----------------------------------------------------------------------
ud_summary <- comparison %>%
  group_by(metric, UD_level) %>%
  dplyr::summarize(boot = list(Hmisc::smean.cl.boot(observed_value)), .groups = "drop") %>%
  mutate(
    mean_val = purrr::map_dbl(boot, ~ .x[["Mean"]]),
    lwr      = purrr::map_dbl(boot, ~ .x[["Lower"]]),
    upr      = purrr::map_dbl(boot, ~ .x[["Upper"]]),
    metric   = factor(metric, levels = names(metric_labels))   # match SES fig panel order
  ) %>%
  select(-boot)

## Plain table for the manuscript / layman's guide:
ud_summary_wide <- ud_summary %>%
  select(metric, UD_level, mean_val, lwr, upr) %>%
  pivot_wider(
    names_from = UD_level,
    values_from = c(mean_val, lwr, upr),
    names_glue = "{.value}_ud{UD_level}"
  )
## ud_summary_wide %>% write.csv("output/tier1_ud_level_summary.csv", row.names = FALSE)

## -----------------------------------------------------------------------
## 2. UD COMPARISON PLOT: one panel per metric, core-use vs. full range
## -----------------------------------------------------------------------
metric_labels_units <- c(
  "area_mn"  = "Mean patch area (ha)",
  "cai_mn"   = "Core area index (%)",
  "clumpy"   = "Aggregation (CLUMPY)",
  "enn_mn"   = "Mean nearest-neighbor distance (m)",
  "pland"    = "Percent thornscrub (%)",
  "shape_mn" = "Mean shape complexity"
)

# Valid theoretical range for the bounded metrics only —
# area_mn and enn_mn have no fixed bound (landscape-dependent), so omitted
bounds_df <- tibble::tribble(
  ~metric,    ~bound,
  "pland",       0,
  "pland",     100,
  "cai_mn",      0,
  "cai_mn",    100,
  "clumpy",     -1,
  "clumpy",      1,
  "cohesion",    0,
  "cohesion",  100,
  "shape_mn",    1
) %>%
  mutate(metric = factor(metric, levels = names(metric_labels)))

plot_ud_dumbbell <- function(df) {
  ggplot(df, aes(x = factor(UD_level, levels = c("0.5", "0.95")), y = mean_val)) +
    geom_hline(data = bounds_df, aes(yintercept = bound),
               linetype = "dotted", color = "grey60") +
    geom_pointrange(aes(ymin = lwr, ymax = upr), size = 1.5) +
    facet_wrap(~metric, scales = "free", ncol = 2, labeller = as_labeller(metric_labels_units)) +
    scale_x_discrete(labels = ud_labels) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 20)
}

plot_ud_dumbbell(ud_summary)
ggsave("output/fig_ud_dumbbell.pdf", plot = plot_ud_dumbbell(ud_summary))