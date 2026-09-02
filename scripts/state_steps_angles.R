library(ggplot2); library(dplyr); library(tidyr)

hmm_model_list <- readRDS("output/objects/hmm_model_list.rds")
m3_cosinor <- hmm_model_list[[4]]

# Extract fitted parameters (mean-covariate transition matrix, but emission params
# are covariate-independent here so these apply regardless of hour)
step_par <- m3_cosinor$mle$step   # matrix: mean, sd per state
angle_par <- m3_cosinor$mle$angle # concentration per state

state_labels <- c("1" = "Local movement", "2" = "Encamped/core", "3" = "Traveling")
state_colors <- c("Encamped/core" = "#0072B2", "Local movement" = "#E69F00", "Traveling" = "#D55E00")

# ---- Step length density overlay ----
# Convert gamma mean/sd to shape/rate for dgamma()
gamma_shape <- (step_par["mean", ] / step_par["sd", ])^2
gamma_rate  <- step_par["mean", ] / (step_par["sd", ]^2)

step_seq <- seq(0.001, quantile(data_adj$step, 0.995, na.rm = TRUE), length.out = 500)

step_dens <- expand.grid(step = step_seq, state = 1:3) %>%
  rowwise() %>%
  mutate(density = dgamma(step, shape = gamma_shape[state], rate = gamma_rate[state])) %>%
  ungroup() %>%
  mutate(state_label = recode(as.character(state), !!!state_labels))

p_step <- ggplot() +
  geom_histogram(data = data_adj %>% filter(step <= quantile(step, 0.995, na.rm = TRUE)), 
                  aes(x = step, y = after_stat(density)),
                  bins = 100, fill = "grey85", color = "grey70") +
  geom_line(data = step_dens, aes(x = step, y = density, color = state_label), linewidth = 1) +
  scale_color_manual(values = state_colors, name = "State") +
  coord_cartesian(xlim = c(0, quantile(data_adj$step, 0.99, na.rm = TRUE))) +
  labs(x = "Step length (km)", y = "Density") +
  theme_bw(base_size = 20)
# ---- Turning angle density overlay ----
angle_seq <- seq(-pi, pi, length.out = 500)

angle_dens <- expand.grid(angle = angle_seq, state = 1:3) %>%
  rowwise() %>%
  mutate(density = CircStats::dvm(angle, mu = 0, kappa = angle_par["concentration", state])) %>%
  ungroup() %>%
  mutate(state_label = recode(as.character(state), !!!state_labels))

p_angle <- ggplot() +
  geom_histogram(data = data_adj, aes(x = angle, y = after_stat(density)),
                  bins = 40, fill = "grey85", color = "grey70") +
  geom_line(data = angle_dens, aes(x = angle, y = density, color = state_label), linewidth = 1) +
  scale_color_manual(values = state_colors, name = "State") +
  scale_x_continuous(breaks = c(-pi, -pi/2, 0, pi/2, pi),
                      labels = c("-π", "-π/2", "0", "π/2", "π")) +
  labs(x = "Turning angle (radians)", y = "Density") +
  theme_bw(base_size = 20)

# ---- Combine ----
library(patchwork)
p_step / p_angle
