library(momentuHMM); library(dplyr); library(lubridate)
library(ggplot2); library(sf)

# ---------------------------------------------------------------
# Setup
# ---------------------------------------------------------------

# Input created in relocation_data_to_ctmm.R
# Regularizes hourly fixes into bursts if longer than 6 hours
# Also only predicts no more than 6 hours worth of fixes
# Only keep bursts if more than 20 relocations
hmm_input <- readRDS("output/hmm_input.rds") %>% rename(BurstID = burst)

# prepData needs: ID (or burst) column, x, y coordinates
hmm_data <- hmm_input %>% st_drop_geometry() %>% 
  select(BurstID, x, y, timestamp) %>% rename(ID = BurstID)   

hmm_data <- as.data.frame(hmm_data)

prep_data <- prepData(
  data = hmm_data,
  type = "UTM",              # projected coordinates, not lon/lat
  coordNames = c("x", "y")
)

prep_data$step_km <- prep_data$step / 1000
prep_data_km <- prep_data %>% rename(step_orig = step) %>% rename(step = step_km)

# prep_data_km$step_adj = step + 10m GPS-error floor offset (in km)
prep_data_km$hour <- hour(hmm_input$timestamp)  # confirm row order/ID+timestamp alignment
stopifnot(nrow(prep_data_km) == length(prep_data_km$hour))

prep_data_km$step_adj <- prep_data_km$step + 0.010  # +10m GPS error floor, km
data_adj <- prep_data_km %>% mutate(step = step_adj)  # build once, reuse everywhere below

dists <- list(step = "gamma", angle = "vm")
formula_hour <- ~cosinor(hour, period = 24)

saveRDS(data_adj, "output/objects/hmm_data_adj.rds")
# ---------------------------------------------------------------
# 2-state NULL
# ---------------------------------------------------------------
stateNames2 <- c("encamped", "traveling")
DM_null2 <- list(step = list(mean = ~1, sd = ~1), angle = list(concentration = ~1))
Par0_null2 <- list(step = c(0.025, 0.31, 0.025, 0.31), angle = c(0.1, 0.8))

m2_null <- fitHMM(
  data = data_adj, nbStates = 2, dist = dists, Par0 = Par0_null2,
  estAngleMean = list(angle = FALSE), stateNames = stateNames2,
  DM = DM_null2, retryFits = 10
)
m2_null$mod$code
print(m2_null)
table(viterbi(m2_null)) / nrow(data_adj)

# ---------------------------------------------------------------
# 2-state COSINOR
# ---------------------------------------------------------------
Par0_cos2 <- getPar0(model = m2_null, DM = DM_null2, formula = formula_hour)

m2_cosinor <- fitHMM(
  data = data_adj, nbStates = 2, dist = dists, Par0 = Par0_cos2$Par,
  estAngleMean = list(angle = FALSE), stateNames = stateNames2,
  DM = DM_null2, beta0 = Par0_cos2$beta, formula = formula_hour, retryFits = 10
)
m2_cosinor$mod$code
print(m2_cosinor)
table(viterbi(m2_cosinor)) / nrow(data_adj)

# ---------------------------------------------------------------
# 3-state NULL
# ---------------------------------------------------------------
stateNames3 <- c("state1", "state2", "state3")
DM_null3 <- list(step = list(mean = ~1, sd = ~1), angle = list(concentration = ~1))
Par0_null3 <- list(
  step  = c(0.013, 0.025, 0.31,
            0.013, 0.025, 0.31),
  angle = c(0.05, 0.2, 0.8)
)

m3_null <- fitHMM(
  data = data_adj, nbStates = 3, dist = dists, Par0 = Par0_null3,
  estAngleMean = list(angle = FALSE), stateNames = stateNames3,
  DM = DM_null3, retryFits = 10
)
m3_null$mod$code
print(m3_null)
table(viterbi(m3_null)) / nrow(data_adj)

# ---------------------------------------------------------------
# 3-state COSINOR
# ---------------------------------------------------------------
Par0_cos3 <- getPar0(model = m3_null, DM = DM_null3, formula = formula_hour)

m3_cosinor <- fitHMM(
  data = data_adj, nbStates = 3, dist = dists, Par0 = Par0_cos3$Par,
  estAngleMean = list(angle = FALSE), stateNames = stateNames3,
  DM = DM_null3, beta0 = Par0_cos3$beta, formula = formula_hour, retryFits = 10
)
m3_cosinor$mod$code
print(m3_cosinor)
table(viterbi(m3_cosinor)) / nrow(data_adj)
table(viterbi(m3_null)) / nrow(data_adj)

plot(m3_cosinor)

hmm_model_list <- list(m2_null, m2_cosinor, m3_null, m3_cosinor)
saveRDS(hmm_model_list, "output/hmm_model_list.rds")

# ---------------------------------------------------------------
# AIC comparison
# ---------------------------------------------------------------
AIC(m2_null, m2_cosinor, m3_null, m3_cosinor)

# ---------------------------------------------------------------
# Compare null model to cosinor
# ---------------------------------------------------------------
viterbi_null2 <- viterbi(m2_null)
viterbi_cos2  <- viterbi(m2_cosinor)
mean(viterbi_null2 != viterbi_cos2)
table(null = viterbi_null2, cosinor = viterbi_cos2)

# ---------------------------------------------------------------
# Pseudo-residual ACF
# ---------------------------------------------------------------
pr_null3 <- pseudoRes(m3_null)
pr_cos3  <- pseudoRes(m3_cosinor)
acf(pr_null3$stepRes[is.finite(pr_null2$stepRes)], lag.max = 100, main = "Null 2-state")
acf(pr_cos3$stepRes[is.finite(pr_cos3$stepRes)], lag.max = 100, main = "Cosinor 3-state")

# ---------------------------------------------------------------
# Per-individual sample size check for 3-state (needs animal ID, not BurstID)
# ---------------------------------------------------------------
animal_map <- hmm_input %>% distinct(BurstID, Deployment_ID)
data_adj <- data_adj %>% left_join(animal_map, by = c("ID" = "BurstID"))
state_by_animal <- data_adj %>% mutate(state3 = viterbi(m3_null)) %>%
  group_by(Deployment_ID) %>% count(state3) %>%
  tidyr::pivot_wider(names_from = state3, values_from = n, values_fill = 0)

print(state_by_animal, n = 32)
