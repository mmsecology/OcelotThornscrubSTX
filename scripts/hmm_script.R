library(momentuHMM); library(dplyr)
library(sf); library(lubridate)
library(ggplot2)

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


# quick check - distribution of step lengths and turning angles
summary(prep_data$step)
summary(prep_data$angle)
hist(prep_data$step, breaks = 100)

# NA steps should roughly equal number of bursts (one NA per burst start)
length(unique(prep_data$ID))
sum(is.na(prep_data$step))

# NA angles should be roughly 2x that (first two points of each burst)
sum(is.na(prep_data$angle))

zero_step_count <- prep_data %>% filter(step == 0) %>% nrow()
zero_step_count

n_bursts <- length(unique(prep_data$ID))
n_bursts  # expect ~2x this contributing to angle NAs from burst starts alone
zero_step_count  # additional angle NAs from zero-length steps

###--- Descriptive stats to inform Par0 ---###
summary(prep_data$step); mean(prep_data$step, na.rm = TRUE); sd(prep_data$step, na.rm = TRUE)
hist(prep_data$step, breaks = 100, xlab = "step length", main = "")
(length(which(prep_data$step == 0)) / nrow(prep_data))  # zero-step proportion
hist(prep_data$angle, breaks = seq(-pi, pi, length = 15), xlab = "angle", main = "")
acf(prep_data$step[!is.na(prep_data$step)], lag.max = 100)

###--- Model setup ---###
stateNames <- c("encamped", "traveling")
dists <- list(step = "gamma", angle = "vm")

# Means/SDs based on the summary above 
Par0_null <- list(step = c(164, 233, 164, 233, 0.01, 0.0005), angle = c(0.1, 0.8))

DM_null_ZeroMass <- list(
  step  = list(mean = ~1, sd = ~1, zeromass = ~1),
  angle = list(concentration = ~1)
)

###--- Fit ---###
tictoc::tic()
m1 <- fitHMM(
  data = prep_data, nbStates = 2, dist = dists, Par0 = Par0_null,
  estAngleMean = list(angle = FALSE), stateNames = stateNames,
  DM = DM_null_ZeroMass, retryFits = 50
)
tictoc::toc()

###--- Fit diagnostics ---###

# Convergence code (should be 1; anything else = check optimizer output)
m1$mod$code

# Was the best fit found more than once across retries? 
m1$mle_run <- m1$mod$minimum
print(m1$mod$minimum)

# Parameter estimates + CIs
print(m1)
ci_est <- CIbeta(m1, alpha = 0.95)
print(ci_est)

# State assignment via Viterbi — check proportions
states <- viterbi(m1)
table(states) / length(states)

# Pseudo-residuals — check for approx N(0,1) and no autocorrelation
pr <- pseudoRes(m1)
qqnorm(pr$stepRes); qqline(pr$stepRes)
qqnorm(pr$angleRes); qqline(pr$angleRes)
acf(pr$stepRes[!is.na(pr$stepRes)], lag.max = 100)
acf(pr$angleRes[!is.na(pr$angleRes)], lag.max = 100)

# Visual check — state-colored step/angle histograms + example tracks
plot(m1, plotCI = TRUE, ask = FALSE)

# Stationary state distribution 
stationary(m1)


workBounds <- list(
  step = matrix(c(
    -5, 10,    # mean_1 (log scale): natural range ~0.007 to ~22,000
    -5, 10,    # mean_2
    -5, 10,    # sd_1
    -5, 10,    # sd_2
    -10, 10,   # zeromass_1 (logit scale): natural range ~0.00005 to ~0.99995
    -10, 10    # zeromass_2
  ), ncol = 2, byrow = TRUE),
  angle = matrix(c(
    -10, 10,   # concentration_1 (log scale)
    -10, 10    # concentration_2
  ), ncol = 2, byrow = TRUE)
)

m1_test <- fitHMM(
  data = prep_data, nbStates = 2, dist = dists, Par0 = Par0_null,
  estAngleMean = list(angle = FALSE), stateNames = stateNames,
  DM = DM_null_ZeroMass, workBounds = workBounds, retryFits = 0
)

print(m1_test)
table(viterbi(m1_test)) / nrow(prep_data)


Par0_null_v2 <- list(step = c(60, 150, 60, 150, 0.005, 0.005), angle = c(0.3, 0.8))

null_fit <- function(data, iter, DM) {
  allm <- list(); parm_tested <- list()
  loglik <- numeric(iter); code <- numeric(iter)
  
  for (i in 1:iter) {
    stepMean0 <- runif(2, min = c(30, 100), max = c(100, 300))
    stepSD0   <- runif(2, min = c(30, 100), max = c(100, 300))
    angleCon0 <- runif(2, min = c(0.1, 0.5), max = c(0.5, 1))
    par0 <- list(step = c(sort(stepMean0), sort(stepSD0), 0.005, 0.005), angle = angleCon0)
    parm_tested[[i]] <- par0
    
    fit <- tryCatch(
      fitHMM(data = data, nbStates = 2, dist = dists, Par0 = par0,
             estAngleMean = list(angle = FALSE), stateNames = stateNames, DM = DM),
      error = function(e) NULL
    )
    allm[[i]] <- fit
    loglik[i] <- if (!is.null(fit)) fit$mod$minimum else NA
    code[i]   <- if (!is.null(fit)) fit$mod$code else NA
  }
  list(models = allm, par0 = parm_tested, loglik = loglik, code = code)
}

result <- null_fit(prep_data, iter = 10, DM = DM_null_ZeroMass)

table(result$code)
sort(result$loglik)
sapply(result$models, function(m) if (!is.null(m)) table(viterbi(m)) / nrow(prep_data) else NA)

Par0_null_v2 <- list(step = c(80, 250, 80, 250, 0.01, 0.005), angle = c(0.1, 0.8))

null_fit <- function(data, iter, DM) {
  allm <- list(); parm_tested <- list()
  loglik <- numeric(iter); code <- numeric(iter)
  
  for (i in 1:iter) {
    m1 <- runif(1, 40, 120)    # encamped mean, below pooled mean/median-ish
    m2 <- runif(1, 180, 320)   # traveling mean, around/above pooled mean
    par0 <- list(step = c(m1, m2, m1, m2, 0.01, 0.005), angle = c(0.1, 0.8))
    parm_tested[[i]] <- par0
    
    fit <- tryCatch(
      fitHMM(data = data, nbStates = 2, dist = dists, Par0 = par0,
             estAngleMean = list(angle = FALSE), stateNames = stateNames, DM = DM),
      error = function(e) NULL
    )
    allm[[i]] <- fit
    loglik[i] <- if (!is.null(fit)) fit$mod$minimum else NA
    code[i]   <- if (!is.null(fit)) fit$mod$code else NA
  }
  list(models = allm, par0 = parm_tested, loglik = loglik, code = code)
}

result <- null_fit(prep_data, iter = 10, DM = DM_null_ZeroMass)
table(result$code)
sort(result$loglik)
sapply(result$models, function(m) if (!is.null(m)) table(viterbi(m)) / nrow(prep_data) else NA)

zero_rows <- prep_data %>% filter(step == 0)
table(zero_rows$ID)  #

set.seed(1)
prep_data$step[prep_data$step == 0 & !is.na(prep_data$step)] <- 
  runif(sum(prep_data$step == 0, na.rm = TRUE), min = 0.1, max = 2)

DM_null_noZero <- list(
  step  = list(mean = ~1, sd = ~1),
  angle = list(concentration = ~1)
)

Par0_null_v3 <- list(step = c(80, 250, 80, 250), angle = c(0.1, 0.8))

m1_test <- fitHMM(
  data = prep_data, nbStates = 2, dist = dists, Par0 = Par0_null_v3,
  estAngleMean = list(angle = FALSE), stateNames = stateNames,
  DM = DM_null_noZero, retryFits = 10
)
table(viterbi(m1_test)) / nrow(prep_data)
print(m1_test)

burst_lengths <- prep_data %>% group_by(ID) %>% summarise(n = n())
summary(burst_lengths$n)
hist(burst_lengths$n, breaks = 30)

# Does every burst actually contain at least some longer/traveling-type steps,
# or do some (many?) bursts consist entirely of short encamped-type steps?
burst_step_summary <- prep_data %>% group_by(ID) %>% 
  summarise(max_step = max(step, na.rm = TRUE), mean_step = mean(step, na.rm = TRUE))
summary(burst_step_summary$max_step)
sum(burst_step_summary$max_step < 100)  # how many bursts never show a "long" step at all?



workBounds <- list(
  step = matrix(c(
    -3, 8,     # mean_1 (log scale): natural range ~0.05 to ~2981
    -3, 8,     # mean_2
    -3, 8,     # sd_1
    -3, 8      # sd_2
  ), ncol = 2, byrow = TRUE),
  angle = matrix(c(
    -5, 5,     # concentration_1
    -5, 5      # concentration_2
  ), ncol = 2, byrow = TRUE)
)

Par0_null_v3 <- list(step = c(80, 250, 80, 250), angle = c(0.1, 0.8))

m1_bounded <- fitHMM(
  data = prep_data, nbStates = 2, dist = dists, Par0 = Par0_null_v3,
  estAngleMean = list(angle = FALSE), stateNames = stateNames,
  DM = DM_null_noZero, workBounds = workBounds,
  retryFits = 10, stepmax = 3
)

print(m1_bounded)
table(viterbi(m1_bounded)) / nrow(prep_data)


one_id <- burst_step_summary %>% arrange(desc(max_step)) %>% slice(5) %>% pull(ID)  # pick one, not the single most extreme outlier
single_data <- prep_data %>% filter(ID == one_id)

m1_single <- fitHMM(
  data = single_data, nbStates = 2, dist = dists, Par0 = Par0_null_v3,
  estAngleMean = list(angle = FALSE), stateNames = stateNames,
  DM = DM_null_noZero, retryFits = 10
)
print(m1_single)
table(viterbi(m1_single)) / nrow(single_data)

# log-scale density is often where step-length bimodality becomes visually obvious
prep_data %>% filter(!is.na(step), step > 0) %>%
  ggplot(aes(x = log(step))) + geom_density()

# or classic histogram check
hist(log(prep_data$step[prep_data$step > 0]), breaks = 100)

dup_check <- prep_data %>% group_by(ID) %>% 
  mutate(dup = x == lag(x) & y == lag(y)) %>% 
  ungroup() %>% summarise(n_dup = sum(dup, na.rm = TRUE))
dup_check

gap_check <- prep_data %>% group_by(ID) %>% 
  arrange(timestamp) %>% 
  mutate(dt = as.numeric(difftime(timestamp, lag(timestamp), units = "hours"))) %>%
  ungroup()
summary(gap_check$dt)  # should be tightly clustered near 1 hour if truly hourly/regularized

prep_data$step_km <- prep_data$step / 1000
prep_data_km <- prep_data %>% rename(step_orig = step) %>% rename(step = step_km)

Par0_km <- list(step = c(0.015, 0.3, 0.015, 0.3), angle = c(0.1, 0.8))

m1_km <- fitHMM(
  data = prep_data_km, nbStates = 2, dist = dists, Par0 = Par0_km,
  estAngleMean = list(angle = FALSE), stateNames = stateNames,
  DM = DM_null_noZero, retryFits = 10
)
print(m1_km)
table(viterbi(m1_km)) / nrow(prep_data_km)


# --- Step-only (already sketched, using the km-rescaled data) ---
m1_step_only <- fitHMM(
  data = prep_data_km, nbStates = 2, dist = list(step = "gamma"),
  Par0 = list(step = Par0_km$step), stateNames = stateNames,
  DM = list(step = DM_null_noZero$step), retryFits = 10
)
print(m1_step_only)
table(viterbi(m1_step_only)) / nrow(prep_data_km)

# --- Angle-only ---
m1_angle_only <- fitHMM(
  data = prep_data_km, nbStates = 2, dist = list(angle = "vm"),
  Par0 = list(angle = c(0.1, 0.8)), estAngleMean = list(angle = FALSE),
  stateNames = stateNames, DM = list(angle = DM_null_noZero$angle), retryFits = 10
)
print(m1_angle_only)
table(viterbi(m1_angle_only)) / nrow(prep_data_km)


Par0_km_zm <- list(step = c(0.015, 0.3, 0.015, 0.3, 0.01, 0.001), angle = c(0.1, 0.8))

m1_km_zm <- fitHMM(
  data = prep_data_km, nbStates = 2, dist = dists, Par0 = Par0_km_zm,
  estAngleMean = list(angle = FALSE), stateNames = stateNames,
  retryFits = 2
)
print(m1_km_zm)
table(viterbi(m1_km_zm)) / nrow(prep_data_km)




## USED WITH KM INSTEAD 

hmm_data_km <- hmm_data %>% mutate(x = x / 1000, y = y / 1000)  # rescale coordinates to km, before prepData

prep_data_km <- prepData(data = hmm_data_km, type = "UTM", coordNames = c("x", "y"))

class(prep_data_km)  # confirm momentuHMMData class is intact
summary(prep_data_km$step)  # sanity check - should now be in km directly

Par0_km_zm <- list(step = c(0.015, 0.3, 0.015, 0.3, 0.01, 0.0005), angle = c(0.1, 0.8))

m1_km_zm <- fitHMM(
  data = prep_data_km, nbStates = 2, dist = dists, Par0 = Par0_km_zm,
  estAngleMean = list(angle = FALSE), stateNames = stateNames,
  DM = DM_null_ZeroMass, retryFits = 2
)
print(m1_km_zm)
table(viterbi(m1_km_zm)) / nrow(prep_data_km)


result <- null_fit(prep_data_km, iter = 10, DM = DM_null_ZeroMass)
table(result$code)
sort(result$loglik)
sapply(result$models, function(m) if (!is.null(m)) table(viterbi(m)) / nrow(prep_data) else NA)

best_idx <- which(result$loglik < 120000)  # should grab the 4 good ones
best_idx

for (i in best_idx) {
  cat("--- Model", i, "---\n")
  print(result$models[[i]])
  print(table(viterbi(result$models[[i]])) / nrow(prep_data_km))
  cat("\n")
}

set.seed(1)
zero_idx <- which(prep_data_km$step == 0 & !is.na(prep_data_km$step))
length(zero_idx)  # 168

# jitter within known ~10m GPS error, converted to km: 0.001-0.01 km
prep_data_km$step[zero_idx] <- runif(length(zero_idx), min = 0.001, max = 0.01)

DM_null_noZero <- list(step = list(mean = ~1, sd = ~1), angle = list(concentration = ~1))
Par0_jit <- list(step = c(0.015, 0.3, 0.015, 0.3), angle = c(0.1, 0.8))

m1_jit <- fitHMM(
  data = prep_data_km, nbStates = 2, dist = dists, Par0 = Par0_jit,
  estAngleMean = list(angle = FALSE), stateNames = stateNames,
  DM = DM_null_noZero, retryFits = 10
)
print(m1_jit)
table(viterbi(m1_jit)) / nrow(prep_data_km)


# how many/what fraction of steps are very small (near GPS-noise scale), not just exact zeros
sum(prep_data_km$step < 0.003, na.rm = TRUE)   # < 3m
sum(prep_data_km$step < 0.01, na.rm = TRUE)    # < 10m
sum(prep_data_km$step < 0.02, na.rm = TRUE)    # < 20m

# and look at just this low tail on its own histogram
hist(prep_data_km$step[prep_data_km$step < 0.05], breaks = 100)

prep_data_km$step_adj <- prep_data_km$step + 0.010  # +10m GPS error floor, km

Par0_adj <- list(step = c(0.025, 0.31, 0.025, 0.31), angle = c(0.1, 0.8))

m1_adj <- fitHMM(
  data = prep_data_km %>% mutate(step = step_adj), nbStates = 2, dist = dists,
  Par0 = Par0_adj, estAngleMean = list(angle = FALSE), stateNames = stateNames,
  DM = DM_null_noZero, retryFits = 10
)
print(m1_adj)
table(viterbi(m1_adj)) / nrow(prep_data_km)



# Convergence code (should be 1; anything else = check optimizer output)
m1_adj$mod$code

# Was the best fit found more than once across retries? 
m1_adj$mle_run <- m1_adj$mod$minimum
print(m1_adj$mod$minimum)

# Parameter estimates + CIs
print(m1_adj)
ci_est <- CIbeta(m1_adj, alpha = 0.95)
print(ci_est)

# State assignment via Viterbi — check proportions
states <- viterbi(m1_adj)
table(states) / length(states)

# Pseudo-residuals — check for approx N(0,1) and no autocorrelation
pr <- pseudoRes(m1_adj)
qqnorm(pr$stepRes); qqline(pr$stepRes)
qqnorm(pr$angleRes); qqline(pr$angleRes)
acf(pr$stepRes[!is.na(pr$stepRes)], lag.max = 100)
acf(pr$angleRes[!is.na(pr$angleRes)], lag.max = 100)

# Visual check — state-colored step/angle histograms + example tracks
plot(m1_adj, plotCI = TRUE, ask = FALSE)

# Stationary state distribution 
stationary(m1_adj)


str(pr)
length(pr$stepRes)
sum(is.na(pr$stepRes))
sum(is.infinite(pr$stepRes))
summary(pr$stepRes)

step_res_clean <- pr$stepRes[is.finite(pr$stepRes)]
angle_res_clean <- pr$angleRes[is.finite(pr$angleRes)]

qqnorm(step_res_clean); qqline(step_res_clean)
qqnorm(angle_res_clean); qqline(angle_res_clean)
acf(step_res_clean, lag.max = 100)
acf(angle_res_clean, lag.max = 100)
step_res_clean <- pr$stepRes[is.finite(pr$stepRes)]
length(step_res_clean)   # sanity check - should be ~131035 (131241 - 205 NA - 1 Inf)
summary(step_res_clean)  # confirm this actually has values, no errors
qqnorm(step_res_clean)

stateNames <- c("stationary", "encamped", "traveling")
dists <- list(step = "gamma", angle = "vm")

Par0_3state <- list(
  step  = c(0.003, 0.015, 0.30,   # means (km): stationary, encamped, traveling
            0.003, 0.015, 0.30),  # sds (Michelot: mean = sd per state)
  angle = c(0.05, 0.2, 0.8)       # concentration: near-uniform, weak, directed
)

DM_3state <- list(
  step  = list(mean = ~1, sd = ~1),
  angle = list(concentration = ~1)
)

m3 <- fitHMM(
  data = prep_data_km, nbStates = 3, dist = dists, Par0 = Par0_3state,
  estAngleMean = list(angle = FALSE), stateNames = stateNames,
  DM = DM_3state, retryFits = 10
)
print(m3)
table(viterbi(m3)) / nrow(prep_data_km)