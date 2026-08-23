library(momentuHMM)

# Input created in relocation_data_to_ctmm.R
# Regularizes hourly fixes into bursts if longer than 6 hours
# Also only predicts no more than 6 hours worth of fixes
# Only keep bursts if more than 20 relocations
hmm_input <- readRDS("output/hmm_input.rds") %>% rename(BurstID = burst)

# prepData needs: ID (or burst) column, x, y coordinates
hmm_data <- hmm_input %>%   # or hmm_input if you didn't switch to calcBurst
  st_drop_geometry() %>%
  select(BurstID, x, y, timestamp) %>%
  rename(ID = BurstID)   # prepData expects the grouping column to be named "ID"

hmm_data <- as.data.frame(hmm_data)

prep_data <- prepData(
  data = hmm_data,
  type = "UTM",              # projected coordinates, not lon/lat
  coordNames = c("x", "y")
)

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
