library(amt); library(dplyr); library(lubridate); library(sf)
library(momentuHMM)

# Read in needed data and outputs
hmm_input <- readRDS("output/objects/hmm_input.rds") %>% rename(BurstID = burst)
data_adj <- readRDS("output/objects/hmm_data_adj.rds")
hmm_model_list <- readRDS("output/objects/hmm_model_list.rds")
m3_cosinor <- hmm_model_list[[4]]

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













