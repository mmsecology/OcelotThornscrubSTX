library(amt); library(dplyr); library(lubridate); library(sf)
library(momentuHMM)

# Read in need data and outputs
hmm_input <- readRDS("output/hmm_input.rds") %>% rename(BurstID = burst)
data_adj <- readRDS("output/objects/hmm_data_adj.rds")
hmm_model_list <- readRDS("output/objects/hmm_model_list.rds")
m3_cosinor <- hmm_model_list[[4]]

# Attach decoded 3-state Viterbi output, real-world coords, and animal ID
data_adj$state <- viterbi(m3_cosinor)
data_adj$state_label <- recode(as.character(data_adj$state),
                                "1" = "local_movement", "2" = "encamped", "3" = "traveling")

animal_map <- hmm_input %>% st_drop_geometry() %>% distinct(BurstID, Deployment_ID)
data_adj <- data_adj %>% left_join(animal_map, by = c("ID" = "BurstID"))

# reattach meters-scale x/y (undo km rescale) + timestamp
data_adj <- data_adj %>%
  left_join(hmm_input %>% st_drop_geometry() %>% select(BurstID, timestamp, x, y),
            by = c("ID" = "BurstID", "timestamp" = "timestamp"))

traveling_only <- data_adj %>% filter(state_label == "traveling")














