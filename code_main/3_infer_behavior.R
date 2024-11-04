#############################################################
# 5_infer_behavior.R
#
# This file compares observed proportion uniqueness to the null
# distribution for each user/taxon to infer behavior type 
# based on the results from script 2
#############################################################


source("revision_thru2023/functions.R")

summaries <- list.files("intermediate/", pattern = "user_summary*.", 
                        full.names = TRUE)
# stop("Update nchar below")
taxa <- substr(summaries, 22, nchar(summaries)-4)

summary_df <- bind_rows(map2_df(summaries, taxa, function(x, y) {
  temp <- read_csv(x)
  temp$taxon <- y
  temp$first_obs = as.Date(temp$first_obs)
  temp
}))

helper <- function(path){
  temp <- read_csv(path)
  temp$first_obs = as.Date(temp$first_obs)
  temp
}

user_values <- bind_rows(lapply(list.files("revision_thru2023/", pattern = 'user_summary_',
                                           full.names = TRUE), helper))
user_dists <- bind_rows(lapply(list.files("revision_thru2023/", pattern = 'userDist',
                                          full.names = TRUE), read_csv)) %>% 
  filter(bias_index == 0)

user_values$inferred_behavior <- "Null"
for (i in 1:nrow(user_values)) {
  this_dist <- user_dists %>% 
    filter(user == user_values$user_login[[i]], taxon == user_values$taxon[[i]])
  
  gt_rate <- mean(user_values$observed_propUnique[i] > this_dist$propUnique)
  lt_rate <- mean(user_values$observed_propUnique[i] < this_dist$propUnique)
  
  # Handle these separately to make sure obs == alt goes to null
  if (lt_rate > 0.975) {
    user_values$inferred_behavior[i] <- "Specialization"
  }
  if (gt_rate > 0.975) {
    user_values$inferred_behavior[i] <- "Novelty"
  }
}

write_csv(user_values, "intermediate/inference_decisions.csv")
