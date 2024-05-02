#############################################################
# 3_calculate_NS_index.R
#
# This file contains the code for calculating the NS index
# based on the results from script 2
#############################################################

library(tidyverse)
library(mixtools)
library(parallel)

# Load in helper functions
source("code_helper/functions.R")

ncores <- 40

### Read in the user-level summary results
summaries <- list.files("intermediate", pattern = "user_summary*.", full.names = TRUE)
# might need to update nchar depending on filepaths:
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

# Read in user results in detail and user histories
user_values <- bind_rows(lapply(list.files("intermediate/", pattern = 'user_summary_',
                                           full.names = TRUE), helper))
user_dists <- bind_rows(lapply(list.files("intermediate/", pattern = 'userDist',
                                          full.names = TRUE), read_csv))

# Loop over all specs, in parallel, and calculate the NS index based on a 
# mixture of normals
if (ncores > 1) {
  
  cl <- makeCluster(ncores)
  
  prep <- clusterEvalQ(cl, {
    library(tidyverse)
    library(mixtools)
  })
  
  prep <- clusterExport(cl, varlist = list("get_bias_index",
                                           "get_index_weights",
                                           "get_density",
                                           "user_dists",
                                           "user_values"))
  
  # apply the function (see functions.R)
  result_list <- parLapply(cl, 1:nrow(user_values), function(i) {
    get_bias_index(user_dists = user_dists,
                   user_values = user_values, 
                   this_user_login = user_values$user_login[i],
                   user_values$taxon[i])
  })
  
  result_df <- bind_rows(result_list) %>% 
    mutate(type = ifelse(
      sign(Q025) == sign(Q975), 
      ifelse(sign(Q025) == -1, "Favoritism", "Novelty"),
      "Null"
    )) 
  
  stopCluster(cl)
  
} else {
  result_df <- map2_df(
    user_values$user_login[user_values$taxon != "All"], 
    user_values$taxon[user_values$taxon != "All"],
    get_bias_index, 
    user_dists = user_dists, 
    user_values = user_values
  ) %>% 
    mutate(type = ifelse(
      sign(Q025) == sign(Q975), 
      ifelse(sign(Q025) == -1, "Favoritism", "Novelty"),
      "Null"
    )) 
}


write_csv(result_df, "intermediate/estimated_user_BIs_allPA.csv")
