library(tidyverse)
library(mixtools)
library(parallel)

source("helper_code/functions.R")

ncores <- 40

summaries <- list.files("intermediate", pattern = "user_summary*.", full.names = TRUE)
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

user_values <- bind_rows(lapply(list.files("intermediate/", pattern = 'user_summary_',
                                           full.names = TRUE), helper))
user_dists <- bind_rows(lapply(list.files("intermediate/", pattern = 'userDist',
                                          full.names = TRUE), read_csv))



## these are way out of the distributions
# plot_user(this_taxon = "Herptiles", this_user = "kyleshikes", user_dists, user_values)
# plot_user(this_taxon = "Herptiles", this_user = "wabbytwax", user_dists, user_values)


# plot_user(this_taxon = "Herptiles", this_user = "dnydick", user_dists, user_values)






## when out of bounds need something
#tryThis = get_bias_index(user_dists, user_values,this_user_login = "kyleshikes", this_taxon = "Herptiles")

#tryThis = get_bias_index(user_dists, user_values,this_user_login = "dnydick", this_taxon = "Herptiles")


## this part is pretty slow
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

result_df <- read.csv("intermediate/estimated_user_BIs_allPA.csv")

result_df %>% 
  count(taxon, type) %>% 
  group_by(taxon) %>% 
  mutate(pct = n / sum(n)) %>% 
  ggplot(aes(type, n, fill = type)) +
  geom_col() + 
  facet_wrap(~taxon, scales = "free_y")

#

ggplot(result_df) +
  geom_histogram(aes(bias_index, fill = type))


result_df %>% 
  left_join(user_info) %>% 
  filter(taxon == "All") %>% 
  ggplot(aes(log(nobsT), bias_index, ymin = Q025, ymax = Q975,
             col = type)) +
  geom_pointrange()

result_df %>% 
  left_join(user_info) %>% 
  # filter(taxon == "All") %>% 
  ggplot(aes(log(nobsT), bias_index, ymin = Q025, ymax = Q975,
             col = type)) +
  geom_pointrange(alpha = 0.8) +
  theme_minimal() +
  facet_wrap(~taxon)

