#############################################################
# validation1_usertypes.R
#
# This file conducts a series of simulations testing 
# our ability to recover true biases depending on the overall
# makeup of the user pool. We do not present these results in
# the main manuscript
#############################################################

library(tidyverse)
library(parallel)

if (!dir.exists("intermediate/validation_sim")) {
  dir.create("intermediate/validation_sim")
}

ncores <- 25
n <- 250 # Number of draws per distribution

source("code_helper/functions.R")
source("code_helper/sim_fns.R")

#### Run some simulations ####

# Simulation with no bias
result_list_no_bias <- run_validation_simulation(
  nusers = 100, 
  nspecies = 500, 
  nhex = 500, 
  nobs_per_user = 200, 
  bias_indexes = rep(0, 100), 
  ncores = ncores,
  sim_name = "no_bias"
)
saveRDS(result_list_no_bias, "intermediate/validation_sim/result_no_bias.RDS")

# Simulation w/ half neutral, half no bias
result_list_specialization_bias <- run_validation_simulation(
  nusers = 100, 
  nspecies = 500, 
  nhex = 500, 
  nobs_per_user = 200, 
  bias_indexes = c(rep(0, 50), rep(-0.5, 50)), 
  ncores = ncores,
  sim_name = "half_specialization"
)
saveRDS(result_list_specialization_bias, "intermediate/validation_sim/result_halfspec.RDS")


result_list_novelty_bias <- run_validation_simulation(
  nusers = 100, 
  nspecies = 500, 
  nhex = 500, 
  nobs_per_user = 200, 
  bias_indexes = c(rep(0, 50), rep(0.5, 50)), 
  ncores = ncores,
  sim_name = "half_novelty"
)
saveRDS(result_list_novelty_bias, "intermediate/validation_sim/result_halfnovel.RDS")



#### Look at the results ####
colors <- c("Favoritism" = "#f0027f",
            "Null" = "#bbbbbb",
            "Novelty" = "#386cb0")

result_list_no_bias <- readRDS("intermediate/validation_sim/result_no_bias.RDS")
result_list_specialization_bias <- readRDS("intermediate/validation_sim/result_halfspec.RDS")
result_list_novelty_bias <- readRDS("intermediate/validation_sim/result_halfnovel.RDS")


result_list_no_bias$result_df %>% 
  ggplot() + 
  geom_histogram(aes(bias_index, fill = type)) +
  scale_fill_manual(values = colors) +
  facet_wrap(~paste0("True bias = ", true_bias_index)) +
  ggtitle("Validation simulation 1")


result_list_specialization_bias$result_df %>% 
  ggplot() + 
  geom_histogram(aes(bias_index, fill = type)) +
  scale_fill_manual(values = colors) +
  facet_wrap(~paste0("True bias = ", true_bias_index)) +
  ggtitle("Validation simulation 2")


result_list_novelty_bias$result_df %>% 
  ggplot() + 
  geom_histogram(aes(bias_index, fill = type)) +
  scale_fill_manual(values = colors) +
  facet_wrap(~paste0("True bias = ", true_bias_index)) +
  ggtitle("Validation simulation 3")
