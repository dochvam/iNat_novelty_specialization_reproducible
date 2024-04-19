#############################################################
# validation2_usertypes.R
#
# This file conducts a series of simulations testing 
# our ability to recover true biases as the number of users
# and number of species varies. This simulation is presented
# in the main manuscript.
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

#### Test: how does changing the size of the species pool impact power to detect bias? ####

# For this simulation, we set a pool of users with bias indexes from -0.5 to 0.5,
# in increments of 0.1. At each level, there are 5 users with a detection history
# of 40 obs, 60 obs, etc. up to 100. We then run simulations with species
# pools of different sizes.

user_profiles <- expand.grid(nobs = c(40, 80, 120, 160, 200),
                             bias_index = seq(from = -0.5, to = 0.5, by = 0.1))
nrep_users <- 20


# Simulation with low number of species (Test that 40 is a small enough det hist)

nspec_vec <- c(100, 200, 500, 2000)

for (i in 1:length(nspec_vec)) {
  writeLines(paste0("Running sim for ", nspec_vec[i], " species"))
  
  this_result_list <- run_validation_simulation(
    nusers = nrow(user_profiles) * nrep_users, 
    nspecies = nspec_vec[i], 
    nhex = 500, 
    nobs_per_user = rep(user_profiles$nobs, each = nrep_users), 
    bias_indexes = rep(user_profiles$bias_index, each = nrep_users), 
    ncores = ncores,
    sim_name = paste0(nspec_vec[i], " species")
  )
  
  saveRDS(this_result_list, 
          paste0("intermediate/validation_sim/result_list_",
                 nspec_vec[i], "spec.RDS"))
}

set.seed(36137)
nspec_vec <- c(100, 200, 500, 2000)
real_dat <- read_csv("intermediate/data_user_species.csv") %>% 
  count(scientific_name) %>% 
  sample_frac()

for (i in 1:length(nspec_vec)) {
  writeLines(paste0("Running sim for ", nspec_vec[i], " species"))
  
  this_result_list <- run_validation_simulation(
    nusers = nrow(user_profiles) * nrep_users, 
    nspecies = nspec_vec[i], 
    nhex = 500, 
    real_det_rates = real_dat$n[1:nspec_vec[i]],
    nobs_per_user = rep(user_profiles$nobs, each = nrep_users), 
    bias_indexes = rep(user_profiles$bias_index, each = nrep_users), 
    ncores = ncores,
    sim_name = paste0(nspec_vec[i], " species")
  )
  
  saveRDS(this_result_list, 
          paste0("intermediate/validation_sim/result_list_",
                 nspec_vec[i], "_realistic_densities_spec.RDS"))
}



#### Look at the results ####
colors <- c("Specialization" = "#f0027f",
            "Null" = "#bbbbbb",
            "Novelty" = "#386cb0")

results_files <- list.files("intermediate/validation_sim", pattern = "_spec.RDS",
                            full.names = TRUE)
# results_files <- list.files("intermediate/validation_sim", pattern = "0spec.RDS",
#                             full.names = TRUE)

results_df <- lapply(results_files, function(x) {
  temp <- readRDS(x)
  temp$result_df
}) %>% 
  bind_rows() %>% 
  mutate(type = recode(type, "Favoritism" = "Specialization"))

  
power_levels <- results_df %>% 
  group_by(sim_name, true_nobs, true_bias_index) %>% 
  summarize(power = mean(type != "Null"))

min_breaks <- power_levels %>% 
  filter(power >= 0.75, true_bias_index < 0) %>% 
  filter(true_bias_index == max(true_bias_index))
max_breaks <- power_levels %>% 
  filter(power >= 0.75, true_bias_index > 0) %>% 
  filter(true_bias_index == min(true_bias_index))

fac_levels_nobs <- paste0(c(40, 80, 120, 160, 200), " obs per user")

sim_name_levels <- paste0(nspec_vec, " species")

results_df %>%
  ggplot() +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(aes(true_bias_index, bias_index, col = type)) +
  scale_color_manual(values = colors) +
  facet_grid(factor(sim_name, levels = sim_name_levels)~factor(paste0(true_nobs, " obs per user"),
                             levels = fac_levels_nobs)) +
  theme_bw() +
  xlab("True bias index") + ylab("Point est. of bias index") 




results_df %>%
  ggplot() +
  geom_point(aes(true_bias_index, bs_se, col = type)) +
  scale_color_manual(values = colors) +
  facet_grid(factor(sim_name, levels = sim_name_levels)~factor(paste0(true_nobs, " obs per user"),
                                                               levels = fac_levels_nobs)) +
  theme_bw() +
  xlab("True bias index") + ylab("SE") +
  geom_vline(aes(xintercept = true_bias_index), data = max_breaks,
             color = colors[3]) +
  geom_vline(aes(xintercept = true_bias_index), data = min_breaks,
             color = colors[1])


results_df %>%
  ggplot() +
  geom_bar(aes(factor(type, levels = names(colors)), fill = type),
           show.legend = F) +
  scale_fill_manual("", values = colors) +
  facet_grid(factor(sim_name, levels = sim_name_levels)~factor(paste0(true_nobs, " obs per user"),
                                                               levels = fac_levels_nobs)) +
  theme_bw() +
  xlab("Estimated bias type") + ylab("Num. individuals") 



good_results_plot <- results_df %>%
  ggplot() +
  geom_bar(aes(true_bias_index,
               group = factor(type, levels = names(colors)), fill = type)) +
  geom_bar(data = results_df, 
           aes(true_bias_index, linewidth = true_bias_index == 0),
           fill = NA, show.legend = F, 
           color = "black") +
  scale_linewidth_manual("", values = c(NA, 0.5)) +
  scale_fill_manual("Inferred bias type", values = colors) +
  facet_grid(factor(sim_name, levels = sim_name_levels)~factor(paste0(true_nobs, " obs per user"),
                                                               levels = fac_levels_nobs)) +
  theme_bw() +
  xlab("True NS index") + ylab("Num. simulated users") 


ggsave(filename = "regroup/figs/sim_results.jpg", plot = good_results_plot,
       width = 8, height = 5)


# Calculate the overall false positive rate and false negative rate per sim, effort
results_df %>% 
  mutate(correct = (true_bias_index == 0 & type == "Null") |
                   (true_bias_index != 0 & type != "Null")) %>% 
  group_by(true_nobs, sim_name, truly_null = true_bias_index == 0) %>% 
  summarize(false_rate = 1 - mean(correct)) %>% 
  mutate(false_type = ifelse(truly_null, "False positive", "False negative")) %>% 
  ggplot() +
  geom_line(aes(true_nobs, false_rate, group = sim_name)) +
  geom_point(aes(true_nobs, false_rate, group = sim_name)) +
  facet_wrap(~false_type)



type_colors <- c(
  "True positive" = "darkgray",
  "True negative" = "lightgray",
  "False negative" = "#2277aa",
  "False positive" = "#991599",
  "Erroneous" = "darkred"
)

results_df %>%
  mutate(result_type = ifelse(
    true_bias_index < 0 & type == "Favoritism" | 
      true_bias_index > 0 & type == "Novelty", 
    "True positive", ifelse(
      true_bias_index == 0 & type == "Null",
      "True negative", ifelse(
        type == "Null", "False negative", ifelse(
          true_bias_index != 0, "Erroneous", "False positive"
        )
    )
  ))) %>%
  ggplot() +
  geom_bar(aes(true_bias_index,
               group = factor(result_type), fill = result_type)) +
  scale_fill_manual("Inferred bias type", values = type_colors) +
  facet_grid(factor(sim_name, levels = sim_name_levels)~factor(paste0(true_nobs, " obs per user"),
                                                               levels = fac_levels_nobs)) +
  theme_bw() +
  xlab("True bias index") + ylab("Num. individuals") 


# 
