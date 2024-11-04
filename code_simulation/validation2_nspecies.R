library(tidyverse)
library(parallel)

library(mixtools)
if (!dir.exists("validation_sim/output")) { 
  dir.create("validation_sim/output")
}

source("validation_sim/sim_fns.R")
source("revision_thru2023/functions.R")

ncores <- 16
nrep_per <- 250

#### Test: how does changing the size of the species pool impact power to detect bias? ####

# For this simulation, we set a pool of users with bias indexes from -0.5 to 0.5,
# in increments of 0.1. At each level, there are 5 users with a detection history
# of 40 obs, 60 obs, etc. up to 100. We then run simulations with species
# pools of different sizes.

user_profiles <- expand.grid(nobs = c(40, 80, 120, 160, 200),
                             bias_index = seq(from = -0.5, to = 0.5, by = 0.1))
nrep_users <- 20



set.seed(36137)
nspec_vec <- c(5, 100, 200, 500, 2000, 25000)
real_dat <- read_csv("revision_thru2023/data_user_species.csv") %>% 
  count(scientific_name, iconic_taxon_name) %>% 
  sample_frac()

for (i in 2:length(nspec_vec)) {
  writeLines(paste0("Running sim for ", nspec_vec[i], " species"))
  
  this_result <- run_validation_simulation(
    nusers = nrow(user_profiles) * nrep_users,
    nspecies = nspec_vec[i],
    nrep_per = nrep_per,
    nhex = 500,
    real_det_rates = sample(real_dat$n, size = nspec_vec[i], replace = TRUE),
    nobs_per_user = rep(user_profiles$nobs, each = nrep_users),
    bias_indexes = rep(user_profiles$bias_index, each = nrep_users),
    ncores = ncores,
    sim_name = paste0(nspec_vec[i], " species")
  )
  
  write_csv(ungroup(this_result), 
            paste0("validation_sim/output/result_list_",
                   nspec_vec[i], "_realistic_densities_spec.csv"))
}



#### Look at the results ####
colors <- c("Specialization" = "#f0027f",
            "Null" = "#bbbbbb",
            "Novelty" = "#386cb0")

results_files <- list.files("validation_sim/output", pattern = "_spec.csv",
                            full.names = TRUE)
# results_files <- list.files("validation_sim/output", pattern = "0spec.RDS",
#                             full.names = TRUE)

results_df <- lapply(results_files, read_csv) %>% 
  bind_rows() %>% 
  mutate(type = recode(type, "Favoritism" = "Specialization"))


power_levels <- results_df %>% 
  group_by(sim_name, true_nobs, true_NS) %>% 
  summarize(power = mean(type != "Null"))

min_breaks <- power_levels %>% 
  filter(power >= 0.75, true_NS < 0) %>% 
  filter(true_NS == max(true_NS))
max_breaks <- power_levels %>% 
  filter(power >= 0.75, true_NS > 0) %>% 
  filter(true_NS == min(true_NS))

fac_levels_nobs <- paste0(c(40, 80, 120, 160, 200), " obs per user")

sim_name_levels <- paste0(nspec_vec, " species")

good_results_plot <- results_df %>%
  ggplot() +
  geom_bar(aes(true_NS,
               group = factor(type, levels = names(colors)), fill = type)) +
  geom_bar(data = results_df, 
           aes(true_NS, linewidth = true_NS == 0),
           fill = NA, show.legend = F, 
           color = "black") +
  scale_linewidth_manual("", values = c(NA, 0.5)) +
  scale_fill_manual("Inferred\nbehavior", values = colors) +
  facet_grid(factor(sim_name, levels = sim_name_levels)~factor(paste0(true_nobs, " obs per user"),
                                                               levels = fac_levels_nobs)) +
  theme_bw() +
  xlab("True NS index") + ylab("Num. simulated users") +
  ggtitle("Simulation 2") + theme(plot.title = element_text(size=18))

ggsave(filename = "revision_thru2023/figs/validation_sim_2.jpg", plot = good_results_plot,
       width = 8, height = 7)
