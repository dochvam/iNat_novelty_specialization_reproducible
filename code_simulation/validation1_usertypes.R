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
  nusers = 200, 
  nspecies = 1000, 
  nhex = 500, 
  nrep_per = nrep_per,
  nobs_per_user = 200, 
  bias_indexes = rep(0, 200), 
  ncores = ncores,
  sim_name = "no_bias"
)
write_csv(result_list_no_bias, "validation_sim/output/result_no_bias.csv")

# Simulation w/ half neutral, half no bias
result_list_specialization_bias <- run_validation_simulation(
  nusers = 200, 
  nspecies = 1000, 
  nhex = 500, 
  nrep_per = nrep_per,
  nobs_per_user = 200, 
  bias_indexes = c(rep(0, 50), rep(-0.5, 150)), 
  ncores = ncores,
  sim_name = "half_specialization"
)
write_csv(result_list_specialization_bias, "validation_sim/output/result_halfspec.csv")


result_list_novelty_bias <- run_validation_simulation(
  nusers = 200, 
  nspecies = 1000, 
  nhex = 500, 
  nrep_per = nrep_per,
  nobs_per_user = 200, 
  bias_indexes = c(rep(0, 50), rep(0.5, 150)), 
  ncores = ncores,
  sim_name = "half_novelty"
)
write_csv(result_list_novelty_bias, "validation_sim/output/result_halfnovel.csv")


#### Look at the results ####
colors <- c("Specialization" = "#f0027f",
            "Null" = "#bbbbbb",
            "Novelty" = "#386cb0")

result_list_no_bias <- read_csv("validation_sim/output/result_no_bias.csv")
result_list_specialization_bias <- read_csv("validation_sim/output/result_halfspec.csv")
result_list_novelty_bias <- read_csv("validation_sim/output/result_halfnovel.csv")

get_legend <- function(myggplot) {
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

p1 <- bind_rows(
  result_list_no_bias %>% mutate(sim = "Simulation 1.1"),
  result_list_specialization_bias %>% mutate(sim = "Simulation 1.2"),
  result_list_novelty_bias %>% mutate(sim = "Simulation 1.3")
)%>% 
  ggplot() + 
  geom_bar(aes(as.factor(true_NS), fill = type), position = "fill") +
  scale_fill_manual("Inferred\nbehavior", values = colors) +
  facet_wrap(~sim, ncol = 1) +
  theme_minimal() +
  xlab("True NS") + ylab("Pct. of users") +
  coord_flip()


ggsave("revision_thru2023/figs/validation_sim_1.jpg", p1,
       width = 5, height = 3, dpi = "retina")
