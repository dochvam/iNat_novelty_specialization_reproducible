#############################################################
# functions.R
#
# This file contains support functions for executing the 
# validation simulations
#############################################################

#### function to run the simulation ####

run_validation_simulation <- function(nusers, nspecies, nhex, nrep_per,
                                      nobs_per_user, bias_indexes, ncores,
                                      sim_name, real_det_rates = NULL) {
  
  if (length(nobs_per_user) == 1) {
    nobs_per_user <- rep(nobs_per_user, nusers)
  } else if (length(nobs_per_user) != nusers) {
    stop("Parameter nobs_per_user must be length 1 or length nusers")
  }
  
  if (!is.null(real_det_rates)) {
    abund_vec <- real_det_rates
  } else {
    abund_vec <- 0.5 + 6 * 1:nspecies / nspecies
  }
  
  # Simulate prevalence of the species
  grid_summary_true <- simulate_prevalence(nspecies,
                                           abundance = abund_vec,
                                           nhex,
                                           zif = 0.5)
  
  # Simulate user grid cells
  user_gridcell_hists <- simulate_gridcell_hists(nusers, nhex, nobs_per_user,
                                                 grid_summary_true)
  
  
  # Generate user data from prevalence and grid cell history, given known bias
  # indexes
  user_data <- generate_user_data(grid_summary_true, 
                                  user_gridcell_hists,
                                  bias_indexes)
  
  user_data$gridcell <- user_data$hex
  user_data$iconic_taxon_name <- "Sim" # Placeholder
  
  user_values <<- user_data %>% 
    group_by(user_login) %>% 
    summarize(observed_propUnique = length(unique(scientific_name)) / n()) %>% 
    mutate(taxon = "Sim", 
           sim_name = sim_name,
           true_NS = bias_indexes,
           sim_name = sim_name,
           true_nobs = nobs_per_user)
  
  
  grid_summary_observed <- make_grid_summary(user_data)
  
  # Now, forget the true prevalence, and use the simulated user histories
  # to do the entire workflow.
  byUserResultsAlt <- vector("list", nusers)
  
  bias_index_grid <- c(0)
  
  npoints_fine <- 0
  
  cl <- makeCluster(ncores)
  
  prep <- clusterEvalQ(cl, {
    library(tidyverse)
    source("code_helper/functions.R")
    NULL
  })
  
  byUserResultsAlt <- parLapply(cl, 
                                1:nusers,
                                runOneUser, 
                                target_users = 1:nusers, 
                                bias_index_grid = bias_index_grid, 
                                grid_summary = grid_summary_observed, 
                                dat = user_data,
                                npoints_fine = npoints_fine,
                                n = nrep_per,
                                taxon = "Sim",
                                write.out = FALSE)
  stopCluster(cl)
  
  user_dists <<- bind_rows(byUserResultsAlt) %>% 
    mutate(taxon = "Sim")
  
  result_df <- user_dists %>%
    left_join(user_values, by = c("user" = "user_login", "taxon")) %>% 
    group_by(user, taxon, true_NS, sim_name, true_nobs, observed_propUnique) %>% 
    summarize(gt_rate = mean(observed_propUnique > propUnique),
              lt_rate = mean(observed_propUnique < propUnique)) %>% 
    mutate(type = ifelse(
      gt_rate > 0.975, 
      "Novelty",
      ifelse(lt_rate > 0.975, "Specialization", "Null")
    ))
  
  return(result_df)
}


simulate_prevalence <- function(nspecies, abundance, nhex, zif) {
  
  if (length(abundance) == 1) {
    abundance <- rep(abundance, nspecies)
  } else if (length(abundance) != nspecies) {
    stop("Abundance must be vector of length nspecies, or scalar")
  }
  
  prevalence_mtx <- matrix(integer(), nrow = nspecies, ncol = nhex)
  
  for (j in 1:nhex) {
    nonzero <- rbinom(nspecies,1,1-zif)
    prevalence_mtx[, j] <- rpois(nspecies, nonzero * abundance)
  }
  
  if (any(colSums(prevalence_mtx) == 0)) {
    warning("Generated cells with no observations")
  }
  
  grid_summary <- list()
  for (i in 1:nhex) {
    grid_summary[[i]] <- data.frame(
      scientific_name = 1:nspecies,
      weight = prevalence_mtx[, i] / sum(prevalence_mtx[, i])
    ) %>% 
      filter(weight > 0)
  }
  
  return(grid_summary)
}


simulate_gridcell_hists <- function(nusers, nhex, nobs_per_user, grid_summary_true) {
  
  user_df_list <- list()
  
  obs_per_hex <- lapply(grid_summary_true, nrow) %>% as.numeric()
  
  for (i in 1:nusers) {
    user_df_list[[i]] <- data.frame(
      user_login = i,
      obs_num = 1:nobs_per_user[i], 
      eventDate = as.Date("2000-01-01") + 1:nobs_per_user[i],
      hex = sample(which(obs_per_hex > 0), size = nobs_per_user[i], replace = T)
    )
  }  
  
  return(bind_rows(user_df_list))
}

generate_user_data <- function(grid_summary_true, user_gridcell_hists,
                               bias_indexes) {
  
  user_data <- user_gridcell_hists %>% 
    mutate(scientific_name = "")
  
  for (i in 1:length(bias_indexes)) {
    this_inds <- which(user_data$user_login == i)
    user_data$scientific_name[this_inds] <- draw_new_samples_spatial(
      grid_summary_true, 
      grid_sampling_history = user_data$hex[this_inds], 
      bias_index = bias_indexes[i]
    )
  }
  
  return(user_data)
}








