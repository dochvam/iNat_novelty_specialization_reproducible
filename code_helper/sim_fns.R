#### function to run the simulation ####

run_validation_simulation <- function(nusers, nspecies, nhex, 
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
  user_gridcell_hists <- simulate_gridcell_hists(nusers, nhex, nobs_per_user)
  
  
  # Generate user data from prevalence and grid cell history, given known bias
  # indexes
  user_data <- generate_user_data(grid_summary_true, 
                                  user_gridcell_hists,
                                  bias_indexes)
  
  user_data$gridcell <- user_data$hex
  user_data$iconic_taxon_name <- "Arachnidae" # Placeholder
  
  user_values <<- user_data %>% 
    group_by(user_login) %>% 
    summarize(observed_propUnique = length(unique(scientific_name)) / n()) %>% 
    mutate(taxon = "Arachnidae",
           sim_name = sim_name)
  
  
  grid_summary_observed <- make_grid_summary(user_data)
  
  # Now, forget the true prevalence, and use the simulated user histories
  # to do the entire workflow.
  byUserResultsAlt <- vector("list", nusers)
  
  bias_index_grid <- c(-0.99, seq(from = -0.8, to = 0.8, by = 0.2), 0.99)
  
  npoints_fine <- 12
  
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
                                n = n,
                                taxon = "Arachnidae",
                                write.out = TRUE)
  
  user_dists <<- bind_rows(byUserResultsAlt) %>% 
    mutate(taxon = "Arachnidae")
  
  prep <- clusterExport(cl, varlist = list("get_bias_index",
                                           "get_index_weights",
                                           "get_density",
                                           "user_dists",
                                           "user_values"))
  
  result_list <- parLapply(cl, 1:nusers, function(i) {
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
    )) %>% 
    mutate(true_bias_index = bias_indexes,
           sim_name = sim_name)
  
  stopCluster(cl)
  
  result_df$true_bias_index <- bias_indexes
  result_df$true_nobs <- nobs_per_user
  result_df$sim_name <- sim_name
  
  return(list(
    result_df = result_df,
    user_values = user_values
  ))
}


simulate_prevalence <- function(nspecies, abundance, nhex, zif) {
  
  if (length(abundance) == 1) {
    abundance <- rep(abundance, nspecies)
  } else if (length(abundance) != nspecies) {
    stop("Abundance must be vector of length nspecies, or scalar")
  }
  
  prevalence_mtx <- matrix(integer(), nrow = nspecies, ncol = nhex)

  
  
  for (i in 1:nspecies) {
    for (j in 1:nhex) {
      nonzero <- rbinom(1,1,1-zif)
      prevalence_mtx[i, j] <- rpois(1, nonzero * abundance[i])
    }
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


simulate_gridcell_hists <- function(nusers, nhex, nobs_per_user) {
  
  user_df_list <- list()

  for (i in 1:nusers) {
    user_df_list[[i]] <- data.frame(
      user_login = i,
      obs_num = 1:nobs_per_user[i], 
      eventDate = as.Date("2000-01-01") + 1:nobs_per_user[i],
      hex = sample(1:nhex, size = nobs_per_user[i], replace = T)
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




