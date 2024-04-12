library(tidyverse)
library(lubridate)
library(pbapply)


ncores <- 40

# # from cleanPA1.R

data_user_species <- read_csv("intermediate/data_user_species.csv")
user_info <- read_csv("intermediate/pa_user_info.csv")

source("helper_code/functions.R")


target_taxa <- c("Insects", "Plants", "Herptiles", "Arachnids", "Birds", "Mammals", "All")
n <- 250 # Number of draws per distribution

# Number of obs. a user needs to be included - all
userCutoff_All <- 100

# Number of obs. a user needs to be included - taxon-specific
userCutoff_Taxon <- 25

## numbers to get us to 25 obs

for (t_ind in 1:length(target_taxa)) {
  thistaxon <- target_taxa[[t_ind]]
  cat("Running model for taxon:", thistaxon)
  
  if (thistaxon == "All") {
    dat_thistaxon <- data_user_species
    user_info <- dat_thistaxon %>% group_by(user_login) %>% 
      summarize(nobsT = n(), first_obs = min(eventDate)) %>% 
      arrange(-nobsT)
    user_sub <- user_info %>% filter(nobsT > userCutoff_All)
    
  } else if (thistaxon == "Birds") {
    dat_thistaxon <- data_user_species %>% 
      filter(iconic_taxon_name == "Aves")
    user_info <- dat_thistaxon %>% group_by(user_login) %>% 
      summarize(nobsT = n(), first_obs = min(eventDate)) %>% 
      arrange(-nobsT)
    user_sub <- user_info %>% filter(nobsT > userCutoff_Taxon)
    
  } else if (thistaxon == "Mammals") {
    dat_thistaxon <- data_user_species %>% 
      filter(iconic_taxon_name == "Mammalia")
    user_info <- dat_thistaxon %>% group_by(user_login) %>% 
      summarize(nobsT = n(), first_obs = min(eventDate)) %>% 
      arrange(-nobsT)
    user_sub <- user_info %>% filter(nobsT > userCutoff_Taxon)
    
  } else if (thistaxon == "Insects") {
    dat_thistaxon <- data_user_species %>% 
      filter(iconic_taxon_name == "Insecta")
    user_info <- dat_thistaxon %>% group_by(user_login) %>% 
      summarize(nobsT = n(), first_obs = min(eventDate)) %>% 
      arrange(-nobsT)
    user_sub <- user_info %>% filter(nobsT > userCutoff_Taxon)
  } else if (thistaxon == "Plants") {
    dat_thistaxon <- data_user_species %>% 
      filter(iconic_taxon_name == "Plantae")
    user_info <- dat_thistaxon %>% group_by(user_login) %>% 
      summarize(nobsT = n(), first_obs = min(eventDate)) %>% 
      arrange(-nobsT)
    user_sub <- user_info %>% filter(nobsT > userCutoff_Taxon)
  } else if (thistaxon == "Herptiles") {
    dat_thistaxon <- data_user_species %>% 
      filter(iconic_taxon_name %in% c("Reptilia", "Amphibia"))
    user_info <- dat_thistaxon %>% group_by(user_login) %>% 
      summarize(nobsT = n(), first_obs = min(eventDate)) %>% 
      arrange(-nobsT)
    user_sub <- user_info %>% filter(nobsT > userCutoff_Taxon)
  } else if (thistaxon == "Arachnids") {
    dat_thistaxon <- data_user_species %>% 
      filter(iconic_taxon_name %in% "Arachnida")
    user_info <- dat_thistaxon %>% group_by(user_login) %>% 
      summarize(nobsT = n(), first_obs = min(eventDate)) %>% 
      arrange(-nobsT)
    user_sub <- user_info %>% filter(nobsT > userCutoff_Taxon)
  } else {
    stop(print0("Target taxon ", thistaxon, " not recognized."))
  }
  
  # Get a list of target users
  target_users <- user_sub$user_login
  
  user_sub$observed_propUnique <- lapply(target_users, function(x) {
    specs_observed <- dat_thistaxon %>% filter(user_login == x) %>% .$scientific_name
    length(unique(specs_observed)) / length(specs_observed)
  }) %>% unlist()
  
  user_sub$observed_timeToFirstDup <- lapply(target_users, function(x) {
    specs_observed <- dat_thistaxon %>% filter(user_login == x) %>% .$scientific_name
    min(which(duplicated(specs_observed)))
  }) %>% unlist()
  
  user_sub <- user_sub %>% 
    mutate(taxon = thistaxon)
  write_csv(user_sub, paste0("intermediate/user_summary_", thistaxon, ".csv"))
  
  # Associate the data with a spatial grid and summarize species counts by grid
  dat <- associate_wgrid(input_dat = dat_thistaxon,
                           hex_short_diameter = 25000)
  dat$eventDate <- as.Date(dat$eventDate)
  
  grid_summary <- make_grid_summary(dat)
  
  
  #### Null and alternative distributions ####
  
  byUserResultsAlt <- vector("list", length(target_users))
  
  bias_index_grid <- c(-0.99, seq(from = -0.8, to = 0.8, by = 0.2), 0.99)

  npoints_fine <- 12
  
  
  if (ncores > 1) {
    #### Iterate over all users in parallel ####
    
    library(parallel)
    cl <- makeCluster(ncores)
    
    prep <- clusterEvalQ(cl, {
      library(tidyverse)
      source("helper_code/functions.R")
      NULL
    })
    
    byUserResultsAlt <- parLapply(cl, 
                                  sample(1:nrow(user_sub), size = nrow(user_sub)),
                                  runOneUser, 
                                  target_users = target_users, 
                                  bias_index_grid = bias_index_grid, 
                                  grid_summary = grid_summary, 
                                  dat = dat,
                                  npoints_fine = npoints_fine,
                                  n = n,
                                  taxon = thistaxon,
                                  write.out = TRUE)
    
    stopCluster(cl)
    # saveRDS(list(target_users = target_users,
    #              bias_index_grid = bias_index_grid,
    #              grid_summary = grid_summary,
    #              dat = dat), "../test.RDS")
    
  } else {
    #### Sequential version ####
    pb <- progress::progress_bar$new(
      total = length(target_users) * (length(bias_index_grid) + npoints_fine)
    )
    
    for (i in 1:length(target_users)) {
      
      byUserResultsAlt[[i]] <- 
        runOneUser(i, target_users, bias_index_grid, grid_summary, dat,
                   npoints_fine, n)
      
    }
  }
   ## 3.4 hours
  
  
  userDistributions <- bind_rows(byUserResultsAlt) %>% 
    mutate(taxon = thistaxon)
  write_csv(userDistributions, paste0("intermediate/userDistributions_", thistaxon, ".csv"))
}
