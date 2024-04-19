#############################################################
# functions.R
#
# This file contains support functions for executing the main
# data analyses
#############################################################

library(sf)
library(mixtools)

# This function puts the input iNaturalist data on a hexagonal grid
associate_wgrid <- function(input_dat, hex_short_diameter) {
  obs_points <- st_as_sf(x = input_dat, coords = c("longitude", "latitude"))
  st_crs(obs_points) <- st_crs("+proj=longlat +datum=WGS84")
  
  obs_points <- st_transform(obs_points, st_crs("+init=epsg:2163 +datum=WGS84"))
  
  bbox <- st_bbox(obs_points)
  bbox[1:2] <- bbox[1:2] - 1.5 * (hex_short_diameter)
  bbox[3:4] <- bbox[3:4] + 1.5 * (hex_short_diameter)
  
  hex_grid <<- st_make_grid(bbox, 
                            cellsize = hex_short_diameter, square = FALSE)
  
  input_dat$gridcell <- as.numeric(unlist(lapply(st_intersects(obs_points, hex_grid),
                                                 function(x) x[[1]])))
  
  return(input_dat)
}


# This function creates a summary list giving the weighted species pool in
# each hex as realized in input_dat
make_grid_summary <- function(input_dat) {
  
  grid_summary_list <- list()
  for (i in 1:max(input_dat$gridcell)) {
    grid_summary_list[[i]] <- input_dat %>% 
      filter(gridcell == i) %>% 
      count(scientific_name, iconic_taxon_name) %>% 
      ungroup() %>% 
      mutate(weight = n / sum(n))
  }
  
  return(grid_summary_list)
}


# This function returns one user's spatial sampling history
# (which grids they sampled in, in order)
get_user_gridhist <- function(input_dat, username) {
  toReturn <- input_dat %>% 
    filter(user_login == username) %>% 
    arrange(eventDate) %>% 
    .$gridcell
  
  return(toReturn)
}


# This function produces the permuted distributions according to each user's
# spatially explicit sampling history, for a given "true" bias index.
perm_dist_spatial <- function(grid_summary, grid_sampling_history, bias_index = 0) {
  #browser()
  
  # Get the overall species pool
  taxon_key <- distinct(do.call(rbind, grid_summary), scientific_name, iconic_taxon_name)
  
  # Make a data frame to store the new sampling history
  new_samples <- data.frame(
    observed_on = 1:length(grid_sampling_history),
    scientific_name = draw_new_samples_spatial(
      grid_summary, grid_sampling_history, bias_index = bias_index
    )
  ) %>% 
    left_join(taxon_key, by = "scientific_name") %>% 
    rename(taxa = iconic_taxon_name)

  # Summarize the new samples
  test <- new_samples %>%
    group_by(taxa, scientific_name) %>%
    arrange(observed_on) %>%
    mutate(species_id_ct = row_number())
  
  # Is each species new?
  test <- test %>% mutate(isNew = (species_id_ct == 1))
  
  ## number unique up to this point, number obs up to this point
  test <- test %>%
    arrange(observed_on) %>%
    mutate(dummy = 1)
  
  test$numUnique <- cumsum(test$isNew)
  
  ## test statistic --> number of species that are unique per number of observations
  
  ## time to first dup
  timeToFirstDup <- test$observed_on[(which(diff(test$numUnique) == 0)[1] + 1)]
  
  byT <- test %>%
    group_by(taxa) %>%
    summarise(count = n()) %>%
    mutate(prop = count / sum(count))
  
  byS <- test %>%
    group_by(scientific_name) %>%
    summarise(count = n()) %>%
    mutate(prop = count / sum(count))
  

  
  
  
  if (nrow(test) >=50){
    numU50 = max(test$numUnique[1:50])
  } else {
    numU50 = NA
  }
  
  # Return the test statistic -- proportion of unique specs.
  max(test$numUnique) / nrow(test)
}

# This function draws new samples from a user's sampling history
draw_new_samples_spatial <- function(grid_summary, 
                                     grid_sampling_history, 
                                     bias_index) {
  
  if (bias_index < -1 | bias_index > 1) {
    stop("bias_index must be between -1 and 1")
  }
  
  new_samples_spec <- character(length(grid_sampling_history))

    
  if (bias_index == 0) {
    # Sample species neutrally, according to availability
    for (i in 1:length(grid_sampling_history)) {
      sampled_index <- sample(size = 1,
                              1:nrow(grid_summary[[grid_sampling_history[i]]]), 
                              prob = grid_summary[[grid_sampling_history[i]]]$weight)
      
      new_samples_spec[i] <- grid_summary[[grid_sampling_history[i]]]$scientific_name[sampled_index]
    }
  } else if (bias_index > 0) {
    # Sample species preferentially for new species
    for (i in 1:length(grid_sampling_history)) {
      
      report_newspec_only <- rbinom(1, 1, bias_index)
      
      if (report_newspec_only == 1) {
        
        this_weights <- grid_summary[[grid_sampling_history[i]]]$weight *
          (1 - (grid_summary[[grid_sampling_history[i]]]$scientific_name %in% new_samples_spec))
        
        if (sum(this_weights) <= 0) {
          this_weights[] <- grid_summary[[grid_sampling_history[i]]]$weight
        }
        
        sampled_index <- sample(size = 1,
                                1:nrow(grid_summary[[grid_sampling_history[i]]]), 
                                prob = this_weights
        )
        
        new_samples_spec[i] <- grid_summary[[grid_sampling_history[i]]]$scientific_name[sampled_index]
        
      } else {
        sampled_index <- sample(size = 1,
                                1:nrow(grid_summary[[grid_sampling_history[i]]]), 
                                prob = grid_summary[[grid_sampling_history[i]]]$weight)
        
        new_samples_spec[i] <- grid_summary[[grid_sampling_history[i]]]$scientific_name[sampled_index]
      }
    }
  } else {
    # Sample species preferentially for those previously observed
    for (i in 1:length(grid_sampling_history)) {
      
      report_oldspec_only <- rbinom(1, 1, -bias_index)
      
      if (report_oldspec_only == 1) {
        
        this_weights <- grid_summary[[grid_sampling_history[i]]]$weight *
          (grid_summary[[grid_sampling_history[i]]]$scientific_name %in% new_samples_spec)
        
        if (sum(this_weights) <= 0) {
          this_weights[] <- grid_summary[[grid_sampling_history[i]]]$weight
        }
        
        sampled_index <- sample(size = 1,
                                1:nrow(grid_summary[[grid_sampling_history[i]]]), 
                                prob = this_weights
        )
        
        new_samples_spec[i] <- grid_summary[[grid_sampling_history[i]]]$scientific_name[sampled_index]
        
      } else {
        sampled_index <- sample(size = 1,
                                1:nrow(grid_summary[[grid_sampling_history[i]]]), 
                                prob = grid_summary[[grid_sampling_history[i]]]$weight)
        
        new_samples_spec[i] <- grid_summary[[grid_sampling_history[i]]]$scientific_name[sampled_index]
      }
    }
  }
  
  new_samples_spec
}


# This function runs the whole exercise to generate the permutations
# of species histories for a given user
runOneUser <- function(i, target_users, bias_index_grid, grid_summary, dat,
                       npoints_fine, n, taxon = NULL, write.out = FALSE) {
  
  this_user <- list()
  user_grid_history <- get_user_gridhist(input_dat = dat, username = target_users[i])
  
  for (j in 1:length(bias_index_grid)) {
    
    this_user[[j]] <- data.frame(
      propUnique = unlist(
        replicate(n, perm_dist_spatial(grid_summary = grid_summary,
                                       grid_sampling_history = user_grid_history,
                                       bias_index = bias_index_grid[j]), simplify = F)
      ),
      bias_index = bias_index_grid[j],
      user = target_users[i]
    )
  }
  
  species <- dat %>% 
    filter(user_login == target_users[i]) %>% 
    .$scientific_name
  
  observed_propUnique <- length(unique(species)) / length(species)
  
  
  byUserResultsTemp <- bind_rows(this_user) %>% 
    group_by(bias_index) %>% 
    summarize(min = min(propUnique), max = max(propUnique))
  
  BI_lower_bound <- byUserResultsTemp$bias_index[max(which(byUserResultsTemp$max < observed_propUnique))]
  BI_upper_bound <- byUserResultsTemp$bias_index[min(which(byUserResultsTemp$min > observed_propUnique))]
  if (is.na(BI_lower_bound) | is.infinite(BI_lower_bound)) BI_lower_bound <- -1
  if (is.na(BI_upper_bound) | is.infinite(BI_upper_bound)) BI_upper_bound <- 1
  
  finer_bias_index_grid <- seq(BI_lower_bound, BI_upper_bound, length.out = npoints_fine)
  finer_bias_index_grid <- finer_bias_index_grid[!finer_bias_index_grid %in% bias_index_grid]
  
  this_user_finer <- list()
  for (j in 1:length(finer_bias_index_grid)) {
    this_user_finer[[j]] <- data.frame(
      propUnique = unlist(
        replicate(n, perm_dist_spatial(grid_summary = grid_summary,
                                       grid_sampling_history = user_grid_history,
                                       bias_index = finer_bias_index_grid[j]), simplify = F)
      ),
      bias_index = finer_bias_index_grid[j],
      user = target_users[i]
    )
  }
  
  output <- bind_rows(this_user, this_user_finer)
  
  if (write.out) {
    write_csv(output, paste0("regroup/intermediate/", taxon, "_", i, ".csv"))
  }
  
  return(bind_rows(this_user, this_user_finer))
}

# Create a plot of one user's observed uniqueness value and the alternative
# distributions against which it is compared
plot_user <- function(user_dists, user_values,
                      this_user_login = NULL, this_taxon = NULL) {
  
  if (is.null(this_user_login) && is.null(this_taxon)) {
    this_comparison <- user_values %>% 
      distinct(user_login, taxon) %>% 
      sample_n(size = 1)
    this_user_login <- this_comparison$user_login
    this_taxon <- this_comparison$taxon
  } else if (is.null(this_user_login)) {
    this_comparison <- user_values %>% 
      distinct(user_login, taxon) %>% 
      filter(taxon == this_taxon) %>% 
      sample_n(size = 1)
    this_user_login <- this_comparison$user_login
  }
  
  this_user_dists <- user_dists %>% filter(user == this_user_login, 
                                           taxon == this_taxon)
  this_user_values <- user_values %>% filter(user_login == this_user_login,
                                             taxon == this_taxon)
  ggplot() +
    geom_density(data = this_user_dists, 
                 aes(propUnique, group = bias_index, fill = bias_index), 
                 bins = 100,
                 alpha = 0.5) +
    scale_fill_gradient2(mid = "gray") +
    geom_vline(data = this_user_values, aes(xintercept=observed_propUnique)) +
    theme_minimal() +
    ylab("") +
    ggtitle(paste0("User ", this_user_login, ", taxon: ", this_taxon,
                   ". nobs = ", this_user_values$nobsT))
}

# Helper function for mixture of normals
get_density <- function(val, distr) {
  test <- density(distr, cut = 0)
  grid <- cbind.data.frame(x = test$x, f.x = test$y)
  
  if (val < min(grid$x) || val > max(grid$x)) {
    return(0)
  }
  
  grid$f.x[which.min(abs(grid$x - val))]
} 

# Get mixture weights from alternative distributions given a true observed value
get_index_weights <- function(val, distr_df) {
  
  weights <- distr_df %>% 
    group_by(bias_index) %>% 
    summarize(density = get_density(val, propUnique)) %>%
    ungroup() %>% 
    mutate(weight = density / sum(density, na.rm = T))
  
  return(weights)
}


# Compute the bias index for a target user, given permutation output
get_bias_index <- function(user_dists, user_values,
                           this_user_login = NULL, this_taxon = NULL) {
  cat(this_user_login, "-", this_taxon, "\n")
  #browser()
  this_user_dists <- user_dists %>% filter(user == this_user_login, 
                                           taxon == this_taxon)
  this_user_values <- user_values %>% filter(user_login == this_user_login,
                                             taxon == this_taxon)
  
  weights <- get_index_weights(val = this_user_values$observed_propUnique, 
                               distr_df = this_user_dists)
  
  if(sum(weights$density) == 0){
    
    ## figure out sign if way out of bounds
    
    assigned_val = ifelse(this_user_values$observed_propUnique> max(this_user_dists$propUnique), 1, -1)
    
    return(data.frame(
      user_login = this_user_login,
      taxon = this_taxon,
      bias_index = assigned_val,
      bs_median = NA,
      bs_mean = NA,
      bs_se = NA,
      Q025 = NA,
      Q975 = NA,
      numNA = NA
    ))
    
  }else{
    
    normal_approx <- this_user_dists %>% 
      group_by(bias_index) %>% 
      summarize(mean = mean(propUnique), sd = sd(propUnique))
    
    new_dist <- rnormmix(500, weights$weight, mu = normal_approx$mean, 
                         sigma = normal_approx$sd) ## will do the normalize for us in future
    
    dist_of_BI <- unlist(lapply(new_dist, function(x) {
      temp <- get_index_weights(x, distr_df = this_user_dists)
      sum(temp$weight * temp$bias_index)
    }))
    
    
    # print(sum(is.na(dist_of_BI)))
    
    return(data.frame(
      user_login = this_user_login,
      taxon = this_taxon,
      bias_index = sum(weights$weight * weights$bias_index),
      bs_median = median(dist_of_BI),
      bs_mean = mean(dist_of_BI, na.rm = TRUE),
      bs_se = sd(dist_of_BI, na.rm = TRUE),
      Q025 = quantile(dist_of_BI, probs = 0.025, na.rm = TRUE),
      Q975 = quantile(dist_of_BI, probs = 0.975, na.rm = TRUE),
      numNA = sum(is.na(dist_of_BI))
    ))
  }
}

