library(tidyverse)
source("regroup/functions.R")


colors <- c("Specialization" = "#f0027f",
            "Null" = "#bbbbbb",
            "Novelty" = "#386cb0")


result_df <- read_csv("regroup/estimated_user_BIs_allPA.csv") %>% 
  mutate(type = recode(type, "Favoritism" = "Specialization")) %>% 
  mutate(type = factor(type, levels = names(colors)))
  

user_info <- lapply(list.files("regroup", pattern = "user_summary_",
                               full.names = TRUE), read_csv) %>% 
  bind_rows()

#### Make some example figs ####

user_dists <- bind_rows(lapply(list.files("regroup/", pattern = 'userDist',
                                          full.names = TRUE), read_csv))

plot_user(this_taxon = "All", this_user = "catspit", user_dists, user_info)


user_dists %>% filter(taxon == "All", 
                      user == "catspit", 
                      abs(bias_index - -0.2) <= 1e-3 | 
                      abs(bias_index - 0.2) <= 1e-3) %>% 
  ggplot() +
  geom_density(aes(propUnique, group = bias_index, fill = bias_index),
               alpha = 0.5) +
  scale_fill_gradient2(mid = "gray") +
  theme_minimal() +
  ylab("Density")


#### 

(taxon_counts_all <- result_df %>% 
  count(taxon, type) %>%
  filter(taxon == "All") %>%
  group_by(taxon) %>% 
  mutate(pct = n / sum(n),
         taxon = paste0(taxon, " (n=", sum(n), ")")) %>% 
  ggplot(aes(type, n, fill = type)) +
  geom_col(show.legend = F) + 
  facet_wrap(~taxon, scales = "free_y") +
  scale_fill_manual("Bias type", values = colors) +
  xlab("") + ylab("Number of users") +
  theme_minimal())
(taxon_counts_taxon <- result_df %>% 
  count(taxon, type) %>%
    filter(taxon != "All") %>% 
  group_by(taxon) %>% 
  mutate(pct = n / sum(n),
         taxon = paste0(taxon, " (n=", sum(n), ")")) %>% 
  ggplot(aes(type, n, fill = type)) +
  geom_col(show.legend = F) + 
  facet_wrap(~taxon, scales = "free_y") +
  scale_fill_manual("Bias type", values = colors) +
  xlab("") + ylab("Number of users") +
  theme_minimal())

lmtx <- matrix(c(NA,2,1,2,NA,2), nrow = 2)

taxon_counts_fig <- gridExtra::arrangeGrob(taxon_counts_all, taxon_counts_taxon, layout_matrix = lmtx,
                                 widths = c(0.5, 1, 0.5), heights = c(1, 1.5))

#

ggplot(result_df) +
  geom_histogram(aes(bias_index, fill = type),
                 position = "identity", alpha = 0.5) +
  scale_fill_manual("Bias type", values = colors) +
  theme_minimal()


pointcloud_fig <- result_df %>% 
  left_join(user_info) %>% 
  ggplot(aes(log(nobsT), bias_index, ymin = Q025, ymax = Q975,
             col = type)) +
  geom_pointrange() +
  facet_wrap(~taxon) +
  xlab("Log number of obs.") +
  ylab("Bias index (95%CI)") +
  scale_color_manual("Bias type", values = colors) + 
  theme_minimal()

alphas <- c(1, 0.5, 1)
names(alphas) <- names(colors)
(se_fig_All <- result_df %>% 
    filter(taxon == "All") %>% 
    left_join(user_info) %>% 
    ggplot(aes(log(nobsT), Q975-Q025, group = taxon, color = type, alpha = type)) +
    geom_point(show.legend = FALSE) +
    # geom_smooth(color = "black", alpha = 0.5,
    #             method = lm, formula = y~ poly(x, 2)) +
    xlab("Log number of obs.") +
    ylab("Width of 95%CI on bias index") +
    scale_color_manual("Bias type", values = colors) +
    scale_alpha_manual("Bias type", values = alphas) +
    facet_wrap(~taxon) +
    theme_minimal())
(se_fig_taxa <- result_df %>% 
  filter(taxon != "All") %>% 
  left_join(user_info) %>% 
  ggplot(aes(log(nobsT), Q975-Q025, group = taxon, color = type, alpha = type)) +
  geom_point() +
  # geom_smooth(color = "black", alpha = 0.5,
  #             method = lm, formula = y~ poly(x, 2)) +
  xlab("Log number of obs.") +
  ylab("Width of 95%CI on bias index") +
  scale_color_manual("Bias type", values = colors) +
  scale_alpha_manual("Bias type", values = alphas) +
  facet_wrap(~taxon) +
  theme_minimal())

lmtx <- matrix(c(NA,2,1,2,NA,2), nrow = 2)

se_fig <- gridExtra::arrangeGrob(se_fig_All, se_fig_taxa, layout_matrix = lmtx,
                                 widths = c(0.5, 1, 0.5), heights = c(1, 1.5))

#### How do users compare across taxa? ####
result_all <- result_df %>% 
  select(user_login, taxon, type) %>% 
  filter(taxon == "All") %>% 
  rename(type_all = type) %>% 
  select(-taxon)

# Calculate the percentages
percentage_df <- result_df %>%
  filter(user_login %in% result_all$user_login, taxon != "All") %>%
  select(user_login, taxon, type) %>% 
  left_join(result_all) %>%
  count(type_all, type) %>%
  mutate(percentage = n / sum(n))
percentage_df_B <- result_df %>%
  filter(user_login %in% result_all$user_login, taxon != "All") %>%
  select(user_login, taxon, type) %>% 
  left_join(result_all) %>%
  count(type_all, type, taxon) %>%
  group_by(taxon) %>% 
  mutate(percentage = n / sum(n))

# Create the heatmap plot
cross_taxon_fig_A <- percentage_df %>%
  ggplot(aes(x = type_all, y = type, fill = percentage)) +
  geom_tile() + # using geom_tile for the heatmap
  geom_text(aes(label = paste0(round(percentage * 100, 1), "%"),
                color = ifelse(percentage > 0.5, "black", "white")),
            show.legend = F) +
  scale_fill_viridis_c("% of users", labels = scales::percent, limits = c(0, max(percentage_df_B$percentage))) + # to display the fill as percentages
  xlab("Overall bias type") +
  ylab("Taxon-specific bias type") +
  theme_minimal() + 
  coord_fixed() +
  scale_color_manual(values = c("black" = "black", "white" = "white")) +
  ggtitle("A. All taxa")


cross_taxon_fig_B <- percentage_df_B %>% 
  ggplot(aes(x = type_all, y = type, fill = percentage)) +
  geom_tile() + # using geom_tile for the heatmap
  geom_text(aes(label = paste0(round(percentage * 100, 1), "%"),
                color = ifelse(percentage > 0.5, "black", "white")),
            show.legend = F, size = 2.8) +
  scale_fill_viridis_c("% of users", labels = scales::percent) + # to display the fill as percentages
  xlab("Overall bias type") +
  ylab("Taxon-specific bias type") +
  theme_minimal() + 
  coord_fixed() +
  facet_wrap(~taxon, nrow = 1) +
  scale_color_manual(values = c("black" = "black", "white" = "white")) +
  ggtitle("B. Separate by taxon")


cross_taxon_fig <- gridExtra::arrangeGrob(grobs  = list(
  cross_taxon_fig_A, cross_taxon_fig_B
), nrow = 2)

#### save some plots ####

ggsave("regroup/figs/fig1_counts.jpg", taxon_counts_fig, width = 7, height = 6)
ggsave("regroup/figs/fig2_ptcloud.jpg", pointcloud_fig, width = 9, height = 6)
ggsave("regroup/figs/fig3_SEs.jpg", se_fig, width = 9, height = 7)
ggsave("regroup/figs/fig4_crossplot.jpg", cross_taxon_fig, width = 12, height = 5.5)

#### Representative user plots (supplement?) ####

user_values <- bind_rows(lapply(list.files("regroup/", pattern = 'user_summary_',
                                           full.names = TRUE), read_csv))
user_dists <- bind_rows(lapply(list.files("regroup/", pattern = 'userDist',
                                          full.names = TRUE), read_csv))


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
                 alpha = 0.5) +
    scale_fill_gradient2(mid = "gray") +
    geom_vline(data = this_user_values, aes(xintercept=observed_propUnique)) +
    theme_minimal() +
    ylab("") +
    ggtitle(paste0("User ", this_user_login, ", taxon: ", this_taxon,
                   ". nobs = ", this_user_values$nobsT))
}

plot_user(this_taxon = "All", this_user = "kyleshikes", 
          user_dists, user_values)

# Figures:
# > Counts by taxon
# > Plots of estimates, SEs w/r/t number of obs.
# > Plot of who we'd need to drop based on licensing?

# How does novelism vs. specialism interact for a user across taxa? 

# Main takeaways:
#  > 


#### Plot the hex grid over PA ####
library(terra)
library(tidyterra)

# Make hex counts
source("regroup/functions.R")
data_user_species <- read_csv("regroup/data_user_species.csv")

dat <- associate_wgrid(input_dat = data_user_species,
                       hex_short_diameter = 25000)


usa_map <- geodata::gadm(country = "USA", level = 1, path = "data/geodata")
PA_map <- usa_map[usa_map$NAME_1 == "Pennsylvania", ]

obs_per_hex <- dat %>% count(gridcell)
specs_per_hex <- dat %>% 
  distinct(gridcell, scientific_name) %>% 
  count(gridcell) %>% 
  rename(nspec = n)

hex_grid_vect <- vect(hex_grid)
hex_grid_vect$gridcell <- 1:nrow(hex_grid_vect)
hex_grid_vect <- left_join(hex_grid_vect, obs_per_hex, by = "gridcell") %>% 
  left_join(specs_per_hex, by = "gridcell") %>% 
  filter(n > 0) %>% 
  terra::project(PA_map)

gridmap <- ggplot() +
  geom_spatvector(data = hex_grid_vect, aes(fill = log(n))) +
  geom_spatvector(data = PA_map, fill = NA, col = "black", linewidth = 0.8) +
  theme_minimal() +
  scale_fill_viridis_c("Log number of \niNaturalist obs.")
gridmap2 <- ggplot() +
  geom_spatvector(data = hex_grid_vect, aes(fill = log(nspec))) +
  geom_spatvector(data = PA_map, fill = NA, col = "black", linewidth = 0.8) +
  theme_minimal() +
  scale_fill_viridis_c("Log number of \nspecies observed")

ggsave("regroup/figs/hexmap.jpg", gridmap, width = 6, height = 3.5)


#### Visualize basics ####

# Start date
user_info %>% 
  filter(taxon == "All") %>% 
  count(year = year(first_obs)) %>% 
  ggplot() + 
  geom_col(aes(as.factor(year), n)) +
  theme_minimal() + xlab("Year") + ylab("Num. users") +
  ggtitle("Year of first observation")


# Num. species per taxon
all_dat <- read_csv("regroup/data_user_species.csv") %>% 
  mutate(iconic_taxon_name = recode(iconic_taxon_name, "Reptilia" = "Herptiles",
                                    "Amphibia" = "Herptiles",
                                    "Arachnida" = "Arachnids", "Aves" = "Birds",
                                    "Insecta" = "Insects", 
                                    "Mammalia" = "Mammals", "Plantae" = "Plants"
                                    
                                    ))

spec_counts <- data.frame(
  taxon = c("All", "Arachnids", "Birds", "Herptiles",
            "Insects", "Mammals", "Plants"),
  nspec = NA
)

for (i in 1:nrow(spec_counts)) {
  this_users <- user_info %>% filter(taxon == spec_counts$taxon[i])
  obs <- all_dat %>% 
    filter(user_login %in% this_users$user_login)
  if (spec_counts$taxon[i] != "All") {
    obs <- obs %>% filter(iconic_taxon_name == spec_counts$taxon[i])
  }
  spec_counts$nspec[i] <- length(unique(obs$scientific_name))
}

  
  
  
spec_counts %>% 
  ggplot() +
  geom_col(aes(iconic_taxon_name, n)) +
  theme_minimal() + xlab("Taxonomic group") + ylab("Num. species") +
  ggtitle("Number of species in each taxonomic group")

#### Visualize top orders of favoritists ####
# Make hex counts
source("regroup/functions.R")
data_user_species <- read_csv("regroup/data_user_species.csv")

favoritists <- result_df %>% 
  filter(taxon == "All", type == "Specialization")


find_most_common <- function(df, x1, x2) {
  df %>%
    group_by(!!sym(x1)) %>%
    count(!!sym(x2)) %>%
    group_by(!!sym(x1)) %>%
    mutate(fraction = n / sum(n)) %>%
    arrange(desc(n)) %>%
    slice(1) %>%
    select(!!sym(x1), most_common = !!sym(x2), fraction)
}

top_species_data <- data_user_species %>% 
  filter(user_login %in% result_df$user_login) %>%
  find_most_common("user_login", "scientific_name") %>% 
  left_join(result_df[, c("user_login", "type"), by = "user_login"])

ggplot(top_species_data) +
  geom_density(aes(nimble::logit(fraction), group = type, col = type))


add_taxon_group <- function(df) {
  df %>% 
    mutate(taxon = case_match(
      iconic_taxon_name,
      "Insecta" ~ "Insects", 
      "Plantae" ~ "Plants", 
      c("Reptilia", "Amphibia") ~ "Herptiles", 
      "Arachnida" ~ "Arachnids", 
      "Aves" ~ "Birds", 
      "Mammalia" ~ "Mammals",
      .default = "Other"
    ))
}


top_taxon_data <- data_user_species %>%
  add_taxon_group() %>% 
  filter(user_login %in% result_df$user_login) %>%
  find_most_common("user_login", "taxon") %>% 
  left_join(result_df[, c("user_login", "type"), by = "user_login"])

ggplot(top_taxon_data) +
  geom_density(aes(fraction, group = type, col = type)) +
  xlab("Pct. of observations which are the top taxon") +
  ylab("Density") +
  theme_minimal()

