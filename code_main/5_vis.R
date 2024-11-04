#############################################################
# 5_vis.R
#
# This file produces a variety of visualizations based on the
# outputs of script 1-4
#############################################################

library(tidyverse)
source("code_helper/functions.R")


colors <- c("Specialization" = "#f0027f",
            "Null" = "#bbbbbb",
            "Novelty" = "#386cb0")


result_df_hypos <- read_csv("intermediate/inference_decisions.csv") %>% 
  mutate(type = recode(inferred_behavior, "Favoritism" = "Specialization")) %>% 
  mutate(type = factor(type, levels = names(colors)))
result_df_values <- read_csv("intermediate/estimated_user_BIs_allPA.csv") %>% 
  select(-type) %>% 
  left_join(result_df_hypos[, c("user_login", "taxon", "type")])


user_info <- lapply(list.files("intermediate", pattern = "user_summary_",
                               full.names = TRUE), read_csv) %>% 
  bind_rows()

#### Make some example figs ####

user_dists <- bind_rows(lapply(list.files("intermediate/", pattern = 'userDist',
                                          full.names = TRUE), read_csv))

plot_user(this_taxon = "All", this_user = "catspit", user_dists, user_info)


#### 

(taxon_counts_all <- result_df_hypos %>% 
    count(taxon, type) %>%
    filter(taxon == "All") %>%
    group_by(taxon) %>% 
    mutate(pct = n / sum(n),
           taxon = paste0(taxon, " (n=", sum(n), ")")) %>% 
    ggplot(aes(type, n, fill = type)) +
    geom_col(show.legend = F) + 
    facet_wrap(~taxon, scales = "free_y") +
    scale_fill_manual("Inferred\nbehavior", values = colors) +
    xlab("") + ylab("Number of users") +
    theme_minimal())
(taxon_counts_taxon <- result_df_hypos %>% 
    count(taxon, type) %>%
    filter(taxon != "All") %>% 
    group_by(taxon) %>% 
    mutate(pct = n / sum(n),
           taxon = paste0(taxon, " (n=", sum(n), ")")) %>% 
    ggplot(aes(type, n, fill = type)) +
    geom_col(show.legend = F) + 
    facet_wrap(~taxon, scales = "free_y") +
    scale_fill_manual("Inferred\nbehavior", values = colors) +
    xlab("") + ylab("Number of users") +
    theme_minimal())

lmtx <- matrix(c(NA,2,1,2,NA,2), nrow = 2)

taxon_counts_fig <- gridExtra::arrangeGrob(taxon_counts_all, taxon_counts_taxon, layout_matrix = lmtx,
                                           widths = c(0.5, 1, 0.5), heights = c(1, 1.5))

#

ggplot(result_df_values) +
  geom_histogram(aes(bias_index, fill = type),
                 position = "identity", alpha = 0.5) +
  scale_fill_manual("Inferred\nbehavior", values = colors) +
  theme_minimal()


pointcloud_fig <- result_df_values %>% 
  left_join(user_info) %>% 
  ggplot(aes(log(nobsT), bias_index, ymin = Q025, ymax = Q975,
             col = type)) +
  geom_pointrange() +
  facet_wrap(~taxon) +
  xlab("Log number of obs.") +
  ylab("Bias index (95%CI)") +
  scale_color_manual("Inferred\nbehavior", values = colors) + 
  theme_minimal()

alphas <- c(1, 0.5, 1)
names(alphas) <- names(colors)
(se_fig_All <- result_df_values %>% 
    filter(taxon == "All") %>% 
    left_join(user_info) %>% 
    ggplot(aes(log(nobsT), Q975-Q025, group = taxon, color = type, alpha = type)) +
    geom_point(show.legend = FALSE) +
    # geom_smooth(color = "black", alpha = 0.5,
    #             method = lm, formula = y~ poly(x, 2)) +
    xlab("Log number of obs.") +
    ylab("Width of 95%CI on bias index") +
    scale_color_manual("Inferred\nbehavior", values = colors) +
    scale_alpha_manual("Inferred\nbehavior", values = alphas) +
    facet_wrap(~taxon) +
    theme_minimal())
(se_fig_taxa <- result_df_values %>% 
    filter(taxon != "All") %>% 
    left_join(user_info) %>% 
    ggplot(aes(log(nobsT), Q975-Q025, group = taxon, color = type, alpha = type)) +
    geom_point() +
    # geom_smooth(color = "black", alpha = 0.5,
    #             method = lm, formula = y~ poly(x, 2)) +
    xlab("Log number of obs.") +
    ylab("Width of 95%CI on bias index") +
    scale_color_manual("Inferred\nbehavior", values = colors) +
    scale_alpha_manual("Inferred\nbehavior", values = alphas) +
    facet_wrap(~taxon) +
    theme_minimal())

lmtx <- matrix(c(NA,2,1,2,NA,2), nrow = 2)

se_fig <- gridExtra::arrangeGrob(se_fig_All, se_fig_taxa, layout_matrix = lmtx,
                                 widths = c(0.5, 1, 0.5), heights = c(1, 1.5))


#### Figure illustrating the NS values we estimated ####


NS_values_figure_all <- left_join(
  result_df_hypos %>% select(user_login, taxon, type),
  result_df_values %>% select(user_login, taxon, bias_index)
) %>% 
  filter(taxon == "All") %>% 
  ggplot() + 
  geom_histogram(aes(bias_index, fill = type),
                 binwidth = 0.1, boundary = 0) +
  geom_vline(xintercept = 0, col = "black") +
  scale_fill_manual("Inferred\nbehavior", values = colors) +
  theme_minimal() + xlab("Point estimate of NS index") + ylab("Num. observers") +
  ggtitle("(A) All taxa")

NS_values_figure_taxa <- left_join(
  result_df_hypos %>% select(user_login, taxon, type),
  result_df_values %>% select(user_login, taxon, bias_index)
) %>% 
  filter(taxon != "All") %>% 
  ggplot() + 
  geom_histogram(aes(bias_index, fill = type),
                 binwidth = 0.1,
                 boundary = 0, show.legend = F) +
  geom_vline(xintercept = 0, col = "black") +
  scale_fill_manual("Inferred\nbehavior", values = colors) +
  facet_wrap(~taxon, scales = "free_y") +
  theme_minimal() + xlab("Point estimate of NS index") + ylab("Num. observers") +
  ggtitle("(B) Taxon-specific")

NS_values <- gridExtra::arrangeGrob(NS_values_figure_all, NS_values_figure_taxa,
                                    ncol = 1, heights = c(0.6, 1))
ggsave("figs/NS_values.jpg", NS_values,
       width= 7, height  = 6)

#### How do users compare across taxa? ####
result_all <- result_df_hypos %>% 
  select(user_login, taxon, type) %>% 
  filter(taxon == "All") %>% 
  rename(type_all = type) %>% 
  select(-taxon)

# Calculate the percentages
percentage_df <- result_df_hypos %>%
  filter(user_login %in% result_all$user_login, taxon != "All") %>%
  select(user_login, taxon, type) %>% 
  left_join(result_all) %>%
  count(type_all, type) %>%
  mutate(percentage = n / sum(n))
percentage_df_B <- result_df_hypos %>%
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
  xlab("Overall behavior") +
  ylab("Taxon-specific behavior") +
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
  xlab("Overall behavior") +
  ylab("Taxon-specific behavior") +
  theme_minimal() + 
  coord_fixed() +
  facet_wrap(~taxon, nrow = 1) +
  scale_color_manual(values = c("black" = "black", "white" = "white")) +
  ggtitle("B. Separate by taxon")


cross_taxon_fig <- gridExtra::arrangeGrob(grobs  = list(
  cross_taxon_fig_A, cross_taxon_fig_B
), nrow = 2)

#### save some plots ####

ggsave("figs/fig1_counts.jpg", taxon_counts_fig, width = 7, height = 6)
ggsave("figs/fig2_ptcloud.jpg", pointcloud_fig, width = 9, height = 6)
ggsave("figs/fig3_SEs.jpg", se_fig, width = 9, height = 7)
ggsave("figs/fig4_crossplot.jpg", cross_taxon_fig, width = 12, height = 5.5)

#### Representative user plots (supplement?) ####

user_values <- bind_rows(lapply(list.files("intermediate/", pattern = 'user_summary_',
                                           full.names = TRUE), read_csv))
user_dists <- bind_rows(lapply(list.files("intermediate/", pattern = 'userDist',
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

plot_user(this_taxon = "All", this_user = "pachogut", 
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
library(tidyverse)
library(tidyterra)

# Make hex counts
source("code_helper/functions.R")
data_user_species <- read_csv("intermediate/data_user_species.csv")

dat <- associate_wgrid(input_dat = data_user_species,
                       hex_short_diameter = 25000)


usa_map <- geodata::gadm(country = "USA", level = 1, path = "data/geodata")
PA_map <- usa_map[usa_map$NAME_1 == "Pennsylvania", ]

obs_per_hex <- dat %>% count(gridcell)
specs_per_hex <- dat %>% 
  distinct(gridcell, scientific_name) %>% 
  count(gridcell) %>% 
  rename(nspec = n)
users_per_hex <- dat %>% 
  distinct(gridcell, user_login) %>% 
  count(gridcell) %>% 
  rename(nuser = n)
max_user_per_hex <- dat %>% 
  count(gridcell, user_login) %>% 
  group_by(gridcell) %>% 
  summarize(prop_max = max(n) / sum(n)) %>% 
  ungroup()

hex_grid_vect <- vect(hex_grid)
hex_grid_vect$gridcell <- 1:nrow(hex_grid_vect)
hex_grid_vect <- left_join(hex_grid_vect, obs_per_hex, by = "gridcell") %>% 
  left_join(specs_per_hex, by = "gridcell") %>% 
  left_join(users_per_hex, by = "gridcell") %>% 
  left_join(max_user_per_hex, by = "gridcell") %>% 
  filter(n > 0) %>% 
  terra::project(PA_map)

gridmap <- ggplot() + 
  geom_spatvector(data = hex_grid_vect, aes(fill = n)) +
  geom_spatvector(data = PA_map, fill = NA, col = "black", linewidth = 0.8) +
  theme_minimal() +
  scale_fill_viridis_c("Number of \niNaturalist obs.", trans = "log",
                       breaks = c(1, 20, 400, 8000))
gridmap2 <- ggplot() +
  geom_spatvector(data = hex_grid_vect, aes(fill = nspec)) +
  geom_spatvector(data = PA_map, fill = NA, col = "black", linewidth = 0.8) +
  theme_minimal() +
  scale_fill_viridis_c("Number of \nspecies observed", trans = "log",
                       breaks = c(1, 8, 160, 2400))
gridmap3 <- ggplot() +
  geom_spatvector(data = hex_grid_vect, aes(fill = nuser)) +
  geom_spatvector(data = PA_map, fill = NA, col = "black", linewidth = 0.8) +
  theme_minimal() +
  scale_fill_viridis_c("Number of \nunique observers", trans = "log",
                       breaks = c(1, 12, 120, 1400))
gridmap4 <- ggplot() +
  geom_spatvector(data = hex_grid_vect, aes(fill = prop_max)) +
  geom_spatvector(data = PA_map, fill = NA, col = "black", linewidth = 0.8) +
  theme_minimal() +
  scale_fill_viridis_c("Pct. obs. from\ntop user")


hex_map <- gridExtra::arrangeGrob(grobs = list(gridmap, gridmap2, gridmap3, gridmap4),
                                  ncol = 1)

ggsave("figs/hexmap.jpg", hex_map, width = 6, height = 9)


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
all_dat <- read_csv("intermediate/data_user_species.csv") %>% 
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




spec_counts_by_taxon <- spec_counts %>% 
  filter(taxon != "All") %>% 
  ggplot() +
  geom_col(aes(taxon, nspec)) +
  theme_minimal() + xlab("Taxonomic group") + ylab("Num. species") +
  ggtitle("Number of species per taxonomic group") +
  scale_y_log10(breaks = c(6, 20, 60, 200, 600, 2000, 6000)) + 
  theme(panel.grid.minor = element_blank())
ggsave("figs/taxon_spec_counts.jpg", spec_counts_by_taxon,
       width = 5, height = 4)

#### Visualize top orders of favoritists ####
# Make hex counts
source("code_helper/functions.R")
data_user_species <- read_csv("intermediate/data_user_species.csv")

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

