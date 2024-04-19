#############################################################
# 1_process_data.R
#
# This file contains the code we used to process the iNaturalist
# data before conducting all analyses
#############################################################


# quality_grade=research&identifications=any&place_id=42&d1=2008-01-01&d2=2018-12-31
library(tidyverse)
library(lubridate)

# Because of the size of the data, we retrieved iNaturalist data from multiple
# downloads, meaning that we had to process multiple input files:

# Load in each data file, apply filters, and concatenate
data1 <- read_csv("data/observations-365569.csv") %>%
  filter(!is.na(scientific_name)) %>%
  filter(quality_grade == "research") %>% 
  filter(scientific_name != "") %>%
  mutate(observed_on = as.Date(observed_on, format = "%Y-%m-%d")) %>% 
  rename(species = scientific_name,
         eventDate = observed_on) %>% 
  mutate(year = lubridate::year(eventDate))  %>%
  filter(!is.na(eventDate))

data2 <- read_csv("data/observations-365582.csv") %>%
  filter(!is.na(scientific_name)) %>%
  filter(quality_grade == "research") %>% 
  filter(scientific_name != "") %>%
  mutate(observed_on = as.Date(observed_on, format = "%Y-%m-%d")) %>% 
  rename(species = scientific_name,
         eventDate = observed_on) %>% 
  mutate(year = lubridate::year(eventDate))  %>%
  filter(!is.na(eventDate))

data3 <- read_csv("data/observations-366081.csv") %>%
  filter(!is.na(scientific_name)) %>%
  filter(quality_grade == "research") %>% 
  filter(scientific_name != "") %>%
  mutate(observed_on = as.Date(observed_on, format = "%Y-%m-%d")) %>% 
  rename(species = scientific_name,
         eventDate = observed_on) %>% 
  mutate(year = lubridate::year(eventDate))  %>%
  filter(!is.na(eventDate))

data <- rbind.data.frame(data1, data2, data3) %>% 
  filter(!duplicated(id))

# Summarize by species
bySpecies <- data %>%
  group_by(species, iconic_taxon_name) %>%
  summarise(species_count = n())
bySpecies <- bySpecies %>% mutate(obs_prev = species_count / nrow(data))

data$eventDate = as.Date(data$eventDate)


# how many total observations for each user?
user_info <- data %>%
  group_by(user_login) %>%
  summarize(
    first_obs = min(eventDate),
    nobsT = n(), .groups = "drop"
  ) %>%
  arrange(desc(nobsT))

# how many total observations for each user (by taxon)?
user_infoT <- data %>%
  group_by(user_login, iconic_taxon_name) %>%
  summarize(
    first_obs = min(eventDate),
    nobsT = n(), .groups = "drop"
  ) %>%
  arrange(desc(nobsT)) 

## label each observation with the number
data_user <- data %>%
  group_by(user_login) %>%
  arrange(eventDate) %>%
  mutate(user_id_ct = row_number()) %>%
  left_join(user_info, by = c("user_login" = "user_login"))

## is duplicate up until now?
data_user_species <- data %>%
  group_by(user_login, species) %>%
  arrange(eventDate) %>%
  mutate(user_species_id_ct = row_number()) %>%
  left_join(user_info, by = c("user_login" = "user_login"))

data_user_species <- data_user_species %>% 
  mutate(isNew = (user_species_id_ct == 1))


## number unique up to this point, number obs up to this point
data_user_species <- data_user_species %>%
  group_by(user_login) %>%
  arrange(eventDate) %>%
  mutate(dummy = 1) %>%
  mutate(numUnique = cumsum(isNew), numObs = cumsum(dummy))


## rest of the species accumulation curve as reference?
species_info <- data_user %>%
  group_by(species) %>%
  summarize(
    first_obs = min(eventDate),
    nobsT = n(), .groups = "drop"
  )

data_overall_species <- data %>%
  group_by(species) %>%
  arrange(eventDate) %>%
  mutate(species_id_ct = row_number())

data_overall_species <- data_overall_species %>% mutate(isNew = (species_id_ct == 1))


## number unique up to this point, number obs up to this point
data_overall_species <- data_overall_species %>%
  arrange(eventDate) %>%
  mutate(dummy = 1)


data_overall_species$numUnique <- cumsum(data_overall_species$isNew)
data_overall_species$numObs <- cumsum(data_overall_species$dummy)

# Select the columns to keep
data_user_species <- data_user_species %>% 
  select(id, eventDate, species,
         user_species_id_ct, first_obs, nobsT,
         latitude, longitude,
         isNew, numUnique, numObs, iconic_taxon_name, 
         scientific_name = species, user_login = user_login)


# Write out data frames for later use
write_csv(data_user_species, "intermediate/data_user_species.csv")
write_csv(user_info, "intermediate/pa_user_info.csv")


