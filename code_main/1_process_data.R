#############################################################
# 1_process_data.R
#
# This file contains the code we used to process the iNaturalist
# data before conducting all analyses
#############################################################
# quality_grade=research&identifications=any&place_id=42&d1=2008-01-01&d2=2018-12-31

if (!dir.exists("intermediate")) dir.create("intermediate")

library(tidyverse)
library(lubridate)

data_files <- list.files("input_data/", full.names = TRUE,
                         pattern = "*.csv.zip$")

all_dat_raw <- lapply(data_files, read_csv, locale=locale(encoding="latin1")) %>% 
  bind_rows()

ggplot(all_dat_raw) +
  geom_histogram(aes(as.Date(observed_on)), binwidth = 30)

all_dat <- all_dat_raw %>% 
  filter(!duplicated(id)) %>% 
  filter(!is.na(scientific_name)) %>%
  filter(grepl(" ", scientific_name)) %>% 
  filter(quality_grade == "research") %>% 
  filter(scientific_name != "") %>%
  filter(!is.na(observed_on)) %>% 
  mutate(observed_on = as.Date(observed_on, format = "%Y-%m-%d")) %>% 
  rename(species = scientific_name,
         eventDate = observed_on) %>% 
  mutate(year = lubridate::year(eventDate))

observers_to_drop <- all_dat_raw %>% 
  distinct(user_login, license) %>% 
  filter(!license %in% c("CC0", "CC-BY", "CC-BY-NC"))

all_dat <- all_dat %>% 
  filter(!user_login %in% observers_to_drop$user_login)



bySpecies <- all_dat %>%
  group_by(species, iconic_taxon_name) %>%
  summarise(species_count = n())
bySpecies <- bySpecies %>% mutate(obs_prev = species_count / nrow(all_dat))

all_dat$eventDate = as.Date(all_dat$eventDate)


## how many total observations?
user_info <- all_dat %>%
  group_by(user_login) %>%
  summarize(
    first_obs = min(eventDate),
    nobsT = n(), .groups = "drop"
  ) %>%
  arrange(desc(nobsT))

user_infoT <- all_dat %>%
  group_by(user_login, iconic_taxon_name) %>%
  summarize(
    first_obs = min(eventDate),
    nobsT = n(), .groups = "drop"
  ) %>%
  arrange(desc(nobsT)) 

## label each observation with the number
data_user <- all_dat %>%
  group_by(user_login) %>%
  arrange(eventDate) %>%
  mutate(user_id_ct = row_number()) %>%
  left_join(user_info, by = c("user_login" = "user_login"))

## is duplicate up until now?
data_user_species <- all_dat %>%
  group_by(user_login, species) %>%
  arrange(eventDate) %>%
  mutate(user_species_id_ct = row_number()) %>%
  left_join(user_info, by = c("user_login" = "user_login")) %>% 
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

data_overall_species <- all_dat %>%
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


data_user_species <- data_user_species %>% 
  select(id, eventDate, species,
         user_species_id_ct, first_obs, nobsT,
         latitude, longitude,
         isNew, numUnique, numObs, iconic_taxon_name, 
         scientific_name = species, user_login = user_login)

write_csv(data_user_species, "intermediate/data_user_species.csv")
write_csv(user_info, "intermediate/pa_user_info.csv")


