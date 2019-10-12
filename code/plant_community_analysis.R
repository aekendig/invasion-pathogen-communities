## Goal: examine relationships between invasion and richness/evenness/density


#### set up ####

# clear everything except final data
rm(list = ls())

# load libraries
library(tidyverse)
library(vegan)
library(GGally)

# import data
plant15.1 <- read_csv("./data/background_plants_raw_transect_2015.csv")
plant16T.1 <- read_csv("./data/background_plants_raw_transect_2016.csv")
plant16C.3 <- read_csv("./data/background_plants_raw_competition_2016.csv")

bg15 <- read_csv("./data/background_plants_transect_2015.csv")
bg16T <- read_csv("./data/background_plants_transect_2016.csv")
bg16C <- read_csv("./data/background_plants_competition_2016.csv")


#### edit data ####

# community datasets
unique(plant15.1$sp) # all are grasses
unique(plant16T.1$sp) # forb and fab not grasses
unique(plant16C.3$sp) # forb and fab not grasses

dat15 <- plant15.1 %>%
  group_by(year, experiment, plot, subplot) %>%
  summarise(host.density = sum(abundance),
            host.richness = specnumber(abundance),
            host.diversity = diversity(abundance)) %>%
  mutate(evenness = host.diversity / log(host.richness)) %>%
  ungroup() %>%
  full_join(bg15)

dat16T <- plant16T.1 %>%
  filter(!(sp %in% c("Forb", "Fabaceae"))) %>%
  group_by(year, experiment, plot, subplot) %>%
  summarise(host.density = sum(abundance),
            host.richness = specnumber(abundance),
            host.diversity = diversity(abundance)) %>%
  mutate(evenness = host.diversity / log(host.richness)) %>%
  ungroup() %>%
  full_join(bg16T)

dat16C <- plant16C.3 %>%
  filter(!(sp %in% c("Forb", "Fabaceae"))) %>%
  group_by(year, experiment, plot, subplot) %>%
  summarise(host.density = sum(abundance),
            host.richness = specnumber(abundance),
            host.diversity = diversity(abundance)) %>%
  mutate(evenness = host.diversity / log(host.richness)) %>%
  ungroup() %>%
  full_join(bg16C)

# species datasets
datw15 <- plant15.1 %>%
  select(-c(species, BGStatus, CalfloraInvasive)) %>%
  spread(key = sp, value = abundance)

datw16T <- plant16T.1 %>%
  filter(!(sp %in% c("Forb", "Fabaceae", "POA"))) %>%
  select(-c(sp, BGStatus, CalfloraInvasive)) %>%
  spread(key = species, value = abundance)

datw16C <- plant16C.3 %>%
  filter(!(sp %in% c("Forb", "Fabaceae"))) %>%
  select(-c(sp, BGStatus, CalfloraInvasive)) %>%
  spread(key = species, value = abundance)


#### correlations ####

dat15 %>%
  select(-c(year:subplot)) %>%
  ggpairs()
# non-native highly correlated with host density, but not much else, slightly positively correlated with richness

dat16T %>%
  select(-c(year:subplot)) %>%
  ggpairs()
# similar to above

dat16C %>%
  select(-c(year:subplot, bg.species, competition.density)) %>%
  ggpairs()
# similar to above, decreased diversity and evenness

datw15 %>%
  select(-c(year:subplot)) %>%
  ggpairs()
# no strong relationships

datw16T %>%
  select(-c(year:subplot)) %>%
  ggpairs()
# no strong relationships

datw16C %>%
  select(-c(year:subplot)) %>%
  ggpairs()
# no strong relationships