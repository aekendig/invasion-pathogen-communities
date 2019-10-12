## Goal: format damage data


#### set up ####

# clear everything except final data
rm(list = ls())

# load libraries
library(tidyverse)

# import data
dam15 <- read_csv("./data/JRBPmidAprPathPrcntDamCulturesData_25Apr2015ERS.csv") # transect 2015
dam16 <- read_csv("./data/JRBP_PathogenDamage_21July2016_Full.csv",
                  col_types = cols(Transect = col_double(), Plot = col_character())) # transect and competition 2016
compPlot <- read_csv("./data/2016ExperimentPlotList_26Apr16_PlotAssignments.csv")


#### edit data ####

# examine data
dam15 # format Transect and Plot, need Species
unique(dam15$Transect)
unique(dam15$Plot)
unique(dam15$Species)
length(unique(dam15$Individual))
nrow(dam15) #each row is a unique entry

dam16
unique(dam16$Transect)
unique(dam16$Plot)
unique(dam16$Species)
length(unique(dam16$Plant.ID))
nrow(dam16) #some don't have plant ID's or they're duplicated
sum(duplicated(dam16$Plant.ID)) # 1001 duplicated

# check on potential duplication
dup16 <- dam16 %>%
  group_by(Experiment, Transect, Plot, Treatment, Species, Plant.ID) %>%
  summarise(n = n()) %>%
  filter(n > 1) %>%
  data.frame()
dup16 # all except one are missing an ID, could be true duplicates

dam16 %>%
  inner_join(select(dup16, -n)) %>%
  rowwise() %>%
  mutate(dam = sum(c(leaf1, leaf2, leaf3, leaf4, leaf5, leaf6), na.rm = T)) %>%
  filter(dam > 0)
# none of them have damage

#From Erin's code (2016 Seed Data Processing_EM.R): mean.dam = take the mean of all damage values for the different leaves on a plant and divide by 100, prop.dam = the proportion of leaves with damage greater than 0, remove NA's (don't count as healthy leaves)

# 2015 data
# calculate mean damage and prop damage
# remove samples that had no damage data
# add columns
# remove unnecessary columns
dam15.1 <- dam15 %>%
  select(Transect, Plot, Species, leaf1, leaf2, leaf3, leaf4, leaf5, leaf6) %>%
  rowwise() %>%
  mutate(mean.dam = mean(c(leaf1, leaf2, leaf3, leaf4, leaf5, leaf6), na.rm = T)/100,
         leaf1f = as.numeric(leaf1 > 0),
         leaf2f = as.numeric(leaf2 > 0),
         leaf3f = as.numeric(leaf3 > 0),
         leaf4f = as.numeric(leaf4 > 0),
         leaf5f = as.numeric(leaf5 > 0),
         leaf6f = as.numeric(leaf6 > 0),
         leaves.dam = sum(c(leaf1f, leaf2f, leaf3f, leaf4f, leaf5f, leaf6f), na.rm = T),
         leaves.tot = sum(!is.na(c(leaf1f, leaf2f, leaf3f, leaf4f, leaf5f, leaf6f))),
         prop.dam = leaves.dam/leaves.tot) %>%
  filter(!is.na(mean.dam)) %>%
  mutate(year = 2015,
         experiment = "transect",
         plot = paste("transect", gsub("[^[:digit:]]", "", Transect), sep = "_"),
         subplot = paste(plot, Plot, sep = "")) %>%
  rename(host = Species) %>%
  ungroup()

dam15plant <- dam15.1 %>%
  select(year, experiment, plot, subplot, host, mean.dam, leaves.dam, leaves.tot, prop.dam)

dam15leaf <- dam15.1 %>%
  select(year, experiment, plot, subplot, host, leaf1f:leaf6f) %>%
  mutate(plant = 1:nrow(dam15.1)) %>%
  gather(key = leaf, value = infected, -c(year:host, plant))
  
# separate 2016 data by experiment
dam16C <- filter(dam16, Experiment == "Competition")
dam16T <- filter(dam16, Experiment == "Transect")

# check treatments
unique(dam16T$Treatment)

dam16C %>%
  group_by(Plot) %>%
  filter(length(unique(Treatment)) > 1) %>%
  summarise(treatments = paste(unique(Treatment), collapse = "_")) %>%
  data.frame() # no duplicate treatments

# 2016 transect data
# calculate mean damage and prop damage
# remove samples that had no damage data
# add columns
# remove unnecessary columns
dam16T.1 <- dam16T %>%
  filter(Treatment == "ambient") %>%
  select(Transect, Plot, Species, leaf1, leaf2, leaf3, leaf4, leaf5, leaf6) %>%
  rowwise() %>%
  mutate(mean.dam = mean(c(leaf1, leaf2, leaf3, leaf4, leaf5, leaf6), na.rm = T)/100,
         leaf1f = as.numeric(leaf1 > 0),
         leaf2f = as.numeric(leaf2 > 0),
         leaf3f = as.numeric(leaf3 > 0),
         leaf4f = as.numeric(leaf4 > 0),
         leaf5f = as.numeric(leaf5 > 0),
         leaf6f = as.numeric(leaf6 > 0),
         leaves.dam = sum(c(leaf1f, leaf2f, leaf3f, leaf4f, leaf5f, leaf6f), na.rm = T),
         leaves.tot = sum(!is.na(c(leaf1f, leaf2f, leaf3f, leaf4f, leaf5f, leaf6f))),
         prop.dam = leaves.dam/leaves.tot) %>%
  filter(!is.na(mean.dam)) %>%
  mutate(year = 2016,
         experiment = "transect",
         plot = paste("transect", gsub("[^[:digit:]]", "", Transect), sep = "_"),
         subplot = paste(plot, Plot, sep = "")) %>%
  rename(host = Species) %>%
  ungroup()

dam16Tplant <- dam16T.1 %>%
  select(year, experiment, plot, subplot, host, mean.dam, leaves.dam, leaves.tot, prop.dam)

dam16Tleaf <- dam16T.1 %>%
  select(year, experiment, plot, subplot, host, leaf1f:leaf6f) %>%
  mutate(plant = 1:nrow(dam16T.1)) %>%
  gather(key = leaf, value = infected, -c(year:host, plant))

# 2016 competition data
# calculate mean damage and prop damage
# remove samples that had no damage data
# add columns
# remove unnecessary columns
dam16C.1 <- dam16C %>%
  filter(Treatment == "ambient") %>%
  select(Plot, BackgroundSpecies, Density, Species, leaf1, leaf2, leaf3, leaf4, leaf5, leaf6) %>%
  rowwise() %>%
  mutate(mean.dam = mean(c(leaf1, leaf2, leaf3, leaf4, leaf5, leaf6), na.rm = T)/100,
         leaf1f = as.numeric(leaf1 > 0),
         leaf2f = as.numeric(leaf2 > 0),
         leaf3f = as.numeric(leaf3 > 0),
         leaf4f = as.numeric(leaf4 > 0),
         leaf5f = as.numeric(leaf5 > 0),
         leaf6f = as.numeric(leaf6 > 0),
         leaves.dam = sum(c(leaf1f, leaf2f, leaf3f, leaf4f, leaf5f, leaf6f), na.rm = T),
         leaves.tot = sum(!is.na(c(leaf1f, leaf2f, leaf3f, leaf4f, leaf5f, leaf6f))),
         prop.dam = leaves.dam/leaves.tot) %>%
  filter(!is.na(mean.dam)) %>%
  mutate(year = 2016,
         experiment = "competition",
         plot = paste("competition", Plot, sep = "_"),
         subplot = plot) %>%
  rename(host = Species, 
         bg.species = BackgroundSpecies,
         competition.density = Density) %>%
  ungroup()

dam16Cplant <- dam16C.1 %>%
  select(year, experiment, plot, subplot, bg.species, competition.density, host, mean.dam, leaves.dam, leaves.tot, prop.dam)

dam16Cleaf <- dam16C.1 %>%
  select(year, experiment, plot, subplot, bg.species, competition.density, host, leaf1f:leaf6f) %>%
  mutate(plant = 1:nrow(dam16C.1)) %>%
  gather(key = leaf, value = infected, -c(year:host, plant))

  
#### outputs ####
write_csv(dam15plant, "./data/damage_plant_transect_2015.csv")
write_csv(dam16Tplant, "./data/damage_plant_transect_2016.csv")
write_csv(dam16Cplant, "./data/damage_plant_competition_2016.csv")

write_csv(dam15leaf, "./data/damage_leaf_transect_2015.csv")
write_csv(dam16Tleaf, "./data/damage_leaf_transect_2016.csv")
write_csv(dam16Cleaf, "./data/damage_leaf_competition_2016.csv")