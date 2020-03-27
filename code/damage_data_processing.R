#### set up ####

# clear everything except final data
rm(list = ls())

# load libraries
library(tidyverse)

# import data
dam15m <- read_csv("./data/JRBPmidMarPathPrcntDamData_5Apr2015ERS.csv")
dam15a1 <- read_csv("./data/JRBPmidAprPathPrcntDamCulturesData_25Apr2015ERS.csv") # transect 2015
dam15a2 <- read_csv("./data/JRBPtransects11to14PathPrcntDam28Apr_7May2015ERS.csv")
dam16 <- read_csv("./data/JRBP_PathogenDamage_21July2016_Full.csv",
                  col_types = cols(Transect = col_double(), Plot = col_character())) # transect and competition 2016
compPlot <- read_csv("./data/2016ExperimentPlotList_26Apr16_PlotAssignments.csv")


#### examine data ####

# 2015 March
dam15m
unique(dam15m$notes)
unique(dam15m$trans)
unique(dam15m$plot) # not sure why some start with P and others with T, but they don't seem to be overlapping
unique(dam15m$spcode)
length(unique(dam15m$indiv))
nrow(dam15m)
# one duplicate individual ID

# check on potential duplication
dup15m <- dam15m %>%
  mutate(dups = duplicated(indiv)) %>%
  filter(dups == T)
filter(dam15m, indiv %in% dup15m$indiv)

# note that it's duplicated
filter(dam15m, indiv %in% dup15m$indiv) %>%
  select(notes) %>%
  data.frame()

# 2015 April 1
dam15a1
unique(dam15a1$notes)
unique(dam15a1$Transect)
unique(dam15a1$Plot)
unique(dam15a1$Species)
length(unique(dam15a1$Individual))
nrow(dam15a1) 
#each row is a unique entry

# check plants mislabelled
dam15a1 %>%
  filter(substr(notes, 1, 13) == "plant was mis") %>%
  select(Species, notes)

# 2015 April 2
dam15a2
unique(dam15a2$notes)
unique(dam15a2$Transect)
unique(dam15a2$Species)
nrow(dam15a2) 
# no ID, check for row duplication

# check roadside plant
filter(dam15a2, notes == "this AB individual was on roadside not along transect")
# 5 plants

# 2016
dam16
unique(dam16$Notes)
unique(dam16$Transect)
unique(dam16$Plot)
unique(dam16$Species)
length(unique(dam16$Plant.ID))
nrow(dam16) #some don't have plant ID's or they're duplicated
sum(duplicated(dam16$Plant.ID)) # 1001 duplicated

# look at notes closer
dam16 %>%
  filter(!is.na(Notes)) %>%
  select(Species, Plant.ID, Notes) %>%
  data.frame()

# one plant assessed as two
filter(dam16, Plant.ID %in% c("X437", "X438")) %>%
  data.frame()
# one is all NA - will be removed

# check duplicates
dup16 <- dam16 %>%
  group_by(Experiment, Transect, Plot, Treatment, Species, Plant.ID) %>%
  summarise(n = n()) %>%
  filter(n > 1) %>%
  data.frame()
# 88 duplicates when other info is included
dup16 # all except one are missing an ID, could be true duplicates

dam16 %>%
  inner_join(select(dup16, -n)) %>%
  select(leaf1, leaf2, leaf3, leaf4, leaf5, leaf6) %>%
  unique()
# none of them have data - will be removed

# look at types of scores assigned
dam15m %>%
  select(leaf1, leaf2, leaf3, leaf4, leaf5, leaf6) %>%
  rbind(select(dam15a1, leaf1, leaf2, leaf3, leaf4, leaf5, leaf6)) %>%
  rbind(select(dam15a2, leaf1, leaf2, leaf3, leaf4, leaf5, leaf6)) %>%
  rbind(select(dam16, leaf1, leaf2, leaf3, leaf4, leaf5, leaf6)) %>%
  gather(key = leaf, value = score) %>%
  select(score) %>%
  unique()

#### edit data ####

#From Erin's code (2016 Seed Data Processing_EM.R): mean.dam = take the mean of all damage values for the different leaves on a plant and divide by 100, prop.dam = the proportion of leaves with damage greater than 0, remove NA's (don't count as healthy leaves)

# function to process plant-level data
plant_fun <- function(dat){
  
  # format damage data
  dat2 <- dat %>%
    select(year, month, experiment, plot, subplot, host, leaf1, leaf2, leaf3, leaf4, leaf5, leaf6) %>%
    rowwise() %>%
    mutate(mean.dam = mean(c(leaf1, leaf2, leaf3, leaf4, leaf5, leaf6), na.rm = T)/100,
           leaf.i.1 = as.numeric(leaf1 > 0),
           leaf.i.2 = as.numeric(leaf2 > 0),
           leaf.i.3 = as.numeric(leaf3 > 0),
           leaf.i.4 = as.numeric(leaf4 > 0),
           leaf.i.5 = as.numeric(leaf5 > 0),
           leaf.i.6 = as.numeric(leaf6 > 0),
           leaf.s.1 = leaf1 / 100,
           leaf.s.2 = leaf2 / 100,
           leaf.s.3 = leaf3 / 100,
           leaf.s.4 = leaf4 / 100,
           leaf.s.5 = leaf5 / 100,
           leaf.s.6 = leaf6 / 100,
           leaves.dam = sum(c(leaf.i.1, leaf.i.2, leaf.i.3, leaf.i.4, leaf.i.5, leaf.i.6), na.rm = T),
           leaves.tot = sum(!is.na(c(leaf.i.1, leaf.i.2, leaf.i.3, leaf.i.4, leaf.i.5, leaf.i.6))),
           prop.dam = leaves.dam/leaves.tot,
           grass.group = ifelse(host %in% c("SP","EG"), "native\nperennial", "non-native\nannual"),
           nonnative = ifelse(grass.group == "non-native\nannual", 1, 0)) %>%
    filter(!is.na(mean.dam)) %>%
    ungroup() %>%
    select(-c(leaf1:leaf6))
  
  # output
  return(dat2)
}

# function to process leaf-level data
leaf_fun <- function(dat){
  
  # make long by leaf
  dat2 <- dat %>%
    select(year, month, experiment, plot, subplot, host, grass.group, nonnative, leaf.s.1:leaf.s.6, leaf.i.1:leaf.i.6) %>%
    mutate(plant = 1:nrow(dat)) %>%
    gather(key, value, -c(year:nonnative, plant)) %>%
    extract(key, c("metric", "leaf"), "(leaf\\..)\\.(.)") %>%
    spread(metric, value) %>%
    rename(infected = leaf.i, surface = leaf.s) %>%
    filter(!is.na(infected))
  
  # output
  return(dat2)
}

# 2015 March
dam15mplant <- dam15m %>%
  mutate(year = 2015,
         month = "March",
         experiment = "transect",
         subplot = paste("transect", substr(plot, 2, 4), sep = "_"),
         plot = paste("transect", gsub("[^[:digit:]]", "", plot), sep = "_")) %>%
  rename(host = spcode) %>%
  plant_fun

dam15mleaf <- leaf_fun(dam15mplant)

# 2015 April 1
dam15a1plant <- dam15a1 %>%
  mutate(year = 2015,
         month = "April",
         experiment = "transect",
         plot = paste("transect", gsub("[^[:digit:]]", "", Transect), sep = "_"),
         subplot = paste(plot, Plot, sep = "")) %>%
  rename(host = Species) %>%
  plant_fun

dam15a1leaf <- leaf_fun(dam15a1plant)

# 2015 April 2
dam15a2plant <- dam15a2 %>%
  filter(is.na(notes) | notes != "this AB individual was on roadside not along transect") %>%
  mutate(year = 2015,
         month = "April",
         experiment = "transect",
         plot = paste("transect", Transect, sep = "_"),
         subplot = plot) %>%
  rename(host = Species) %>%
  plant_fun

dam15a2leaf <- leaf_fun(dam15a2plant)

# 2016 transect
dam16Tplant <- dam16 %>%
  filter(Experiment == "Transect" & Treatment == "ambient") %>%
  mutate(year = 2016,
         month = NA_character_,
         experiment = "transect",
         plot = paste("transect", gsub("[^[:digit:]]", "", Transect), sep = "_"),
         subplot = paste(plot, Plot, sep = "")) %>%
  rename(host = Species) %>%
  plant_fun

dam16Tleaf <- leaf_fun(dam16Tplant)

# 2016 competition
dam16Cplant <- dam16 %>%
  filter(Experiment == "Competition" & Treatment == "ambient") %>%
  select(Plot, BackgroundSpecies, Density, Species, leaf1, leaf2, leaf3, leaf4, leaf5, leaf6) %>%
  rowwise() %>%
  mutate(mean.dam = mean(c(leaf1, leaf2, leaf3, leaf4, leaf5, leaf6), na.rm = T)/100,
         leaf.i.1 = as.numeric(leaf1 > 0),
         leaf.i.2 = as.numeric(leaf2 > 0),
         leaf.i.3 = as.numeric(leaf3 > 0),
         leaf.i.4 = as.numeric(leaf4 > 0),
         leaf.i.5 = as.numeric(leaf5 > 0),
         leaf.i.6 = as.numeric(leaf6 > 0),
         leaf.s.1 = leaf1 / 100,
         leaf.s.2 = leaf2 / 100,
         leaf.s.3 = leaf3 / 100,
         leaf.s.4 = leaf4 / 100,
         leaf.s.5 = leaf5 / 100,
         leaf.s.6 = leaf6 / 100,
         leaves.dam = sum(c(leaf.i.1, leaf.i.2, leaf.i.3, leaf.i.4, leaf.i.5, leaf.i.6), na.rm = T),
         leaves.tot = sum(!is.na(c(leaf.i.1, leaf.i.2, leaf.i.3, leaf.i.4, leaf.i.5, leaf.i.6))),
         prop.dam = leaves.dam/leaves.tot,
         grass.group = ifelse(Species %in% c("SP","EG"), "native\nperennial", "non-native\nannual"),
         nonnative = ifelse(grass.group == "non-native\nannual", 1, 0),
         year = 2016,
         experiment = "competition",
         plot = paste("competition", Plot, sep = "_"),
         subplot = plot) %>%
  rename(host = Species, 
         bg.species = BackgroundSpecies,
         competition.density = Density) %>%
  ungroup() %>%
  filter(!is.na(mean.dam)) %>%
  select(year, experiment, plot, subplot, bg.species, competition.density, host, grass.group, nonnative, mean.dam, leaf.i.1:prop.dam)

dam16Cleaf <- dam16Cplant %>%
  select(year:leaf.s.6) %>%
  mutate(plant = 1:nrow(dam16Cplant)) %>%
  gather(key, value, -c(year:nonnative, plant)) %>%
  extract(key, c("metric", "leaf"), "(leaf\\..)\\.(.)") %>%
  spread(metric, value) %>%
  rename(infected = leaf.i, surface = leaf.s) %>%
  filter(!is.na(infected))

#### combine data ####

dam15aplant = rbind(dam15a1plant, dam15a2plant)
dam15aleaf = rbind(dam15a1leaf, dam15a2leaf)


#### outputs ####

write_csv(dam15mplant, "./data/damage_plant_transect_2015_march.csv")
write_csv(dam15aplant, "./data/damage_plant_transect_2015_april.csv")
write_csv(dam16Tplant, "./data/damage_plant_transect_2016.csv")
write_csv(dam16Cplant, "./data/damage_plant_competition_2016.csv")

write_csv(dam15mleaf, "./data/damage_leaf_transect_2015_march.csv")
write_csv(dam15aleaf, "./data/damage_leaf_transect_2015_april.csv")
write_csv(dam16Tleaf, "./data/damage_leaf_transect_2016.csv")
write_csv(dam16Cleaf, "./data/damage_leaf_competition_2016.csv")