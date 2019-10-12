## Goal: format background plant data


#### set up ####

# clear everything except final data
rm(list = ls())

# load libraries
library(tidyverse)

# import data
plant15 <- read_csv("./data/JRBP_plant_count_combined_2015-04-10.csv")
plant16T <- read_csv("./data/transect_seeds_composition_damage.csv")
plant16C <- read_csv("./data/competition_plot_seeds_composition_damage.csv")
compPlot <- read_csv("./data/2016ExperimentPlotList_26Apr16_PlotAssignments.csv")
bgPlants <- read_csv("./data/BgPlantSizeStatus_JepsonCalFlora_072718.csv")
dat <- read_csv("./data/fungal_pathogens_2015_2017.csv") # ignore warnings if not using sentinel bin
dato <- read.csv("../Data/Data File 3 - Full Dataset.csv")


#### edit data ####

# need it to match dat
dat
unique(dat$plot)
unique(dat$subplot)
colnames(dat)
unique(dat$bg.species)
unique(dat$competition.density)

# examine data
plant15
unique(plant15$Transect) # transect only
# transect number = plot, Plot = subplot
plant16T
unique(plant16T$Transect)
unique(plant16T$Plot) # subplot
plant16C
unique(plant16C$Plot)
unique(plant16C$Treatment)
unique(plant16C$BG.spp)

# competition data structure
plant16C %>%
  group_by(Plot) %>%
  filter(length(unique(Treatment)) > 1 | length(unique(BG.spp)) > 1 | length(unique(Density)) > 1) %>%
  summarise(treatments = paste(unique(Treatment), collapse = "_"),
            BG = paste(unique(BG.spp), collapse = "_"),
            density = paste(unique(Density), collapse = "_")) %>%
  data.frame() # multiple treatments per plot number, only one bg sp and density, 46 plots

# original fungal data (before ambient removed)
dato %>%
  filter(!is.na(isolate.id) & 
           soil.type == "non-serpentine" & 
           (plant.tissue == "foliar" | is.na(plant.tissue)) &
           !is.na(host) & !(host %in% c("?", "SSP", "BRAC", "BRMA"))) %>%
  filter(experiment == "competition") %>%
  group_by(competition.plot) %>%
  filter(length(unique(environmental.treatment)) > 1 | length(unique(bg.species)) > 1 | length(unique(competition.density)) > 1) %>%
  summarise(treatments = paste(unique(environmental.treatment), collapse = "_"),
            BG = paste(unique(bg.species), collapse = "_"),
            density = paste(unique(competition.density), collapse = "_")) %>%
  data.frame() # no issues here

# check if composition differs
plant16C %>%
  select(Plot, Treatment, EG.adult:GAST) %>%
  unique()  %>%
  group_by(Plot, Treatment) %>%
  rowwise() %>%
  mutate(tot = sum(c(EG.adult, EG.seed, SP.adult, SP.seed, BH, BD, Avena, FEMY, BRADIS, Forb, Fabaceae, GAST))) %>%
  group_by(Plot) %>%
  filter(length(unique(Treatment)) > 1) %>%
  select(Plot, Treatment, tot) %>%
  spread(key = Treatment, value = tot) %>%
  data.frame() # they're the same

# look at one example
filter(plant16C, Plot == 9) %>%
  select(Plot, Treatment, EG.adult:GAST) %>%
  unique() # exactly the same
  
# data structure
plant16T %>%
  group_by(Plot) %>%
  filter(length(unique(Treatment)) > 1) %>%
  summarise(treatments = paste(unique(Treatment), collapse = "_")) %>%
  data.frame() # multiple treatments per plot letter

# check if composition differs
plant16T %>%
  select(Transect, Plot, Treatment, EG.adult:POA) %>%
  unique()  %>%
  group_by(Transect, Plot, Treatment) %>%
  rowwise() %>%
  mutate(tot = sum(c(EG.adult, SP.adult, BH, BD, Avena, Forb, Fabaceae, BRADIS, FEMY, POA))) %>%
  group_by(Transect, Plot) %>%
  filter(length(unique(Treatment)) > 1) %>%
  select(Transect, Plot, Treatment, tot) %>%
  spread(key = Treatment, value = tot) %>%
  data.frame() # these are different

# edit 2015 transect data
# add columns
# remove unnecessary columns
# make long
# combine with native/non-native info
plant15.1 <- plant15 %>%
  mutate(year = 2015,
         experiment = "transect",
         plot = paste("transect", gsub("[^[:digit:]]", "", Transect), sep = "_"),
         subplot = paste(plot, Plot, sep = "")) %>%
  select(-c(Transect:Plot)) %>%
  gather(key = "sp", value = "abundance", -c(year:subplot)) %>%
  left_join(select(bgPlants, abbreviation.2015, species, BGStatus, CalfloraInvasive) %>%
              rename(sp = abbreviation.2015))

# make sure all plants have a status
filter(plant15.1, is.na(BGStatus)) # yes

# edit 2016 transect data
# only keep ambient treatments
# remove individual data
# add columns
# remove unnecessary columns
# make long
# combine with native/non-native info
plant16T.1 <- plant16T %>% 
  filter(Treatment == "ambient") %>%
  select(Transect,Plot,EG.adult:POA) %>%
  unique()  %>%
  mutate(year = 2016,
         experiment = "transect",
         plot = paste("transect", Transect, sep = "_"),
         subplot = paste(plot, Plot, sep = "")) %>%
  select(-c(Transect, Plot)) %>%
  gather(key = "sp", value = "abundance", -c(year:subplot)) %>%
  left_join(select(bgPlants, abbreviation.2016, species, BGStatus, CalfloraInvasive) %>%
              rename(sp = abbreviation.2016))

# make sure all have a species
sum(is.na(plant16T.1$species))

# edit 2016 competition data
# remove individual data
# add plot set-up info
plant16C.1 <- plant16C %>%
  select(Plot, Treatment, BG.spp, Density, EG.adult:GAST) %>%
  unique()  %>%
  left_join(compPlot %>%
              rename(Plot = Number, PlantingType = Plot.type.code, Treatment2 = Treatment) %>%
              select(Plot, Treatment2, PlantingType))

filter(plant16C.1, Treatment != Treatment2 | (is.na(Treatment) & !is.na(Treatment2))) # 50 plots, makes sense because 4 of the plots have 3 treatments

plant16C.1 %>%
  group_by(Plot) %>%
  filter(length(unique(Treatment)) > 1 | Treatment != Treatment2) %>%
  summarise(treatments = paste(unique(Treatment), collapse = "_"),
            treatments2 = paste(unique(Treatment2), collapse = "_")) %>%
  data.frame() # all have replacements, none are NA

# only keep ambient treatments
# add columns
# combine seedling and adult
# remove unnecessary columns
# make long
# combine with native/non-native info
plant16C.2 <- plant16C.1 %>%
  select(-Treatment) %>%
  unique() %>%
  filter(Treatment2 == "ambient") %>%
  select(Plot, BG.spp, Density, EG.adult:GAST) %>%
  unique()  %>%
  mutate(year = 2016,
         experiment = "competition",
         plot = paste("competition", Plot, sep = "_"),
         subplot = plot,
         EG.adult = EG.adult + EG.seed,
         SP.adult = SP.adult + SP.seed) %>%
  select(-c(Plot, EG.seed, SP.seed)) %>%
  gather(key = "sp", value = "abundance", -c(BG.spp, Density, year:subplot)) %>%
  left_join(select(bgPlants, abbreviation.2016, species, BGStatus, CalfloraInvasive) %>%
              rename(sp = abbreviation.2016)) %>% 
  rename(bg.species = BG.spp,
         competition.density = Density)

# make sure plots match data
plant16C.2 %>% 
  select(year, experiment, plot, bg.species, competition.density) %>%
  mutate(check = "yes") %>%
  right_join(dat) %>%
  filter(year == 2016 & experiment == "competition" & is.na(check)) # merges

# make sure all have a species
sum(is.na(plant16C.2$species))

# remove columns
plant16C.3 <- select(plant16C.2, -c(bg.species, competition.density))


#### group background by origin ####

# status and invasive
plant15.1 %>%
  group_by(BGStatus, CalfloraInvasive) %>%
  summarise(tot = sum(abundance))

plant16T.1 %>%
  group_by(BGStatus, CalfloraInvasive) %>%
  summarise(tot = sum(abundance))

filter(plant16T.1, is.na(BGStatus)) %>% select(sp) %>% unique() # remove forb and fab, quantify poa (unknown grasses)

plant16T.1 %>%
  filter(!(sp %in% c("Forb", "Fabaceae"))) %>%
  group_by(BGStatus, CalfloraInvasive) %>%
  summarise(tot = sum(abundance)) # 160 individuals POA

plant16C.3 %>%
  group_by(BGStatus, CalfloraInvasive) %>%
  summarise(tot = sum(abundance))

filter(plant16C.3, is.na(BGStatus)) %>% select(sp) %>% unique() # remove forb and fab

# remove non-grass and combine datasets
bgdat <- full_join(plant15.1, plant16T.1) %>%
  full_join(plant16C.3) %>%
  filter(!(sp %in% c("Forb", "Fabaceae")))

# table of species
bgsptab <- bgdat %>%
  mutate(invasive = ifelse(CalfloraInvasive == 1, "yes", "no")) %>%
  group_by(experiment, year, species, BGStatus, invasive) %>%
  summarise(abundance = round(sum(abundance))) %>%
  ungroup() %>%
  rename(origin = BGStatus) %>%
  select(experiment, year, species, abundance, origin, invasive)

filter(bgsptab, species == "unidentified grass")

# summarise by origin
bgsum15 <- plant15.1 %>%
  filter(!is.na(BGStatus) & !(sp %in% c("Forb", "Fabaceae"))) %>%
  group_by(year, experiment, plot, subplot, BGStatus) %>%
  summarise(abundance = sum(abundance)) %>%
  ungroup() %>%
  spread(key = BGStatus, value = abundance) %>%
  rename(native.density = native, nonnative.density = nonnative)

bgsum16T <- plant16T.1 %>%
  filter(!is.na(BGStatus) & !(sp %in% c("Forb", "Fabaceae"))) %>%
  group_by(year, experiment, plot, subplot, BGStatus) %>%
  summarise(abundance = sum(abundance)) %>%
  ungroup() %>%
  spread(key = BGStatus, value = abundance) %>%
  rename(native.density = native, nonnative.density = nonnative)

bgsum16C <- plant16C.2 %>%
  filter(!is.na(BGStatus) & !(sp %in% c("Forb", "Fabaceae"))) %>%
  group_by(year, experiment, plot, subplot, bg.species, competition.density, BGStatus) %>%
  summarise(abundance = sum(abundance)) %>%
  ungroup() %>%
  spread(key = BGStatus, value = abundance) %>%
  rename(native.density = native, nonnative.density = nonnative)

# merge with fungal data
bgsum <- bgdat %>%
  filter(!is.na(BGStatus) & !(sp %in% c("Forb", "Fabaceae"))) %>%
  group_by(year, experiment, plot, subplot, BGStatus) %>%
  summarise(abundance = sum(abundance)) %>%
  ungroup() %>%
  spread(key = BGStatus, value = abundance) %>%
  rename(native.density = native, nonnative.density = nonnative)

dat2 <- dat %>%
  left_join(bgsum)
nrow(dat)
nrow(dat2)

# make sure all fungal samples have background info
dat2 %>%
  filter((year == 2015 & experiment == "transect") | (year == 2016 & experiment == "transect") | (year == 2016 & experiment == "competition")) %>%
  filter(is.na(native.density) | is.na(nonnative.density)) %>%
  select(year, experiment, plot, subplot) %>%
  unique()
# transects 11 - 15 missing in 2015 - okay


#### outputs ####
write_csv(bgsptab, "./output/background_species_table.csv")
write_csv(bgsum15, "./data/background_plants_transect_2015.csv")
write_csv(bgsum16T, "./data/background_plants_transect_2016.csv")
write_csv(bgsum16C, "./data/background_plants_competition_2016.csv")
write_csv(plant15.1, "./data/background_plants_raw_transect_2015.csv")
write_csv(plant16T.1, "./data/background_plants_raw_transect_2016.csv")
write_csv(plant16C.3, "./data/background_plants_raw_competition_2016.csv")
