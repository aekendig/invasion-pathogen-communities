## Goal: select data to be used in analyses


#### set up ####

# clear everything except final data
rm(list = ls())

# load libraries
library(tidyverse)

# import data
dat <- read.csv("../Data/Data File 3 - Full Dataset.csv") # see metadata file for full information


#### edit data ####

# look at data
colnames(dat)
str(dat)

# see which experimental treatments should be used
unique(dat$host) # will eventually select
filter(dat, is.na(host)) # one is an Avena smut sample, the others are from sentinel sites
filter(dat, host == "?") # species not recorded
unique(dat$experiment) # all
unique(dat$year)
unique(dat$environmental.treatment) # ambient or NA, removes water, fungicide, and inoculation
unique(dat$soil.type) # non-serpentine
unique(dat$sentinel.site) 
unique(subset(dat,substr(sentinel.site,1,4)=="SERP")$soil.type) # serpentine will be removed
unique(dat$site.type) # PERENNIAL, ANNUAL, NA (not SERPENTINE)
unique(filter(dat, site.type == "SERPENTINE")$soil.type) # all serpentine
unique(dat$transect.density) # all, there are two types of C
unique(dat$gce.treatment) # ID numbers, info should be in environmental treatment
unique(dat$plant.tissue) # foliar or NA
unique(dat$bg.species) # all
unique(dat$competition.type) # all, lambda plots not included in dataset
unique(dat$jef.site) # all
unique(dat$fungicide.treated) # FALSE or NA
unique(subset(dat,fungicide.treated=="TRUE")$environmental.treatment) # fungicide will be removed

# fix transect density type-o
# add in an origin column
# scale down year
# remove data without pathogen isolate ID or from experimental treatments
# remove SP grown in serpentine
# remove BRAC and BRMA - only sampled once
dat2 <- dat %>%
  mutate(transect.density = recode(transect.density, "C " = "C"),
         grass.status = ifelse(host %in% c("SP","EG"), "native", "non-native"),
         native = ifelse(grass.status == "native", 1, 0),
         year.s = year - 2014,
         year.f = paste("year", year.s, sep = " ")) %>%
  filter(!is.na(isolate.id) & 
           (environmental.treatment == "ambient" | is.na(environmental.treatment)) & 
           soil.type == "non-serpentine" & 
           (plant.tissue == "foliar" | is.na(plant.tissue)) &
           !is.na(host) & !(host %in% c("?", "SSP", "BRAC", "BRMA"))) %>%
  as_tibble()


#### save data ####

write_csv(dat2, "./data/fungal_pathogens_2015_2017.csv")