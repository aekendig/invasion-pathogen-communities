##### info ####

# file: pathogen-community-analysis.R
# author: Amy Kendig
# date last edited: 4/26/19
# goal: characterize pathogen communities on native and invasive grasses


#### set up ####

# clear everything except final data
rm(list = ls())

# load libraries
library(tidyverse)

# import data
d <- read.csv("./data/fungal-pathogens-2015-2017.csv")


#### edit data ####

# see which experimental treatments should be avoided
unique(d$environmental.treatment) # ambient or NA
unique(d$soil.type) # non-serpentine
unique(d$sentinel.site) 
unique(subset(d,substr(sentinel.site,1,4)=="SERP")$soil.type) # serpentine will be removed
unique(d$site.type)
unique(d$transect.density) # A-E with two types of C
unique(d$gce.treatment) # ID numbers, info should be in environmental treatment
unique(d$plant.tissue) # foliar or NA
unique(d$bg.species)
unique(d$competition.type) # lambda plots not included in dataset
unique(d$jef.site)
unique(d$fungicide.treated) # FALSE
unique(subset(d,fungicide.treated=="TRUE")$environmental.treatment) # fungicide will be removed

# fix transect density type-o
# add in an origin column
# remove data without pathogen isolate ID or from experimental treatments
d2 <- d %>%
  mutate(transect.density = recode(transect.density, "C " = "C"),
         grass.status = ifelse(host %in% c("SP","EG"), "native", "non-native")) %>%
  filter(!is.na(isolate.id) & 
           (environmental.treatment == "ambient" | is.na(environmental.treatment)) & 
           soil.type == "non-serpentine" & 
           (plant.tissue == "foliar" | is.na(plant.tissue)))
