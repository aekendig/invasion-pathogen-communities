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
dat2 <- dat %>%
  mutate(transect.density = recode(transect.density, "C " = "C"),
         grass.status = ifelse(host %in% c("SP","EG"), "native", "non-native"),
         native = ifelse(grass.status == "native", 1, 0),
         year.s = year - 2015) %>%
  filter(!is.na(isolate.id) & 
           (environmental.treatment == "ambient" | is.na(environmental.treatment)) & 
           soil.type == "non-serpentine" & 
           (plant.tissue == "foliar" | is.na(plant.tissue)) &
           !is.na(host) & host != "?") %>%
  as_tibble()


#### sample sizes ####

# check that all isolates are unique
sum(is.na(dat2$isolate.id))
length(dat2$isolate.id) - length(unique(dat2$isolate.id))

# plots within experiments
dat2 %>%
  group_by(experiment) %>%
  summarise(plots = length(unique(plot)),
            first = unique(plot)[1]) # all have plots

# look at what's in plots
dat2 %>%
  group_by(year, experiment, plot) %>%
  summarise(origins = length(unique(grass.status)), samples = length(isolate.id)) %>%
  data.frame() # JEF and competition only have one origin per plot

# look at species within experiments
dat2 %>%
  group_by(year, experiment, host) %>%
  summarise(samples = length(isolate.id)) %>%
  data.frame() # JEF only has SP, remove to help make native and non-native numbers more similar

# group by community type
dat2 %>%
  filter(experiment != "JEF transect") %>%
  group_by(year, experiment, plot, host) %>%
  summarise(isolates = length(isolate.id)) %>%
  ggplot(aes(x = isolates, fill = host)) + 
  geom_histogram(binwidth = 1) + 
  geom_vline(xintercept = 4, color = "red", linetype = "dashed") # a lot are less than 4, which is minimum for Chao

# counts by origin
dat2 %>%
  filter(experiment != "JEF transect") %>%
  group_by(year, experiment, plot, host, grass.status) %>%
  summarise(isolates = length(isolate.id)) %>%
  filter(isolates > 3) %>%
  group_by(grass.status) %>%
  summarise(communities = length(grass.status), samples = sum(isolates))

# group by community type, origin instead of host
dat2 %>%
  filter(experiment != "JEF transect") %>%
  group_by(year, experiment, plot, grass.status) %>%
  summarise(isolates = length(isolate.id)) %>%
  ggplot(aes(x = isolates, fill = grass.status)) + 
  geom_histogram(binwidth = 1) + 
  geom_vline(xintercept = 4, color = "red", linetype = "dashed") 

# counts by origin, origin instead of host
dat2 %>%
  filter(experiment != "JEF transect") %>%
  group_by(year, experiment, plot, grass.status) %>%
  summarise(isolates = length(isolate.id)) %>%
  filter(isolates > 3) %>%
  group_by(grass.status) %>%
  summarise(communities = length(grass.status), samples = sum(isolates)) # many more samples are retained and sample numbers are similar, similar when rare communities are included too

# see if one host type is disproportionately represented in one experiment
dat2 %>%
  filter(experiment != "JEF transect") %>%
  group_by(year, experiment, grass.status) %>%
  summarise(isolates = length(isolate.id)) %>%
  data.frame() # yes, but not extreme

# rarefaction curves with host species
dat2 %>%
  filter(experiment != "JEF transect") %>%
  group_by(year, experiment, plot, host, grass.status) %>%
  summarise(isolates = length(isolate.id), richness = length(unique(otu.id))) %>%
  ggplot(aes(x = isolates, y = richness, color = host, shape = grass.status)) +
  geom_point(size = 2) +
  geom_line()

# rarefaction curves with host species
dat2 %>%
  filter(experiment != "JEF transect") %>%
  group_by(year, experiment, plot, grass.status) %>%
  summarise(isolates = length(isolate.id), richness = length(unique(otu.id))) %>%
  ggplot(aes(x = isolates, y = richness, color = grass.status)) +
  geom_point(size = 2) +
  geom_line() # more similar than host species and appear to asymptote


#### edit data ####

# remove JEF experiment
dat3 <- dat2 %>%
  filter(experiment != "JEF transect")


#### save data ####

write_csv(dat3, "./data/fungal_pathogens_2015_2017.csv")
