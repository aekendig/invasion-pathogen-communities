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
library(SpadeR)
library(brms)
library(tidybayes)

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
# scale down year
# remove data without pathogen isolate ID or from experimental treatments
d2 <- d %>%
  mutate(transect.density = recode(transect.density, "C " = "C"),
         grass.status = ifelse(host %in% c("SP","EG"), "native", "non-native"),
         native = ifelse(grass.status == "native", 1, 0),
         year.s = year - 2015) %>%
  filter(!is.na(isolate.id) & 
           (environmental.treatment == "ambient" | is.na(environmental.treatment)) & 
           soil.type == "non-serpentine" & 
           (plant.tissue == "foliar" | is.na(plant.tissue)))


#### estimated richness (iChao) ####

# hosts of interest
hostList = c("AB","BD","BH","EG","SP")

# dataset (at least four samples)
dc <- d2 %>%
  filter(host %in% hostList) %>%
  group_by(year, year.s, experiment, host, grass.status, native, otu.id) %>%
  summarise(abundance = length(isolate.id)) %>%
  group_by(year, year.s, experiment, host, grass.status, native) %>%
  mutate(isolates = sum(abundance)) %>% 
  filter(isolates > 3) %>% 
  summarise(otus = length(otu.id), 
            isolates = sum(abundance), 
            chao = ChaoSpecies(abundance, datatype="abundance")$Species_table[5,"Estimate"])

# data-only figure
dc %>%
  ggplot(aes(x = grass.status, y = chao, shape = host)) +
  geom_point() + 
  facet_wrap(~year)

# models
mc1 <- brm(data = dc, family = gaussian,
           chao ~ native * year.s + (1|experiment) + (1|host),
           prior <- c(prior(normal(0, 100), class = Intercept),
                      prior(normal(0, 10), class = b),
                      prior(cauchy(0, 1), class = sd)),
           iter = 5000, warmup = 1000, chains = 3, cores = 2,
           control = list(adapt_delta = 0.9))
summary(mc1)
marginal_effects(mc1)
plot(mc1)

# subset models
mc2 <- update(mc1, .~. -native:year.s)
summary(mc2)
mc3 <- update(mc2, .~. -native)
plot(mc3)
mc4 <- update(mc2, .~. -year.s)
plot(mc4)
mc5 <- update(mc4, .~. -native)

# model comparison
lc = list(loo(mc1, reloo = T), loo(mc2, reloo = T), loo(mc3, reloo = T), loo(mc4, reloo = T), loo(mc5, reloo = T))
lc
loo_model_weights(lc)
