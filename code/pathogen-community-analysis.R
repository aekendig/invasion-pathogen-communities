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

# first full model
mc1 <- brm(data = dc, family = gaussian,
           chao ~ native * year.s + (1|experiment) + (1|host),
           prior <- c(prior(normal(0, 100), class = Intercept),
                      prior(normal(0, 10), class = b),
                      prior(cauchy(0, 1), class = sd)),
           iter = 6000, warmup = 1000, chains = 1, cores = 1,
           control = list(adapt_delta = 0.99))
summary(mc1)
pp_check(mc1, nsamples = 100) # values go below zero making the model distributions wider than the observed

# second full model - log-transform richness
dc$log_chao = log(dc$chao)
mc1b <- update(mc1, formula. = log_chao ~ native * year.s + (1|experiment) + (1|host), newdata = dc)
pp_check(mc1b, nsamples = 100) # closer than before
summary(mc1b)

# run full model with more chains
mc1c <- update(mc1b, chains = 3, cores = 2,
               control = list(adapt_delta = 0.9999, max_treedepth = 15))
pp_check(mc1c, nsamples = 100)
summary(mc1c)
plot(mc1c) # convergence among chains
plot(marginal_effects(mc1c), points = T) # richness is higher for native plants, increasing over time

# subset models
mc2 <- update(mc1c, .~. -native:year.s)
mc3 <- update(mc2, .~. -native)
mc4 <- update(mc2, .~. -year.s)
mc5 <- update(mc4, .~. -native)

# model comparison
lc = list(loo(mc1c, reloo = T), loo(mc2, reloo = T), loo(mc3, reloo = T), loo(mc4, reloo = T), loo(mc5, reloo = T))
lc # check that k values are low and that p_loo values are near number of parameters
loo_model_weights(lc) # full model receives most of the weight followed by intercept-only model

# save models
save(mc1c, file = "./output/pathogen-community-analysis-chao-full.rda")
save(mc2, file = "./output/pathogen-community-analysis-chao-main.rda")
save(mc3, file = "./output/pathogen-community-analysis-chao-year.rda")
save(mc4, file = "./output/pathogen-community-analysis-chao-origin.rda")
save(mc5, file = "./output/pathogen-community-analysis-chao-intercept.rda")


#### observed richness ####

# dataset  - use same as for iChao

# data-only figure
dc %>%
  mutate(richness = otus/isolates) %>%
  ggplot(aes(x = grass.status, y = richness, shape = host)) +
  geom_point() + 
  facet_wrap(~year)

# check mean and variance
mean(dc$otus)
var(dc$otus) # almost double

# first full model
mr1 <- brm(data = dc, family = negbinomial(),
           otus ~ offset(isolates) + native * year.s + (1|experiment) + (1|host),
           prior <- c(prior(normal(0, 100), class = Intercept),
                      prior(normal(0, 10), class = b),
                      prior(cauchy(0, 1), class = sd)),
           iter = 6000, warmup = 1000, chains = 1, cores = 1,
           control = list(adapt_delta = 0.999))
summary(mr1)
#### start here - predicted values very large, making this show up weird ####
mr1_yrep <- posterior_predict(mr1)
pp_check(mr1, nsamples = 100)
# ppc_loo_pit_overlay

# run full model with more chains
mr1b <- update(mr1, chains = 3, cores = 2)
pp_check(mr1b)
summary(mr1b)
plot(mc1c) # convergence among chains
plot(marginal_effects(mc1c), points = T) # richness is higher for native plants, increasing over time

# subset models
mc2 <- update(mc1c, .~. -native:year.s)
mc3 <- update(mc2, .~. -native)
mc4 <- update(mc2, .~. -year.s)
mc5 <- update(mc4, .~. -native)

# model comparison
lc = list(loo(mc1c, reloo = T), loo(mc2, reloo = T), loo(mc3, reloo = T), loo(mc4, reloo = T), loo(mc5, reloo = T))
lc # check that k values are low and that p_loo values are near number of parameters
loo_model_weights(lc) # full model receives most of the weight followed by intercept-only model
