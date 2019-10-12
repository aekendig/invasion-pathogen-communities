## Goal: analyze pathogen community richness, evenness, and diversity


#### set up ####

# clear everything except final data
rm(list = ls())

# load libraries
library(tidyverse)
library(vegan)
library(SpadeR)
library(glmmTMB)

# import data
dat <- read.csv("./data/fungal_pathogens_2015_2017.csv")


#### functions ####

# function to format data for Chao estimator
chaodat_fun <- function(df, origin, min.iso, max.iso){
  
  dout <- df %>%
    mutate(comm = paste(host, plot, year, sep = "_")) %>%
    filter(grass.status == origin) %>%
    group_by(comm, otu.id) %>%
    summarise(abundance = length(isolate.id)) %>%
    group_by(comm) %>%
    mutate(isolates = sum(abundance)) %>%
    filter(isolates > min.iso & isolates < max.iso) %>%
    ungroup() %>%
    select(-isolates) %>%
    spread(comm, abundance) %>%
    select(-otu.id)%>%
    data.frame()
  
  dout[is.na(dout)] = 0
  dout <- as.matrix(dout)
  
  return(dout)
}


#### edit data ####

# add sampling intensity and richness
rdat <- dat %>%
  group_by(year, year.f, experiment, plot, host, grass.status, otu.id) %>%
  summarise(abundance = length(isolate.id)) %>%
  group_by(year, year.f, experiment, plot, host, grass.status) %>%
  summarise(isolates = sum(abundance), 
            richness = length(otu.id),
            diversity = diversity(abundance),
            evenness = diversity(abundance)/log(specnumber(abundance))) %>%
  ungroup()
# checked that length(otu.id) = specnumber(abundance)

# look at sampling for origin types
rdat %>%
  ggplot(aes(x = isolates, y = richness, color = grass.status)) +
  geom_point(size = 2, position = position_jitter(width = 0.1, height = 0.1), alpha = 0.5) 

# counts by origin
rdat %>%
  group_by(grass.status) %>%
  summarise(communities = length(grass.status), samples = sum(isolates)) 
# more non-native communities were sampled, but native communities had higher sampling intensity

# correct by sampling intensity
rdat %>%
  ggplot(aes(x = isolates, y = log(richness/isolates), color = grass.status)) +
  geom_point(size = 2, position = position_jitter(width = 0.1, height = 0.1), alpha = 0.5) 
# log(isolates) will cause undetermined values for isolates = 1

# counts by origin, singletons removed
rdat %>%
  filter(isolates > 1) %>%
  group_by(grass.status) %>%
  summarise(communities = length(grass.status), samples = sum(isolates))
# increased the difference in sample number, but decreased the difference in community number

# subset data
rdat1 <- rdat %>% filter(isolates > 1)
rdat2 <- rdat %>% filter(isolates > 3) # same dataset used in composition analysis figure


#### correlations ####

# diversity metrics
rdat1 %>% ggplot(aes(richness, diversity)) + geom_point() # highly correlated
rdat1 %>% ggplot(aes(richness, evenness)) + geom_point() # not correlated
rdat1 %>% ggplot(aes(log(richness/isolates), evenness)) + geom_point() # ?
rdat1 %>% ggplot(aes(evenness, diversity)) + geom_point() # ?

cor.test(~log(richness/isolates) + evenness, rdat1) # yes
cor.test(~evenness + diversity, rdat1) # no

# random effects
rdat1 %>%
  select(host, year, plot) %>%
  unique() %>%
  data.frame()

rdat1 %>%
  group_by(host, plot) %>%
  summarise(nyears = length(unique(year))) %>%
  filter(nyears > 1)


#### richness models ####

# variance vs. mean
var(rdat1$richness)
mean(rdat1$richness) # similar in magnitude

var(rdat2$richness)
mean(rdat2$richness) # similar in magnitude

# distributions
rdat1 %>% ggplot(aes(richness)) + geom_histogram()
rdat2 %>% ggplot(aes(richness)) + geom_histogram()
rdat1 %>% ggplot(aes(log(richness/isolates))) + geom_histogram()
rdat2 %>% ggplot(aes(log(richness/isolates))) + geom_histogram()

# note: models did not converge using glmer or brm functions

# corrected richness, communities greater than 1
rmod1 <- glmmTMB(richness ~ offset(log(isolates)) + grass.status + (1|host) + (1|year.f) + (1|plot), data = rdat1, family = poisson)
summary(rmod1)
# no sig diff by origin, variance explained by plot, host, year

# raw richness, communities greater than 1
rmod1b <- glmmTMB(richness ~ grass.status + (1|host) + (1|year.f) + (1|plot), data = rdat1, family = poisson)
summary(rmod1b)
# non-native lower, variance explained by plot, host, year

# corrected richness, communities greater than 3
rmod2 <- glmmTMB(richness ~ offset(log(isolates)) + grass.status + (1|host) + (1|year.f) + (1|plot), data = rdat2, family = poisson)
# model does not converge

# raw richness, communities greater than 3
rmod2b <- glmmTMB(richness ~ grass.status + (1|host) + (1|year.f) + (1|plot), data = rdat2, family = poisson)
# model does not converge


#### evenness models ####

# check values
filter(rdat1, is.na(evenness)) %>% select(richness,isolates) %>% unique() # richness = 1
filter(rdat1, evenness == 1) %>% select(richness,isolates) %>% unique()

# maximum value other than one (more difficult to model y = 1, learn gamlss or zoib packages to use beta regressions with y = 1 and random effects)
filter(rdat1, evenness < 0.9999999) %>% select(evenness) %>% max(na.rm = T) # 0.9823888

# remove NA values and adjust values of 1
edat1 <- rdat1 %>% filter(!is.na(evenness)) %>% mutate(evenness2 = evenness - 1e-7)
edat2 <- rdat2 %>% filter(!is.na(evenness)) %>% mutate(evenness2 = evenness - 1e-7)

# distributions
edat1 %>% ggplot(aes(evenness)) + geom_histogram()
edat1 %>% ggplot(aes(evenness2)) + geom_histogram()
edat2 %>% ggplot(aes(evenness)) + geom_histogram()
edat2 %>% ggplot(aes(evenness2)) + geom_histogram()

# evenness, communities greater than 1
emod1 <- glmmTMB(evenness2 ~ grass.status + (1|host) + (1|year.f) + (1|plot), data = edat1, family = beta_family)
summary(emod1) # grass status not sig, host and plot, followed by year

# evenness, communities greater than 3
emod2 <- glmmTMB(evenness2 ~ grass.status + (1|host) + (1|year.f) + (1|plot), data = edat2, family = beta_family)
summary(emod2) # grass status barely sig (non-native higher), host and plot, followed by year


#### diversity models ####

# distributions
rdat1 %>% ggplot(aes(diversity)) + geom_histogram()
rdat2 %>% ggplot(aes(diversity)) + geom_histogram()

# diversity, communities greater than 1
dmod1 <- glmmTMB(diversity ~ grass.status + (1|host) + (1|year.f) + (1|plot), data = rdat1, family = gaussian)
summary(dmod1) # grass status not sig, plot, host, year

# diversity, communities greater than 3
dmod2 <- glmmTMB(diversity ~ grass.status + (1|host) + (1|year.f) + (1|plot), data = rdat2, family = gaussian)
summary(dmod2) # grass status sig (non-native lower), plot, followed by host and year

#### variance components - see bookmarked brms page ####
#### start here ####

#### Chao estimated richness ####
 
cdatn <- chaodat_fun(dat, "native", 0, 10) # richness = 58 (in 134 communities)
ChaoSpecies(cdatn, datatype = "abundance") # 89 (69, 125)

cdati <- chaodat_fun(dat, "non-native", 0, 10) # richness = 53 (in 197 communities)
ChaoSpecies(cdati, datatype = "abundance") # 131 (99, 187)

cdatn2 <- chaodat_fun(dat, "native", 0, 100) # richness = 65 (in 143 communities)
ChaoSpecies(cdatn2, datatype = "abundance") # 94 (79, 117)


#### diversity indices ####

Diversity(cdatn, datatype = "abundance") # 70 (53, 88)
Diversity(cdatn2, datatype = "abundance") # 73 (55, 92)
Diversity(cdati, datatype = "abundance") # 105 (76, 133)
