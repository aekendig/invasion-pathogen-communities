## Goal: analyze pathogen community richness, evenness, and diversity


#### set up ####

# clear everything except final data
rm(list = ls())

# load libraries
library(tidyverse)
library(vegan)
library(lme4)
library(SpadeR)

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
            evenness = diversity(abundance)/log(specnumber(abundance)))

# look at sampling for origin types
rdat %>%
  ggplot(aes(x = isolates, y = richness, color = grass.status)) +
  geom_point(size = 2) +
  geom_line() # cut-off top 8

# counts by origin
rdat %>%
  filter(isolates > 3 & isolates < 10) %>%
  group_by(grass.status) %>%
  summarise(communities = length(grass.status), samples = sum(isolates)) # similar when minimum is 3, less similar when it's 0 (native lower)

# subset data
rdat1 <- rdat %>% filter(isolates > 0 & isolates < 10)
rdat2 <- rdat %>% filter(isolates > 3 & isolates < 10)
rdat3 <- rdat %>% filter(isolates > 3) # same dataset used in composition analysis


#### richness models ####

# look at variance vs. mean
var(rdat1$richness)
mean(rdat1$richness) # similar in magnitude

var(rdat2$richness)
mean(rdat2$richness) # mean is a little larger

var(rdat3$richness)
mean(rdat3$richness)  # similar in magnitude

# richness, small communities included
rmod1 <- glmer(richness ~ grass.status + (1|host) + (1|year.f) + (1|plot), data = rdat1, family = poisson)
summary(rmod1)
# no sig effect of status (slightly lower for invasive), plot explains the most variance

# richness, small communities excluded
rmod2 <- glmer(richness ~ grass.status + (1|host) + (1|year.f) + (1|plot), data = rdat2, family = poisson)
summary(rmod2)
# singular fit

# richness, small communities excluded, large communities included
rmod3 <- glmer(richness ~ grass.status + (1|host) + (1|year.f) + (1|plot), data = rdat3, family = poisson)
summary(rmod3)
# singular fit


#### evenness models ####

# look at values
rdat1 %>% ggplot(aes(evenness)) + geom_histogram()
rdat1 %>% ggplot(aes(log(evenness))) + geom_histogram()


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
