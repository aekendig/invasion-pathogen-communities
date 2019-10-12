## Goal: determine whether plant species within the same origin group have more similar pathogen communities than between groups


#### set up ####

# clear everything except final data
rm(list = ls())

# load libraries
library(tidyverse)
library(vegan)

# import data
dat <- read.csv("./data/fungal_pathogens_2015_2017.csv")

# function to create wide data
datw_fun <- function(df, min.iso) {
  dout <- df %>%
    group_by(year, year.f, host, grass.status, otu.id) %>%
    summarise(abundance = length(isolate.id)) %>%
    group_by(year, year.f, host, grass.status) %>%
    mutate(isolates = sum(abundance)) %>%
    filter(isolates > min.iso) %>%
    ungroup() %>%
    spread(otu.id, abundance)
  
  return(dout)
}

# function to create community matrix
cdat_fun <- function(datw){
  cdat <- datw %>%
    select(-c(year:isolates)) %>%
    data.frame()
  
  cdat[is.na(cdat)] = 0
  cdat <- as.matrix(cdat)
  
  return(cdat)
}


#### edit data ####

# remove JEF to make numbers between native and non-native more comparable
dat1 <- dat %>%
  filter(experiment != "JEF transect")

# wide dataset
(datw <- datw_fun(dat1, 0))
# create one that excludes BRAD (much lower sampling than others)
(datw1 <- datw_fun(dat1, 4))

# community matrix
cdat1 <- cdat_fun(datw1)

# environmental matrix
edat1 <- datw1 %>% select(c(year:grass.status))


#### PERMANOVA ####
pmod1 <- adonis(cdat1 ~ grass.status + host + year.f, data = edat1, method="chao")
pmod1$aov.tab
# all three are sig, status is the least important


#### NMDS ####

# models
nmds1 <- metaMDS(cdat1, distance = "chao", halfchange = F, expand = F, trymax = 100)
nmds1b <- metaMDS(cdat1, distance = "bray", halfchange = F, expand = F, trymax = 100)

# goodness of fit
gof1 <- goodness(nmds1)
gof1b <- goodness(nmds1b)

# extract data scores from nmds
ndat1 <- as.data.frame(scores(nmds1)) %>%
  cbind(edat1) %>%
  cbind(gof1)
head(ndat1) 

ndat1b <- as.data.frame(scores(nmds1b)) %>%
  cbind(edat1) %>%
  cbind(gof1b)
head(ndat1b) 

# visualize
ndat1 %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = host, shape = grass.status, size = year.f))

ndat1b %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = host, shape = grass.status, size = year.f))
