## Goal: analyze pathogen community composition


#### set up ####

# clear everything except final data
rm(list = ls())

# load libraries
library(tidyverse)
library(vegan)

# import data
dat <- read.csv("./data/fungal_pathogens_2015_2017.csv")


#### edit data ####

# make data wide
dat2 <- dat %>%
  group_by(year, experiment, plot, host, grass.status, otu.id) %>%
  summarise(abundance = length(isolate.id)) %>%
  ungroup() %>%
  spread(otu.id, abundance)

# community matrix 
cdat <- dat2 %>%
  select(-c(year:grass.status)) %>%
  data.frame()

cdat[is.na(cdat)] = 0
cdat <- as.matrix(cdat)

# environmental matrix
edat <- dat2 %>%
  select(c(year:grass.status))


#### analyze data ####

# PERMANOVA
pmod <- adonis(cdat ~ grass.status + host + year + plot, data = edat, method="chao")
pmod$aov.tab

# NMDS - start here, nmds not converging, probably because of small communities
nmds <- metaMDS(cdat, distance = "chao", halfchange = F, expand = F, trymax = 100)
nmds
