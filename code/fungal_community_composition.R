## Goal: analyze pathogen community composition


#### set up ####

# clear everything except final data
rm(list = ls())

# load libraries
library(tidyverse)
library(vegan)
library(bipartite)

# import data
dat <- read.csv("./data/fungal_pathogens_2015_2017.csv")


#### functions ####

# function to create wide data
datw_fun <- function(df, min.iso) {
  dout <- df %>%
    group_by(year, year.f, experiment, plot, host, grass.status, otu.id) %>%
    summarise(abundance = length(isolate.id)) %>%
    group_by(year, experiment, plot, host, grass.status) %>%
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

# function to extract ellipses
# edited version of : https://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo
veganCovEllipse <- function (df){
  scale <- 1
  npoints <- 100
  cov <- with(df, cov.wt(cbind(NMDS1, NMDS2), wt=rep(1/length(NMDS1), length(NMDS1)))$cov)
  center <- with(df, c(mean(NMDS1), mean(NMDS2)))
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}


#### edit data ####

# wide datasets
datw1 <- datw_fun(dat, 0)
datw2 <- datw_fun(dat, 1)
datw3 <- datw_fun(dat, 2)
datw4 <- datw_fun(dat, 3)

# grass status replication for datasets
datw1 %>% group_by(grass.status) %>% summarise(n = length(plot))
datw2 %>% group_by(grass.status) %>% summarise(n = length(plot))
datw3 %>% group_by(grass.status) %>% summarise(n = length(plot))
datw4 %>% group_by(grass.status) %>% summarise(n = length(plot))

# community matrices
cdat1 <- cdat_fun(datw1)
cdat2 <- cdat_fun(datw2)
cdat3 <- cdat_fun(datw3)
cdat4 <- cdat_fun(datw4)

# environmental matrices
edat1 <- datw1 %>% select(c(year:grass.status))
edat2 <- datw2 %>% select(c(year:grass.status))
edat3 <- datw3 %>% select(c(year:grass.status))
edat4 <- datw4 %>% select(c(year:grass.status))


#### PERMANOVA ####

pmod1 <- adonis(cdat1 ~ grass.status + host + year.f + plot, data = edat1, method="chao")
pmod1$aov.tab

pmod2 <- adonis(cdat2 ~ grass.status + host + year.f + plot, data = edat2, method="chao")
pmod2$aov.tab

pmod3 <- adonis(cdat3 ~ grass.status + host + year.f + plot, data = edat3, method="chao")
pmod3$aov.tab

pmod4 <- adonis(cdat4 ~ grass.status + host + year.f + plot, data = edat4, method="chao")
pmod4$aov.tab


#### NMDS ####
# nmds1 <- metaMDS(cdat1, distance = "chao", halfchange = F, expand = F, trymax = 100) # nmds can't converge with cdat using chao or bray

nmds2 <- metaMDS(cdat2, distance = "chao", halfchange = F, expand = F, trymax = 100)

nmds3 <- metaMDS(cdat3, distance = "bray", halfchange = F, expand = F, trymax = 300) # doesn't converge with chao, doesn't consistently converge

nmds4 <- metaMDS(cdat4, distance = "chao", halfchange = F, expand = F, trymax = 300) # doesn't consistently converge

# goodness of fit
gof2 <- goodness(nmds2)
gof3 <- goodness(nmds3)
gof4 <- goodness(nmds4)

# extract data scores from nmds
ndat2 <- as.data.frame(scores(nmds2)) %>%
  cbind(edat2) %>%
  cbind(gof2)
head(ndat2) 

ndat3 <- as.data.frame(scores(nmds3)) %>%
  cbind(edat3) %>%
  cbind(gof3)
head(ndat3) 

ndat4 <- as.data.frame(scores(nmds4)) %>%
  cbind(edat4) %>%
  cbind(gof4)
head(ndat4) 

# extract ellipses
# edited version of: https://stackoverflow.com/questions/56646820/functional-programming-use-broom-nest-tidy-unnest-and-map-within-a-function
ell4 <- ndat4 %>%
  group_by(grass.status) %>%
  nest() %>%
  mutate(out = map(data, veganCovEllipse),
         tidied = map(out, as_tibble)) %>%
  unnest(tidied, .drop = T)


#### visualize NMDS ####

# communities with at least 2
ndat2 %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = host)) # difficult to see most of the communities because of one very different one

filter(ndat2, NMDS1 > 0.8) # FP in a transect
nrow(filter(ndat2, host == "FP")) # 4 FP communities

ndat2 %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = grass.status)) +
  xlim(-0.01, 0.01)

# communities with at least 3
ndat3 %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = host)) # difficult to see most of the communities because of one very different one

filter(ndat3, NMDS1 > 0.8)  # PA in a transect
nrow(filter(ndat2, host == "PA")) # 5 PA communities

# communities with at least 4
ndat4 %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = plot, shape = grass.status), size = gof4*75) +
  geom_path(data = ell4, aes(linetype = grass.status), size=1) +  
  scale_color_viridis_d(guide = F)
  

#### visualize bipartite plot - didn't finish this ####

# create web
web1 <- dat %>%
  group_by(year, experiment, plot, host, grass.status, otu.id) %>%
  summarise(abundance = length(isolate.id)) %>%
  group_by(year, experiment, plot, host, grass.status) %>%
  mutate(isolates = sum(abundance),
         site = paste(plot, host, year, sep = "_")) %>%
  filter(isolates > 1)

# create plot
plotweb(web1,method = "normal", text.rot = 90, col.high = "white", col.low = "white", col.interaction = "black", labsize = 0.75) 

