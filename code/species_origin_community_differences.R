## Goal: determine whether plant species within the same origin group have more similar pathogen communities than between groups

# notes: this is an updated version of origin_community_differences.R


#### set up ####

# clear everything except final data
rm(list = ls())

# load libraries
library(tidyverse)
library(vegan)
library(nationalparkcolors)
library(sjPlot)

# figure settings
pal <- c(park_palette("Everglades")[1:5], "black", "gray", "chartreuse")
axisText = 12
axisTitle = 14
legendText = 12
legendTitle = 12

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

# function to extract ellipses
# edited version of : https://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplot
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


#### output taxonomy table ####

# taxonomy data
tdat <- dat1 %>%
  group_by(year, year.f, otu.id, taxonomy) %>%
  summarise(abundance = length(isolate.id),
            host.species = paste(unique(host), collapse = ", "),
            host.origin = paste(unique(grass.status), collapse = ", ")) %>%
  ungroup() %>%
  group_by(year, year.f) %>%
  mutate(rank = rank(-abundance, ties.method = "random")) %>%
  ungroup()

# look at all fungi
tdat %>%
  ggplot(aes(x = rank, y = abundance)) +
  geom_point(size = 2) +
  facet_wrap(~year)

# looks like a natural break at different abundance levels
breakdat <- tibble(year = c(2015, 2016, 2017), rank.break = c(5.5, 5.5, 4.5))

# add in breaks
tdat %>%
  ggplot(aes(x = rank, y = abundance)) +
  geom_point(size = 2) +
  geom_vline(data = breakdat, aes(xintercept = rank.break), linetype = "dashed") +
  facet_wrap(~year)

# extract top fungal species
fdat <- tdat %>%
  full_join(breakdat) %>%
  filter(rank < rank.break) %>%
  select(-rank.break)

# number of fungal species
length(unique(fdat$otu.id)) # 7 (one only occurs on native species)

# save data
write_csv(fdat, "./output/taxonomy_species_origin_data.csv")


#### PERMANOVA ####
pmod1 <- adonis(cdat1 ~ grass.status + host + year.f, data = edat1, method="chao")
pmod1$aov.tab
# all three are sig, status is the least important

tab_df(round(pmod1$aov.tab, 3))


#### NMDS ####

# models
nmds1 <- metaMDS(cdat1, distance = "chao", halfchange = F, expand = F, trymax = 100) # error: results may be meaningless with non-integer data in method “chao” (every value is an integer)
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

# extract ellipses
# edited version of: https://stackoverflow.com/questions/56646820/functional-programming-use-broom-nest-tidy-unnest-and-map-within-a-function
ell1 <- ndat1 %>%
  group_by(grass.status) %>%
  nest() %>%
  mutate(out = map(data, veganCovEllipse),
         tidied = map(out, as_tibble)) %>%
  unnest(tidied)

# extract means for labels
lab1 <- ndat1 %>%
  group_by(grass.status) %>%
  summarise(NMDS1 = mean(NMDS1),
            NMDS2 = mean(NMDS2)) %>%
  mutate(grass.status = recode(grass.status, "non-native" = "invasive"))


#### visualize ####

# Chao nmds
pdf("./output/species_origin_community_differences.pdf", width = 5.5, height = 4)
ndat1 %>%
  mutate(host = recode(host, AB = "A. barbata", AF = "A. fatua", BD = "B. diandrus", BH = "B. hordeaceus", EG = "E. glaucus", FP = "F. perennis", PA = "P. aquatica", SP = "S. pulchra"),
         grass.status = recode(grass.status, "non-native" = "invasive")) %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(fill = host, shape = grass.status, group = year.f), size = 3) +
  geom_path(data = ell1, aes(linetype = grass.status), size=1) +  
  geom_text(data = lab1, aes(label = grass.status)) +
  scale_fill_manual(values = pal, name = "Host species") +
  scale_shape_manual(values = c(22, 25), name = "Origin") +
  guides(fill=guide_legend(override.aes=list(shape=c(rep(25, 4), 22, 25, 25, 22))), linetype = "none") +
  theme_bw() +
  theme(axis.text = element_text(size = axisText, color="black"),
        axis.title = element_text(size = axisTitle),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size=legendText),
        legend.title = element_text(size=legendTitle))
dev.off()

# bray-curtis nmds
ndat1b %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = host, shape = grass.status, size = year.f))

# fungal taxonomy
pdf("./output/species_origin_taxonomy_rank.pdf", width = 6, height = 3.5)
tdat %>%
  ggplot(aes(x = rank, y = abundance)) +
  geom_point(size = 1.5) +
  geom_vline(data = breakdat, aes(xintercept = rank.break), linetype = "dashed", size = 0.2) +
  geom_text(data = fdat, aes(label = taxonomy), hjust = -0.05, vjust = -0.05, size = 1.5) +
  facet_wrap(~year, scales = "free")  +  
  xlab("Rank") +
  ylab("Abundance") +
  theme_bw() +
  theme(axis.text = element_text(size = 10, color="black"),
        axis.title = element_text(size = 12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 10, color="black"))
dev.off()
