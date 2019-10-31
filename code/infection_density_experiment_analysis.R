## Goal: look at infection prevalence across density gradient

# notes: updated of infection_experiment_analysis.R


#### set up ####

# clear everything except final data
rm(list = ls())

# load libraries
library(tidyverse)
library(glmmTMB)
library(DHARMa) # plot glmmTMB
library(MuMIn) # dredge

# import data
dat <- read_csv("./data/fungal_pathogens_2015_2017.csv")
bg15T <- read_csv("./data/background_plants_transect_2015.csv")
bg16T <- read_csv("./data/background_plants_transect_2016.csv")
bgC <- read_csv("./data/background_plants_competition_2016.csv")
fun <- read_csv("./output/taxonomy_species_origin_data.csv")


#### edit data ####

# list of fungi to analyze
select(fun, taxonomy, otu.id) %>% unique()

# combine transect background
bgT <- full_join(bg15T, bg16T)

# transect data
datT <- dat %>%
  inner_join(bgT) %>%
  mutate(nonnative = ifelse(host %in% c("EG", "SP"), 0, 1),
         origin = ifelse(nonnative == 1, "non-native", "native"),
         total.density = native.density + nonnative.density,
         natdens.s = scale(native.density)[,1],
         nondens.s = scale(nonnative.density)[,1],
         nonnative.rel = nonnative.density / total.density,
         ainf = ifelse(otu.id == 1, 1, 0),
         rpro = ifelse(otu.id == 3, 1, 0),
         pcha = ifelse(otu.id == 4, 1, 0),
         plol = ifelse(otu.id == 5, 1, 0),
         ptri = ifelse(otu.id == 8, 1, 0),
         dres = ifelse(otu.id == 2, 1, 0),
         pave = ifelse(otu.id == 7, 1, 0))

# check for ones that weren't added
subplotT <- filter(bgT, !(subplot %in% datT$subplot)) %>% select(subplot) %>% unique()
filter(dat, subplot %in% subplotT) 
filter(dat, subplot == "transect_7A") # checked the ones that weren't B and D

# competition data
datC <- dat %>%
  inner_join(bgC) %>%
  mutate(nonnative = ifelse(host %in% c("EG", "SP"), 0, 1),
         origin = ifelse(nonnative == 1, "non-native", "native"),
         total.density = native.density + nonnative.density,
         natdens.s = scale(native.density)[,1],
         nondens.s = scale(nonnative.density)[,1],
         nonnative.rel = nonnative.density / total.density,
         ainf = ifelse(otu.id == 1, 1, 0),
         rpro = ifelse(otu.id == 3, 1, 0),
         pcha = ifelse(otu.id == 4, 1, 0),
         plol = ifelse(otu.id == 5, 1, 0),
         ptri = ifelse(otu.id == 8, 1, 0),
         dres = ifelse(otu.id == 2, 1, 0),
         pave = ifelse(otu.id == 7, 1, 0))

# check ones that weren't added
subplotC <- filter(bgC, !(subplot %in% datC$subplot)) %>% select(subplot) %>% unique()
filter(dat, subplot %in% subplotC) 


#### ainf ####

# check for presence
sum(datT$ainf) # 33
sum(datC$ainf) # 30

# transect, absolute abundance
aiamodT <- glmmTMB(ainf ~ nonnative * (natdens.s + nondens.s) + (1|subplot), data = datT, family = binomial)
summary(aiamodT) # non-native origin increases
plot(simulateResiduals(aiamodT))
aiamodTa <- model.avg(get.models(dredge(aiamodT), subset = cumsum(weight) <= .95))
summary(aiamodTa) # origin, 2/3 by nat and nonnat dens

# transect, relative abundance
airmodT <- glmmTMB(ainf ~ nonnative * nonnative.rel + (1|subplot), data = datT, family = binomial)
summary(airmodT) # non-native origin increases
plot(simulateResiduals(airmodT))
airmodTa <- model.avg(get.models(dredge(airmodT), subset = cumsum(weight) <= .95))
summary(airmodTa) # 1/2 by rel. ab

# competition, absolute abundance
aiamodC <- glmmTMB(ainf ~ nonnative * (natdens.s + nondens.s) + (1|subplot), data = datC, family = binomial)
summary(aiamodC) # none
plot(simulateResiduals(aiamodC))
aiamodCa <- model.avg(get.models(dredge(aiamodC), subset = cumsum(weight) <= .95))
summary(aiamodCa) # main effects ~ 1/3

# competition, relative abundance
airmodC <- glmmTMB(ainf ~ nonnative * nonnative.rel + (1|subplot), data = datC, family = binomial)
summary(airmodC) # none
plot(simulateResiduals(airmodC))
airmodCa <- model.avg(get.models(dredge(airmodC), subset = cumsum(weight) <= .95))
summary(airmodCa) # main ~ 1/4


#### rpro ####

# check for presence
sum(datT$rpro) # 18
sum(datC$rpro) # 1

# transect, absolute abundance
rpamodT <- glmmTMB(rpro ~ nonnative * (natdens.s + nondens.s) + (1|subplot), data = datT, family = binomial)
summary(rpamodT) # none
plot(simulateResiduals(rpamodT))
rpamodTa <- model.avg(get.models(dredge(rpamodT), subset = cumsum(weight) <= .95))
summary(rpamodTa) # native density is the strongest predictor, followed closely by non-native

# transect, relative abundance
rprmodT <- glmmTMB(rpro ~ nonnative * nonnative.rel + (1|subplot), data = datT, family = binomial)
summary(rprmodT) # none
plot(simulateResiduals(rprmodT))
rprmodTa <- model.avg(get.models(dredge(rprmodT), subset = cumsum(weight) <= .95))
summary(rprmodTa) # relative abundance and origin tied


#### pcha ####

# check for presence
sum(datT$pcha) # 11
sum(datC$pcha) # 8

# transect, absolute abundance
pcamodT <- glmmTMB(pcha ~ nonnative * (natdens.s + nondens.s) + (1|subplot), data = datT, family = binomial)
summary(pcamodT) # output is unreliable - all z-values are 0

# transect, relative abundance
pcrmodT <- glmmTMB(pcha ~ nonnative * nonnative.rel + (1|subplot), data = datT, family = binomial)
summary(pcrmodT) #  output is unreliable - all z-values are 0

# competition, absolute abundance
pcamodC <- glmmTMB(pcha ~ nonnative * (natdens.s + nondens.s) + (1|subplot), data = datC, family = binomial)
summary(pcamodC) #  output is unreliable - all z-values are 0

# competition, relative abundance
pcrmodC <- glmmTMB(pcha ~ nonnative * nonnative.rel + (1|subplot), data = datC, family = binomial) #  model convergence error


#### plol ####

# check for presence
sum(datT$plol) # 24
sum(datC$plol) # 7

#### start here

# 2015 transect
plmod15 <- glmmTMB(plol ~ nonnative * (natdens.s + nondens.s) + (1|subplot), data = dat15, family = binomial)
summary(plmod15) # origin:native and main effects
plot(simulateResiduals(plmod15))
plmod15a <- model.avg(get.models(dredge(plmod15), subset = cumsum(weight) <= .95))
summary(plmod15a) # strong effect of interaction

# 2016 transect
plmod16T <- glmmTMB(plol ~ nonnative * (natdens.s + nondens.s) + (1|subplot), data = dat16T, family = binomial) # convergence error

# 2016 competition
plmod16C <- glmmTMB(plol ~ nonnative * (natdens.s + nondens.s) + (1|subplot), data = dat16C, family = binomial)
summary(plmod16C) # none
plot(simulateResiduals(plmod16C))
plmod16Ca <- model.avg(get.models(dredge(plmod16C), subset = cumsum(weight) <= .95)) # convergence errors


#### dres by density ####

# check for presence
sum(dat15$dres) # 6
sum(dat16T$dres) # 2
sum(dat16C$dres) # 66

# 2015 transect
drmod15 <- glmmTMB(dres ~ nonnative * (natdens.s + nondens.s) + (1|subplot), data = dat15, family = binomial)
summary(drmod15) # non dens
plot(simulateResiduals(drmod15))
drmod15a <- model.avg(get.models(dredge(drmod15), subset = cumsum(weight) <= .95))
summary(drmod15a) # non-native density and origin

# 2016 transect
drmod16T <- glmmTMB(dres ~ nonnative * (natdens.s + nondens.s) + (1|subplot), data = dat16T, family = binomial) # convergence error

# 2016 competition
drmod16C <- glmmTMB(dres ~ nonnative * (natdens.s + nondens.s) + (1|subplot), data = dat16C, family = binomial)
summary(drmod16C) # none
plot(simulateResiduals(drmod16C))
drmod16Ca <- model.avg(get.models(dredge(drmod16C), subset = cumsum(weight) <= .95))
summary(drmod16Ca) # density of each and origin all about the same


#### outputs ####

write_csv(dat15, "./output/infect_density_transect15_data.csv")
write_csv(dat16T, "./output/infect_density_transect16_data.csv")
write_csv(dat16C, "./output/infect_density_competition_data.csv")

save(rpmod15a, file = "./output/infection_density_rpro_transect15_avg_model.rda")
save(rpmod16Ta, file = "./output/infection_density_rpro_transect16_avg_model.rda")
save(plmod15a, file = "./output/infection_density_plol_transect15_avg_model.rda")
save(pcmod15a, file = "./output/infection_density_pcha_transect15_avg_model.rda")
save(aimod15a, file = "./output/infection_density_ainf_transect15_avg_model.rda")
save(aimod16Ta, file = "./output/infection_density_ainf_transect16_avg_model.rda")
save(aimod16Ca, file = "./output/infection_density_ainf_competition_avg_model.rda")
save(drmod15a, file = "./output/infection_density_dres_transect15_avg_model.rda")
save(drmod16Ca, file = "./output/infection_density_dres_competition_avg_model.rda")