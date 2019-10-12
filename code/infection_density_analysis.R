## Goal: look at infection prevalence across density gradient


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
bg15 <- read_csv("./data/background_plants_transect_2015.csv")
bg16T <- read_csv("./data/background_plants_transect_2016.csv")
bg16C <- read_csv("./data/background_plants_competition_2016.csv")


#### edit data ####
dat15 <- dat %>%
  inner_join(bg15) %>%
  mutate(nonnative = ifelse(host %in% c("EG", "SP"), 0, 1),
         origin = ifelse(nonnative == 1, "non-native", "native"),
         total.density = native.density + nonnative.density,
         natdens.s = scale(native.density)[,1],
         nondens.s = scale(nonnative.density)[,1],
         rpro = ifelse(taxonomy == "Ramularia cf. proteae", 1, 0),
         plol = ifelse(taxonomy == "Pyrenophora cf. lolii", 1, 0),
         pcha = ifelse(taxonomy == "Pyrenophora cf. chaetomioides", 1, 0),
         ainf = ifelse(taxonomy == "Alternaria cf. infectoria", 1, 0),
         dres = ifelse(taxonomy == "Drechslera sp.", 1, 0))

dat16T <- dat %>%
  inner_join(bg16T) %>%
  mutate(nonnative = ifelse(host %in% c("EG", "SP"), 0, 1),
         origin = ifelse(nonnative == 1, "non-native", "native"),
         total.density = native.density + nonnative.density,
         natdens.s = scale(native.density)[,1],
         nondens.s = scale(nonnative.density)[,1],
         rpro = ifelse(taxonomy == "Ramularia cf. proteae", 1, 0),
         plol = ifelse(taxonomy == "Pyrenophora cf. lolii", 1, 0),
         pcha = ifelse(taxonomy == "Pyrenophora cf. chaetomioides", 1, 0),
         ainf = ifelse(taxonomy == "Alternaria cf. infectoria", 1, 0),
         dres = ifelse(taxonomy == "Drechslera sp.", 1, 0))

dat16C <- dat %>%
  inner_join(bg16C) %>%
  mutate(nonnative = ifelse(host %in% c("EG", "SP"), 0, 1),
         origin = ifelse(nonnative == 1, "non-native", "native"),
         total.density = native.density + nonnative.density,
         natdens.s = scale(native.density)[,1],
         nondens.s = scale(nonnative.density)[,1],
         rpro = ifelse(taxonomy == "Ramularia cf. proteae", 1, 0),
         plol = ifelse(taxonomy == "Pyrenophora cf. lolii", 1, 0),
         pcha = ifelse(taxonomy == "Pyrenophora cf. chaetomioides", 1, 0),
         ainf = ifelse(taxonomy == "Alternaria cf. infectoria", 1, 0),
         dres = ifelse(taxonomy == "Drechslera sp.", 1, 0))


#### rpro by density ####

# check for presence
sum(dat15$rpro) # 15
sum(dat16T$rpro) # 3
sum(dat16C$rpro) # 1

# 2015 transect
rpmod15 <- glmmTMB(rpro ~ nonnative * (natdens.s + nondens.s) + (1|subplot), data = dat15, family = binomial)
summary(rpmod15) # none
plot(simulateResiduals(rpmod15))
rpmod15a <- model.avg(get.models(dredge(rpmod15), subset = cumsum(weight) <= .95))
summary(rpmod15a) # native density is the strongest predictor

# 2016 transect
rpmod16T <- glmmTMB(rpro ~ nonnative * (natdens.s + nondens.s) + (1|subplot), data = dat16T, family = binomial)
summary(rpmod16T) # none
plot(simulateResiduals(rpmod16T))
rpmod16Ta <- model.avg(get.models(dredge(rpmod16T), subset = cumsum(weight) <= .95))
summary(rpmod16Ta) # main effects are even


#### plol by density ####

# check for presence
sum(dat15$plol) # 21
sum(dat16T$plol) # 4
sum(dat16C$plol) # 4

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


#### pcha by density ####

# check for presence
sum(dat15$pcha) # 10
sum(dat16T$pcha) # 1
sum(dat16C$pcha) # 9

# 2015 transect
pcmod15 <- glmmTMB(pcha ~ nonnative * (natdens.s + nondens.s) + (1|subplot), data = dat15, family = binomial)
summary(pcmod15) # none
plot(simulateResiduals(pcmod15))
pcmod15a <- model.avg(get.models(dredge(pcmod15), subset = cumsum(weight) <= .95))
summary(pcmod15a) # origin most important, each density somewhat important

# 2016 competition
pcmod16C <- glmmTMB(pcha ~ nonnative * (natdens.s + nondens.s) + (1|subplot), data = dat16C, family = binomial) # convergence error


#### ainf by density ####

# check for presence
sum(dat15$ainf) # 16
sum(dat16T$ainf) # 17
sum(dat16C$ainf) # 26

# 2015 transect
aimod15 <- glmmTMB(ainf ~ nonnative * (natdens.s + nondens.s) + (1|subplot), data = dat15, family = binomial)
summary(aimod15) # origin
plot(simulateResiduals(aimod15))
aimod15a <- model.avg(get.models(dredge(aimod15), subset = cumsum(weight) <= .95))
summary(aimod15a) # strong effect of origin (higher on non-native)

# 2016 transect
aimod16T <- glmmTMB(ainf ~ nonnative * (natdens.s + nondens.s) + (1|subplot), data = dat16T, family = binomial)
summary(aimod16T) # none
plot(simulateResiduals(aimod16T))
aimod16Ta <- model.avg(get.models(dredge(aimod16T), subset = cumsum(weight) <= .95))
summary(aimod16Ta) # origin most influential

# 2016 competition
aimod16C <- glmmTMB(ainf ~ nonnative * (natdens.s + nondens.s) + (1|subplot), data = dat16C, family = binomial)
summary(aimod16C) # none
plot(simulateResiduals(aimod16C))
aimod16Ca <- model.avg(get.models(dredge(aimod16C), subset = cumsum(weight) <= .95))
summary(aimod16Ca) # density of each and origin all about a third


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