## Goal: relationship between plant density and damage

# notes: update of damage_density_analysis.R with transect experiment combined across years


#### set up ####

# clear everything except final data
rm(list = ls())

# load libraries
library(tidyverse)
library(glmmTMB)
library(lme4)
library(DHARMa) # plot glmmTMB
library(MuMIn) # dredge

# import data
bg15 <- read_csv("./data/background_plants_transect_2015.csv")
bg16T <- read_csv("./data/background_plants_transect_2016.csv")
bg16C <- read_csv("./data/background_plants_competition_2016.csv")

dam15plant <- read_csv("./data/damage_plant_transect_2015.csv")
dam16Tplant <- read_csv("./data/damage_plant_transect_2016.csv")
dam16Cplant <- read_csv("./data/damage_plant_competition_2016.csv")

dam15leaf <- read_csv("./data/damage_leaf_transect_2015.csv")
dam16Tleaf <- read_csv("./data/damage_leaf_transect_2016.csv")
dam16Cleaf <- read_csv("./data/damage_leaf_competition_2016.csv")


#### edit data ####

# 2015 transect
bg15
dam15plant
dam15leaf

dat15plant <- dam15plant %>%
  filter(host != "PA") %>%
  left_join(bg15) %>%
  mutate(nonnative = ifelse(host %in% c("EG", "SP"), 0, 1),
         origin = ifelse(nonnative == 1, "non-native", "native"),
         total.density = native.density + nonnative.density + other.density)

dat15leaf <- dam15leaf %>%
  filter(host != "PA") %>%
  left_join(bg15) %>%
  mutate(nonnative = ifelse(host %in% c("EG", "SP"), 0, 1),
         origin = ifelse(nonnative == 1, "non-native", "native"),
         total.density = native.density + nonnative.density + other.density) %>%
  filter(!is.na(infected))

filter(dat15plant, is.na(native.density) | is.na(nonnative.density)) # all plots have info
filter(dat15leaf, is.na(native.density) | is.na(nonnative.density)) # all plots have info

# 2016 transect
bg16T
dam16Tplant
dam16Tleaf

dat16Tplant <- dam16Tplant %>%
  left_join(bg16T) %>%
  mutate(nonnative = ifelse(host %in% c("EG", "SP"), 0, 1),
         origin = ifelse(nonnative == 1, "non-native", "native"),
         total.density = native.density + nonnative.density + other.density)

dat16Tleaf <- dam16Tleaf %>%
  left_join(bg16T) %>%
  mutate(nonnative = ifelse(host %in% c("EG", "SP"), 0, 1),
         origin = ifelse(nonnative == 1, "non-native", "native"),
         total.density = native.density + nonnative.density + other.density) %>%
  filter(!is.na(infected))

filter(dat16Tplant, is.na(native.density) | is.na(nonnative.density)) # all plots have info
filter(dat16Tleaf, is.na(native.density) | is.na(nonnative.density)) # all plots have info

# combine 2015 and 2016 transect
datTplant <- full_join(dat15plant, dat16Tplant) %>%
  mutate(natdens.s = scale(native.density)[,1],
         nondens.s = scale(nonnative.density)[,1],
         othdens.s = scale(other.density)[,1],
         nonnative.rel = nonnative.density / total.density,
         leaves.healthy = leaves.tot - leaves.dam)

datTleaf <- full_join(dat15leaf, dat16Tleaf) %>%
  mutate(natdens.s = scale(native.density)[,1],
         nondens.s = scale(nonnative.density)[,1],
         othdens.s = scale(other.density)[,1],
         nonnative.rel = nonnative.density / total.density)

# 2016 competition
bg16C
dam16Cplant
dam16Cleaf

datCplant <- dam16Cplant %>%
  left_join(bg16C %>% mutate(bg.species = recode(bg.species, "EG adults" = "EG adult"))) %>%
  mutate(nonnative = ifelse(host %in% c("EG", "SP"), 0, 1),
         origin = ifelse(nonnative == 1, "non-native", "native"),
         total.density = native.density + nonnative.density,
         natdens.s = scale(native.density)[,1],
         nondens.s = scale(nonnative.density)[,1],
         nonnative.rel = nonnative.density / total.density,
         leaves.healthy = leaves.tot - leaves.dam)

datCleaf <- dam16Cleaf %>%
  left_join(bg16C %>% mutate(bg.species = recode(bg.species, "EG adults" = "EG adult"))) %>%
  mutate(nonnative = ifelse(host %in% c("EG", "SP"), 0, 1),
         origin = ifelse(nonnative == 1, "non-native", "native"),
         total.density = native.density + nonnative.density,
         natdens.s = scale(native.density)[,1],
         nondens.s = scale(nonnative.density)[,1],
         nonnative.rel = nonnative.density / total.density) %>%
  filter(!is.na(infected))

filter(datCplant, is.na(native.density) | is.na(nonnative.density)) # all plots have info
filter(datCleaf, is.na(native.density) | is.na(nonnative.density)) # all plots have info

# all plants
datplant <- full_join(datTplant, datCplant)

datleaf <- full_join(datTleaf, datCleaf)

# sample sizes
nrow(datplant)
unique(datplant$host)

# plot-scale dataset
plotdat <- datplant %>%
  dplyr::select(year, experiment, plot, subplot, bg.species, competition.density, native.density, nonnative.density, other.density, total.density, natdens.s, nondens.s, othdens.s, nonnative.rel) %>%
  unique() %>%
  mutate(year.f = factor(year))


#### visualize ####

# native density and damage
datplant %>%
  ggplot(aes(x = log(native.density+1), y = mean.dam, color = origin)) +
  geom_point(position = position_dodge(0.1), alpha = 0.5) +
  facet_wrap(~experiment, scales = "free")

datplant %>%
  ggplot(aes(x = natdens.s, y = mean.dam, color = origin)) +
  geom_point(position = position_dodge(0.1), alpha = 0.5) +
  facet_wrap(~experiment, scales = "free")

datplant %>%
  ggplot(aes(x = log(native.density+1), y = prop.dam, color = origin)) +
  geom_point(position = position_dodge(0.1), alpha = 0.5) +
  facet_wrap(~experiment, scales = "free")

datplant %>%
  ggplot(aes(x = natdens.s, y = prop.dam, color = origin)) +
  geom_point(position = position_dodge(0.1), alpha = 0.5) +
  facet_wrap(~experiment, scales = "free")

# non-native density and damage
datplant %>%
  ggplot(aes(x = log(nonnative.density+1), y = mean.dam, color = origin)) +
  geom_point(position = position_dodge(0.1), alpha = 0.5) +
  facet_wrap(~experiment, scales = "free")

datplant %>%
  ggplot(aes(x = nondens.s, y = mean.dam, color = origin)) +
  geom_point(position = position_dodge(0.1), alpha = 0.5) +
  facet_wrap(~experiment, scales = "free")

datplant %>%
  ggplot(aes(x = log(nonnative.density+1), y = prop.dam, color = origin)) +
  geom_point(position = position_dodge(0.1), alpha = 0.5) +
  facet_wrap(~experiment, scales = "free")

datplant %>%
  ggplot(aes(x = nondens.s, y = prop.dam, color = origin)) +
  geom_point(position = position_dodge(0.1), alpha = 0.5) +
  facet_wrap(~experiment, scales = "free")


#### correlations ####

# native/non-native
plotdat %>%
  filter(experiment == "transect") %>%
  cor.test(~ native.density + nonnative.density, data = .) # negative cor, but not strong

plotdat %>%
  filter(experiment == "transect") %>%
  cor.test(~ other.density + nonnative.density, data = .) # no cor

plotdat %>%
  filter(experiment == "transect") %>%
  cor.test(~ native.density + other.density, data = .) # no cor

plotdat %>%
  filter(experiment == "competition") %>%
  cor.test(~ native.density + nonnative.density, data = .) # no cor

# density/relative
plotdat %>%
  cor.test(~ nondens.s + nonnative.rel, data = .) # highly sig cor

plotdat %>%
  filter(experiment == "transect") %>%
  cor.test(~ nondens.s + nonnative.rel, data = .) # highly sig cor

plotdat %>%
  filter(experiment == "competition") %>%
  cor.test(~ nondens.s + nonnative.rel, data = .) # highly sig cor


#### prop.dam by density ####

# transect, absolute abundance
pamodT <- glmmTMB(infected ~ nonnative * (natdens.s + nondens.s + othdens.s) + (1|year/subplot/plant), data = datTleaf, family = binomial)
summary(pamodT)
plot(simulateResiduals(pamodT))
pamodTa <- model.avg(get.models(dredge(pamodT), subset = cumsum(weight) <= .95))
summary(pamodTa)

# alternate form of model
pamodTb <- glmer(cbind(leaves.dam, leaves.healthy) ~ nonnative * (natdens.s + nondens.s + othdens.s) + (1|year/subplot), data = datTplant, family = binomial) # failed to converge

# transect, relative abundance
prmodT <- glmmTMB(infected ~ nonnative * nonnative.rel + (1|year/subplot/plant), data = datTleaf, family = binomial)
summary(prmodT)
plot(simulateResiduals(prmodT))
prmodTa <- model.avg(get.models(dredge(prmodT), subset = cumsum(weight) <= .95))
summary(prmodTa)

# competition, absolute abundance
pamodC <- glmmTMB(infected ~ nonnative * (natdens.s + nondens.s) + (1|subplot/plant), data = datCleaf, family = binomial)
summary(pamodC)
plot(simulateResiduals(pamodC))
pamodCa <- model.avg(get.models(dredge(pamodC), subset = cumsum(weight) <= .95))
summary(pamodCa)

# competition, relative abundance
prmodC <- glmmTMB(infected ~ nonnative * nonnative.rel + (1|subplot/plant), data = datCleaf, family = binomial)
summary(prmodC)
plot(simulateResiduals(prmodC))
prmodCa <- model.avg(get.models(dredge(prmodC), subset = cumsum(weight) <= .95))
summary(prmodCa)


#### mean.dam by density ####

# subset data
mdatT <- filter(datTplant, mean.dam > 0)
mdatC <- filter(datCplant, mean.dam > 0)

# transect, absolute abundance
mamodT <- glmmTMB(mean.dam ~ nonnative * (natdens.s + nondens.s + othdens.s) + (1|year/subplot), data = mdatT, family = beta_family)
summary(mamodT)
plot(simulateResiduals(mamodT)) # significant deviation of expected from observed
mamodTa <- model.avg(get.models(dredge(mamodT), subset = cumsum(weight) <= .95)) # warnings: value out of range in 'lgamma', safe to ignore: https://mran.microsoft.com/snapshot/2017-11-04/web/packages/glmmTMB/vignettes/troubleshooting.html
summary(mamodTa)

# transect, relative abundance
mrmodT <- glmmTMB(mean.dam ~ nonnative * nonnative.rel + (1|year/subplot), data = mdatT, family = beta_family)
summary(mrmodT)
plot(simulateResiduals(mrmodT))
mrmodTa <- model.avg(get.models(dredge(mrmodT), subset = cumsum(weight) <= .95))
summary(mrmodTa)

# competition, absolute abundance
mamodC <- glmmTMB(mean.dam ~ nonnative * (natdens.s + nondens.s) + (1|subplot), data = mdatC, family = beta_family)
summary(mamodC)
plot(simulateResiduals(mamodC))
mamodCa <- model.avg(get.models(dredge(mamodC), subset = cumsum(weight) <= .95))
summary(mamodCa)

# competition, relative abundance
mrmodC <- glmmTMB(mean.dam ~ nonnative * nonnative.rel + (1|subplot), data = mdatC, family = beta_family)
summary(mrmodC)
plot(simulateResiduals(mrmodC))
mrmodCa <- model.avg(get.models(dredge(mrmodC), subset = cumsum(weight) <= .95))
summary(mrmodCa)


#### outputs ####
# write_csv(plotdat, "./output/plot_scale_density_data.csv") # exported this in damage_density_analysis

write_csv(mdatT, "./output/damage_density_experiment_meandam_transect_data.csv")
write_csv(mdatC, "./output/damage_density_experiment_meandam_competition_data.csv")

write_csv(datTleaf, "./output/damage_density_experiment_propdam_transect_data.csv")
write_csv(datCleaf, "./output/damage_density_experiment_propdam_competition_data.csv")

write_csv(datTplant, "./output/damage_density_experiment_plant_transect_data.csv")
write_csv(datCplant, "./output/damage_density_experiment_plant_competition_data.csv")

save(mamodTa, file = "./output/damage_density_experiment_meandam_absolute_transect_avg_amodel.rda")
save(mamodCa, file = "./output/damage_density_experiment_meandam_absolute_competition_avg_amodel.rda")

save(pamodTa, file = "./output/damage_density_experiment_propdam_absolute_transect_avg_amodel.rda")
save(pamodCa, file = "./output/damage_density_experiment_propdam_absolute_competition_avg_amodel.rda")