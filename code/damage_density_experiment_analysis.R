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
  left_join(bg15) %>%
  mutate(nonnative = ifelse(host %in% c("EG", "SP"), 0, 1),
         origin = ifelse(nonnative == 1, "non-native", "native"),
         total.density = native.density + nonnative.density)

dat15leaf <- dam15leaf %>%
  left_join(bg15) %>%
  mutate(nonnative = ifelse(host %in% c("EG", "SP"), 0, 1),
         origin = ifelse(nonnative == 1, "non-native", "native"),
         total.density = native.density + nonnative.density) %>%
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
         total.density = native.density + nonnative.density)

dat16Tleaf <- dam16Tleaf %>%
  left_join(bg16T) %>%
  mutate(nonnative = ifelse(host %in% c("EG", "SP"), 0, 1),
         origin = ifelse(nonnative == 1, "non-native", "native"),
         total.density = native.density + nonnative.density) %>%
  filter(!is.na(infected))

filter(dat16Tplant, is.na(native.density) | is.na(nonnative.density)) # all plots have info
filter(dat16Tleaf, is.na(native.density) | is.na(nonnative.density)) # all plots have info

# combine 2015 and 2016 transect
datTplant <- full_join(dat15plant, dat16Tplant) %>%
  mutate(natdens.s = scale(native.density)[,1],
         nondens.s = scale(nonnative.density)[,1])

datTleaf <- full_join(dat15leaf, dat16Tleaf) %>%
  mutate(natdens.s = scale(native.density)[,1],
         nondens.s = scale(nonnative.density)[,1])

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
         nondens.s = scale(nonnative.density)[,1])

datCleaf <- dam16Cleaf %>%
  left_join(bg16C %>% mutate(bg.species = recode(bg.species, "EG adults" = "EG adult"))) %>%
  mutate(nonnative = ifelse(host %in% c("EG", "SP"), 0, 1),
         origin = ifelse(nonnative == 1, "non-native", "native"),
         total.density = native.density + nonnative.density,
         natdens.s = scale(native.density)[,1],
         nondens.s = scale(nonnative.density)[,1]) %>%
  filter(!is.na(infected))

filter(dat16Cplant, is.na(native.density) | is.na(nonnative.density)) # all plots have info
filter(dat16Cleaf, is.na(native.density) | is.na(nonnative.density)) # all plots have info

# all plants
datplant <- full_join(datTplant, datCplant)

datleaf <- full_join(datTleaf, datCleaf)

# plot-scale dataset
plotdat <- datplant %>%
  dplyr::select(year, experiment, plot, subplot, bg.species, competition.density, native.density, nonnative.density, total.density, natdens.s, nondens.s) %>%
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
  cor.test(~ native.density + nonnative.density, data = .) # no cor

plotdat %>%
  filter(experiment == "competition") %>%
  cor.test(~ native.density + nonnative.density, data = .) # no cor

plotdat %>%
  filter(experiment == "transect") %>%
  cor.test(~ natdens.s + nondens.s, data = .) # no cor

plotdat %>%
  filter(experiment == "competition") %>%
  cor.test(~ natdens.s + nondens.s, data = .) # no cor


#### prop.dam by density ####

# transect
pmodT <- glmmTMB(infected ~ nonnative * (natdens.s + nondens.s) + (1|year/subplot/plant), data = datTleaf, family = binomial)
summary(pmodT) # non-native have lower
plot(simulateResiduals(pmodT))
pmodTa <- model.avg(get.models(dredge(pmodT), subset = cumsum(weight) <= .95))
summary(pmodTa) # origin is the most influential, native and non-native density are similar (both increase)

# alternate form of model
pmodTb <- glmer(cbind(leaves.dam, leaves.tot) ~ nonnative * (natdens.s + nondens.s) + (1|year/subplot), data = datTplant, family = binomial) # singular fit

# competition
pmodC <- glmmTMB(infected ~ nonnative * (natdens.s + nondens.s) + (1|subplot/plant), data = datCleaf, family = binomial)
summary(pmodC) # non-native have lower
plot(simulateResiduals(pmodC))
pmodCa <- model.avg(get.models(dredge(pmodC), subset = cumsum(weight) <= .95))
summary(pmodCa) # origin is the most influential, native and non-native are  a third as important (both decrease)


#### mean.dam by density ####

# subset data
mdatT <- filter(datTplant, mean.dam > 0)
mdatC <- filter(datCplant, mean.dam > 0)

# transect
mmodT <- glmmTMB(mean.dam ~ nonnative * (natdens.s + nondens.s) + (1|year/subplot), data = mdatT, family = beta_family)
summary(mmodT) # no sig effects
plot(simulateResiduals(mmodT))
mmodTa <- model.avg(get.models(dredge(mmodT), subset = cumsum(weight) <= .95))
summary(mmodTa) # origin is the most influential, followed closely by non-native and then native density (both increase)

# competition
mmodC <- glmmTMB(mean.dam ~ nonnative * (natdens.s + nondens.s) + (1|subplot), data = mdatC, family = beta_family)
summary(mmodC) # non-native have lower
plot(simulateResiduals(mmodC))
mmodCa <- model.avg(get.models(dredge(mmodC), subset = cumsum(weight) <= .95))
summary(mmodCa) # origin the most important, followed by non-native and native density (both decrease)


#### outputs ####
# write_csv(plotdat, "./output/plot_scale_density_data.csv") # same file as in damage_density_analysis

write_csv(mdatT, "./output/damage_density_experiment_meandam_transect16_data.csv")
write_csv(mdatC, "./output/damage_density_experiment_meandam_competition_data.csv")

write_csv(datTleaf, "./output/damage_density_experiment_propdam_transect_data.csv")
write_csv(datCleaf, "./output/damage_density_experiment_propdam_competition_data.csv")

write_csv(datTplant, "./output/damage_density_experiment_plant_transect_data.csv")
write_csv(datCplant, "./output/damage_density_experiment_plant_competition_data.csv")

save(mmodT, file = "./output/damage_density_experiment_meandam_transect16_model.rda")
save(mmodC, file = "./output/damage_density_experiment_meandam_competition_model.rda")

save(pmodT, file = "./output/damage_density_experiment_propdam_transect16_model.rda")
save(pmodC, file = "./output/damage_density_experiment_propdam_competition_model.rda")

save(mmodTa, file = "./output/damage_density_experiment_meandam_transect16_avg_model.rda")
save(mmodCa, file = "./output/damage_density_experiment_meandam_competition_avg_model.rda")

save(pmodTa, file = "./output/damage_density_experiment_propdam_transect16_avg_model.rda")
save(pmodCa, file = "./output/damage_density_experiment_propdam_competition_avg_model.rda")
