## Goal: relationship between plant density and damage

# updated version is damage_density_experiment_analysis.R


#### set up ####

# clear everything except final data
rm(list = ls())

# load libraries
library(tidyverse)
library(glmmTMB)
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

# 2015
bg15
dam15plant
dam15leaf

dat15plant <- dam15plant %>%
  left_join(bg15) %>%
  mutate(nonnative = ifelse(host %in% c("EG", "SP"), 0, 1),
         origin = ifelse(nonnative == 1, "non-native", "native"),
         total.density = native.density + nonnative.density,
         natdens.s = scale(native.density)[,1],
         nondens.s = scale(nonnative.density)[,1])

dat15leaf <- dam15leaf %>%
  left_join(bg15) %>%
  mutate(nonnative = ifelse(host %in% c("EG", "SP"), 0, 1),
         origin = ifelse(nonnative == 1, "non-native", "native"),
         total.density = native.density + nonnative.density,
         natdens.s = scale(native.density)[,1],
         nondens.s = scale(nonnative.density)[,1]) %>%
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
         total.density = native.density + nonnative.density,
         natdens.s = scale(native.density)[,1],
         nondens.s = scale(nonnative.density)[,1])

dat16Tleaf <- dam16Tleaf %>%
  left_join(bg16T) %>%
  mutate(nonnative = ifelse(host %in% c("EG", "SP"), 0, 1),
         origin = ifelse(nonnative == 1, "non-native", "native"),
         total.density = native.density + nonnative.density,
         natdens.s = scale(native.density)[,1],
         nondens.s = scale(nonnative.density)[,1]) %>%
  filter(!is.na(infected))

filter(dat16Tplant, is.na(native.density) | is.na(nonnative.density)) # all plots have info
filter(dat16Tleaf, is.na(native.density) | is.na(nonnative.density)) # all plots have info

# 2016 competition
bg16C
dam16Cplant
dam16Cleaf

dat16Cplant <- dam16Cplant %>%
  left_join(bg16C) %>%
  mutate(nonnative = ifelse(host %in% c("EG", "SP"), 0, 1),
         origin = ifelse(nonnative == 1, "non-native", "native"),
         total.density = native.density + nonnative.density,
         natdens.s = scale(native.density)[,1],
         nondens.s = scale(nonnative.density)[,1])

dat16Cleaf <- dam16Cleaf %>%
  left_join(bg16C) %>%
  mutate(nonnative = ifelse(host %in% c("EG", "SP"), 0, 1),
         origin = ifelse(nonnative == 1, "non-native", "native"),
         total.density = native.density + nonnative.density,
         natdens.s = scale(native.density)[,1],
         nondens.s = scale(nonnative.density)[,1]) %>%
  filter(!is.na(infected))

filter(dat16Cplant, is.na(native.density) | is.na(nonnative.density)) # all plots have info
filter(dat16Cleaf, is.na(native.density) | is.na(nonnative.density)) # all plots have info

# all plants
datplant <- full_join(dat15plant, dat16Tplant) %>%
  full_join(dat16Cplant)

datleaf <- full_join(dat15leaf, dat16Tleaf) %>%
  full_join(dat16Cleaf)

# plot-scale dataset
plotdat <- datplant %>%
  dplyr::select(year, experiment, plot, subplot, bg.species, competition.density, native.density, nonnative.density, total.density, natdens.s, nondens.s) %>%
  unique() %>%
  mutate(year.f = factor(year))


#### visualize ####

# relationship between densities
plotdat %>%
  ggplot(aes(x = native.density, y = nonnative.density, color = experiment, shape = year.f)) +
  geom_point()

plotdat %>%
  filter(nonnative.density <= max(native.density)) %>%
  ggplot(aes(x = native.density, y = nonnative.density, color = experiment, shape = year.f)) +
  geom_point() # not the best way to down-sample - come back to this

plotdat %>%
  ggplot(aes(x = log(native.density+1), y = log(nonnative.density+1), color = experiment, shape = year.f)) +
  geom_point()

plotdat %>%
  ggplot(aes(x = natdens.s, y = nondens.s, color = experiment, shape = year.f)) +
  geom_point()

plotdat %>%
  ggplot(aes(x = log(native.density+1), y = log(nonnative.density+1))) +
  geom_point() + 
  facet_wrap(year~experiment, scales = "free")

plotdat %>%
  ggplot(aes(x = natdens.s, y = nondens.s)) +
  geom_point() + 
  facet_wrap(year~experiment, scales = "free")

# distribution of density
plotdat%>%
  ggplot(aes(x = native.density)) +
  geom_histogram()

plotdat%>%
  ggplot(aes(x = nonnative.density)) +
  geom_histogram() # want to make this plot look like native, find dist of native and then sample from non-native

# native density and damage
datplant %>%
  ggplot(aes(x = log(native.density+1), y = mean.dam, color = origin)) +
  geom_point(position = position_dodge(0.1), alpha = 0.5) +
  facet_wrap(year~experiment, scales = "free")

datplant %>%
  ggplot(aes(x = natdens.s, y = mean.dam, color = origin)) +
  geom_point(position = position_dodge(0.1), alpha = 0.5) +
  facet_wrap(year~experiment, scales = "free")

datplant %>%
  ggplot(aes(x = log(native.density+1), y = prop.dam, color = origin)) +
  geom_point(position = position_dodge(0.1), alpha = 0.5) +
  facet_wrap(year~experiment, scales = "free")

datplant %>%
  ggplot(aes(x = natdens.s, y = prop.dam, color = origin)) +
  geom_point(position = position_dodge(0.1), alpha = 0.5) +
  facet_wrap(year~experiment, scales = "free")

# non-native density and damage
datplant %>%
  ggplot(aes(x = log(nonnative.density+1), y = mean.dam, color = origin)) +
  geom_point(position = position_dodge(0.1), alpha = 0.5) +
  facet_wrap(year~experiment, scales = "free")

datplant %>%
  ggplot(aes(x = nondens.s, y = mean.dam, color = origin)) +
  geom_point(position = position_dodge(0.1), alpha = 0.5) +
  facet_wrap(year~experiment, scales = "free")

datplant %>%
  ggplot(aes(x = log(nonnative.density+1), y = prop.dam, color = origin)) +
  geom_point(position = position_dodge(0.1), alpha = 0.5) +
  facet_wrap(year~experiment, scales = "free")

datplant %>%
  ggplot(aes(x = nondens.s, y = prop.dam, color = origin)) +
  geom_point(position = position_dodge(0.1), alpha = 0.5) +
  facet_wrap(year~experiment, scales = "free")


#### correlations ####

# native/non-native
plotdat %>%
  filter(year == 2015 & experiment == "transect") %>%
  cor.test(~ native.density + nonnative.density, data = .) # no cor

plotdat %>%
  filter(year == 2016 & experiment == "transect") %>%
  cor.test(~ native.density + nonnative.density, data = .) # marginally sig, -0.4

plotdat %>%
  filter(year == 2016 & experiment == "competition") %>%
  cor.test(~ native.density + nonnative.density, data = .) # no cor

plotdat %>%
  filter(year == 2015 & experiment == "transect") %>%
  cor.test(~ natdens.s + nondens.s, data = .) # no cor

plotdat %>%
  filter(year == 2016 & experiment == "transect") %>%
  cor.test(~ natdens.s + nondens.s, data = .) # marginally sig, -0.4

plotdat %>%
  filter(year == 2016 & experiment == "competition") %>%
  cor.test(~ natdens.s + nondens.s, data = .) # no cor

# native/total
plotdat %>%
  filter(year == 2015 & experiment == "transect") %>%
  cor.test(~ native.density + total.density, data = .) # no cor

plotdat %>%
  filter(year == 2016 & experiment == "transect") %>%
  cor.test(~ native.density + total.density, data = .) # marginally sig, -0.4

plotdat %>%
  filter(year == 2016 & experiment == "competition") %>%
  cor.test(~ native.density + total.density, data = .) # no cor

# total/non-native
plotdat %>%
  filter(year == 2015 & experiment == "transect") %>%
  cor.test(~ total.density + nonnative.density, data = .) # highly cor

plotdat %>%
  filter(year == 2016 & experiment == "transect") %>%
  cor.test(~ total.density + nonnative.density, data = .) # highly cor

plotdat %>%
  filter(year == 2016 & experiment == "competition") %>%
  cor.test(~ total.density + nonnative.density, data = .) # highly cor


#### prop.dam by density ####

# 2015 transect - not sure if this formula is right - look up weights examples
pmod15 <- glmmTMB(infected ~ nonnative * (natdens.s + nondens.s) + (1|subplot/plant), data = dat15leaf, family = binomial)
summary(pmod15) # non-native have lower
plot(simulateResiduals(pmod15))
pmod15a <- model.avg(get.models(dredge(pmod15), subset = cumsum(weight) <= .95))
summary(pmod15a) # origin is the most influential, native and non-native density are similar (both increase)

# 2016 transect
pmod16T <- glmmTMB(infected ~ nonnative * (natdens.s + nondens.s) + (1|subplot/plant), data = dat16Tleaf, family = binomial)
summary(pmod16T) # non-native have lower
plot(simulateResiduals(pmod16T))
pmod16Ta <- model.avg(get.models(dredge(pmod16T), subset = cumsum(weight) <= .95))
summary(pmod16Ta) # origin is the most influential, native is half as important, then non-native (both increase)

# 2016 competition
pmod16C <- glmmTMB(infected ~ nonnative * (natdens.s + nondens.s) + (1|subplot/plant), data = dat16Cleaf, family = binomial)
summary(pmod16C) # non-native have lower
plot(simulateResiduals(pmod16C))
pmod16Ca <- model.avg(get.models(dredge(pmod16C), subset = cumsum(weight) <= .95))
summary(pmod16Ca) # origin is the most influential, native and non-native are a third of importance (both decrease)


#### mean.dam by density ####

# subset data
mdat15 <- filter(dat15plant, mean.dam > 0)
mdat16T <- filter(dat16Tplant, mean.dam > 0)
mdat16C <- filter(dat16Cplant, mean.dam > 0)

# 2015 transect
mmod15 <- glmmTMB(mean.dam ~ nonnative * (natdens.s + nondens.s) + (1|subplot), data = mdat15, family = beta_family)
summary(mmod15) # no sig effects
plot(simulateResiduals(mmod15))
mmod15a <- model.avg(get.models(dredge(mmod15), subset = cumsum(weight) <= .95))
summary(mmod15a) # native density is the most influential, followed by non-native (both increase) then origin

# 2016 transect
mmod16T <- glmmTMB(mean.dam ~ nonnative * (natdens.s + nondens.s) + (1|subplot), data = mdat16T, family = beta_family)
summary(mmod16T) # no sig effects
plot(simulateResiduals(mmod16T))
mmod16Ta <- model.avg(get.models(dredge(mmod16T), subset = cumsum(weight) <= .95))
summary(mmod16Ta) # origin is the most influential, but still not very, followed by native and non-native density (native increase, non-native decrease)

# 2016 competition
mmod16C <- glmmTMB(mean.dam ~ nonnative * (natdens.s + nondens.s) + (1|subplot), data = mdat16C, family = beta_family)
summary(mmod16C) # non-native have lower
plot(simulateResiduals(mmod16C))
mmod16Ca <- model.avg(get.models(dredge(mmod16C), subset = cumsum(weight) <= .95))
summary(mmod16Ca) # origin the most important, followed by non-native and native density (both decrease)


#### outputs ####
write_csv(plotdat, "./output/plot_scale_density_data.csv")

write_csv(mdat15, "./output/damage_density_meandam_transect15_data.csv")
write_csv(mdat16T, "./output/damage_density_meandam_transect16_data.csv")
write_csv(mdat16C, "./output/damage_density_meandam_competition_data.csv")

write_csv(dat15leaf, "./output/damage_density_propdam_transect15_data.csv")
write_csv(dat16Tleaf, "./output/damage_density_propdam_transect16_data.csv")
write_csv(dat16Cleaf, "./output/damage_density_propdam_competition_data.csv")

write_csv(dat15plant, "./output/damage_density_plant_transect15_data.csv")
write_csv(dat16Tplant, "./output/damage_density_plant_transect16_data.csv")
write_csv(dat16Cplant, "./output/damage_density_plant_competition_data.csv")

save(mmod15, file = "./output/damage_density_meandam_transect15_model.rda")
save(mmod16T, file = "./output/damage_density_meandam_transect16_model.rda")
save(mmod16C, file = "./output/damage_density_meandam_competition_model.rda")

save(pmod15, file = "./output/damage_density_propdam_transect15_model.rda")
save(pmod16T, file = "./output/damage_density_propdam_transect16_model.rda")
save(pmod16C, file = "./output/damage_density_propdam_competition_model.rda")

save(mmod15a, file = "./output/damage_density_meandam_transect15_avg_model.rda")
save(mmod16Ta, file = "./output/damage_density_meandam_transect16_avg_model.rda")
save(mmod16Ca, file = "./output/damage_density_meandam_competition_avg_model.rda")

save(pmod15a, file = "./output/damage_density_propdam_transect15_avg_model.rda")
save(pmod16Ta, file = "./output/damage_density_propdam_transect16_avg_model.rda")
save(pmod16Ca, file = "./output/damage_density_propdam_competition_avg_model.rda")
