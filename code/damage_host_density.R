## Goal: How do native perennial and non-native annual grass densities affect disease severity?


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

dam15mleaf <- read_csv("./data/damage_leaf_transect_2015_march.csv")
dam15aleaf <- read_csv("./data/damage_leaf_transect_2015_april.csv")
dam16Tleaf <- read_csv("./data/damage_leaf_transect_2016.csv")
dam16Cleaf <- read_csv("./data/damage_leaf_competition_2016.csv")

dam15mplant <- read_csv("./data/damage_plant_transect_2015_march.csv")
dam15aplant <- read_csv("./data/damage_plant_transect_2015_april.csv")
dam16Tplant <- read_csv("./data/damage_plant_transect_2016.csv")
dam16Cplant <- read_csv("./data/damage_plant_competition_2016.csv")

# function to transform data to account for 0's and 1's
transform01 <- function(x) {
  (x * (length(x) - 1) + 0.5) / (length(x))
}

backtransform01 <- function(x) {
  (x * length(x) - 0.5) / (length(x) - 1)
}  


#### edit data ####

# 2015 transect
bg15
dam15mplant %>%
  select(year, experiment, subplot)
unique(dam15mplant$host)
dam15aplant %>%
  select(year, experiment, subplot)
unique(dam15aplant$host)

# add total density
bg15.2 <- bg15 %>%
  mutate(total.density = native.density + nonnative.density + other.density)

dat15mplant <- dam15mplant %>%
  filter(host != "PA") %>%
  inner_join(bg15.2) 

dat15aplant <- dam15aplant %>%
  filter(host != "PA") %>%
  inner_join(bg15.2)

dat15mleaf <- dam15mleaf %>%
  filter(host != "PA") %>%
  inner_join(bg15.2) %>%
  filter(!is.na(infected))

dat15aleaf <- dam15aleaf %>%
  filter(host != "PA") %>%
  inner_join(bg15.2) %>%
  filter(!is.na(infected))

filter(dat15mplant, is.na(native.density) | is.na(nonnative.density)) 
filter(dat15aplant, is.na(native.density) | is.na(nonnative.density)) %>%
  select(subplot) %>%
  unique()
filter(dat15mleaf, is.na(native.density) | is.na(nonnative.density))
filter(dat15aleaf, is.na(native.density) | is.na(nonnative.density))%>%
  select(subplot) %>%
  unique()

# 2016 transect
bg16T
dam16Tplant %>%
  select(year, experiment, subplot)
unique(dam16Tplant$host)

bg16T.2 <- bg16T %>%
  mutate(total.density = native.density + nonnative.density + other.density)

dat16Tplant <- dam16Tplant %>%
  inner_join(bg16T.2)

dat16Tleaf <- dam16Tleaf %>%
  inner_join(bg16T.2) %>%
  filter(!is.na(infected))

filter(dat16Tplant, is.na(native.density) | is.na(nonnative.density)) # all plots have info
filter(dat16Tleaf, is.na(native.density) | is.na(nonnative.density)) # all plots have info

# combine 2015 and 2016 transect
datTmplant <- full_join(dat15mplant, dat16Tplant %>%
                          mutate(month = NA_character_)) %>%
  mutate(mean.dam.scaled = transform01(mean.dam),
         natdens.s = scale(native.density)[,1],
         nondens.s = scale(nonnative.density)[,1],
         othdens.s = scale(other.density)[,1],
         nonnative.rel = nonnative.density / total.density,
         leaves.healthy = leaves.tot - leaves.dam)

datTaplant <- full_join(dat15aplant, dat16Tplant %>%
                          mutate(month = NA_character_)) %>%
  mutate(mean.dam.scaled = transform01(mean.dam),
         natdens.s = scale(native.density)[,1],
         nondens.s = scale(nonnative.density)[,1],
         othdens.s = scale(other.density)[,1],
         nonnative.rel = nonnative.density / total.density,
         leaves.healthy = leaves.tot - leaves.dam)

datTmleaf <- full_join(dat15mleaf, dat16Tleaf %>%
                         mutate(month = NA_character_)) %>%
  mutate(surface.scaled = transform01(surface),
         natdens.s = scale(native.density)[,1],
         nondens.s = scale(nonnative.density)[,1],
         othdens.s = scale(other.density)[,1],
         nonnative.rel = nonnative.density / total.density)

datTaleaf <- full_join(dat15aleaf, dat16Tleaf %>%
                         mutate(month = NA_character_)) %>%
  mutate(surface.scaled = transform01(surface),
         natdens.s = scale(native.density)[,1],
         nondens.s = scale(nonnative.density)[,1],
         othdens.s = scale(other.density)[,1],
         nonnative.rel = nonnative.density / total.density)

# 2016 competition
bg16C
dam16Cplant %>%
  select(year, experiment, subplot, bg.species, competition.density)
unique(dam16Cplant$host)
unique(dam16Cplant$bg.species)

datCplant <- dam16Cplant %>%
  inner_join(bg16C %>% mutate(bg.species = recode(bg.species, "EG adults" = "EG adult"))) %>%
  mutate(mean.dam.scaled = transform01(mean.dam),
         total.density = native.density + nonnative.density,
         natdens.s = scale(native.density)[,1],
         nondens.s = scale(nonnative.density)[,1],
         nonnative.rel = nonnative.density / total.density,
         leaves.healthy = leaves.tot - leaves.dam)

datCleaf <- dam16Cleaf %>%
  inner_join(bg16C %>% mutate(bg.species = recode(bg.species, "EG adults" = "EG adult"))) %>%
  mutate(surface.scaled = transform01(surface),
         total.density = native.density + nonnative.density,
         natdens.s = scale(native.density)[,1],
         nondens.s = scale(nonnative.density)[,1],
         nonnative.rel = nonnative.density / total.density) %>%
  filter(!is.na(infected))

filter(datCplant, is.na(native.density) | is.na(nonnative.density)) # all plots have info
filter(datCleaf, is.na(native.density) | is.na(nonnative.density)) # all plots have info

# create datasets without zeros
datTmleaf2 <- filter(datTmleaf, surface > 0)
datTaleaf2 <- filter(datTaleaf, surface > 0)
datCleaf2 <- filter(datCleaf, surface > 0)
# all plants
dataplant <- full_join(datTaplant, datCplant)
datmplant <- full_join(datTmplant, datCplant)

# sample sizes
nrow(datmplant)
nrow(dataplant)

# plot-scale dataset
plotmdat <- datmplant %>%
  select(year, experiment, subplot, bg.species, competition.density, native.density:nonnative.rel) %>%
  unique() %>%
  mutate(year.f = factor(year))

plotadat <- dataplant %>%
  select(year, experiment, subplot, bg.species, competition.density, native.density:nonnative.rel) %>%
  unique() %>%
  mutate(year.f = factor(year))


#### visualize ####

# native density and damage
datmplant %>%
  ggplot(aes(x = native.density, y = mean.dam.scaled, color = grass.group)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "glm") +
  facet_wrap(~experiment, scales = "free")

dataplant %>%
  ggplot(aes(x = native.density, y = mean.dam.scaled, color = grass.group)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "glm") +
  facet_wrap(~experiment, scales = "free")

datmplant %>%
  ggplot(aes(x = natdens.s, y = prop.dam, color = grass.group)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "glm") +
  facet_wrap(~experiment, scales = "free")

dataplant %>%
  ggplot(aes(x = natdens.s, y = prop.dam, color = grass.group)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "glm") +
  facet_wrap(~experiment, scales = "free")

# non-native density and damage
datmplant %>%
  ggplot(aes(x = nonnative.density, y = mean.dam.scaled, color = grass.group)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "glm") +
  facet_wrap(~experiment, scales = "free")

dataplant %>%
  ggplot(aes(x = nonnative.density, y = mean.dam.scaled, color = grass.group)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "glm") +
  facet_wrap(~experiment, scales = "free")

datmplant %>%
  ggplot(aes(x = nondens.s, y = prop.dam, color = grass.group)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "glm") +
  facet_wrap(~experiment, scales = "free")

dataplant %>%
  ggplot(aes(x = nondens.s, y = prop.dam, color = grass.group)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "glm") +
  facet_wrap(~experiment, scales = "free")


#### correlations ####

# transect March
plotmdat %>%
  filter(experiment == "transect") %>%
  cor.test(~ natdens.s + nondens.s, data = .)

plotmdat %>%
  filter(experiment == "transect") %>%
  cor.test(~ natdens.s + othdens.s, data = .)

plotmdat %>%
  filter(experiment == "transect") %>%
  cor.test(~ nondens.s + othdens.s, data = .)

# transect April
plotadat %>%
  filter(experiment == "transect") %>%
  cor.test(~ natdens.s + nondens.s, data = .)
# weakly negative

plotadat %>%
  filter(experiment == "transect") %>%
  cor.test(~ natdens.s + othdens.s, data = .)

plotadat %>%
  filter(experiment == "transect") %>%
  cor.test(~ nondens.s + othdens.s, data = .)

# Competition
plotmdat %>%
  filter(experiment == "competition") %>%
  cor.test(~ natdens.s + nondens.s, data = .)


#### infected by density ####

# transect March
imodTM <- glmmTMB(infected ~ nonnative * (natdens.s + nondens.s + othdens.s) + (1|year) + (1|subplot/plant), data = datTmleaf, family = binomial)
summary(imodTM)
plot(simulateResiduals(imodTM))
# sig deviation
imodTMa <- model.avg(get.models(dredge(imodTM), subset = cumsum(weight) <= .95))
summary(imodTMa)

# transect April
imodTA <- glmmTMB(infected ~ nonnative * (natdens.s + nondens.s + othdens.s) + (1|year) + (1|subplot/plant), data = datTaleaf, family = binomial)
summary(imodTA)
plot(simulateResiduals(imodTA))
# sig deviation
imodTAa <- model.avg(get.models(dredge(imodTA), subset = cumsum(weight) <= .95))
summary(imodTAa)

# competition
imodC <- glmmTMB(infected ~ nonnative * (natdens.s + nondens.s) + (1|subplot/plant), data = datCleaf, family = binomial)
summary(imodC)
plot(simulateResiduals(imodC))
imodCa <- model.avg(get.models(dredge(imodC), subset = cumsum(weight) <= .95))
summary(imodCa)


#### surface by density ####

# transect March
smodTM <- glmmTMB(surface ~ nonnative * (natdens.s + nondens.s + othdens.s) + (1|year) + (1|subplot/plant), data = datTmleaf2, family = beta_family)
summary(smodTM)
plot(simulateResiduals(smodTM))
# sig deviation
smodTMa <- model.avg(get.models(dredge(smodTM), subset = cumsum(weight) <= .95))
# multiple model convergence errors
summary(smodTMa)
# no estimates

# transect March, remove other density to simplify model
smodTM2 <- glmmTMB(surface ~ nonnative * (natdens.s + nondens.s) + (1|year) + (1|subplot/plant), data = datTmleaf2, family = beta_family)
summary(smodTM2)
plot(simulateResiduals(smodTM2))
# sig deviation
smodTM2a <- model.avg(get.models(dredge(smodTM2), subset = cumsum(weight) <= .95))
# multiple model convergence errors

# transect April
smodTA <- glmmTMB(surface ~ nonnative * (natdens.s + nondens.s + othdens.s) + (1|year) + (1|subplot/plant), data = datTaleaf2, family = beta_family)
# model convergence issue
summary(smodTA)

# transect April, remove other density to simplify model
smodTA2 <- glmmTMB(surface ~ nonnative * (natdens.s + nondens.s) + (1|year) + (1|subplot/plant), data = datTaleaf2, family = beta_family)
# model convergence issue
summary(smodTA2)

# transect April, remove interactions
smodTA3 <- glmmTMB(surface ~ nonnative + natdens.s + nondens.s + (1|year) + (1|subplot/plant), data = datTaleaf2, family = beta_family)
# model convergence issue
summary(smodTA3)

# year has a small variance
filter(datTaleaf2, year == 2016) %>%
  nrow()

# transect April, remove year random effect
smodTA4 <- glmmTMB(surface ~ nonnative * (natdens.s + nondens.s + othdens.s) + (1|subplot/plant), data = datTaleaf2, family = beta_family)
# model convergence issue
summary(smodTA4)
plot(simulateResiduals(smodTA4))
# sig deviation
smodTAa <- model.avg(get.models(dredge(smodTA4), subset = cumsum(weight) <= .95))
# multiple model convergence errors
summary(smodTAa)

# competition
smodC <- glmmTMB(surface ~ nonnative * (natdens.s + nondens.s) + (1|subplot/plant), data = datCleaf2, family = beta_family)
summary(smodC)
plot(simulateResiduals(smodC))
# sig deviation
smodCa <- model.avg(get.models(dredge(smodC), subset = cumsum(weight) <= .95))
summary(smodCa)


#### surface scaled by density ####

# transect March
scmodTM <- glmmTMB(surface.scaled ~ nonnative * (natdens.s + nondens.s + othdens.s) + (1|year) + (1|subplot/plant), data = datTmleaf, family = beta_family)
summary(scmodTM)
plot(simulateResiduals(scmodTM))
# sig deviation
scmodTMa <- model.avg(get.models(dredge(scmodTM), subset = cumsum(weight) <= .95))
# couldn't converge

# transect March, remove other density to simplify model
scmodTM2 <- glmmTMB(surface.scaled ~ nonnative * (natdens.s + nondens.s) + (1|year) + (1|subplot/plant), data = datTmleaf, family = beta_family)
summary(scmodTM2)
plot(simulateResiduals(scmodTM2))
# sig deviation
scmodTM2a <- model.avg(get.models(dredge(scmodTM2), subset = cumsum(weight) <= .95))
summary(scmodTM2a)
# all NA's

# check model fit
datTmleaf %>%
  mutate(pred = predict(scmodTM, type = "response", newdata = datTmleaf)) %>%
  ggplot(aes(surface.scaled, pred)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0)

# transect April
scmodTA <- glmmTMB(surface.scaled ~ nonnative * (natdens.s + nondens.s + othdens.s) + (1|year) + (1|subplot/plant), data = datTaleaf, family = beta_family)
summary(scmodTA)
# all NA's

# transect April, remove other density to simplify model
scmodTA2 <- glmmTMB(surface.scaled ~ nonnative * (natdens.s + nondens.s) + (1|year) + (1|subplot/plant), data = datTaleaf, family = beta_family)
summary(scmodTA2)
scmodTAa <- model.avg(get.models(dredge(scmodTA2), subset = cumsum(weight) <= .95))
summary(scmodTAa)
# a bunch of NA's

# check model fit
datTaleaf %>%
  mutate(pred = predict(scmodTA2, type = "response", newdata = datTaleaf)) %>%
  ggplot(aes(surface.scaled, pred)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0)

# competition
scmodC <- glmmTMB(surface.scaled ~ nonnative * (natdens.s + nondens.s) + (1|subplot/plant), data = datCleaf, family = beta_family)
summary(scmodC)
plot(simulateResiduals(scmodC))
# sig deviation
scmodCa <- model.avg(get.models(dredge(scmodC), subset = cumsum(weight) <= .95))
summary(scmodCa)
# all NA's

# check model fit
datCleaf %>%
  mutate(pred = predict(scmodC, type = "response", newdata = datCleaf)) %>%
  ggplot(aes(surface.scaled, pred)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0)


#### mean damage per plant scaled by density ####

# transect March
scpmodTM <- glmmTMB(mean.dam.scaled ~ nonnative * (natdens.s + nondens.s + othdens.s) + (1|year) + (1|subplot), data = datTmplant, family = beta_family)
summary(scpmodTM)
# estimates similar to leaf-level model
plot(simulateResiduals(scpmodTM))
# sig deviation
scpmodTMa <- model.avg(get.models(dredge(scpmodTM), subset = cumsum(weight) <= .95))
summary(scpmodTMa)
# couldn't converge

# check model fit
datTmplant %>%
  mutate(pred = predict(scpmodTM, type = "response", newdata = datTmplant)) %>%
  ggplot(aes(mean.dam.scaled, pred)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0)

# transect April
scpmodTA <- glmmTMB(mean.dam.scaled ~ nonnative * (natdens.s + nondens.s + othdens.s) + (1|year) + (1|subplot), data = datTaplant, family = beta_family)
summary(scpmodTA)
plot(simulateResiduals(scpmodTA))
# sig deviation

# check model fit
datTaplant %>%
  mutate(pred = predict(scpmodTA, type = "response", newdata = datTaplant)) %>%
  ggplot(aes(mean.dam.scaled, pred)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0)

# competition
scpmodC <- glmmTMB(mean.dam.scaled ~ nonnative * (natdens.s + nondens.s) + (1|subplot), data = datCplant, family = beta_family)
summary(scpmodC)
# estimates similar to leaf-level model
plot(simulateResiduals(scpmodC))
# sig deviation

# check model fit
datCplant %>%
  mutate(pred = predict(scpmodC, type = "response", newdata = datCplant)) %>%
  ggplot(aes(mean.dam.scaled, pred)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0)

#### outputs ####
write_csv(plotmdat, "./output/plot_scale_damage_march_data.csv") 
write_csv(plotadat, "./output/plot_scale_damage_april_data.csv") 

write_csv(datTaplant, "./output/damage_host_density_plant_april_transect_data.csv")
write_csv(datTmplant, "./output/damage_host_density_plant_march_transect_data.csv")
write_csv(datCplant, "./output/damage_host_density_plant_competition_data.csv")

write_csv(datTaleaf, "./output/damage_host_density_leaf_april_transect_data.csv")
write_csv(datTmleaf, "./output/damage_host_density_leaf_march_transect_data.csv")
write_csv(datCleaf, "./output/damage_host_density_leaf_competition_data.csv")

save(imodTMa, file = "./output/damage_host_density_infection_march_transect_model.rda")
save(imodTAa, file = "./output/damage_host_density_infection_april_transect_model.rda")
save(imodCa, file = "./output/damage_host_density_infection_competition_model.rda")

save(smodTM, file = "./output/damage_host_density_surface_march_transect_model.rda")
save(smodTAa, file = "./output/damage_host_density_surface_april_transect_model.rda")
save(smodCa, file = "./output/damage_host_density_surface_competition_model.rda")

save(scmodTM, file = "./output/damage_host_density_surface_scaled_march_transect_model.rda")
save(scmodTA2, file = "./output/damage_host_density_surface_scaled_april_transect_model.rda")
save(scmodC, file = "./output/damage_host_density_surface_scaled_competition_model.rda")