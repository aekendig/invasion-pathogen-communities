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

# sample sizes
nrow(datC) + nrow(datT)
unique(datC$host)
unique(datT$host)


#### plot density correlations ####

datC %>%
  select(plot, subplot, natdens.s, nondens.s) %>%
  unique() %>%
  cor.test(~ natdens.s + nondens.s, data = .) # not correlated

datC %>%
  select(plot, subplot, native.density, nonnative.density) %>%
  unique() %>%
  cor.test(~ log(native.density + 1) + log(nonnative.density + 1), data = .) 


#### ainf ####

# check for presence
sum(datT$ainf) # 33
sum(datC$ainf) # 30

# transect, absolute abundance
aiamodT <- glmmTMB(ainf ~ nonnative * (natdens.s + nondens.s) + (1|year/subplot), data = datT, family = binomial)
summary(aiamodT) # non-native origin increases
plot(simulateResiduals(aiamodT))
aiamodTa <- model.avg(get.models(dredge(aiamodT), subset = cumsum(weight) <= .95))
summary(aiamodTa) # origin, 2/3 by nat and nonnat dens

# transect, relative abundance
airmodT <- glmmTMB(ainf ~ nonnative * nonnative.rel + (1|year/subplot), data = datT, family = binomial)
summary(airmodT) # none
plot(simulateResiduals(airmodT))
airmodTa <- model.avg(get.models(dredge(airmodT), subset = cumsum(weight) <= .95))
summary(airmodTa) # 1/3 by rel. ab, origin sig

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

# check host type
filter(datT, rpro == 1) %>% select(year, nonnative, natdens.s, nondens.s) %>% unique() # lots of combinations, but only two are from 2016

# transect, absolute abundance
rpamodT <- glmmTMB(rpro ~ nonnative * (natdens.s + nondens.s) + (1|year/subplot), data = datT, family = binomial)
summary(rpamodT) # none
plot(simulateResiduals(rpamodT))
rpamodTa <- model.avg(get.models(dredge(rpamodT), subset = cumsum(weight) <= .95))
summary(rpamodTa) # covergence issue (didn't exist when year was left out)
rpamodTb <- glmmTMB(rpro ~ nonnative + natdens.s + nondens.s + (1|year/subplot), data = datT, family = binomial) # convergence issue

# transect, relative abundance
rprmodT <- glmmTMB(rpro ~ nonnative * nonnative.rel + (1|year/subplot), data = datT, family = binomial) # convergence issue (didn't exist when year was left out)
rprmodTb <- glmmTMB(rpro ~ nonnative + nonnative.rel + (1|year/subplot), data = datT, family = binomial) 
summary(rprmodTb) # none
plot(simulateResiduals(rprmodTb))
rprmodTa <- model.avg(get.models(dredge(rprmodTb), subset = cumsum(weight) <= .95)) # convergence issue


#### pcha ####

# check for presence
sum(datT$pcha) # 11
sum(datC$pcha) # 8

# nonnative only?
filter(datT, pcha == 1) %>% select(origin) %>% unique() # yes
filter(datC, pcha == 1) %>% select(origin) %>% unique() # yes

# subset data
nondatT <- filter(datT, origin == "non-native")
nondatC <- filter(datC, origin == "non-native")

# transect, absolute abundance
pcamodT <- glmmTMB(pcha ~ nonnative * (natdens.s + nondens.s) + (1|year/subplot), data = datT, family = binomial)
summary(pcamodT) # output is unreliable - all z-values are 0
pcamodTb <- glmmTMB(pcha ~ natdens.s + nondens.s + (1|year/subplot), data = nondatT, family = binomial)
summary(pcamodTb) # none
plot(simulateResiduals(pcamodTb))
pcamodTa <- model.avg(get.models(dredge(pcamodTb), subset = cumsum(weight) <= .95))
summary(pcamodTa) # nondens 1/4, natdens close

# transect, relative abundance
pcrmodTb <- glmmTMB(pcha ~ nonnative.rel + (1|year/subplot), data = nondatT, family = binomial)
summary(pcrmodTb) # none, value for rel kind of high (-2)
plot(simulateResiduals(pcrmodTb))

# competition, absolute abundance
pcamodCb <- glmmTMB(pcha ~ natdens.s + nondens.s + (1|subplot), data = nondatC, family = binomial)
summary(pcamodCb) # none
plot(simulateResiduals(pcamodCb))
pcamodCa <- model.avg(get.models(dredge(pcamodCb), subset = cumsum(weight) <= .95))
summary(pcamodCa) # nondens = 1

# competition, relative abundance
pcrmodCb <- glmmTMB(pcha ~ nonnative.rel + (1|subplot), data = nondatC, family = binomial)
summary(pcrmodCb) # none
plot(simulateResiduals(pcrmodCb))

# species-specific?
nondatT %>%
  group_by(host) %>%
  summarise(pc = sum(pcha),
            n = length(pcha),
            pc.prop = pc/n)

nondatC %>%
  group_by(host) %>%
  summarise(pc = sum(pcha),
            n = length(pcha),
            pc.prop = pc/n)


#### plol ####

# check for presence
sum(datT$plol) # 24
sum(datC$plol) # 7

# transect, absolute abundance
plamodT <- glmmTMB(plol ~ nonnative * (natdens.s + nondens.s) + (1|year/subplot), data = datT, family = binomial)
summary(plamodT) # none
plot(simulateResiduals(plamodT))
plamodTa <- model.avg(get.models(dredge(plamodT), subset = cumsum(weight) <= .95))
summary(plamodTa) # origin, then native density

# transect, relative abundance
plrmodT <- glmmTMB(plol ~ nonnative * nonnative.rel + (1|year/subplot), data = datT, family = binomial)
summary(plrmodT) # none
plot(simulateResiduals(plrmodT))
plrmodTa <- model.avg(get.models(dredge(plrmodT), subset = cumsum(weight) <= .95))
summary(plrmodTa) # origin

# competition, absolute abundance
plamodC <- glmmTMB(plol ~ nonnative * (natdens.s + nondens.s) + (1|subplot), data = datC, family = binomial)
# model convergence error
filter(datC, plol == 1) %>% select(nonnative, natdens.s, nondens.s) # unclear what the issue is
plamodCb <- glmmTMB(plol ~ nonnative + natdens.s + nondens.s + (1|subplot), data = datC, family = binomial)
summary(plamodCb) # natdens
plot(simulateResiduals(plamodCb))
plamodCa <- model.avg(get.models(dredge(plamodCb), subset = cumsum(weight) <= .95))
summary(plamodCa) # natdens

# competition, relative abundance
plrmodC <- glmmTMB(plol ~ nonnative * nonnative.rel + (1|subplot), data = datC, family = binomial)
summary(plrmodC) # none
plot(simulateResiduals(plrmodC))
plrmodCa <- model.avg(get.models(dredge(plrmodC), subset = cumsum(weight) <= .95))
summary(plrmodCa) # main effects 0.22


#### ptri ####

# check for presence
sum(datT$ptri) # 20
sum(datC$ptri) # 3

# native only?
filter(datT, ptri == 1) %>% select(origin) %>% unique() # yes
filter(datC, ptri == 1) %>% select(origin) %>% unique() # yes

# subset data
natdatT <- filter(datT, origin == "native")
natdatC <- filter(datC, origin == "native")

# transect, absolute abundance
ptamodT <- glmmTMB(ptri ~ nonnative * (natdens.s + nondens.s) + (1|year/subplot), data = datT, family = binomial)
summary(ptamodT) # none, z values are zero
ptamodTb <- glmmTMB(ptri ~ natdens.s + nondens.s + (1|year/subplot), data = natdatT, family = binomial)
summary(ptamodTb) # none
plot(simulateResiduals(ptamodTb))
ptamodTa <- model.avg(get.models(dredge(ptamodTb), subset = cumsum(weight) <= .95))
summary(ptamodTa) # each density ~ 1/3

# transect, relative abundance
ptrmodTb <- glmmTMB(ptri ~ nonnative.rel + (1|year/subplot), data = natdatT, family = binomial)
summary(ptrmodTb) # none
plot(simulateResiduals(ptrmodTb))

# competition, absolute abundance
ptamodCb <- glmmTMB(ptri ~ natdens.s + nondens.s + (1|subplot), data = natdatC, family = binomial)
summary(ptamodCb) # none
plot(simulateResiduals(ptamodCb))
ptamodCa <- model.avg(get.models(dredge(ptamodCb), subset = cumsum(weight) <= .95))
summary(ptamodCa) # each density ~ 1/5

# competition, relative abundance
ptrmodCb <- glmmTMB(ptri ~ nonnative.rel + (1|subplot), data = natdatC, family = binomial)
summary(ptrmodCb) # none
plot(simulateResiduals(ptrmodCb))

# check host specificity
natdatT %>%
  group_by(host) %>%
  summarise(pt = sum(ptri),
            n = length(ptri),
            pt.prop = pt/n)

natdatC %>%
  group_by(host) %>%
  summarise(pt = sum(ptri),
            n = length(ptri),
            pt.prop = pt/n)


#### dres ####

# check for presence
sum(datT$dres) # 8
sum(datC$dres) # 68

# check year difference
datT %>% 
  group_by(year, origin) %>% 
  summarise(dr = sum(dres),
            n = length(dres),
            dr.prop = dr/n) # only in natives in 2016, but in both in 2015

# transect, absolute abundance
dramodT <- glmmTMB(dres ~ nonnative * (natdens.s + nondens.s) + (1|year/subplot), data = datT, family = binomial)
summary(dramodT) # origin:natdens, origin, nondens - only working with a few points to get interaction
plot(simulateResiduals(dramodT))
dramodTa <- model.avg(get.models(dredge(dramodT), subset = cumsum(weight) <= .95)) # convergence issue
dramodTb <- glmmTMB(dres ~ nonnative + natdens.s + nondens.s + (1|year/subplot), data = datT, family = binomial)
summary(dramodTb) # none
plot(simulateResiduals(dramodTb))
dramodTa <- model.avg(get.models(dredge(dramodTb), subset = cumsum(weight) <= .95)) 
summary(dramodTa) # main effects ~ 1/3

# transect, relative abundance
drrmodT <- glmmTMB(dres ~ nonnative * nonnative.rel + (1|year/subplot), data = datT, family = binomial)
# convergence error
filter(datT, dres == 1) %>% select(nonnative.rel) %>% unique() # all very high nonnative rel abu.
drrmodTb <- glmmTMB(dres ~ nonnative + nonnative.rel + (1|year/subplot), data = datT, family = binomial)
summary(drrmodTb) # none
plot(simulateResiduals(drrmodTb))
drrmodTa <- model.avg(get.models(dredge(drrmodTb), subset = cumsum(weight) <= .95))
summary(drrmodTa) # nonnative rel 1/3

# competition, absolute abundance
dramodC <- glmmTMB(dres ~ nonnative * (natdens.s + nondens.s) + (1|subplot), data = datC, family = binomial)
summary(dramodC) # nonnative-natdens int
plot(simulateResiduals(dramodC))
dramodCa <- model.avg(get.models(dredge(dramodC), subset = cumsum(weight) <= .95))
summary(dramodCa) # nat dens and nondens highly impt, followed by origin

# competition, relative abundance
drrmodC <- glmmTMB(dres ~ nonnative * nonnative.rel + (1|subplot), data = datC, family = binomial)
summary(drrmodC) # none
filter(datC, dres == 1) %>% select(nonnative.rel) %>% unique() # pretty high, but not exclusively
plot(simulateResiduals(drrmodC))
drrmodCa <- model.avg(get.models(dredge(drrmodC), subset = cumsum(weight) <= .95))
summary(drrmodCa) # origin > rel


#### pave ####

# check for presence
sum(datT$pave) # 12
sum(datC$pave) # 24

# transect, absolute abundance
paamodT <- glmmTMB(pave ~ nonnative * (natdens.s + nondens.s) + (1|year/subplot), data = datT, family = binomial)
summary(paamodT) # none
plot(simulateResiduals(paamodT))
paamodTa <- model.avg(get.models(dredge(paamodT), subset = cumsum(weight) <= .95)) # convergence
paamodTb <- glmmTMB(pave ~ nonnative + natdens.s + nondens.s + (1|year/subplot), data = datT, family = binomial)
summary(paamodTb) # origin
paamodTa <- model.avg(get.models(dredge(paamodTb), subset = cumsum(weight) <= .95))
summary(paamodTa) # origin

# transect, relative abundance
parmodT <- glmmTMB(pave ~ nonnative * nonnative.rel + (1|year/subplot), data = datT, family = binomial)
summary(parmodT) # none
plot(simulateResiduals(parmodT))
parmodTa <- model.avg(get.models(dredge(parmodT), subset = cumsum(weight) <= .95))
summary(parmodTa) # origin

# competition, absolute abundance
paamodC <- glmmTMB(pave ~ nonnative * (natdens.s + nondens.s) + (1|subplot), data = datC, family = binomial)
summary(paamodC) # none
plot(simulateResiduals(paamodC))
paamodCa <- model.avg(get.models(dredge(paamodC), subset = cumsum(weight) <= .95))
summary(paamodCa) # origin

# competition, relative abundance
parmodC <- glmmTMB(pave ~ nonnative * nonnative.rel + (1|subplot), data = datC, family = binomial)
summary(parmodC) # none
plot(simulateResiduals(parmodC))
parmodCa <- model.avg(get.models(dredge(parmodC), subset = cumsum(weight) <= .95))
summary(parmodCa) # origin


#### outputs ####

write_csv(datT, "./output/infect_density_experiment_transect_data.csv")
write_csv(datC, "./output/infect_density_experiment_competition_data.csv")

save(aiamodTa, file = "./output/infection_density_experiment_ainf_absolute_transect_avg_model.rda")
save(airmodTa, file = "./output/infection_density_experiment_ainf_relative_transect_avg_model.rda")
save(aiamodCa, file = "./output/infection_density_experiment_ainf_absolute_competition_avg_model.rda")
save(airmodCa, file = "./output/infection_density_experiment_ainf_relative_competition_avg_model.rda")

save(pcamodTa, file = "./output/infection_density_experiment_pcha_absolute_transect_avg_model.rda")
save(pcrmodTb, file = "./output/infection_density_experiment_pcha_relative_transect_avg_model.rda")
save(pcamodCa, file = "./output/infection_density_experiment_pcha_absolute_competition_avg_model.rda")
save(pcrmodCb, file = "./output/infection_density_experiment_pcha_relative_competition_avg_model.rda")

save(plamodTa, file = "./output/infection_density_experiment_plol_absolute_transect_avg_model.rda")
save(plrmodTa, file = "./output/infection_density_experiment_plol_relative_transect_avg_model.rda")
save(plamodCa, file = "./output/infection_density_experiment_plol_absolute_competition_avg_model.rda")
save(plrmodCa, file = "./output/infection_density_experiment_plol_relative_competition_avg_model.rda")

save(ptamodTa, file = "./output/infection_density_experiment_ptri_absolute_transect_avg_model.rda")
save(ptrmodTb, file = "./output/infection_density_experiment_ptri_relative_transect_avg_model.rda")
save(ptamodCa, file = "./output/infection_density_experiment_ptri_absolute_competition_avg_model.rda")
save(ptrmodCb, file = "./output/infection_density_experiment_ptri_relative_competition_avg_model.rda")

save(dramodTa, file = "./output/infection_density_experiment_dres_absolute_transect_avg_model.rda")
save(drrmodTa, file = "./output/infection_density_experiment_dres_relative_transect_avg_model.rda")
save(dramodCa, file = "./output/infection_density_experiment_dres_absolute_competition_avg_model.rda")
save(drrmodCa, file = "./output/infection_density_experiment_dres_relative_competition_avg_model.rda")

save(paamodTa, file = "./output/infection_density_experiment_pave_absolute_transect_avg_model.rda")
save(parmodTa, file = "./output/infection_density_experiment_pave_relative_transect_avg_model.rda")
save(paamodCa, file = "./output/infection_density_experiment_pave_absolute_competition_avg_model.rda")
save(parmodCa, file = "./output/infection_density_experiment_pave_relative_competition_avg_model.rda")
