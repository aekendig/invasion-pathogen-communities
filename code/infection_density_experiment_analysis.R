## Goal: look at infection prevalence across density gradient

# notes: updated version of infection_experiment_analysis.R


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

# edit fungal data
dat1 <- dat %>%
  filter(experiment != "JEF transect" & 
           host %in% c("AB", "AF", "BRAD", "BD", "BH", "EG", "SP")) %>%
  mutate(grass_group = recode(grass.status, native = "native\nperennial", "non-native" = "non-native\nannual"),
         host_sp = recode(host, AB = "Avena barbata", AF = "Avena fatua", BRAD = "Brachypodium distachyonm", BD = "Bromus diandrus", BH = "Bromus hordeaceus", EG = "Elymus glaucus", SP = "Stipa pulchra"),
         host_sp = fct_relevel(host_sp, "Avena barbata", "Avena fatua", "Brachypodium distachyonm"))

# list of fungi to analyze
select(fun, pathogen, otu.id) %>% unique()

# combine transect background
bgT <- full_join(bg15T, bg16T)

# transect data
datT <- dat1 %>%
  inner_join(bgT) %>%
  mutate(nonnative = ifelse(host %in% c("EG", "SP"), 0, 1),
         total.density = native.density + nonnative.density + other.density,
         natdens.s = scale(native.density)[,1],
         nondens.s = scale(nonnative.density)[,1],
         othdens.s = scale(other.density)[,1],
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
datC <- dat1 %>%
  inner_join(bgC) %>%
  mutate(nonnative = ifelse(host %in% c("EG", "SP"), 0, 1),
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


#### sample sizes ####

# transect plots
datT %>%
  group_by(year) %>%
  summarise(n = length(unique(subplot)))

# transect isolates
datT %>%
  group_by(year, grass_group) %>%
  summarise(n = length(unique(isolate.id)))

# competition plots
datC %>%
  group_by(bg.species, competition.density, competition.type) %>%
  summarise(n = length(unique(subplot)))

# competition isolates
datC %>%
  group_by(grass_group) %>%
  summarise(n = length(unique(isolate.id)))


#### ainf ####

# check for presence
sum(datT$ainf) # 32
sum(datC$ainf) # 30

# transect, absolute abundance
aiamodT <- glmmTMB(ainf ~ nonnative * (natdens.s + nondens.s + othdens.s) + (1|year/subplot), data = datT, family = binomial)
summary(aiamodT) 
plot(simulateResiduals(aiamodT))
aiamodTa <- model.avg(get.models(dredge(aiamodT), subset = cumsum(weight) <= .95))
summary(aiamodTa) 

# transect, relative abundance
airmodT <- glmmTMB(ainf ~ nonnative * nonnative.rel + (1|year/subplot), data = datT, family = binomial)
summary(airmodT) 
plot(simulateResiduals(airmodT))
airmodTa <- model.avg(get.models(dredge(airmodT), subset = cumsum(weight) <= .95))
summary(airmodTa) 

# competition, absolute abundance
aiamodC <- glmmTMB(ainf ~ nonnative * (natdens.s + nondens.s) + (1|subplot), data = datC, family = binomial)
summary(aiamodC) 
plot(simulateResiduals(aiamodC))
aiamodCa <- model.avg(get.models(dredge(aiamodC), subset = cumsum(weight) <= .95))
summary(aiamodCa) 

# competition, relative abundance
airmodC <- glmmTMB(ainf ~ nonnative * nonnative.rel + (1|subplot), data = datC, family = binomial)
summary(airmodC) 
plot(simulateResiduals(airmodC))
airmodCa <- model.avg(get.models(dredge(airmodC), subset = cumsum(weight) <= .95))
summary(airmodCa) 


#### rpro ####

# check for presence
sum(datT$rpro) # 16
sum(datC$rpro) # 1

# check host type
filter(datT, rpro == 1) %>% select(year, nonnative, natdens.s, nondens.s, othdens.s) %>% unique() # lots of combinations, but only two are from 2016

# subset data if needed
datT15 <- filter(datT, year == 2015)

# transect, absolute abundance
rpamodT <- glmmTMB(rpro ~ nonnative * (natdens.s + nondens.s + othdens.s) + (1|year/subplot), data = datT, family = binomial) # model convergence error
rpamodT2 <- glmmTMB(rpro ~ nonnative * (natdens.s + nondens.s) + othdens.s + (1|year/subplot), data = datT, family = binomial)
summary(rpamodT2) 
plot(simulateResiduals(rpamodT2))
rpamodTa <- model.avg(get.models(dredge(rpamodT2), subset = cumsum(weight) <= .95))
summary(rpamodTa)

# transect, relative abundance
rprmodT <- glmmTMB(rpro ~ nonnative * nonnative.rel + (1|year/subplot), data = datT, family = binomial)
summary(rprmodT) 
plot(simulateResiduals(rprmodT))


#### pcha ####

# check for presence
sum(datT$pcha) # 11
sum(datC$pcha) # 8

# nonnative only?
filter(datT, pcha == 1) %>% select(grass_group) %>% unique() # yes
filter(datC, pcha == 1) %>% select(grass_group) %>% unique() # yes

# subset data
nondatT <- filter(datT, grass_group == "non-native\nannual")
nondatC <- filter(datC, grass_group == "non-native\nannual")

# transect, absolute abundance
pcamodT <- glmmTMB(pcha ~ natdens.s + nondens.s + othdens.s + (1|year/subplot), data = nondatT, family = binomial)
summary(pcamodT)
plot(simulateResiduals(pcamodT))
pcamodTa <- model.avg(get.models(dredge(pcamodT), subset = cumsum(weight) <= .95))
summary(pcamodTa)

# transect, relative abundance
pcrmodT <- glmmTMB(pcha ~ nonnative.rel + (1|year/subplot), data = nondatT, family = binomial)
summary(pcrmodT) 
plot(simulateResiduals(pcrmodT))

# competition, absolute abundance
pcamodC <- glmmTMB(pcha ~ natdens.s + nondens.s + (1|subplot), data = nondatC, family = binomial)
summary(pcamodC)
plot(simulateResiduals(pcamodC))
pcamodCa <- model.avg(get.models(dredge(pcamodC), subset = cumsum(weight) <= .95))
summary(pcamodCa) 

# competition, relative abundance
pcrmodC <- glmmTMB(pcha ~ nonnative.rel + (1|subplot), data = nondatC, family = binomial)
summary(pcrmodC) 
plot(simulateResiduals(pcrmodC))

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
plamodT <- glmmTMB(plol ~ nonnative * (natdens.s + nondens.s + othdens.s) + (1|year/subplot), data = datT, family = binomial) # convergence issue
plamodT2 <- glmmTMB(plol ~ nonnative * (natdens.s + nondens.s) + othdens.s + (1|year/subplot), data = datT, family = binomial) 
summary(plamodT2)
plot(simulateResiduals(plamodT2))
plamodTa <- model.avg(get.models(dredge(plamodT2), subset = cumsum(weight) <= .95))
summary(plamodTa)

# transect, relative abundance
plrmodT <- glmmTMB(plol ~ nonnative * nonnative.rel + (1|year/subplot), data = datT, family = binomial)
summary(plrmodT)
plot(simulateResiduals(plrmodT))

# competition, absolute abundance
plamodC <- glmmTMB(plol ~ nonnative * (natdens.s + nondens.s) + (1|subplot), data = datC, family = binomial) # model convergence error
filter(datC, plol == 1) %>% select(nonnative, natdens.s, nondens.s) # unclear what the issue is
plamodC2 <- glmmTMB(plol ~ nonnative + natdens.s + nondens.s + (1|subplot), data = datC, family = binomial)
summary(plamodC2)
plot(simulateResiduals(plamodC2))
plamodCa <- model.avg(get.models(dredge(plamodC2), subset = cumsum(weight) <= .95))
summary(plamodCa)

# competition, relative abundance
plrmodC <- glmmTMB(plol ~ nonnative * nonnative.rel + (1|subplot), data = datC, family = binomial)
summary(plrmodC)
plot(simulateResiduals(plrmodC))


#### ptri ####

# check for presence
sum(datT$ptri) # 20
sum(datC$ptri) # 3

# native only?
filter(datT, ptri == 1) %>% select(grass_group) %>% unique() # yes
filter(datC, ptri == 1) %>% select(grass_group) %>% unique() # yes

# subset data
natdatT <- filter(datT, grass_group == "native\nperennial")
natdatC <- filter(datC, grass_group == "native\nperennial")

# transect, absolute abundance
ptamodT <- glmmTMB(ptri ~ natdens.s + nondens.s + othdens.s + (1|year/subplot), data = natdatT, family = binomial)
summary(ptamodT) 
plot(simulateResiduals(ptamodT))
ptamodTa <- model.avg(get.models(dredge(ptamodT), subset = cumsum(weight) <= .95))
summary(ptamodTa)

# transect, relative abundance
ptrmodT <- glmmTMB(ptri ~ nonnative.rel + (1|year/subplot), data = natdatT, family = binomial)
summary(ptrmodT) 
plot(simulateResiduals(ptrmodT))

# competition, absolute abundance
ptamodC <- glmmTMB(ptri ~ natdens.s + nondens.s + (1|subplot), data = natdatC, family = binomial)
summary(ptamodC)
plot(simulateResiduals(ptamodC))
ptamodCa <- model.avg(get.models(dredge(ptamodC), subset = cumsum(weight) <= .95))
summary(ptamodCa) 

# competition, relative abundance
ptrmodC <- glmmTMB(ptri ~ nonnative.rel + (1|subplot), data = natdatC, family = binomial)
summary(ptrmodC)
plot(simulateResiduals(ptrmodC))

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
sum(datT$dres) # 7
sum(datC$dres) # 68

# check year difference
datT %>% 
  group_by(year, grass_group) %>% 
  summarise(dr = sum(dres),
            n = length(dres),
            dr.prop = dr/n) # only in natives in 2016, but in both in 2015

# transect, absolute abundance
dramodT <- glmmTMB(dres ~ nonnative * (natdens.s + nondens.s + othdens.s) + (1|year/subplot), data = datT, family = binomial) # convergence problem
dramodT2 <- glmmTMB(dres ~ nonnative * (natdens.s + nondens.s) + othdens.s + (1|year/subplot), data = datT, family = binomial) 
summary(dramodT2)
plot(simulateResiduals(dramodT2))
dramodTa <- model.avg(get.models(dredge(dramodT2), subset = cumsum(weight) <= .95)) # convergence issue
dramodT3 <- glmmTMB(dres ~ nonnative + natdens.s + nondens.s + othdens.s + (1|year/subplot), data = datT, family = binomial) 
summary(dramodT3)
plot(simulateResiduals(dramodT3))
dramodTa <- model.avg(get.models(dredge(dramodT3), subset = cumsum(weight) <= .95))
summary(dramodTa)

# transect, relative abundance
drrmodT <- glmmTMB(dres ~ nonnative * nonnative.rel + (1|year/subplot), data = datT, family = binomial)
# convergence error
filter(datT, dres == 1) %>% select(nonnative.rel) %>% unique() # all very high nonnative rel abu.
drrmodT2 <- glmmTMB(dres ~ nonnative + nonnative.rel + (1|year/subplot), data = datT, family = binomial) # still a convergence problem

# competition, absolute abundance
dramodC <- glmmTMB(dres ~ nonnative * (natdens.s + nondens.s) + (1|subplot), data = datC, family = binomial)
summary(dramodC)
plot(simulateResiduals(dramodC))
dramodCa <- model.avg(get.models(dredge(dramodC), subset = cumsum(weight) <= .95))
summary(dramodCa) 

# competition, relative abundance
drrmodC <- glmmTMB(dres ~ nonnative * nonnative.rel + (1|subplot), data = datC, family = binomial)
summary(drrmodC) 
plot(simulateResiduals(drrmodC))


#### pave ####

# check for presence
sum(datT$pave) # 12
sum(datC$pave) # 24

# transect, absolute abundance
paamodT <- glmmTMB(pave ~ nonnative * (natdens.s + nondens.s + othdens.s) + (1|year/subplot), data = datT, family = binomial) # convergence issue
paamodT2 <- glmmTMB(pave ~ nonnative * (natdens.s + nondens.s) + othdens.s + (1|year/subplot), data = datT, family = binomial) # convergence issue
paamodT3 <- glmmTMB(pave ~ nonnative + natdens.s + nondens.s + (1|year/subplot), data = datT, family = binomial)
summary(paamodT3) 
plot(simulateResiduals(paamodT3))
paamodTa <- model.avg(get.models(dredge(paamodT3), subset = cumsum(weight) <= .95))
summary(paamodTa) 

# transect, relative abundance
parmodT <- glmmTMB(pave ~ nonnative * nonnative.rel + (1|year/subplot), data = datT, family = binomial)
summary(parmodT)
plot(simulateResiduals(parmodT))

# competition, absolute abundance
paamodC <- glmmTMB(pave ~ nonnative * (natdens.s + nondens.s) + (1|subplot), data = datC, family = binomial)
summary(paamodC)
plot(simulateResiduals(paamodC))
paamodCa <- model.avg(get.models(dredge(paamodC), subset = cumsum(weight) <= .95))
summary(paamodCa)

# competition, relative abundance
parmodC <- glmmTMB(pave ~ nonnative * nonnative.rel + (1|subplot), data = datC, family = binomial)
summary(parmodC)
plot(simulateResiduals(parmodC))


#### outputs ####

write_csv(datT, "./output/infect_density_experiment_transect_data.csv")
write_csv(datC, "./output/infect_density_experiment_competition_data.csv")

save(rpamodTa, file = "./output/infection_density_experiment_rpro_absolute_transect_avg_model.rda")
save(aiamodTa, file = "./output/infection_density_experiment_ainf_absolute_transect_avg_model.rda")
save(aiamodCa, file = "./output/infection_density_experiment_ainf_absolute_competition_avg_model.rda")
save(pcamodTa, file = "./output/infection_density_experiment_pcha_absolute_transect_avg_model.rda")
save(pcamodCa, file = "./output/infection_density_experiment_pcha_absolute_competition_avg_model.rda")
save(plamodTa, file = "./output/infection_density_experiment_plol_absolute_transect_avg_model.rda")
save(plamodCa, file = "./output/infection_density_experiment_plol_absolute_competition_avg_model.rda")
save(ptamodTa, file = "./output/infection_density_experiment_ptri_absolute_transect_avg_model.rda")
save(ptamodCa, file = "./output/infection_density_experiment_ptri_absolute_competition_avg_model.rda")
save(dramodTa, file = "./output/infection_density_experiment_dres_absolute_transect_avg_model.rda")
save(dramodCa, file = "./output/infection_density_experiment_dres_absolute_competition_avg_model.rda")
save(paamodTa, file = "./output/infection_density_experiment_pave_absolute_transect_avg_model.rda")
save(paamodCa, file = "./output/infection_density_experiment_pave_absolute_competition_avg_model.rda")
