## Goal: look at infection prevalence in year 2 based on density in year 1


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
fun <- read_csv("./output/taxonomy_species_origin_data.csv")


#### edit data ####

# list of fungi to analyze
select(fun, taxonomy, otu.id) %>% unique()

# combine transect background
bgT <- full_join(bg15T, bg16T)

# transect data - note that this has fewer columns than infection_density_experiment_analysis
datT <- dat %>%
  inner_join(bgT) %>%
  mutate(nonnative = ifelse(host %in% c("EG", "SP"), 0, 1),
         origin = ifelse(nonnative == 1, "non-native", "native"),
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

# years of data for each plot
plot_sub <- datT %>%
  group_by(plot, subplot) %>%
  summarise(nyears = length(unique(year)),
            years = as.character(list(unique(year)))) %>%
  filter(nyears > 1)
# six plots were sampled twice

# subset for isolates
datTfun <- datT %>%
  filter(subplot %in% plot_sub$subplot & year == 2016) %>%
  select(plot, subplot, host, origin, nonnative, otu.id, taxonomy, ainf, rpro, pcha, plol, ptri, dres, pave)
# 39 isolates

# subset for density
datTden <- bgT %>%
  filter(subplot %in% plot_sub$subplot & year == 2015) %>%
  select(plot:nonnative.density)
# 6 plots, mostly low native density

# combine datasets
datTsub <- full_join(datTfun, datTden) %>%
  mutate(natdens.s = scale(native.density)[,1],
         nondens.s = scale(nonnative.density)[,1])

# long version
datTlong <- datTsub %>%
  gather(key = pathogen, value = infection, -c(plot:taxonomy, native.density:nondens.s))

# sample sizes
datTlong %>%
  select(pathogen, infection) %>%
  group_by(pathogen) %>%
  summarise(n = sum(infection))


#### visualize ####

ggplot(datTlong, aes(x = native.density, y = infection, color = pathogen)) +
  stat_summary(geom = "point", fun.y = mean, position = position_dodge(0.3)) +
  stat_summary(geom = "line", fun.y = mean, position = position_dodge(0.3)) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, width = 0.1, position = position_dodge(0.3)) +
  facet_wrap(~origin)

ggplot(datTlong, aes(x = nonnative.density, y = infection, color = pathogen)) +
  stat_summary(geom = "point", fun.y = mean, position = position_dodge(0.3)) +
  stat_summary(geom = "line", fun.y = mean, position = position_dodge(0.3)) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, width = 0.1, position = position_dodge(0.3)) +
  facet_wrap(~origin)


#### models ####

# ainf
aiamodY <- glmmTMB(ainf ~ nonnative * (natdens.s + nondens.s) + (1|subplot), data = datTsub, family = binomial)
summary(aiamodY)
plot(simulateResiduals(aiamodY))
aiamodYa <- model.avg(get.models(dredge(aiamodY), subset = cumsum(weight) <= .95))
summary(aiamodYa) # origin and non-native density

# rpro
rpamodY <- glmmTMB(rpro ~ nonnative * (natdens.s + nondens.s) + (1|subplot), data = datTsub, family = binomial)
summary(rpamodY) # fit is bad

# pave
paamodY <- glmmTMB(pave ~ nonnative * (natdens.s + nondens.s) + (1|subplot), data = datTsub, family = binomial)
# won't converge
