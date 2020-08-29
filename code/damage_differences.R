## Goal: How does disease severity differ between native perennial and non-native annual grasses?


#### set up ####

# clear all existing data
rm(list=ls())

# load libraries
library(tidyverse)
library(cowplot)
library(glmmTMB)
library(DHARMa) # plot glmmTMB
library(MuMIn) # dredge
library(brms)

# import data

# leaf-scale data
dam15mleaf <- read_csv("./data/damage_leaf_transect_2015_march.csv")
dam15aleaf <- read_csv("./data/damage_leaf_transect_2015_april.csv")
dam16Tleaf <- read_csv("./data/damage_leaf_transect_2016.csv")
dam16Cleaf <- read_csv("./data/damage_leaf_competition_2016.csv")

# plant-scale data
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


#### edit raw data ####

# competition plots sampled
dam16Cplant %>%
  group_by(bg.species, competition.density) %>%
  count() %>%
  data.frame()

# combine experiments (March transect)
ldatM <- dam16Cleaf %>%
  select(-c(bg.species, competition.density)) %>%
  full_join(dam15mleaf  %>%
              filter(host != "PA")) %>%
  full_join(dam16Tleaf %>%
              mutate(month = NA_character_)) %>%
  mutate(surface.scaled = transform01(surface))

# filter for positive values
ldatM2 <- ldatM %>%
  filter(surface > 0)

# random effects
ldatM %>%
  group_by(year, experiment, nonnative) %>%
  summarise(subplot = length(unique(subplot)),
            plants = length(unique(plant)),
            leaves = mean(leaf))

ldatM2 %>%
  group_by(year, experiment, nonnative) %>%
  summarise(subplot = length(unique(subplot)),
            plants = length(unique(plant)),
            leaves = mean(leaf))

# combine experiments (April transect)
ldatA <- dam16Cleaf %>%
  select(-c(bg.species, competition.density)) %>%
  full_join(dam15aleaf %>%
              filter(host != "PA")) %>%
  full_join(dam16Tleaf %>%
              mutate(month = NA_character_)) %>%
  mutate(surface.scaled = transform01(surface))

# filter for positive values
ldatA2 <- ldatA %>%
  filter(surface > 0)

# random effects
ldatA %>%
  group_by(year, experiment, nonnative) %>%
  summarise(subplot = length(unique(subplot)),
            plants = length(unique(plant)),
            leaves = mean(leaf))

ldatA2 %>%
  group_by(year, experiment, nonnative) %>%
  summarise(subplot = length(unique(subplot)),
            plants = length(unique(plant)),
            leaves = mean(leaf))


#### sample sizes ####
dam15mplant %>%
  filter(host != "PA") %>% 
  group_by(grass.group) %>%
  summarise(n = n())

dam15mplant %>%
  filter(host != "PA") %>% 
  select(plot) %>%
  unique()

dam15aplant %>%
  filter(host != "PA") %>% 
  group_by(grass.group) %>%
  summarise(n = n())

dam15aplant %>%
  filter(host != "PA") %>% 
  select(plot) %>%
  unique()

dam16Tplant %>% 
  group_by(grass.group) %>%
  summarise(n = n())

dam16Tplant %>% 
  select(plot) %>%
  unique()

dam16Cplant %>% 
  group_by(grass.group) %>%
  summarise(n = n())

unique(dam15mplant$host)
unique(dam15aplant$host)
unique(dam16Tplant$host)
unique(dam16Cplant$host)

unique(dam15mplant$subplot)
unique(dam15aplant$subplot)
unique(dam16Tplant$subplot)
unique(dam16Cplant$subplot)


#### visualize ####

# March surface
ggplot(ldatM, aes(x = surface)) +
  geom_histogram()

ggplot(ldatM, aes(surface, surface.scaled)) +
  geom_point()

# April surface
ggplot(ldatA, aes(x = surface)) +
  geom_histogram()

ggplot(ldatA, aes(surface, surface.scaled)) +
  geom_point()


#### statistical models ####

# March infection
imodM <- glmmTMB(infected ~ nonnative + (1|year) + (1|experiment/subplot/plant), data = ldatM, family = binomial)
summary(imodM)
plot(simulateResiduals(imodM))

# April infection
imodA <- glmmTMB(infected ~ nonnative  + (1|year) + (1|experiment/subplot/plant), data = ldatA, family = binomial)
summary(imodA)
plot(simulateResiduals(imodA))

# March surface
smodM <- glmmTMB(surface ~ nonnative + (1|year) + (1|experiment/subplot/plant), data = ldatM2, family = beta_family)
# convergence error
summary(smodM)
# variance of experiment is small, remove - accounted for with subplot
smodM2 <- glmmTMB(surface ~ nonnative + (1|year) + (1|subplot/plant), data = ldatM2, family = beta_family)
summary(smodM2)
plot(simulateResiduals(smodM2))
# sig deviation

# April surface
smodA <- glmmTMB(surface ~ nonnative + (1|year) + (1|experiment/subplot/plant), data = ldatA2, family = beta_family)
# convergence error
summary(smodA)
smodA2 <- glmmTMB(surface ~ nonnative + (1|year) + (1|subplot/plant), data = ldatA2, family = beta_family)
summary(smodA2)
plot(simulateResiduals(smodA2))
# sig deviation

# March zero-inflated
# zmodM <- brm(data = ldatM, family = zero_inflated_beta,
#              surface ~ nonnative + (1|year) + (1|experiment/subplot/plant), zi ~ nonnative,
#              prior = c(prior(normal(0, 10), class = b)),
#              iter = 6000, warmup = 1000, chains = 1)
# suuuuper slow

# March scaled
scmodM <- glmmTMB(surface.scaled ~ nonnative + (1|year) + (1|experiment/subplot/plant), data = ldatM, family = beta_family)
summary(scmodM)
# AIC values are missing - convergence error?
# variance of experiment is small, remove - accounted for with subplot
scmodM2 <- glmmTMB(surface.scaled ~ nonnative + (1|year) + (1|subplot/plant), data = ldatM, family = beta_family)
summary(scmodM2)
plot(simulateResiduals(scmodM2))

# April surface
scmodA <- glmmTMB(surface.scaled ~ nonnative + (1|year) + (1|experiment/subplot/plant), data = ldatA, family = beta_family)
summary(scmodA)
# no error, but should be comparable to March data
scmodA2 <- glmmTMB(surface.scaled ~ nonnative + (1|year) + (1|subplot/plant), data = ldatA, family = beta_family)
summary(scmodA2)
# doesn't have AIC value
plot(simulateResiduals(scmodA))


#### extract model values ####

# prediction function
pred_fun <- function(dat, mod, sev.type){
  
  # predict mean
  dat$pred <- predict(mod, newdata = dat, re.form = NA, type = "response")
  
  # predict SE
  dat$pred.se <- predict(mod, newdata = dat, re.form = NA, se.fit = T, type = "response")$se.fit
  
  # simplify data
  dat_out <- dat %>%
    select(grass.group, pred, pred.se) %>%
    unique() %>%
    mutate(sev.type = sev.type)
  
  # export
  return(dat_out)
}

# evaluate models
idatMpred <- pred_fun(ldatM, imodM, "leaves")
idatApred <- pred_fun(ldatA, imodA, "leaves")
sdatMpred <- pred_fun(ldatM2, smodM2, "surface")
sdatApred <- pred_fun(ldatA2, smodA2, "surface")

# combine
fdatM <- rbind(idatMpred, sdatMpred)
fdatA <- rbind(idatApred, sdatApred)


#### figure ####

# labels
labs <- tibble(let = LETTERS[1:2], 
               sev.type = c("leaves", "surface"), 
               grass.group = rep("native\nperennial", 2))

# figure settings
axisText = 10
axisTitle = 12
legendText = 10
legendTitle = 10

# figure
figM <- ggplot(fdatM, aes(x = grass.group, y = pred, fill = grass.group)) +
  geom_errorbar(aes(ymin = pred - pred.se, ymax = pred + pred.se), width = 0.1) +
  geom_point(size = 3, shape = 21) +
  scale_fill_manual(values = c("black", "white"), guide = F) +
  facet_wrap(~sev.type, ncol = 1, scales = "free") + 
  geom_text(data = labs,
            mapping = aes(x = -Inf, y = Inf, label = let),
            hjust = -0.3, vjust = 1.5,
            fontface = "bold", size = 4) +
  theme_bw() +
  theme(axis.text = element_text(size = axisText, color="black"),
        axis.title = element_text(size = axisTitle),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size=legendText),
        legend.title = element_text(size=legendTitle),
        strip.background = element_blank(),
        strip.text = element_text(size = axisText)) +
  xlab("Host group") +
  ylab("Disease severity") +
  coord_cartesian(ylim = c(0, 1))
figM

figA <- figM %+%
  fdatA
figA

# save
pdf("./output/figure2_damage_march.pdf", width = 3.15, height = 5.5)
figM
dev.off()

pdf("./output/figureS4_damage_april.pdf", width = 3.15, height = 5.5)
figA
dev.off()


#### percent change ####

# leaves
leaf_nat_m <- exp(fixef(imodM)$cond[1]) / (1 + exp(fixef(imodM)$cond[1]))
leaf_non_m <- exp(fixef(imodM)$cond[1] + fixef(imodM)$cond[2]) / (1 + exp(fixef(imodM)$cond[1] + fixef(imodM)$cond[2]))
(leaf_nat_m - leaf_non_m) / leaf_non_m

leaf_nat_a <- exp(fixef(imodA)$cond[1]) / (1 + exp(fixef(imodA)$cond[1]))
leaf_non_a <- exp(fixef(imodA)$cond[1] + fixef(imodA)$cond[2]) / (1 + exp(fixef(imodA)$cond[1] + fixef(imodA)$cond[2]))
(leaf_nat_a - leaf_non_a) / leaf_non_a

# surface
surf_nat_m <- exp(fixef(smodM2)$cond[1]) / (1 + exp(fixef(smodM2)$cond[1]))
surf_non_m <- exp(fixef(smodM2)$cond[1] + fixef(smodM2)$cond[2]) / (1 + exp(fixef(smodM2)$cond[1] + fixef(smodM2)$cond[2]))
(surf_nat_m - surf_non_m) / surf_non_m

surf_nat_a <- exp(fixef(smodA2)$cond[1]) / (1 + exp(fixef(smodA2)$cond[1]))
surf_non_a <- exp(fixef(smodA2)$cond[1] + fixef(smodA2)$cond[2]) / (1 + exp(fixef(smodA2)$cond[1] + fixef(smodA2)$cond[2]))
(surf_nat_a - surf_non_a) / surf_non_a


#### average change ####

# predicted 
scdatMpred <- pred_fun(ldatM, scmodM2, "surface scaled")
(0.0152 - 0.0112)/0.0152 # 0.26
ldatM %>%
  group_by(grass.group) %>%
  summarise(severity = mean(surface))
(0.0136 - 0.0113)/0.0136 # 0.17
scdatApred <- pred_fun(ldatA, scmodA2, "surface scaled")
(0.0317 - 0.0256)/0.0317 # 0.19
ldatA %>%
  group_by(grass.group) %>%
  summarise(severity = mean(surface))
(0.0202 - 0.0317)/0.0202 # -0.6
