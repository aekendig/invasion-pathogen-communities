## Goal: manuscrpit figures


#### set up ####

# clear all existing data
rm(list=ls())

# load libraries
library(tidyverse)
library(cowplot)

# plotting parameters
axisText = 9
axisTitle = 11
legendText = 8
legendTitle = 8

# import data
rdat <- read_csv("./output/richness_origin_data.csv")
sdat <- read_csv("./output/diversity_origin_data.csv")
tdat <- read_csv("./output/taxonomy_origin_data.csv")

idat15 <- read_csv("./output/infect_density_transect15_data.csv") %>%
  mutate(exp.type = "observational (2015)")
idat16T <- read_csv("./output/infect_density_transect16_data.csv") %>%
  mutate(exp.type = "observational (2016)")
idat16C <- read_csv("./output/infect_density_competition_data.csv") %>%
  mutate(exp.type = "manipulated (2016)")

plotdat <- read_csv("./output/plot_scale_density_data.csv")

mdat15 <- read_csv("./output/damage_density_meandam_transect15_data.csv") %>%
  mutate(exp.type = "observational (2015)")
mdat16T <- read_csv("./output/damage_density_meandam_transect16_data.csv") %>%
  mutate(exp.type = "observational (2016)")
mdat16C <- read_csv("./output/damage_density_meandam_competition_data.csv") %>%
  mutate(exp.type = "manipulated (2016)")

dat15leaf <- read_csv("./output/damage_density_propdam_transect15_data.csv") %>%
  mutate(exp.type = "observational (2015)")
dat16Tleaf <- read_csv("./output/damage_density_propdam_transect16_data.csv") %>%
  mutate(exp.type = "observational (2016)")
dat16Cleaf <- read_csv("./output/damage_density_propdam_competition_data.csv") %>%
  mutate(exp.type = "manipulated (2016)")

dat15plant <- read_csv("./output/damage_density_plant_transect15_data.csv") %>%
  mutate(exp.type = "observational (2015)")
dat16Tplant <- read_csv("./output/damage_density_plant_transect16_data.csv") %>%
  mutate(exp.type = "observational (2016)")
dat16Cplant <- read_csv("./output/damage_density_plant_competition_data.csv") %>%
  mutate(exp.type = "manipulated (2016)")


load("./output/infection_density_rpro_transect15_avg_model.rda")
load("./output/infection_density_rpro_transect16_avg_model.rda")
load("./output/infection_density_plol_transect15_avg_model.rda")
load("./output/infection_density_pcha_transect15_avg_model.rda")
load("./output/infection_density_ainf_transect15_avg_model.rda")
load("./output/infection_density_ainf_transect16_avg_model.rda")
load("./output/infection_density_ainf_competition_avg_model.rda")
load("./output/infection_density_dres_transect15_avg_model.rda")
load("./output/infection_density_dres_competition_avg_model.rda")

load("./output/damage_density_meandam_transect15_avg_model.rda")
load("./output/damage_density_meandam_transect16_avg_model.rda")
load("./output/damage_density_meandam_competition_avg_model.rda")

load("./output/damage_density_propdam_transect15_avg_model.rda")
load("./output/damage_density_propdam_transect16_avg_model.rda")
load("./output/damage_density_propdam_competition_avg_model.rda")

# functions
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

natpredfun <- function(dat, mod){
  out <- tibble(native.density = rep(dat$native.density, 2),
                natdens.s = rep(dat$natdens.s, 2),
                nonnative = rep(c(0, 1), each = nrow(dat)),
                origin = rep(c("native", "non-native"), each = nrow(dat)),
                nondens.s = mean(dat$nondens.s),
                experiment = unique(dat$experiment),
                year = unique(dat$year),
                exp.type = unique(dat$exp.type)) %>%
  mutate(pred = logit2prob(mod$coefficients["full", "cond((Int))"] + 
                                 mod$coefficients["full", "cond(nonnative)"] * nonnative + 
                                 mod$coefficients["full", "cond(natdens.s)"] * natdens.s + 
                                 mod$coefficients["full", "cond(nondens.s)"] * nondens.s + 
                                 mod$coefficients["full", "cond(nondens.s:nonnative)"] * nonnative * nondens.s +
                                 mod$coefficients["full", "cond(natdens.s:nonnative)"] * nonnative * natdens.s),
         density = native.density,
         pred.type = "native")
  
  return(out)
}

natpredfun2 <- function(dat, mod){
  out <- tibble(native.density = rep(dat$native.density, 2),
                natdens.s = rep(dat$natdens.s, 2),
                nonnative = rep(c(0, 1), each = nrow(dat)),
                origin = rep(c("native", "non-native"), each = nrow(dat)),
                nondens.s = mean(dat$nondens.s),
                experiment = unique(dat$experiment),
                year = unique(dat$year),
                exp.type = unique(dat$exp.type)) %>%
    mutate(pred = logit2prob(mod$coefficients["full", "cond((Int))"] + 
                               mod$coefficients["full", "cond(nonnative)"] * nonnative + 
                               mod$coefficients["full", "cond(natdens.s)"] * natdens.s + 
                               mod$coefficients["full", "cond(nondens.s)"] * nondens.s + 
                               mod$coefficients["full", "cond(natdens.s:nonnative)"] * nonnative * natdens.s),
           density = native.density,
           pred.type = "native")
  
  return(out)
}

natpredfun3 <- function(dat, mod){
  out <- tibble(native.density = rep(dat$native.density, 2),
                natdens.s = rep(dat$natdens.s, 2),
                nonnative = rep(c(0, 1), each = nrow(dat)),
                origin = rep(c("native", "non-native"), each = nrow(dat)),
                nondens.s = mean(dat$nondens.s),
                experiment = unique(dat$experiment),
                year = unique(dat$year),
                exp.type = unique(dat$exp.type)) %>%
    mutate(pred = logit2prob(mod$coefficients["full", "cond((Int))"] + 
                               mod$coefficients["full", "cond(nonnative)"] * nonnative + 
                               mod$coefficients["full", "cond(natdens.s)"] * natdens.s + 
                               mod$coefficients["full", "cond(nondens.s)"] * nondens.s + 
                               mod$coefficients["full", "cond(nondens.s:nonnative)"] * nonnative * nondens.s),
           density = native.density,
           pred.type = "native")
  
  return(out)
}

nonpredfun <- function(dat, mod){
  out <- tibble(nonnative.density = rep(dat$nonnative.density, 2),
                nondens.s = rep(dat$nondens.s, 2),
                nonnative = rep(c(0, 1), each = nrow(dat)),
                origin = rep(c("native", "non-native"), each = nrow(dat)),
                natdens.s = mean(dat$natdens.s),
                experiment = unique(dat$experiment),
                year = unique(dat$year),
                exp.type = unique(dat$exp.type)) %>%
    mutate(pred = logit2prob(mod$coefficients["full", "cond((Int))"] + 
                               mod$coefficients["full", "cond(nonnative)"] * nonnative + 
                               mod$coefficients["full", "cond(natdens.s)"] * natdens.s + 
                               mod$coefficients["full", "cond(nondens.s)"] * nondens.s + 
                               mod$coefficients["full", "cond(nondens.s:nonnative)"] * nonnative * nondens.s +
                               mod$coefficients["full", "cond(natdens.s:nonnative)"] * nonnative * natdens.s),
           density = nonnative.density,
           pred.type = "nonnative")
  
  return(out)
}

nonpredfun2 <- function(dat, mod){
  out <- tibble(nonnative.density = rep(dat$nonnative.density, 2),
                nondens.s = rep(dat$nondens.s, 2),
                nonnative = rep(c(0, 1), each = nrow(dat)),
                origin = rep(c("native", "non-native"), each = nrow(dat)),
                natdens.s = mean(dat$natdens.s),
                experiment = unique(dat$experiment),
                year = unique(dat$year),
                exp.type = unique(dat$exp.type)) %>%
    mutate(pred = logit2prob(mod$coefficients["full", "cond((Int))"] + 
                               mod$coefficients["full", "cond(nonnative)"] * nonnative + 
                               mod$coefficients["full", "cond(natdens.s)"] * natdens.s + 
                               mod$coefficients["full", "cond(nondens.s)"] * nondens.s + 
                               mod$coefficients["full", "cond(natdens.s:nonnative)"] * nonnative * natdens.s),
           density = nonnative.density,
           pred.type = "nonnative")
  
  return(out)
}

nonpredfun3 <- function(dat, mod){
  out <- tibble(nonnative.density = rep(dat$nonnative.density, 2),
                nondens.s = rep(dat$nondens.s, 2),
                nonnative = rep(c(0, 1), each = nrow(dat)),
                origin = rep(c("native", "non-native"), each = nrow(dat)),
                natdens.s = mean(dat$natdens.s),
                experiment = unique(dat$experiment),
                year = unique(dat$year),
                exp.type = unique(dat$exp.type)) %>%
    mutate(pred = logit2prob(mod$coefficients["full", "cond((Int))"] + 
                               mod$coefficients["full", "cond(nonnative)"] * nonnative + 
                               mod$coefficients["full", "cond(natdens.s)"] * natdens.s + 
                               mod$coefficients["full", "cond(nondens.s)"] * nondens.s + 
                               mod$coefficients["full", "cond(nondens.s:nonnative)"] * nonnative * nondens.s),
           density = nonnative.density,
           pred.type = "nonnative")
  
  return(out)
}


#### model predictions ####

# simplified datasets
inf15 <- idat15 %>%
  select(experiment, year, exp.type, subplot, native.density, natdens.s, nonnative.density, nondens.s) %>%
  unique()

inf16T <- idat16T %>%
  select(experiment, year, exp.type, subplot, native.density, natdens.s, nonnative.density, nondens.s) %>%
  unique()

inf16C <- idat16C %>%
  select(experiment, year, exp.type, subplot, native.density, natdens.s, nonnative.density, nondens.s) %>%
  unique()

propdam15 <- dat15leaf %>%
  select(experiment, year, exp.type, subplot, native.density, natdens.s, nonnative.density, nondens.s) %>%
  unique()

propdam16T <- dat16Tleaf %>%
  select(experiment, year, exp.type, subplot, native.density, natdens.s, nonnative.density, nondens.s) %>%
  unique()

propdam16C <- dat16Cleaf %>%
  select(experiment, year, exp.type, subplot, native.density, natdens.s, nonnative.density, nondens.s) %>%
  unique()

meandam15 <- mdat15 %>%
  select(experiment, year, exp.type, subplot, native.density, natdens.s, nonnative.density, nondens.s) %>%
  unique()

meandam16T <- mdat16T %>%
  select(experiment, year, exp.type, subplot, native.density, natdens.s, nonnative.density, nondens.s) %>%
  unique()

meandam16C <- mdat16C %>%
  select(experiment, year, exp.type, subplot, native.density, natdens.s, nonnative.density, nondens.s) %>%
  unique()

# predictions
# rpnatpred15 <- natpredfun(dat = inf15, mod = rpmod15a) # no nondens.s:nonnative
rpnatpred15 <- natpredfun2(dat = inf15, mod = rpmod15a) %>% mutate(taxonomy = "Ramularia cf. proteae")
# rpnatpred16T <- natpredfun(dat = inf16T, mod = rpmod16Ta) # no natdens.s:nonnative
rpnatpred16T <- natpredfun3(dat = inf16T, mod = rpmod16Ta) %>% mutate(taxonomy = "Ramularia cf. proteae")
# plnatpred15 <- natpredfun(dat = inf15, mod = plmod15a) # no nondens.s:nonnative
plnatpred15 <- natpredfun2(dat = inf15, mod = plmod15a) %>% mutate(taxonomy = "Pyrenophora cf. lolii")
pcnatpred15 <- natpredfun(dat = inf15, mod = pcmod15a) %>% mutate(taxonomy = "Pyrenophora cf. chaetomioides")
ainatpred15 <- natpredfun(dat = inf15, mod = aimod15a) %>% mutate(taxonomy = "Alternaria cf. infectoria")
# ainatpred16T <- natpredfun(dat = inf16T, mod = aimod16Ta) # no nondens.s:nonnative
ainatpred16T <- natpredfun2(dat = inf16T, mod = aimod16Ta) %>% mutate(taxonomy = "Alternaria cf. infectoria")
ainatpred16C <- natpredfun(dat = inf16C, mod = aimod16Ca) %>% mutate(taxonomy = "Alternaria cf. infectoria")
drnatpred15 <- natpredfun(dat = inf15, mod = drmod15a) %>% mutate(taxonomy = "Drechslera sp.")
drnatpred16C <- natpredfun(dat = inf16C, mod = drmod16Ca) %>% mutate(taxonomy = "Drechslera sp.")

rpnonpred15 <- nonpredfun2(dat = inf15, mod = rpmod15a) %>% mutate(taxonomy = "Ramularia cf. proteae")
rpnonpred16T <- nonpredfun3(dat = inf16T, mod = rpmod16Ta) %>% mutate(taxonomy = "Ramularia cf. proteae")
plnonpred15 <- nonpredfun2(dat = inf15, mod = plmod15a) %>% mutate(taxonomy = "Pyrenophora cf. lolii")
pcnonpred15 <- nonpredfun(dat = inf15, mod = pcmod15a) %>% mutate(taxonomy = "Pyrenophora cf. chaetomioides")
ainonpred15 <- nonpredfun(dat = inf15, mod = aimod15a) %>% mutate(taxonomy = "Alternaria cf. infectoria")
ainonpred16T <- nonpredfun2(dat = inf16T, mod = aimod16Ta) %>% mutate(taxonomy = "Alternaria cf. infectoria")
ainonpred16C <- nonpredfun(dat = inf16C, mod = aimod16Ca) %>% mutate(taxonomy = "Alternaria cf. infectoria")
drnonpred15 <- nonpredfun(dat = inf15, mod = drmod15a) %>% mutate(taxonomy = "Drechslera sp.")
drnonpred16C <- nonpredfun(dat = inf16C, mod = drmod16Ca) %>% mutate(taxonomy = "Drechslera sp.")

natpropdampred15 <- natpredfun(dat = propdam15, mod = pmod15a)
natpropdampred16T <- natpredfun(dat = propdam16T, mod = pmod16Ta)
natpropdampred16C <- natpredfun(dat = propdam16C, mod = pmod16Ca)

nonpropdampred15 <- nonpredfun(dat = propdam15, mod = pmod15a)
nonpropdampred16T <- nonpredfun(dat = propdam16T, mod = pmod16Ta)
nonpropdampred16C <- nonpredfun(dat = propdam16C, mod = pmod16Ca)

# natmeandampred15  <- natpredfun(dat = meandam15, mod = mmod15a) # no natdens.s:nonnative
natmeandampred15  <- natpredfun3(dat = meandam15, mod = mmod15a) 
natmeandampred16T <- natpredfun(dat = meandam16T, mod = mmod16Ta)
natmeandampred16C <- natpredfun(dat = meandam16C, mod = mmod16Ca)

nonmeandampred15  <- nonpredfun3(dat = meandam15, mod = mmod15a)
nonmeandampred16T <- nonpredfun(dat = meandam16T, mod = mmod16Ta)
nonmeandampred16C <- nonpredfun(dat = meandam16C, mod = mmod16Ca)

# below doesn't work because re.form isn't implemented yet
# mutate(pred = predict(model.avg(get.models(pmod16Cd, subset = cumsum(weight) <= .95)), type = "link", backtransform = T, re.form = NA, newdata = .)) 

# combine data
plantdat <- dat15plant %>%
  mutate(pred.type = "native",
         density = native.density) %>%
  full_join(dat15plant %>%
              mutate(pred.type = "nonnative",
                     density = nonnative.density)) %>%
  full_join(dat16Tplant %>%
              mutate(pred.type = "native",
                     density = native.density)) %>%
  full_join(dat16Tplant %>%
              mutate(pred.type = "nonnative",
                     density = nonnative.density)) %>%
  full_join(dat16Cplant %>%
              mutate(pred.type = "native",
                     density = native.density)) %>%
  full_join(dat16Cplant %>%
              mutate(pred.type = "nonnative",
                     density = nonnative.density)) %>%
  mutate(exp.type = factor(exp.type, levels = c("observational (2015)", "observational (2016)", "manipulated (2016)")))

meandat <- mdat15 %>%
  mutate(pred.type = "native",
         density = native.density) %>%
  full_join(mdat15 %>%
              mutate(pred.type = "nonnative",
                     density = nonnative.density)) %>%
  full_join(mdat16T %>%
              mutate(pred.type = "native",
                     density = native.density)) %>%
  full_join(mdat16T %>%
              mutate(pred.type = "nonnative",
                     density = nonnative.density)) %>%
  full_join(mdat16C %>%
              mutate(pred.type = "native",
                     density = native.density)) %>%
  full_join(mdat16C %>%
              mutate(pred.type = "nonnative",
                     density = nonnative.density)) %>%
  mutate(exp.type = factor(exp.type, levels = c("observational (2015)", "observational (2016)", "manipulated (2016)")))

# combine prediction data
infnatpred <- full_join(rpnatpred15, rpnatpred16T) %>%
  full_join(plnatpred15) %>%
  full_join(pcnatpred15) %>%
  full_join(ainatpred15) %>%
  full_join(ainatpred16T) %>%
  full_join(ainatpred16C) %>%
  full_join(drnatpred15) %>%
  full_join(drnatpred16C) %>%
  mutate(tax.exp = paste(taxonomy, exp.type, sep = "\n") %>%
           factor(levels = c("Alternaria cf. infectoria\nobservational (2015)", "Alternaria cf. infectoria\nobservational (2016)", "Alternaria cf. infectoria\nmanipulated (2016)", "Drechslera sp.\nobservational (2015)", "Drechslera sp.\nmanipulated (2016)", "Ramularia cf. proteae\nobservational (2015)", "Ramularia cf. proteae\nobservational (2016)", "Pyrenophora cf. chaetomioides\nobservational (2015)", "Pyrenophora cf. lolii\nobservational (2015)")))

infnonpred <- full_join(rpnonpred15, rpnonpred16T) %>%
  full_join(plnonpred15) %>%
  full_join(pcnonpred15) %>%
  full_join(ainonpred15) %>%
  full_join(ainonpred16T) %>%
  full_join(ainonpred16C) %>%
  full_join(drnonpred15) %>%
  full_join(drnonpred16C) %>%
  mutate(tax.exp = paste(taxonomy, exp.type, sep = "\n") %>%
           factor(levels = c("Alternaria cf. infectoria\nobservational (2015)", "Alternaria cf. infectoria\nobservational (2016)", "Alternaria cf. infectoria\nmanipulated (2016)", "Drechslera sp.\nobservational (2015)", "Drechslera sp.\nmanipulated (2016)", "Ramularia cf. proteae\nobservational (2015)", "Ramularia cf. proteae\nobservational (2016)", "Pyrenophora cf. chaetomioides\nobservational (2015)", "Pyrenophora cf. lolii\nobservational (2015)")))

propdampred <- full_join(natpropdampred15, natpropdampred16T) %>%
  full_join(natpropdampred16C) %>%
  full_join(nonpropdampred15) %>%
  full_join(nonpropdampred16T) %>%
  full_join(nonpropdampred16C) %>%
  mutate(exp.type = factor(exp.type, levels = c("observational (2015)", "observational (2016)", "manipulated (2016)")))

meandampred <- full_join(natmeandampred15, natmeandampred16T) %>%
  full_join(natmeandampred16C) %>%
  full_join(nonmeandampred15) %>%
  full_join(nonmeandampred16T) %>%
  full_join(nonmeandampred16C) %>%
  mutate(exp.type = factor(exp.type, levels = c("observational (2015)", "observational (2016)", "manipulated (2016)")))

# combine infection data
idat <- idat15 %>%
  select(year, experiment, exp.type, plot, subplot, native.density, nonnative.density, origin, rpro, plol, pcha, ainf, dres) %>%
  gather(key = "taxonomy", value = "infection", -(year:origin)) %>%
  mutate(taxonomy = recode(taxonomy, rpro = "Ramularia cf. proteae", plol = "Pyrenophora cf. lolii", pcha = "Pyrenophora cf. chaetomioides", ainf = "Alternaria cf. infectoria", dres = "Drechslera sp.")) %>%
  full_join(select(idat16T, year, experiment, exp.type, plot, subplot, native.density, nonnative.density, origin, rpro, plol, pcha, ainf, dres) %>%
              gather(key = "taxonomy", value = "infection", -(year:origin)) %>%
              mutate(taxonomy = recode(taxonomy, rpro = "Ramularia cf. proteae", plol = "Pyrenophora cf. lolii", pcha = "Pyrenophora cf. chaetomioides", ainf = "Alternaria cf. infectoria", dres = "Drechslera sp."))) %>%
  full_join(select(idat16C, year, experiment, exp.type, plot, subplot, native.density, nonnative.density, origin, rpro, plol, pcha, ainf, dres) %>%
              gather(key = "taxonomy", value = "infection", -(year:origin)) %>%
              mutate(taxonomy = recode(taxonomy, rpro = "Ramularia cf. proteae", plol = "Pyrenophora cf. lolii", pcha = "Pyrenophora cf. chaetomioides", ainf = "Alternaria cf. infectoria", dres = "Drechslera sp."))) %>%
  mutate(tax.exp = paste(taxonomy, exp.type, sep = "\n")) %>%
  filter(tax.exp %in% infnatpred$tax.exp) %>%
  mutate(tax.exp = factor(tax.exp, levels = c("Alternaria cf. infectoria\nobservational (2015)", "Alternaria cf. infectoria\nobservational (2016)", "Alternaria cf. infectoria\nmanipulated (2016)", "Drechslera sp.\nobservational (2015)", "Drechslera sp.\nmanipulated (2016)", "Ramularia cf. proteae\nobservational (2015)", "Ramularia cf. proteae\nobservational (2016)", "Pyrenophora cf. chaetomioides\nobservational (2015)", "Pyrenophora cf. lolii\nobservational (2015)")))


#### communities by origin ####

richfig <- rdat %>%
  ggplot(aes(x = as.factor(year), y = richness, group = grass.status, color = grass.status, shape = type)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, position = position_dodge(0.5), color = "black") +
  geom_point(aes(y = richness), size = 2, position = position_dodge(0.5)) +
  scale_shape_manual(values = c(19, 22)) +
  scale_y_continuous(trans = "log10") +
  ylab(expression(paste(Log[10], "(richness)", sep = ""))) +
  theme_bw() +
  theme(axis.text = element_text(size = axisText, color="black"),
        axis.title.y = element_text(size = axisTitle),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size=legendText),
        legend.title = element_text(size=legendTitle),
        strip.background = element_blank(),
        strip.text = element_text(size = axisText))

divfig <- sdat %>%
  ggplot(aes(x = as.factor(year), y = diversity, group = grass.status, color = grass.status)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, position = position_dodge(0.5), color = "black") +
  geom_point(size = 2, position = position_dodge(0.5), shape = 19) +
  ylab("Shannon diversity") +
  theme_bw() +
  theme(axis.text = element_text(size = axisText, color="black"),
        axis.title.y = element_text(size = axisTitle),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size=legendText),
        legend.title = element_text(size=legendTitle),
        strip.background = element_blank(),
        strip.text = element_text(size = axisText),
        legend.position = "none")

taxfig <- tdat %>%
  ggplot(aes(x = rank, y = abundance, color = grass.status)) +
  geom_point() +
  geom_line() +
  geom_text(aes(label = taxonomy), size = 2, hjust = -0.1) +
  facet_wrap(~as.factor(year), strip.position = "bottom", scales = "free_y") +
  ylab("Isolate abundance") +
  theme_bw() +
  theme(axis.text.y = element_text(size = axisText, color="black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = axisTitle),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size=legendText),
        legend.title = element_text(size=legendTitle),
        strip.background = element_blank(),
        strip.text = element_text(size = axisText),
        strip.placement = "outside",
        legend.position = "none")

richleg <- get_legend(richfig + theme(legend.direction = "horizontal", legend.position = "bottom"))

comfig <- plot_grid(richfig + theme(legend.position = "none"), divfig, 
                    nrow = 1, rel_widths = c(1, 0.95))

pdf("./output/origin_community_figure.pdf", width = 6, height = 6)
plot_grid(comfig, taxfig, richleg,
          nrow = 3, rel_heights = c(1, 1, 0.1), scale = 0.95)
dev.off()


#### taxonomy table ####

taxtab <- tdat %>%
  mutate(grass.status = recode(grass.status, "non-native" = "nonnative")) %>%
  gather(variable, value, -(year:taxonomy)) %>%
  unite(temp, grass.status, variable) %>%
  spread(temp, value) %>%
  filter(native_rank < 3 | nonnative_rank < 3)

write_csv(taxtab, "./output/taxonomy_table.csv")


#### native density infection ####

inatfig <- idat %>%
  ggplot(aes(x = native.density, y = infection, color = origin)) +
  geom_point(alpha = 0.7) +
  geom_line(data = infnatpred, aes(y = pred)) +
  facet_wrap(~ tax.exp, nrow = 3, scales = "free_x") +
  xlab(expression(paste("Native plant density (", m^-1, ")", sep = ""))) +
  ylab("Infection prevalence") +
  theme_bw() +
  theme(axis.text = element_text(size = axisText, color="black"),
        axis.title = element_text(size = axisTitle),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size=legendText),
        legend.title = element_text(size=legendTitle),
        strip.background = element_blank(),
        strip.text = element_text(size = axisText),
        legend.position = "bottom", 
        legend.direction = "horizontal")
  
pdf("./output/infection_native_density_figure.pdf", width = 6, height = 6)
inatfig
dev.off()


#### non-native density infection ####

inonfig <- idat %>%
  ggplot(aes(x = nonnative.density, y = infection, color = origin)) +
  geom_point(alpha = 0.7) +
  geom_line(data = infnonpred, aes(y = pred)) +
  facet_wrap(~ tax.exp, nrow = 3, scales = "free_x") +
  xlab(expression(paste("Non-native plant density (", m^-1, ")", sep = ""))) +
  ylab("Infection prevalence") +
  theme_bw() +
  theme(axis.text = element_text(size = axisText, color="black"),
        axis.title = element_text(size = axisTitle),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size=legendText),
        legend.title = element_text(size=legendTitle),
        strip.background = element_blank(),
        strip.text = element_text(size = axisText),
        legend.position = "bottom", 
        legend.direction = "horizontal")

pdf("./output/infection_nonnative_density_figure.pdf", width = 6, height = 6)
inonfig
dev.off()


#### density correlation ####

corfig <- plotdat %>%
  mutate(exp.type = case_when(experiment == "competition" ~ "manipulated (2016)",
                              experiment == "transect" & year == 2015 ~ "observational (2015)",
                              experiment == "transect" & year == 2016 ~ "observational (2016)") %>% factor(levels = c("observational (2015)", "observational (2016)", "manipulated (2016)"))) %>%
  ggplot(aes(x = natdens.s, y = nondens.s)) +
  geom_point(alpha = 0.7) + 
  facet_wrap(~ exp.type, scales = "free") +
  xlab(expression(paste("Scaled native density (", m^-1, ")", sep = ""))) +
  ylab(expression(paste("Scaled non-native density (", m^-1, ")", sep = ""))) +
  theme_bw() +
  theme(axis.text = element_text(size = axisText, color="black"),
        axis.title = element_text(size = axisTitle),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size=legendText),
        legend.title = element_text(size=legendTitle),
        strip.background = element_blank(),
        strip.text = element_text(size = axisText))

pdf("./output/density_correlation_figure.pdf", width = 6, height = 2.5)
corfig
dev.off()
  

#### prop leaves damaged by origin ####

propfignat <- plantdat %>%
  filter(pred.type == "native") %>%
  ggplot(aes(x = density, y = prop.dam, color = origin)) +
  geom_point(alpha = 0.7) + 
  geom_line(data = filter(propdampred, pred.type == "native"), aes(y = pred)) +
  xlab(expression(paste("Native plant density (", m^-1, ")", sep = ""))) +
  facet_wrap(~exp.type, scales = "free_x") +
  theme_bw() +
  theme(axis.text = element_text(size = axisText, color="black"),
        axis.title.x = element_text(size = axisTitle),
        axis.title.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size=legendText),
        legend.title = element_text(size=legendTitle),
        strip.background = element_blank(),
        strip.text = element_text(size = axisTitle))

propfignon <- plantdat %>%
  filter(pred.type == "nonnative") %>%
  ggplot(aes(x = density, y = prop.dam, color = origin)) +
  geom_point(alpha = 0.7) + 
  geom_line(data = filter(propdampred, pred.type == "nonnative"), aes(y = pred)) +
  xlab(expression(paste("Non-native plant density (", m^-1, ")", sep = ""))) +
  facet_wrap(~exp.type, scales = "free_x") +
  theme_bw() +
  theme(axis.text = element_text(size = axisText, color="black"),
        axis.title.x = element_text(size = axisTitle),
        axis.title.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size=legendText),
        legend.title = element_text(size=legendTitle),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")

propleg <- get_legend(propfignat + theme(legend.direction = "horizontal", legend.position = "bottom"))

pdf("./output/density_propdam_figure.pdf", width = 6, height = 5)
plot_grid(propfignat + theme(legend.position = "none"), propfignon, propleg, 
          nrow = 3, rel_heights = c(1, 1, 0.1), scale = 0.95) +
  draw_label("Proportion of leaves damaged", x =  0, y = 0.5, angle = 90, vjust = 1.2, size = axisTitle)
dev.off()


#### mean damage by origin ####
meanfignat <- meandat %>%
  filter(pred.type == "native") %>%
  ggplot(aes(x = density, y = mean.dam, color = origin)) +
  geom_point(alpha = 0.7) + 
  geom_line(data = filter(meandampred, pred.type == "native"), aes(y = pred)) +
  xlab(expression(paste("Native plant density (", m^-1, ")", sep = ""))) +
  facet_wrap(~exp.type, scales = "free") +
  theme_bw() +
  theme(axis.text = element_text(size = axisText, color="black"),
        axis.title.x = element_text(size = axisTitle),
        axis.title.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size=legendText),
        legend.title = element_text(size=legendTitle),
        strip.background = element_blank(),
        strip.text = element_text(size = axisTitle))

meanfignon <- meandat %>%
  filter(pred.type == "nonnative") %>%
  ggplot(aes(x = density, y = mean.dam, color = origin)) +
  geom_point(alpha = 0.7) + 
  geom_line(data = filter(meandampred, pred.type == "nonnative"), aes(y = pred)) +
  xlab(expression(paste("Non-native plant density (", m^-1, ")", sep = ""))) +
  facet_wrap(~exp.type, scales = "free") +
  theme_bw() +
  theme(axis.text = element_text(size = axisText, color="black"),
        axis.title.x = element_text(size = axisTitle),
        axis.title.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size=legendText),
        legend.title = element_text(size=legendTitle),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")

meanleg <- get_legend(meanfignat + theme(legend.direction = "horizontal", legend.position = "bottom"))

pdf("./output/density_meandam_figure.pdf", width = 6, height = 5)
plot_grid(meanfignat + theme(legend.position = "none"), meanfignon, meanleg, 
          nrow = 3, rel_heights = c(1, 1, 0.1), scale = 0.95) +
  draw_label("Severity of damaged leaves", x =  0, y = 0.5, angle = 90, vjust = 1.2, size = axisTitle)
dev.off()
