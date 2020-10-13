#### set up ####

# clear all existing data
rm(list=ls())

# load libraries
library(sjPlot)
library(tidyverse)
library(cowplot)
library(grid)
library(gridExtra)
library(lemon)
library(MuMIn)

# figure settings
col_pal <- c("#000000", "#56B4E9", "#CC79A7", "#009E73", "#D55E00", "#0072B2", "#F0E442")

axisText = 10
axisTitle = 12
legendText = 10
legendTitle = 10

ggtheme <- theme_bw() +
  theme(axis.text = element_text(size = axisText, color="black"),
        axis.title = element_text(size = axisTitle),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size=legendText),
        legend.title = element_text(size=legendTitle),
        strip.background = element_blank(),
        strip.text = element_text(size = axisText)) 

# import data
idatT <- read_csv("./output/infect_host_density_transect_data.csv")
idatC <- read_csv("./output/infect_host_density_competition_data.csv")

pdatTm <- read_csv("./output/damage_host_density_plant_march_transect_data.csv")
pdatTa <- read_csv("./output/damage_host_density_plant_april_transect_data.csv")
pdatC <- read_csv("./output/damage_host_density_plant_competition_data.csv")

ldatTm <- read_csv("./output/damage_host_density_leaf_march_transect_data.csv")
ldatTa <- read_csv("./output/damage_host_density_leaf_april_transect_data.csv")
ldatC <- read_csv("./output/damage_host_density_leaf_competition_data.csv")

fun <- read_csv("./output/fungal_taxonomy_rank.csv")

# import models
load("./output/infection_density_experiment_rpro_absolute_transect_avg_model.rda")
load("./output/infection_density_experiment_ainf_absolute_transect_avg_model.rda")
load("./output/infection_density_experiment_pcha_absolute_transect_avg_model.rda")
load("./output/infection_density_experiment_plol_absolute_transect_avg_model.rda")
load("./output/infection_density_experiment_ptri_absolute_transect_avg_model.rda")
load("./output/infection_density_experiment_dres_absolute_transect_avg_model.rda")
load("./output/infection_density_experiment_pave_absolute_transect_avg_model.rda")

load("./output/infection_density_experiment_ainf_absolute_competition_avg_model.rda")
load("./output/infection_density_experiment_pcha_absolute_competition_avg_model.rda")
load("./output/infection_density_experiment_plol_absolute_competition_avg_model.rda")
load("./output/infection_density_experiment_ptri_absolute_competition_avg_model.rda")
load("./output/infection_density_experiment_dres_absolute_competition_avg_model.rda")
load("./output/infection_density_experiment_pave_absolute_competition_avg_model.rda")

load("./output/damage_host_density_infection_march_transect_model.rda")
load("./output/damage_host_density_infection_april_transect_model.rda")
load("./output/damage_host_density_infection_competition_model.rda")
load("./output/damage_host_density_surface_march_transect_model.rda")
load("./output/damage_host_density_surface_april_transect_model.rda")
load("./output/damage_host_density_surface_competition_model.rda")
load("./output/damage_host_density_surface_scaled_march_transect_model.rda")
load("./output/damage_host_density_surface_scaled_april_transect_model.rda")
load("./output/damage_host_density_surface_scaled_competition_model.rda")


#### table of pathogens ####

# combine possible host species across years
fun2 <- fun %>%
  select(pathogen, otu.id, host.num, amb.host.num) %>%
  rename(taxonomy = pathogen) %>%
  unique()

# otu's and code names
otus <- tibble(pathogen = c("ainf", "dres", "pave", "pcha", "plol", "ptri", "rpro"),
               otu.id = c(1, 2, 7, 4, 5, 8, 3),
               path.abb = c("A. inf.", "Pyr. sp.", "P. ave.", "P. cha.", "P. lol.", "P. tri.", "R. pro.")) %>%
  left_join(fun2)

# edit data
prevdat <- idatT %>%
  select(experiment, grass.group, ainf, pcha, plol, ptri, dres, pave, rpro) %>%
  full_join(idatC %>% select(experiment, grass.group, ainf, pcha, plol, ptri, dres, pave, rpro)) %>%
  gather(key = pathogen, value = present, -c(experiment, grass.group)) %>%
  left_join(otus) %>%
  mutate(exp.type = case_when(experiment == "transect" ~ "Observational",
                              experiment == "competition" ~ "Manipulated") %>% factor(levels = c("Observational", "Manipulated"))) 

# make table
fun_tab <- prevdat %>%
  group_by(taxonomy, host.num, amb.host.num, grass.group) %>%
  summarise(otus = sum(present),
            prop = round(otus/length(present), 2)) %>%
  mutate(out = paste(otus, prop, sep = " (")) %>%
  select(-c(otus, prop)) %>%
  spread(key = grass.group, value = out) 

tab_df(fun_tab)


#### sample sizes ####

# transect isolates
idatT %>%
  group_by(year, grass.group) %>%
  summarise(n = length(unique(isolate.id)))

# competition isolates
idatC %>%
  group_by(grass.group) %>%
  summarise(n = length(unique(isolate.id)))

# transect hosts isolates
unique(idatT$host)

# competition hosts isoaltes
unique(idatC$host)

# transect severity
pdatTm %>%
  group_by(year, grass.group) %>%
  summarise(n = length(mean.dam))

pdatTa %>%
  group_by(year, grass.group) %>%
  summarise(n = length(mean.dam))

# competition severity
pdatC %>%
  group_by(grass.group) %>%
  summarise(n = length(mean.dam))

# transect hosts severity
unique(pdatTm$host)
unique(pdatTa$host)
unique(pdatC$host)

# competition plots
idatC %>%
  group_by(bg.species, competition.density, competition.type) %>%
  summarise(n = length(unique(subplot)))

pdatC %>%
  group_by(bg.species, competition.density) %>%
  summarise(n = length(unique(subplot))) %>%
  data.frame()

# transect design
idatT %>%
  select(plot, subplot) %>%
  full_join(pdatTm %>% select(plot, subplot)) %>%
  full_join(pdatTa %>% select(plot, subplot)) %>%
  group_by(plot) %>%
  summarise(n = length(unique(subplot)))

# transect plots
idatT %>%
  select(year, subplot) %>%
  full_join(pdatTm %>% select(year, subplot)) %>%
  full_join(pdatTa %>% select(year, subplot)) %>%
  group_by(year) %>%
  summarise(n = length(unique(subplot)))

idatT %>%
  group_by(year) %>%
  summarise(n = length(unique(subplot)))

pdatTm %>%
  group_by(year) %>%
  summarise(n = length(unique(subplot)))

pdatTa %>%
  group_by(year) %>%
  summarise(n = length(unique(subplot)))

# pathogens in group vs. rarer pathogens
idatT %>%
  select(grass.status, otu.id) %>%
  rbind(idatC %>% select(grass.status, otu.id)) %>%
  mutate(common = case_when(otu.id %in% c(1, 4, 5, 8, 2, 7, 3) ~ 1,
                            TRUE ~ 0)) %>%
  group_by(grass.status) %>%
  summarise(com = sum(common),
            tot = length(common)) %>%
  mutate(prop = com/tot)


#### individual density effects on pathogens ####

# standard deviation of each density
(natsdT <- sd(idatT$native.density))
(nonsdT <- sd(idatT$nonnative.density))
(natsdC <- sd(idatC$native.density))
(nonsdC <- sd(idatC$nonnative.density))

# mean of each density
(natmuT <- mean(idatT$native.density))
(nonmuT <- mean(idatT$nonnative.density))
(natmuC <- mean(idatC$native.density))
(nonmuC <- mean(idatC$nonnative.density))

# number of perennial individuals to increase by
ind_p <- 50

# number of annual individuals
ind_a <- 5000

# density table
inddens <- tibble(natdens.s = c(-natmuT/natsdT, (ind_p-natmuT)/natsdT, -natmuT/natsdT, -natmuC/natsdC, (ind_p-natmuC)/natsdC, -natmuC/natsdC),
                nondens.s = c(-nonmuT/nonsdT, -nonmuT/nonsdT,(ind_a-nonmuT)/nonsdT, -nonmuC/nonsdC, -nonmuC/nonsdC, (ind_a-nonmuC)/nonsdC),
                native.density = c(0, ind_p, 0, 0, ind_p, 0),
                nonnative.density = c(0, 0, ind_a, 0, 0, ind_a),
                exp.type = c(rep("Observational", 3), rep("Manipulated", 3)) %>% 
                  factor(levels = c("Observational", "Manipulated")),
                increase = rep(c("none", "nat", "non"), 2)) %>%
  mutate(othdens.s = 0,
         year = NA,
         subplot = NA) %>%
  merge(tibble(nonnative = c(0, 1)), all = T) %>%
  as_tibble()
inddens

# check values
filter(idatC, native.density == min(native.density)) %>% select(native.density, natdens.s) %>% unique()
filter(idatT, native.density == min(native.density)) %>% select(native.density, natdens.s) %>% unique()
filter(idatC, nonnative.density == min(nonnative.density)) %>% select(nonnative.density, nondens.s) %>% unique()
filter(idatT, nonnative.density == min(nonnative.density)) %>% select(nonnative.density, nondens.s) %>% unique() 
# seem similar

# model prediction function
fmod_pred_fun <- function(dat) {
  asim <- dat %>%
    merge(otus, all = T) %>%
    mutate(prev = case_when(pathogen == "rpro" & exp.type == "Observational" ~ 
                              predict(rpamodTa, newdata = ., re.form = NA, type = "response"),
                            pathogen == "ainf" & exp.type == "Observational" ~ 
                              predict(aiamodTa, newdata = ., re.form = NA, type = "response"),
                            pathogen == "pcha" & exp.type == "Observational" ~ 
                              predict(pcamodTa, newdata = ., re.form = NA, type = "response"),
                            pathogen == "plol" & exp.type == "Observational" ~ 
                              predict(plamodTa, newdata = ., re.form = NA, type = "response"),
                            pathogen == "ptri" & exp.type == "Observational" ~ 
                              predict(ptamodTa, newdata = ., re.form = NA, type = "response"),
                            pathogen == "dres" & exp.type == "Observational" ~ 
                              predict(dramodTa, newdata = ., re.form = NA, type = "response"),
                            pathogen == "pave" & exp.type == "Observational" ~ 
                              predict(paamodTa, newdata = ., re.form = NA, type = "response"),
                            pathogen == "ainf" & exp.type == "Manipulated" ~ 
                              predict(aiamodCa, newdata = ., re.form = NA, type = "response"),
                            pathogen == "pcha" & exp.type == "Manipulated" ~ 
                              predict(pcamodCa, newdata = ., re.form = NA, type = "response"),
                            pathogen == "plol" & exp.type == "Manipulated" ~ 
                              predict(plamodCa, newdata = ., re.form = NA, type = "response"),
                            pathogen == "ptri" & exp.type == "Manipulated" ~ 
                              predict(ptamodCa, newdata = ., re.form = NA, type = "response"),
                            pathogen == "dres" & exp.type == "Manipulated" ~ 
                              predict(dramodCa, newdata = ., re.form = NA, type = "response"),
                            pathogen == "pave" & exp.type == "Manipulated" ~ 
                              predict(paamodCa, newdata = ., re.form = NA, type = "response"),
                            TRUE ~ NA_real_),
           prev.se = case_when(pathogen == "rpro" & exp.type == "Observational" ~ 
                                 predict(rpamodTa, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit,
                               pathogen == "ainf" & exp.type == "Observational" ~ 
                                 predict(aiamodTa, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit,
                               pathogen == "pcha" & exp.type == "Observational" ~ 
                                 predict(pcamodTa, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit,
                               pathogen == "plol" & exp.type == "Observational" ~ 
                                 predict(plamodTa, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit,
                               pathogen == "ptri" & exp.type == "Observational" ~ 
                                 predict(ptamodTa, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit,
                               pathogen == "dres" & exp.type == "Observational" ~ 
                                 predict(dramodTa, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit,
                               pathogen == "pave" & exp.type == "Observational" ~ 
                                 predict(paamodTa, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit,
                               pathogen == "ainf" & exp.type == "Manipulated" ~ 
                                 predict(aiamodCa, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit,
                               pathogen == "pcha" & exp.type == "Manipulated" ~ 
                                 predict(pcamodCa, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit,
                               pathogen == "plol" & exp.type == "Manipulated" ~ 
                                 predict(plamodCa, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit,
                               pathogen == "ptri" & exp.type == "Manipulated" ~ 
                                 predict(ptamodCa, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit,
                               pathogen == "dres" & exp.type == "Manipulated" ~ 
                                 predict(dramodCa, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit,
                               pathogen == "pave" & exp.type == "Manipulated" ~ 
                                 predict(paamodCa, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit,
                               TRUE ~ NA_real_)) %>%
    as_tibble() %>%
    filter(!(pathogen == "pcha" & nonnative == 0) & !(pathogen == "ptri" & nonnative == 1) & !(pathogen == "rpro" & exp.type == "Manipulated")) %>%
    mutate(grass.group = case_when(nonnative == 0 ~ "native\nperennial",
                                   nonnative == 1 ~ "non-native\nannual"))
  
  return(asim)
}

# predict prevalence values
aind <- fmod_pred_fun(inddens) 

# select zero values
aind0 <- filter(aind, increase == "none") %>%
  rename(prev0 = prev, prev0.se = prev.se) %>%
  select(exp.type, pathogen, path.abb, grass.group, prev0, prev0.se)

# split by density and calculate change in prevalence, remove specialists
aindnat <- filter(aind, increase == "nat") %>% 
  select(exp.type, pathogen, path.abb, grass.group, prev, prev.se) %>%
  full_join(aind0) %>%
  mutate(prev.change = prev - prev0,
         prev.change.se = sqrt(prev.se^2 + prev0.se^2))

aindnon <- filter(aind, increase == "non") %>% 
  select(exp.type, pathogen, path.abb, grass.group, prev, prev.se) %>%
  full_join(aind0) %>%
  mutate(prev.change = prev - prev0,
         prev.change.se = sqrt(prev.se^2 + prev0.se^2))

# labels
indlabels <- tibble(dens.type = rep(c("native perennial", "non-native annual"), each = 2), exp.type = c("Observational", "Manipulated", "Observational", "Manipulated"), x = rep(0.5, 4), y = rep(1, 4), labels = c("A", "B", "C", "D")) %>%
  mutate(exp.type = factor(exp.type, levels = c("Observational", "Manipulated")))

# figure
changeplotnat <- ggplot(aindnat, aes(x = grass.group, y = prev.change)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_errorbar(aes(ymin = prev.change - prev.change.se, ymax = prev.change + prev.change.se, shape = path.abb), color = "black", position = position_dodge(width = 0.6), width = 0, show.legend = F) +
  geom_point(aes(fill = path.abb, shape = path.abb), position = position_dodge(width = 0.6), size = 3) +
  geom_text(data = filter(indlabels, dens.type == "native perennial"), aes(x = x, y = y, label = labels), size = 4, fontface = "bold") +
  facet_wrap(~ exp.type, strip.position = "right", nrow = 2) +
  ylab(expression(paste("Response to 50 native perennial grasses ", m^-2, sep = ""))) +
  scale_fill_manual(values = col_pal, name = "Pathogen") +
  scale_color_manual(values = col_pal, name = "Pathogen") +
  scale_shape_manual(values = rep(c(21, 24), each = 4), name = "Pathogen") +
  ggtheme +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        strip.text = element_blank()) +
  ylim(-1, 1.05)

changeplotnon <- ggplot(aindnon, aes(x = grass.group, y = prev.change)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_errorbar(aes(ymin = prev.change - prev.change.se, ymax = prev.change + prev.change.se, shape = path.abb), color = "black", position = position_dodge(width = 0.6), width = 0, show.legend = F) +
  geom_point(aes(fill = path.abb, shape = path.abb), position = position_dodge(width = 0.6), size = 3) +
  geom_text(data = filter(indlabels, dens.type == "non-native annual"), aes(x = x, y = y, label = labels), size = 4, fontface = "bold") +
  facet_wrap(~ exp.type, strip.position = "right", nrow = 2) +
  ylab(expression(paste("Response to 5000 non-native annual grasses ", m^-2, sep = ""))) +
  scale_fill_manual(values = col_pal, name = "Pathogen") +
  scale_color_manual(values = col_pal, name = "Pathogen") +
  scale_shape_manual(values = rep(c(21, 24), each = 4), name = "Pathogen") +
  ggtheme +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = axisTitle, hjust = 0.7),
        legend.text = element_text(size = 9, face = "italic"),
        legend.title = element_text(size = 10)) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  ylim(-1, 1.05)

pdf("./output/figure2_prevalence_change_density.pdf", width = 7, height = 4)
cowplot::plot_grid(changeplotnat, changeplotnon, rel_widths = c(0.7, 1))
dev.off()


#### simulated density effects on pathogens ####

# raw data of focal pathogens
idatFnat <- idatT %>%
  select(isolate.id, experiment, grass.group, native.density, ainf, pcha, plol, ptri, dres, pave, rpro) %>%
  full_join(idatC %>% select(isolate.id, experiment, grass.group, native.density, ainf, pcha, plol, ptri, dres, pave, rpro)) %>%
  group_by(experiment, grass.group) %>%
  mutate(native.density.bins = cut_interval(native.density, n = 5) %>%
           as.character()) %>%
  group_by(experiment, grass.group, native.density.bins) %>%
  mutate(min_interval = parse_number(strsplit(native.density.bins, ",")[[1]])[1],
         max_interval = parse_number(strsplit(native.density.bins, ",")[[1]])[2],
         density = (max_interval + min_interval) / 2) %>%
  ungroup() %>%
  gather(key = "pathogen", value = "present", -c(isolate.id:native.density, native.density.bins:density)) %>%
  group_by(experiment, grass.group, native.density.bins, density, pathogen) %>%
  summarise(prev = mean(present),
            prev.se = sd(present) / sqrt(length(present))) %>%
  ungroup() %>%
  left_join(otus) %>%
  mutate(exp.type = case_when(experiment == "transect" ~ "Observational",
                              experiment == "competition" ~ "Manipulated") %>% factor(levels = c("Observational", "Manipulated")))

idatFnon <- idatT %>%
  select(isolate.id, experiment, grass.group, nonnative.density, ainf, pcha, plol, ptri, dres, pave, rpro) %>%
  full_join(idatC %>% select(isolate.id, experiment, grass.group, nonnative.density, ainf, pcha, plol, ptri, dres, pave, rpro)) %>%
  group_by(experiment, grass.group) %>%
  mutate(nonnative.density.bins = cut_interval(nonnative.density, n = 5) %>%
           as.character()) %>%
  group_by(experiment, grass.group, nonnative.density.bins) %>%
  mutate(min_interval = parse_number(strsplit(nonnative.density.bins, ",")[[1]])[1],
         max_interval = parse_number(strsplit(nonnative.density.bins, ",")[[1]])[2],
         density = (max_interval + min_interval) / 2) %>%
  ungroup() %>%
  gather(key = "pathogen", value = "present", -c(isolate.id:nonnative.density, nonnative.density.bins:density)) %>%
  group_by(experiment, grass.group, nonnative.density.bins, density, pathogen) %>%
  summarise(prev = mean(present),
            prev.se = sd(present) / sqrt(length(present))) %>%
  ungroup() %>%
  left_join(otus) %>%
  mutate(exp.type = case_when(experiment == "transect" ~ "Observational",
                              experiment == "competition" ~ "Manipulated") %>% factor(levels = c("Observational", "Manipulated")))

# idat by grass.status
idatT.nat = filter(idatT, grass.group == "native\nperennial")
idatT.non = filter(idatT, grass.group == "non-native\nannual")
idatC.nat = filter(idatC, grass.group == "native\nperennial")
idatC.non = filter(idatC, grass.group == "non-native\nannual")

# density table
idens <- tibble(nonnative = rep(c(0, 1), each = 200),
                exp.type = rep(rep(c("Observational", "Manipulated"), each = 100), 2),
                natdens.s = c(seq(min(idatT.nat$natdens.s), max(idatT.nat$natdens.s), length.out = 100), 
                              seq(min(idatC.nat$natdens.s), max(idatC.nat$natdens.s), length.out = 100),
                              seq(min(idatT.non$natdens.s), max(idatT.non$natdens.s), length.out = 100), 
                              seq(min(idatC.non$natdens.s), max(idatC.non$natdens.s), length.out = 100)),
                nondens.s = c(seq(min(idatT.nat$nondens.s), max(idatT.nat$nondens.s), length.out = 100), 
                              seq(min(idatC.nat$nondens.s), max(idatC.nat$nondens.s), length.out = 100),
                              seq(min(idatT.non$nondens.s), max(idatT.non$nondens.s), length.out = 100), 
                              seq(min(idatC.non$nondens.s), max(idatC.non$nondens.s), length.out = 100)),
                native.density = c(seq(min(idatT.nat$native.density), max(idatT.nat$native.density), length.out = 100), 
                                   seq(min(idatC.nat$native.density), max(idatC.nat$native.density), length.out = 100),
                                   seq(min(idatT.non$native.density), max(idatT.non$native.density), length.out = 100), 
                                   seq(min(idatC.non$native.density), max(idatC.non$native.density), length.out = 100)),
                nonnative.density = c(seq(min(idatT.nat$nonnative.density), max(idatT.nat$nonnative.density), length.out = 100), 
                                      seq(min(idatC.nat$nonnative.density), max(idatC.nat$nonnative.density), length.out = 100),
                                      seq(min(idatT.non$nonnative.density), max(idatT.non$nonnative.density), length.out = 100), 
                                      seq(min(idatC.non$nonnative.density), max(idatC.non$nonnative.density), length.out = 100))) %>%
  mutate(othdens.s = 0,
         year = NA,
         subplot = NA)

# predict by density
asimnat <- idens %>%
  mutate(nondens.s = 0,
         nonnative.density = NA) %>%
  rename(density = native.density) %>%
  fmod_pred_fun(.)
  
asimnon <- idens %>%
  mutate(natdens.s = 0,
         native.density = NA) %>%
  rename(density = nonnative.density) %>%
  fmod_pred_fun(.)

# manuscript figure function
dens_fig_fun <- function(rawdat, preddat, stat, exp, dodge_width){
  
  rawdat2 <- filter(rawdat, grass.group == stat & exp.type == exp)
  preddat2 <- filter(preddat, grass.group == stat & exp.type == exp)
  
  out_plot <- ggplot(rawdat2, aes(x = density, y = prev, ymin = prev - prev.se, ymax = prev + prev.se, color = path.abb, fill = path.abb, linetype = path.abb, shape = path.abb)) +
    geom_ribbon(data = preddat2, color = NA, alpha = 0.1) +
    geom_line(data = preddat2) +
    geom_errorbar(width = 0.1, position = position_dodge(dodge_width)) +
    geom_point(size = 2, position = position_dodge(dodge_width)) +
    scale_color_manual(values = col_pal, name = "Pathogen") +
    scale_fill_manual(values = col_pal, name = "Pathogen") +
    scale_linetype_manual(values = rep(c("solid", "dashed"), each = 4), name = "Pathogen") +
    scale_shape_manual(values = rep(c(19, 17), each = 4), name = "Pathogen") +
    #scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(-0.05, 1)) +
    theme_bw() +
    theme(axis.text = element_text(size = axisText, color="black"),
          axis.title = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = axisText),
          legend.position = "none")
  
  return(out_plot)
}

# save plots
obs.nat.natdens.plot <- dens_fig_fun(idatFnat, asimnat, "native\nperennial", "Observational", 0.3)
obs.non.natdens.plot <- dens_fig_fun(idatFnat, asimnat, "non-native\nannual", "Observational", 0.3)
obs.nat.nondens.plot <- dens_fig_fun(idatFnon, asimnon, "native\nperennial", "Observational", 30)
obs.non.nondens.plot <- dens_fig_fun(idatFnon, asimnon, "non-native\nannual", "Observational", 30)
man.nat.natdens.plot <- dens_fig_fun(idatFnat, asimnat, "native\nperennial", "Manipulated", 3)
man.non.natdens.plot <- dens_fig_fun(idatFnat, asimnat, "non-native\nannual", "Manipulated", 0.3)
man.nat.nondens.plot <- dens_fig_fun(idatFnon, asimnon, "native\nperennial", "Manipulated", 3)
man.non.nondens.plot <- dens_fig_fun(idatFnon, asimnon, "non-native\nannual", "Manipulated", 30)

# max values
filter(idatFnat, grass.group == "native\nperennial" & exp.type == "Manipulated") %>% select(density) %>% max()
filter(idatFnat, grass.group == "non-native\nannual" & exp.type == "Manipulated") %>% select(density) %>% max()
filter(idatFnon, grass.group == "native\nperennial" & exp.type == "Manipulated") %>% select(density) %>% max()
filter(idatFnon, grass.group == "non-native\nannual" & exp.type == "Manipulated") %>% select(density) %>% max()

# plot settings
ymarg = 20
pathlist1 <- c("ainf", "pave", "pcha")
topLabel = 10
sideLabel = 3

# plot for legend
leg.dens.plot1 <- get_legend(
  ggplot(filter(idatFnat, pathogen %in% pathlist1), aes(x = density, y = prev, color = path.abb, linetype = path.abb, shape = path.abb)) +
    geom_point(size = 2) +
    geom_line(data = filter(asimnat, pathogen %in% pathlist1)) +
    scale_color_manual(values = col_pal, name = "Pathogen") +
    scale_linetype_manual(values = rep(c("solid", "dashed"), each = 4), name = "Pathogen") +
    scale_shape_manual(values = rep(c(19, 17), each = 4), name = "Pathogen") +
    ggtheme +
    theme(legend.text = element_text(size = 9, face = "italic"),
          legend.title = element_text(size = 10),
          legend.key = element_blank(),
          legend.key.size = unit(3, "lines"),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.box.margin = margin(-10, -10, -10, -10),
          legend.margin = margin(-10, -10, -10, -10)) +
    guides(color = guide_legend(nrow = 1, title.position = "left"), linetype = guide_legend(nrow = 1, title.position = "left"), shape = guide_legend(nrow = 1, title.position = "left")))

leg.dens.plot2 <- get_legend(
  ggplot(filter(idatFnat, !(pathogen %in% pathlist1)), aes(x = density, y = prev, color = path.abb, linetype = path.abb, shape = path.abb)) +
    geom_point(size = 2) +
    geom_line(data = filter(asimnat, !(pathogen %in% pathlist1))) +
    scale_color_manual(values = col_pal[4:7], name = "Pathogen") +
    scale_linetype_manual(values = rep(c("solid", "dashed"), each = 4)[4:7], name = "Pathogen") +
    scale_shape_manual(values = rep(c(19, 17), each = 4)[4:7], name = "Pathogen") +
    ggtheme +
    theme(legend.text = element_text(size = 9, face = "italic"),
          legend.title = element_blank(),
          legend.key = element_blank(),
          legend.key.size = unit(3, "lines"),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.box.margin = margin(-10, -10, -10, -10),
          legend.margin = margin(-10, -10, -10, -10)) +
    guides(color = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)))


# paste plots
obs.natdens.plot <- cowplot::plot_grid(obs.nat.natdens.plot + 
                                         ggtitle("Observational") + 
                                         annotate(geom = "text", label = "(a)", x = 0, y = 1, color = "black", size = 4, fontface = "bold") +
                                         theme(plot.title = element_text(size = topLabel, hjust = 0.5)), 
                                  obs.non.natdens.plot + 
                                    annotate(geom = "text", label = "(c)", x = 0, y = 1, color = "black", size = 4, fontface = "bold"), 
                                  ncol = 1)

man.natdens.plot <- cowplot::plot_grid(man.nat.natdens.plot + 
                                    ggtitle("Manipulated") + 
                                      annotate(geom = "text", label = "(b)", x = 0, y = 1, color = "black", size = 4, fontface = "bold") + 
                                   annotate(geom = "text", label = "Native perennial hosts", x = 418.8, y = 0.5, angle = 270, color = "black", size = sideLabel) +
                                    coord_cartesian(xlim = c(0, 379), clip = 'off') + 
                                     theme(plot.title = element_text(size = topLabel, hjust = 0.5),
                                           plot.margin = unit(c(5.5, ymarg, 5.5, 5.5), "pt"),
                                           axis.text.y = element_blank()), 
                                  man.non.natdens.plot + 
                                    annotate(geom = "text", label = "(d)", x = 0, y = 1, color = "black", size = 4, fontface = "bold") + 
                                    annotate(geom = "text", label = "Non-native annual hosts", x = 20, y = 0.5, angle =270, color = "black", size = sideLabel) +
                                    coord_cartesian(xlim = c(0, 18), clip = 'off') + 
                                    theme(plot.margin = unit(c(5.5, ymarg, 5.5, 5.5), "pt"),
                                          axis.text.y = element_blank()), 
                                  ncol = 1)

natdens.plot <- cowplot::plot_grid(obs.natdens.plot, 
                                   man.natdens.plot, 
                                   nrow = 1) 

natdens.plot.fin <- ggdraw(add_sub(natdens.plot, expression(paste("Native perennial grass density (", m^-1, ")", sep = "")), vpadding = grid::unit(0,"lines"), size = axisTitle))

obs.nondens.plot <- cowplot::plot_grid(obs.nat.nondens.plot +
                                         annotate(geom = "text", label = "(e)", x = 0, y = 1, color = "black", size = 4, fontface = "bold"),
                                       obs.non.nondens.plot +
                                         annotate(geom = "text", label = "(g)", x = 0, y = 1, color = "black", size = 4, fontface = "bold"), 
                                       ncol = 1)

man.nondens.plot <- cowplot::plot_grid(man.nat.nondens.plot + 
                                         annotate(geom = "text", label = "(f)", x = 0, y = 1, color = "black", size = 4, fontface = "bold") + 
                                         annotate(geom = "text", label = "Native hosts", x = 250.9, y = 0.5, angle = 270, color = "black", size = sideLabel) +
                                         coord_cartesian(xlim = c(0, 227), clip = 'off') + 
                                         theme(plot.margin = unit(c(5.5, ymarg, 5.5, 5.5), "pt"),
                                               axis.text.y = element_blank()), 
                                       man.non.nondens.plot + 
                                         theme(plot.margin = unit(c(5.5, ymarg, 5.5, 5.5), "pt"),
                                               axis.text.y = element_blank()) + 
                                         annotate(geom = "text", label = "(h)", x = 0, y = 1, color = "black", size = 4, fontface = "bold") + 
                                         annotate(geom = "text", label = "Non-native annual hosts", x = 2848.3, y = 0.5, angle = 270, color = "black", size = sideLabel) +
                                         coord_cartesian(xlim = c(0, 2577), clip = 'off'), 
                                       ncol = 1)

nondens.plot <- cowplot::plot_grid(obs.nondens.plot, man.nondens.plot, nrow = 1) 

nondens.plot.fin <- ggdraw(add_sub(nondens.plot, expression(paste("Non-native annual grass density (", m^-1, ")", sep = "")), vpadding = grid::unit(0,"lines"), size = axisTitle))

dens.plot <- cowplot::plot_grid(natdens.plot.fin, nondens.plot.fin, nrow = 2)
leg.dens.plot <- cowplot::plot_grid(leg.dens.plot1, leg.dens.plot2, align = "v", nrow = 2) 
dens.plot.fin <- cowplot::plot_grid(dens.plot, leg.dens.plot, rel_heights = c(0.8, 0.07), rel_widths = c(1, 0.08), align = "v", nrow = 2)

pdf("./output/figureS5_pathogen_relative_abundance_density.pdf", width = 7, height = 10.5)
grid.arrange(arrangeGrob(dens.plot.fin, left = textGrob("Pathogen relative abundance", gp=gpar(fontsize = 13), rot = 90)))
dev.off()


#### simulated density effects on damage ####

# ran the following code twice: once with March data and models and once with April, saved output separately

# change month
# pdatT <- pdatTm
pdatT <- pdatTa

# use plant-scale data because it's the same dataset as that used for mean.dam, but it includes the zeros
pdatT.nat <- filter(pdatT, grass.group == "native\nperennial")
pdatT.non <- filter(pdatT, grass.group == "non-native\nannual")
pdatC.nat <- filter(pdatC, grass.group == "native\nperennial")
pdatC.non <- filter(pdatC, grass.group == "non-native\nannual")

# density table
ddens <- tibble(nonnative = rep(c(0, 1), each = 200),
                exp.type = rep(rep(c("Observational", "Manipulated"), each = 100), 2),
                natdens.s = c(seq(min(pdatT.nat$natdens.s), max(pdatT.nat$natdens.s), length.out = 100), 
                              seq(min(pdatC.nat$natdens.s), max(pdatC.nat$natdens.s), length.out = 100),
                              seq(min(pdatT.non$natdens.s), max(pdatT.non$natdens.s), length.out = 100), 
                              seq(min(pdatC.non$natdens.s), max(pdatC.non$natdens.s), length.out = 100)),
                nondens.s = c(seq(min(pdatT.nat$nondens.s), max(pdatT.nat$nondens.s), length.out = 100), 
                              seq(min(pdatC.nat$nondens.s), max(pdatC.nat$nondens.s), length.out = 100),
                              seq(min(pdatT.non$nondens.s), max(pdatT.non$nondens.s), length.out = 100), 
                              seq(min(pdatC.non$nondens.s), max(pdatC.non$nondens.s), length.out = 100)),
                native.density = c(seq(min(pdatT.nat$native.density), max(pdatT.nat$native.density), length.out = 100), 
                                   seq(min(pdatC.nat$native.density), max(pdatC.nat$native.density), length.out = 100),
                                   seq(min(pdatT.non$native.density), max(pdatT.non$native.density), length.out = 100), 
                                   seq(min(pdatC.non$native.density), max(pdatC.non$native.density), length.out = 100)),
                nonnative.density = c(seq(min(pdatT.nat$nonnative.density), max(pdatT.nat$nonnative.density), length.out = 100), 
                                      seq(min(pdatC.nat$nonnative.density), max(pdatC.nat$nonnative.density), length.out = 100),
                                      seq(min(pdatT.non$nonnative.density), max(pdatT.non$nonnative.density), length.out = 100), 
                                      seq(min(pdatC.non$nonnative.density), max(pdatC.non$nonnative.density), length.out = 100))) %>%
  mutate(othdens.s = 0,
         year = NA,
         subplot = NA,
         plant = NA)

# change models depending on March or April
# smodT <- smodTM
# imodT <- imodTMa
# scmodT <- scmodTM

smodT <- smodTAa
imodT <- imodTAa
scmodT <- scmodTA2

# function for predicting damage from models
dmod_pred_fun <- function(dat) {
  dsim <- dat %>%
    merge(tibble(dam.type = c("surface", "leaves", "scaled surface")), all = T) %>%
    mutate(dam = case_when(dam.type == "surface" & exp.type == "Observational" ~ 
                             predict(smodT, newdata = ., re.form = NA, type = "response"),
                           dam.type == "surface" & exp.type == "Manipulated" ~ 
                             predict(smodCa, newdata = ., re.form = NA, type = "response"),
                           dam.type == "leaves" & exp.type == "Observational" ~ 
                             predict(imodT, newdata = ., re.form = NA, type = "response"),
                           dam.type == "leaves" & exp.type == "Manipulated" ~ 
                             predict(imodCa, newdata = ., re.form = NA, type = "response"),
                           dam.type == "scaled surface" & exp.type == "Observational" ~ 
                             predict(scmodT, newdata = ., re.form = NA, type = "response"),
                           dam.type == "scaled surface" & exp.type == "Manipulated" ~ 
                             predict(scmodC, newdata = ., re.form = NA, type = "response"),
                           TRUE ~ NA_real_),
           dam.se = case_when(dam.type == "surface" & exp.type == "Observational" ~ 
                                 predict(smodT, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit,
                               dam.type == "surface" & exp.type == "Manipulated" ~ 
                                 predict(smodCa, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit,
                               dam.type == "leaves" & exp.type == "Observational" ~ 
                                 predict(imodT, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit,
                               dam.type == "leaves" & exp.type == "Manipulated" ~ 
                                 predict(imodCa, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit,
                              dam.type == "scaled surface" & exp.type == "Observational" ~ 
                                predict(scmodT, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit,
                              dam.type == "scaled surface" & exp.type == "Manipulated" ~ 
                                predict(scmodC, newdata = ., re.form = NA, type = "response", se.fit = T)$se.fit,
                               TRUE ~ NA_real_)) %>%
    as_tibble() %>%
    mutate(grass.group = case_when(nonnative == 0 ~ "native\nperennial",
                                   nonnative == 1 ~ "non-native\nannual"),
           host.group = recode(grass.group, "native\nperennial" = "native perennial", "non-native\nannual" = "non-native annual"),
           exp.type = factor(exp.type, levels = c("Observational", "Manipulated")))
  
  return(dsim)
}


# predict by density
dsimnat <- ddens %>%
  mutate(nondens.s = 0,
         nonnative.density = NA) %>%
  dmod_pred_fun(.) %>%
  rename(density = native.density,
         dens.s = natdens.s) 

dsimnon <- ddens %>%
  mutate(natdens.s = 0,
         native.density = NA) %>%
  dmod_pred_fun(.) %>%
  rename(density = nonnative.density,
         dens.s = nondens.s)

# change month data
# ldatT <- ldatTm
ldatT <- ldatTa

# raw data: surface area
surfdat <- ldatT %>%
  filter(surface > 0) %>%
  select(experiment, grass.group, native.density, nonnative.density, surface) %>%
  full_join(ldatC %>%
              filter(surface > 0) %>%
              select(experiment, grass.group, native.density, nonnative.density, surface)) %>%
  rename(dam = surface) %>%
  mutate(dam.type = "surface")

# raw data: leaves
leafdat <- pdatT %>%
  select(experiment, grass.group, native.density, nonnative.density, prop.dam) %>%
  full_join(pdatC %>%
              select(experiment, grass.group, native.density, nonnative.density, prop.dam)) %>%
  rename(dam = prop.dam) %>%
  mutate(dam.type = "leaves")

# raw data: scaled surface area
scaledsurfdat <- ldatT %>%
  select(experiment, grass.group, native.density, natdens.s, nonnative.density, nondens.s, surface.scaled) %>%
  full_join(ldatC %>%
              select(experiment, grass.group, native.density, natdens.s, nonnative.density, nondens.s, surface.scaled)) %>%
  rename(dam = surface.scaled) %>%
  mutate(dam.type = "scaled surface") %>%
  mutate(exp.type = case_when(experiment == "transect" ~ "Observational",
                              experiment == "competition" ~ "Manipulated") %>%
           factor(levels = c("Observational", "Manipulated")))

# combine raw data
leafsurfdat <- rbind(leafdat, surfdat) %>%
  mutate(exp.type = case_when(experiment == "transect" ~ "Observational",
                              experiment == "competition" ~ "Manipulated") %>%
           factor(levels = c("Observational", "Manipulated")))

# split by density and calculate averages
rawdamdatnat <- leafsurfdat %>%
  select(-nonnative.density) %>% 
  group_by(exp.type) %>%
  mutate(native.density.bins = cut_interval(native.density, n = 5) %>%
           as.character()) %>%
  group_by(exp.type, native.density.bins) %>%
  mutate(min_interval = parse_number(strsplit(native.density.bins, ",")[[1]])[1],
         max_interval = parse_number(strsplit(native.density.bins, ",")[[1]])[2],
         density = (max_interval + min_interval) / 2) %>%
  ungroup() %>%
  mutate(host.group = recode(grass.group, "native\nperennial" = "native perennial", "non-native\nannual" = "non-native annual")) %>%
  group_by(exp.type, grass.group, host.group, native.density.bins, density, dam.type) %>%
  summarise(dam.se = sd(dam) / sqrt(length(dam)),
            dam = mean(dam)) %>%
  ungroup()

rawdamdatnon <- leafsurfdat %>%
  select(-native.density) %>% 
  group_by(exp.type) %>%
  mutate(nonnative.density.bins = cut_interval(nonnative.density, n = 5) %>%
           as.character()) %>%
  group_by(exp.type, nonnative.density.bins) %>%
  mutate(min_interval = parse_number(strsplit(nonnative.density.bins, ",")[[1]])[1],
         max_interval = parse_number(strsplit(nonnative.density.bins, ",")[[1]])[2],
         density = (max_interval + min_interval) / 2) %>%
  ungroup() %>%
  mutate(host.group = recode(grass.group, "native\nperennial" = "native perennial", "non-native\nannual" = "non-native annual")) %>%
  group_by(exp.type, grass.group, host.group, nonnative.density.bins, density, dam.type) %>%
  summarise(dam.se = sd(dam) / sqrt(length(dam)),
            dam = mean(dam)) %>%
  ungroup()

# repeat for scaled surface data
rawdamdatnatscaled <- scaledsurfdat %>%
  mutate(host.group = recode(grass.group, "native\nperennial" = "native perennial", "non-native\nannual" = "non-native annual")) %>%
  select(-c(nondens.s, nonnative.density)) %>% 
  group_by(exp.type, host.group) %>%
  mutate(native.density.bins = cut_width(native.density, 1) %>%
           as.character()) %>%
  group_by(exp.type, host.group, native.density.bins) %>%
  summarise(min_interval = parse_number(strsplit(native.density.bins, ",")[[1]])[1],
         max_interval = parse_number(strsplit(native.density.bins, ",")[[1]])[2],
         density = (max_interval + min_interval) / 2,
         dam.se = sd(dam) / sqrt(length(dam)),
         dam = mean(dam)) %>%
  ungroup()

rawdamdatnonscaled <- scaledsurfdat %>%
  mutate(host.group = recode(grass.group, "native\nperennial" = "native perennial", "non-native\nannual" = "non-native annual")) %>%
  select(-c(natdens.s, native.density)) %>% 
  group_by(exp.type, host.group) %>%
  mutate(nonnative.density.bins = cut_width(nonnative.density, 1) %>%
           as.character()) %>%
  group_by(exp.type, host.group, nonnative.density.bins) %>%
  summarise(min_interval = parse_number(strsplit(nonnative.density.bins, ",")[[1]])[1],
         max_interval = parse_number(strsplit(nonnative.density.bins, ",")[[1]])[2],
         density = (max_interval + min_interval) / 2,
         dam.se = sd(dam) / sqrt(length(dam)),
         dam = mean(dam)) %>%
  ungroup()

# # labels
# damlabels <- tibble(dens.type = rep(c("native perennial", "non-native annual"), each = 2), 
#                     exp.type = c("Observational", "Manipulated", "Observational", "Manipulated"), 
#                     density = rep(0.5, 4), dam = rep(1.1, 4), 
#                     labels = c("A", "B", "C", "D")) %>%
#   mutate(exp.type = factor(exp.type, levels = c("Observational", "Manipulated")),
#          host.group = "native perennial",
#          dam.type = "leaves")

# # native density plot 
# natdamplot <- ggplot(filter(dsimnat, dam.type != "scaled surface"), aes(x = density, y = dam)) +
#   geom_ribbon(alpha = 0.1, color = NA, aes(ymin = dam - dam.se, ymax = dam + dam.se, fill = dam.type, size = host.group)) +
#   geom_line(aes(linetype = host.group, color = dam.type)) +
#   geom_errorbar(data = rawdamdatnat, width = 0.1, position = position_dodge(0.7), aes(ymin = dam - dam.se, ymax = dam + dam.se, fill = host.group, color = dam.type)) +
#   geom_point(data = rawdamdatnat, size = 2, position = position_dodge(0.7), aes(shape = host.group, fill = host.group, color = dam.type)) +
#   geom_text(data = filter(damlabels, dens.type == "native perennial"), aes(label = labels), size = 4, fontface = "bold") +
#   facet_rep_wrap(~ exp.type, scales = "free_x") +
#   scale_linetype_manual(values = c("solid", "dashed"), name = "Host group") +
#   scale_shape_manual(values = c(19, 21), name = "Host group") +
#   scale_color_manual(values = col_pal[1:2], name = "Metric") +
#   scale_fill_manual(values = c("black", "black", "white", "white"), name = "Metric") +
#   xlab(expression(paste("Native perennial grass density (individuals ", m^-2, ")", sep = ""))) +
#   ylab("Disease severity") +
#   ggtheme +
#   theme(legend.position = "none")
# 
# # non-native density plot  
# nondamplot <- ggplot(filter(dsimnon, dam.type != "scaled surface"), aes(x = density, y = dam)) +
#   geom_ribbon(alpha = 0.1, color = NA, aes(ymin = dam - dam.se, ymax = dam + dam.se, fill = dam.type, size = host.group), show.legend = F) +
#   geom_line(aes(linetype = host.group, color = dam.type)) +
#   geom_errorbar(data = rawdamdatnon, width = 0.1, position = position_dodge(0.7), aes(ymin = dam - dam.se, ymax = dam + dam.se, fill = host.group, color = dam.type), show.legend = F) +
#   geom_point(data = rawdamdatnon, size = 2, position = position_dodge(0.7), aes(shape = host.group, fill = host.group, color = dam.type)) +
#   geom_text(data = filter(damlabels, dens.type == "non-native annual"), aes(label = labels), size = 4, fontface = "bold", show.legend = F) +
#   facet_rep_wrap(~ exp.type, scales = "free_x") +
#   scale_linetype_manual(values = c("solid", "dashed"), name = "Host group") +
#   scale_shape_manual(values = c(19, 21), name = "Host group") +
#   scale_color_manual(values = col_pal[1:2], name = "Metric") +
#   scale_fill_manual(values = c("black", "black", "white", "white"), name = "Metric", guide = F) +
#   xlab(expression(paste("Non-native annual grass density (individuals ", m^-2, ")", sep = ""))) +
#   ylab("Disease severity") +
#   ggtheme +
#   theme(strip.text = element_blank(),
#         legend.text = element_text(size = 9),
#         legend.title = element_text(size = 10),
#         legend.key.size = unit(2.5, "lines"),
#         legend.position = "bottom",
#         legend.box = "vertical",
#         legend.margin = margin(0, 0, 0, 0),
#         legend.spacing.y = unit(-0.5, "cm"),
#         legend.box.margin = margin(-10, -10, -10, -10),
#         plot.margin = margin(3, 3, 3, 3)) +
#   guides(shape = guide_legend(override.aes = list(alpha = 1)), 
#          color = guide_legend(override.aes = list(alpha = 1)))
# 
# # combine plots
# dam.plot <- cowplot::plot_grid(natdamplot, nondamplot, nrow = 2, rel_heights = c(0.95, 1))
# note: The ribbons are lost on the surface lines and there are errors about size, but the plots are accurate. You can't see the ribbons on the surface lines even when they're shaded.
# dam.plot

# March labels
# damlabelsscaled <- tibble(dens.type = rep(c("native perennial", "non-native annual"), each = 2), 
#                           exp.type = c("Observational", "Manipulated", "Observational", "Manipulated"), 
#                           density = rep(0.5, 4), dam = c(0.045, 0.045, 0.075, 0.075), 
#                           labels = c("(a)", "(b)", "(c)", "(d)")) %>%
#   mutate(exp.type = factor(exp.type, levels = c("Observational", "Manipulated")),
#          host.group = "native perennial")

# April labels
damlabelsscaled <- tibble(dens.type = rep(c("native perennial", "non-native annual"), each = 2), 
                          exp.type = c("Observational", "Manipulated", "Observational", "Manipulated"), 
                          density = c(0.4, 4, 4, 4), 
                          dam = c(0.11, 0.11, 0.21, 0.21), 
                          labels = c("(a)", "(b)", "(c)", "(d)")) %>%
  mutate(exp.type = factor(exp.type, levels = c("Observational", "Manipulated")),
         host.group = "native perennial")

# native density plot scaled
natdamplotscaled <- ggplot(filter(dsimnat, dam.type == "scaled surface"), aes(x = density, y = dam)) +
  geom_ribbon(color = NA, fill = "black", aes(ymin = dam - dam.se, ymax = dam + dam.se, alpha = host.group)) +
  geom_line(color = "black", aes(linetype = host.group)) +
  #geom_errorbar(data = rawdamdatnatscaled, width = 0, size = 0.3, aes(ymin = dam - dam.se, ymax = dam + dam.se, color = host.group)) +
  geom_point(data = rawdamdatnatscaled, size = 2, color = "black", alpha = 0.8, aes(shape = host.group, fill = host.group)) +
  geom_text(data = filter(damlabelsscaled, dens.type == "native perennial"), aes(label = labels), size = 4, fontface = "bold") +
  facet_rep_wrap(~ exp.type, scales = "free_x") +
  scale_alpha_manual(values = c(0.1, 0.1), guide = F) +
  scale_color_manual(values = c("black", "black"), guide = F) +
  scale_linetype_manual(values = c("solid", "dashed"), name = "Host group") +
  scale_shape_manual(values = c(19, 21), name = "Host group") +
  scale_fill_manual(values = c("black", "white"), name = "Host group") +
  xlab(expression(paste("Native perennial grass density (individuals ", m^-2, ")", sep = ""))) +
  ylab("Scaled disease severity") +
  ggtheme +
  theme(legend.position = "none")

# non-native density plot
nondamplotscaled <- ggplot(filter(dsimnon, dam.type == "scaled surface"), aes(x = density, y = dam)) +
  geom_ribbon(color = NA, fill = "black", aes(ymin = dam - dam.se, ymax = dam + dam.se, alpha = host.group)) +
  geom_line(color = "black", aes(linetype = host.group)) +
  #geom_errorbar(data = rawdamdatnonscaled, width = 0, size = 0.3, show.legend = F, aes(ymin = dam - dam.se, ymax = dam + dam.se, color = host.group)) +
  geom_point(data = rawdamdatnonscaled, size = 2, color = "black", alpha = 0.8, aes(shape = host.group, fill = host.group)) +
  geom_text(data = filter(damlabelsscaled, dens.type == "non-native annual"), aes(label = labels), size = 4, fontface = "bold") +
  facet_rep_wrap(~ exp.type, scales = "free_x") +
  scale_alpha_manual(values = c(0.1, 0.1), guide = F) +
  scale_color_manual(values = c("black", "black"), guide = F) +
  scale_linetype_manual(values = c("solid", "dashed"), name = "Host group") +
  scale_shape_manual(values = c(19, 21), name = "Host group") +
  scale_fill_manual(values = c("black", "white"), name = "Host group", guide = F) +
  xlab(expression(paste("Non-native annual grass density (individuals ", m^-2, ")", sep = ""))) +
  ylab("Scaled disease severity") +
  ggtheme +
  theme(strip.text = element_blank(),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10),
        legend.key.size = unit(2.5, "lines"),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.y = unit(-0.5, "cm"),
        legend.box.margin = margin(-10, -10, -10, -10),
        plot.margin = margin(3, 3, 3, 3)) +
  guides(shape = guide_legend(override.aes = list(alpha = 1)), 
         color = guide_legend(override.aes = list(alpha = 1)))

# combine plots
dam.plot.scaled <- cowplot::plot_grid(natdamplotscaled, nondamplotscaled, nrow = 2, rel_heights = c(0.95, 1))

# pdf("./output/figure3_pathogen_damage_scaled_density_march.pdf", width = 5, height = 6)
# dam.plot.scaled
# dev.off()

pdf("./output/figureS6_pathogen_damage_scaled_density_april.pdf", width = 6, height = 7)
dam.plot.scaled
dev.off()


#### model summary tables ####

# print summaries and paste full average coefficient tables into Excel (couldn't find a way to do it autmatically)
summary(aiamodTa)
summary(dramodTa)
summary(paamodTa)
summary(pcamodTa)
summary(plamodTa)
summary(ptamodTa)
summary(rpamodTa)

summary(aiamodCa)
summary(dramodCa)
summary(paamodCa)
summary(pcamodCa)
summary(plamodCa)
summary(ptamodCa)

summary(smodTM)
summary(imodTMa)
summary(smodCa)
summary(imodCa)

summary(scmodTM)
summary(scmodTA2)
summary(scmodC)

importance(aiamodTa)
importance(dramodTa)
importance(paamodTa)
importance(pcamodTa)
importance(plamodTa)
importance(ptamodTa)
importance(rpamodTa)

importance(aiamodCa)
importance(dramodCa)
importance(paamodCa)
importance(pcamodCa)
importance(plamodCa)
importance(ptamodCa)

importance(imodTMa)
importance(smodCa)
importance(imodCa)


#### density plots ####

# Note that we can't merge the data using scaled values from each analysis because the same plot may have a different value in the damage and infection experiment. New scaled values are created.

# simplify datasets
iplotsT <- idatT %>%
  select(experiment, year, plot, subplot, nonnative.density, native.density, total.density, nonnative.rel) %>%
  unique() %>%
  mutate(infection = 1)

dplotsT <- pdatTm %>%
  select(experiment, year, plot, subplot, nonnative.density, native.density, total.density, nonnative.rel) %>%
  unique() %>%
  mutate(damage = 1)

iplotsC <- idatC %>%
  select(experiment, year, plot, subplot, nonnative.density, native.density, total.density, nonnative.rel) %>%
  unique() %>%
  mutate(infection = 1)

dplotsC <- pdatC %>%
  select(experiment, year, plot, subplot, nonnative.density, native.density, total.density, nonnative.rel) %>%
  unique() %>%
  mutate(damage = 1)

# merge data
plotsT <- full_join(iplotsT, dplotsT) %>%
  mutate(natdens.s = scale(native.density)[,1],
         nondens.s = scale(nonnative.density)[,1])

plotsC <- full_join(iplotsC, dplotsC) %>%
  mutate(natdens.s = scale(native.density)[,1],
         nondens.s = scale(nonnative.density)[,1])

plots <- full_join(plotsC, plotsT) %>%
  mutate(data.type = case_when(damage == 1 & infection == 1 ~ "both",
                               damage == 1 & is.na(infection) ~ "severity only",
                               is.na(damage) & infection == 1 ~ "isolates only") %>%
           factor(levels = c("isolates only", "severity only", "both")),
         native.rel = native.density / total.density,
         year.f = as.factor(year),
         exp.type = case_when(experiment == "transect" ~ "Observational",
                              experiment == "competition" ~ "Manipulated") %>% factor(levels = c("Observational", "Manipulated")))

# labels
denslabels <- tibble(exp.type = c("Observational", "Manipulated", "Observational", "Manipulated"), x = c(0, 0, -0.01, -0.01), y = c(1670, 2720, 1.05, 1.05), plot.type = c("density", "density", "relative", "relative"), labels = c("(a)", "(b)", "(c)", "(d)")) %>%
  mutate(exp.type = factor(exp.type, levels = c("Observational", "Manipulated")))

# figures for manuscript
densplot <- ggplot(plots, aes(native.density, nonnative.density)) +
  geom_point(alpha = 0.6, aes(fill = data.type, shape = year.f)) +
  geom_text(data = filter(denslabels, plot.type == "density"), aes(x = x, y = y, label = labels), size = 4, fontface = "bold") +
  facet_wrap(~exp.type, scales = "free") +
  scale_fill_manual(values = c("white", "black", "red")) +
  scale_shape_manual(values = c(21, 24)) +
  xlab(expression(paste("Native perennial grass density (", m^-2, ")", sep = ""))) +
  ylab(expression(paste("Non-native annual grass density (", m^-2, ")", sep = ""))) +
  ggtheme +
  guides(fill = "none", shape = "none")

relplot2 <- ggplot(plots, aes(native.rel, nonnative.rel)) +
  geom_point(aes(fill = data.type, shape = year.f), alpha = 0.6,  position = position_jitter(0.01, 0.01, seed = 2)) +
  geom_text(data = filter(denslabels, plot.type == "relative"), aes(x = x, y = y, label = labels), size = 4, fontface = "bold") +
  facet_wrap(~exp.type, scales = "free") +
  scale_fill_manual(values = c("white", "black", "red"), name = "Data type") +  
  scale_shape_manual(values = c(21, 24), name = "Year") +
  xlab("Native perennial grass relative abundance") +
  ylab("Non-native annual grass\nrelative abundance") +
  ggtheme +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal") 

leg2 <- get_legend(relplot2 + theme(legend.background = element_blank(), legend.box.background = element_rect(color = "black")) + guides(fill = guide_legend(override.aes = list(shape=21))))

pdf("./output/figureS3_plant_community_density.pdf", width = 7, height = 7)
cowplot::plot_grid(densplot, relplot2 + theme(legend.position = "none"), leg2,
                   nrow = 3,
                   rel_heights = c(1, 1, 0.15))
dev.off()


#### numbers in text ####

# A inf in native vs. non-native, observational
ainf_obs_nat <- exp(aiamodTa$coefficients["full",1]) / (1 + exp(aiamodTa$coefficients["full",1]))
ainf_obs_non <- exp(aiamodTa$coefficients["full",1] + aiamodTa$coefficients["full",2]) / (1 + exp(aiamodTa$coefficients["full",1] + aiamodTa$coefficients["full",2]))
(ainf_obs_non - ainf_obs_nat) / ainf_obs_nat

# Drec with native dens, manipulated
drec_man_natdens_0 <- exp(dramodCa$coefficients["full",1]) / (1 + exp(dramodCa$coefficients["full",1]))
drec_man_natdens_1 <- exp(dramodCa$coefficients["full",1] + dramodCa$coefficients["full",2]) / (1 + exp(dramodCa$coefficients["full",1] + dramodCa$coefficients["full",2]))
(drec_man_natdens_1 - drec_man_natdens_0) / drec_man_natdens_0

# Drec with non-native dens, manipulated
drec_man_nondens_0 <- exp(dramodCa$coefficients["full",1]) / (1 + exp(dramodCa$coefficients["full",1]))
drec_man_nondens_1 <- exp(dramodCa$coefficients["full",1] + dramodCa$coefficients["full",3]) / (1 + exp(dramodCa$coefficients["full",1] + dramodCa$coefficients["full",3]))
(drec_man_nondens_1 - drec_man_nondens_0) / drec_man_nondens_0

# Drec on non with native dens, manipulated
drec_man_non_natdens_0 <- exp(dramodCa$coefficients["full",1] + dramodCa$coefficients["full",4]) / (1 + exp(dramodCa$coefficients["full",1] + dramodCa$coefficients["full",4]))
drec_man_non_natdens_1 <- exp(dramodCa$coefficients["full",1] + dramodCa$coefficients["full",2] + dramodCa$coefficients["full",4] + dramodCa$coefficients["full",5]) / (1 + exp(dramodCa$coefficients["full",1] + dramodCa$coefficients["full",2] + dramodCa$coefficients["full",4] + dramodCa$coefficients["full",5]))
(drec_man_non_natdens_1 - drec_man_non_natdens_0) / drec_man_non_natdens_0

# P ave with non-native dens, manipulated
pave_man_nondens_0 <- exp(paamodCa$coefficients["full",1]) / (1 + exp(paamodCa$coefficients["full",1]))
pave_man_nondens_1 <- exp(paamodCa$coefficients["full",1] + paamodCa$coefficients["full",3]) / (1 + exp(paamodCa$coefficients["full",1] + paamodCa$coefficients["full",3]))
(pave_man_nondens_1 - pave_man_nondens_0) / pave_man_nondens_0

# P cha with non-native dens, manipulated
pcha_man_nondens_0 <- exp(pcamodCa$coefficients["full",1]) / (1 + exp(pcamodCa$coefficients["full",1]))
pcha_man_nondens_1 <- exp(pcamodCa$coefficients["full",1] + pcamodCa$coefficients["full",3]) / (1 + exp(pcamodCa$coefficients["full",1] + pcamodCa$coefficients["full",3]))
(pcha_man_nondens_1 - pcha_man_nondens_0) / pcha_man_nondens_0

# P cha with native dens, manipulated
pcha_man_natdens_1 <- exp(pcamodCa$coefficients["full",1] + pcamodCa$coefficients["full",2]) / (1 + exp(pcamodCa$coefficients["full",1] + pcamodCa$coefficients["full",2]))
(pcha_man_natdens_1 - pcha_man_nondens_0) / pcha_man_nondens_0

# P lol with native dens, observational
plol_obs_natdens_0 <- exp(plamodTa$coefficients["full",1]) / (1 + exp(plamodTa$coefficients["full",1]))
plol_obs_natdens_1 <- exp(plamodTa$coefficients["full",1] + plamodTa$coefficients["full",2]) / (1 + exp(plamodTa$coefficients["full",1] + plamodTa$coefficients["full",2]))
(plol_obs_natdens_1 - plol_obs_natdens_0) / plol_obs_natdens_0

# P lol with native dens, manipulated
plol_man_natdens_0 <- exp(plamodCa$coefficients["full",1]) / (1 + exp(plamodCa$coefficients["full",1]))
plol_man_natdens_1 <- exp(plamodCa$coefficients["full",1] + plamodCa$coefficients["full",2]) / (1 + exp(plamodCa$coefficients["full",1] + plamodCa$coefficients["full",2]))
(plol_man_natdens_1 - plol_man_natdens_0) / plol_man_natdens_0

# R pro with native dens, observational
rpro_man_natdens_0 <- exp(rpamodTa$coefficients["full",1]) / (1 + exp(rpamodTa$coefficients["full",1]))
rpro_man_natdens_1 <- exp(rpamodTa$coefficients["full",1] + rpamodTa$coefficients["full",2]) / (1 + exp(rpamodTa$coefficients["full",1] + rpamodTa$coefficients["full",2]))
(rpro_man_natdens_1 - rpro_man_natdens_0) / rpro_man_natdens_0

# proportion of leaves, observational
prop_obs_nat <- exp(imodTMa$coefficients["full",1]) / (1 + exp(imodTMa$coefficients["full",1]))
prop_obs_non <- exp(imodTMa$coefficients["full",1] + imodTMa$coefficients["full",4]) / (1 + exp(imodTMa$coefficients["full",1] + imodTMa$coefficients["full",4]))
(prop_obs_nat - prop_obs_non) / prop_obs_non

# proportion of leaves, manipulated
prop_man_nat <- exp(imodCa$coefficients["full",1]) / (1 + exp(imodCa$coefficients["full",1]))
prop_man_non <- exp(imodCa$coefficients["full",1] + imodCa$coefficients["full",2]) / (1 + exp(imodCa$coefficients["full",1] + imodCa$coefficients["full",2]))
(prop_man_nat - prop_man_non) / prop_man_non

# proportion with native dens, observational
prop_obs_natdens_0 <- exp(imodTMa$coefficients["full",1]) / (1 + exp(imodTMa$coefficients["full",1]))
prop_obs_natdens_1 <- exp(imodTMa$coefficients["full",1] + imodTMa$coefficients["full",2]) / (1 + exp(imodTMa$coefficients["full",1] + imodTMa$coefficients["full",2]))
(prop_obs_natdens_1 - prop_obs_natdens_0) / prop_obs_natdens_0

# proportion with non-native dens, observational
prop_obs_nondens_0 <- exp(imodTMa$coefficients["full",1]) / (1 + exp(imodTMa$coefficients["full",1]))
prop_obs_nondens_1 <- exp(imodTMa$coefficients["full",1] + imodTMa$coefficients["full",3]) / (1 + exp(imodTMa$coefficients["full",1] + imodTMa$coefficients["full",3]))
(prop_obs_nondens_1 - prop_obs_nondens_0) / prop_obs_nondens_0
