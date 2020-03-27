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
col_pal <- c("#000000", "#56B4E9", "#009E73", "#F0E442", "#000000", "#56B4E9", "#009E73")
col_pal2 <- c("#000000", "#56B4E9", "#009E73", "#F0E442", "gray50", "#81ffff", "#00edad")

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

# function to convert logit to probability
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

# function for extracting coefficients from pathogen models
path_coef_fun <- function(mod, pathogen, exp.type){
  out <- tibble(coef = names(mod$coefficients["full",]), value = mod$coefficients["full",]) %>% 
    mutate(pathogen = pathogen,
           exp.type = exp.type)
}

# table for converting coefficient names
coef_name_conv <- tibble(coef_names = c("cond((Int))", "cond(nonnative)", "cond(natdens.s)", "cond(nondens.s)", "cond(othdens.s)", "cond(natdens.s:nonnative)", "cond(nondens.s:nonnative)", "cond(othdens.s:nonnative)"),
                         orig_names = c("(Intercept)", "nonnative", "natdens.s", "nondens.s", "othdens.s", "nonnative:natdens.s", "nonnative:nondens.s", "nonnative:othdens.s"))

# function for extracting coefficients from damage models
dam_coef_fun <- function(mod, dam.type, exp.type, avg){
  
  if(avg == T){
    coef_names = names(mod$coefficients["full",])
    coef_values = mod$coefficients["full",]
  }

  else{
    out_names = tibble(orig_names = rownames(summary(mod)$coefficients[[1]])) %>%
      left_join(coef_name_conv)
    coef_names = out_names$coef_names
    coef_values = as.numeric(summary(mod)$coefficients[[1]][ ,1])
  }
  
  out <- tibble(coef = coef_names, 
                value = coef_values) %>% 
    mutate(dam.type = dam.type,
           exp.type = exp.type)
  
  return(out)
}


#### table of pathogens ####

# combine possible host species across years
fun2 <- fun %>%
  group_by(pathogen, otu.id, host.num) %>%
  summarise(observed.host.species = toString(host.species)) %>%
  ungroup() %>%
  mutate(observed.host.species = case_when(otu.id == 1 ~ "Ab, Af, Bd, Bh, Eg, Sp",
                                           otu.id == 2 ~ "Ab, Bd, Bh, Eg, Sp",
                                           otu.id == 7 ~ "Bd, Bh, Eg, Sp",
                                           otu.id == 4 ~ "Ab, Af, Bd, Bh, Eg, Sp",
                                           otu.id == 5 ~ "Ab, Af, Bd, Bh, Eg, Sp",
                                           otu.id == 8 ~ "Eg, Sp",
                                           otu.id == 3 ~ "Ab, Bd, Bh, Eg, Sp")) %>% 
  rename(taxonomy = pathogen)

# otu's and code names
otus <- tibble(pathogen = c("ainf", "dres", "pave", "pcha", "plol", "ptri", "rpro"),
               otu.id = c(1, 2, 7, 4, 5, 8, 3),
               path.abb = c("A. inf.", "Drec.", "P. ave.", "P. cha.", "P. lol.", "P. tri.", "R. pro.")) %>%
  left_join(select(fun2, otu.id, taxonomy, observed.host.species, host.num))

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
  group_by(taxonomy, observed.host.species, host.num, grass.group) %>%
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
                exp.type = c(rep("Observational", 3), rep("Manipulated", 3)),
                increase = rep(c("none", "nat", "non"), 2))
inddens

# check values
filter(idatC, native.density == min(native.density)) %>% select(native.density, natdens.s) %>% unique()
filter(idatT, native.density == min(native.density)) %>% select(native.density, natdens.s) %>% unique()
filter(idatC, nonnative.density == min(nonnative.density)) %>% select(nonnative.density, nondens.s) %>% unique()
filter(idatT, nonnative.density == min(nonnative.density)) %>% select(nonnative.density, nondens.s) %>% unique() # seem similar

# extract coefficients for absolute abundance
aind <- rbind(path_coef_fun(rpamodTa, "rpro", "Observational"), path_coef_fun(aiamodTa, "ainf", "Observational"), path_coef_fun(pcamodTa, "pcha", "Observational"), path_coef_fun(plamodTa, "plol", "Observational"), path_coef_fun(ptamodTa, "ptri", "Observational"), path_coef_fun(dramodTa, "dres", "Observational"), path_coef_fun(paamodTa, "pave", "Observational"), path_coef_fun(aiamodCa, "ainf", "Manipulated"), path_coef_fun(pcamodCa, "pcha", "Manipulated"), path_coef_fun(plamodCa, "plol", "Manipulated"), path_coef_fun(ptamodCa, "ptri", "Manipulated"), path_coef_fun(dramodCa, "dres", "Manipulated"), path_coef_fun(paamodCa, "pave", "Manipulated")) %>%
  spread(key = coef, value = value) %>%
  replace(is.na(.), 0) %>%
  full_join(inddens) %>%
  mutate(nat_prev = logit2prob(`cond((Int))` + `cond(natdens.s)` * natdens.s + `cond(nondens.s)` * nondens.s),
         non_prev = logit2prob(`cond((Int))`+ `cond(nonnative)` + `cond(natdens.s)` * natdens.s + `cond(nondens.s)` * nondens.s + `cond(natdens.s:nonnative)` * natdens.s + `cond(nondens.s:nonnative)` * nondens.s)) %>%
  select(exp.type, pathogen, increase, native.density, nonnative.density, nat_prev, non_prev) %>%
  gather(key = host_status, value = prev, -c(pathogen, exp.type, increase, native.density, nonnative.density)) %>%
  mutate(status = recode(host_status, nat_prev = "Native perennial\nhosts", non_prev = "Non-native\nannual hosts") %>% factor(levels = c("Native perennial\nhosts", "Non-native\nannual hosts")),
         exp.type = factor(exp.type, levels = c("Observational", "Manipulated"))) %>%
  left_join(otus)

# split by density and calculate change in prevalence, remove specialists
aindnat <- filter(aind, increase != "non") %>% 
  select(-c(nonnative.density, increase)) %>%
  spread(key = native.density, value = prev) %>%
  mutate(prev_change = `50` - `0`) %>%
  filter(!(pathogen == "pcha" & status == "Native perennial\nhosts") & !(pathogen == "ptri" & status == "Non-native\nannual hosts"))

aindnon <- filter(aind, increase != "nat") %>% 
  select(-c(native.density, increase)) %>%
  spread(key = nonnative.density, value = prev) %>%
  mutate(prev_change = `5000` - `0`) %>%
  filter(!(pathogen == "pcha" & status == "Native perennial\nhosts") & !(pathogen == "ptri" & status == "Non-native\nannual hosts"))

# labels
indlabels <- tibble(dens.type = rep(c("native perennial", "non-native annual"), each = 2), exp.type = c("Observational", "Manipulated", "Observational", "Manipulated"), x = rep(0.5, 4), y = rep(0.97, 4), labels = c("A", "B", "C", "D")) %>%
  mutate(exp.type = factor(exp.type, levels = c("Observational", "Manipulated")))

# figure
changeplotnat <- ggplot(aindnat, aes(x = status, y = prev_change)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_point(aes(fill = path.abb, shape = path.abb), position = position_jitterdodge(dodge.width = 0.1, jitter.width = 0.1, seed = 8), size = 3) +
  geom_text(data = filter(indlabels, dens.type == "native perennial"), aes(x = x, y = y, label = labels), size = 4, fontface = "bold") +
  facet_wrap(~ exp.type, strip.position = "right", nrow = 2) +
  ylab(expression(paste("Response to 50 native perennial grasses ", m^-2, sep = ""))) +
  scale_fill_manual(values = col_pal2, name = "Pathogen") +
  scale_shape_manual(values = rep(c(21, 24), each = 4), name = "Pathogen") +
  ggtheme +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        strip.text = element_blank()) +
  ylim(-0.8, 1)

changeplotnon <- ggplot(aindnon, aes(x = status, y = prev_change)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_point(aes(fill = path.abb, shape = path.abb), position = position_jitterdodge(dodge.width = 0.1, jitter.width = 0.1, seed = 8), size = 3) +
  geom_text(data = filter(indlabels, dens.type == "non-native annual"), aes(x = x, y = y, label = labels), size = 4, fontface = "bold") +
  facet_wrap(~ exp.type, strip.position = "right", nrow = 2) +
  ylab(expression(paste("Response to 5000 non-native annual grasses ", m^-2, sep = ""))) +
  scale_fill_manual(values = col_pal2, name = "Pathogen") +
  scale_shape_manual(values = rep(c(21, 24), each = 4), name = "Pathogen") +
  ggtheme +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = axisTitle, hjust = 0.7)) +
  ylim(-0.8, 1)

pdf("./output/figure4_prevalence_change_density.pdf", width = 7, height = 4)
cowplot::plot_grid(changeplotnat, changeplotnon, rel_widths = c(0.7, 1))
dev.off()


#### simulated density effects on pathogens ####

# raw data of focal pathogens
idatFnat <- idatT %>%
  select(isolate.id, experiment, grass.group, native.density, ainf, pcha, plol, ptri, dres, pave, rpro) %>%
  full_join(idatC %>% select(isolate.id, experiment, grass.group, native.density, ainf, pcha, plol, ptri, dres, pave, rpro)) %>%
  group_by(experiment, grass.group, native.density) %>%
  summarise(ainf = mean(ainf), pcha = mean(pcha), plol = mean(plol), ptri = mean(ptri), dres = mean(dres), pave = mean(pave), rpro = mean(rpro)) %>%
  ungroup() %>%
  gather(key = pathogen, value = prev, -c(experiment:native.density)) %>%
  left_join(otus) %>%
  mutate(exp.type = case_when(experiment == "transect" ~ "Observational",
                              experiment == "competition" ~ "Manipulated") %>% factor(levels = c("Observational", "Manipulated"))) %>%
  rename(density = native.density)

idatFnon <- idatT %>%
  select(isolate.id, experiment, grass.group, nonnative.density, ainf, pcha, plol, ptri, dres, pave, rpro) %>%
  full_join(idatC %>% select(isolate.id, experiment, grass.group, nonnative.density, ainf, pcha, plol, ptri, dres, pave, rpro)) %>%
  group_by(experiment, grass.group, nonnative.density) %>%
  summarise(ainf = mean(ainf), pcha = mean(pcha), plol = mean(plol), ptri = mean(ptri), dres = mean(dres), pave = mean(pave), rpro = mean(rpro)) %>%
  ungroup() %>%
  gather(key = pathogen, value = prev, -c(experiment:nonnative.density)) %>%
  left_join(otus) %>%
  mutate(exp.type = case_when(experiment == "transect" ~ "Observational",
                              experiment == "competition" ~ "Manipulated") %>% factor(levels = c("Observational", "Manipulated"))) %>%
  rename(density = nonnative.density)

# idat by grass.status
idatT.nat = filter(idatT, grass.group == "native\nperennial")
idatT.non = filter(idatT, grass.group == "non-native\nannual")
idatC.nat = filter(idatC, grass.group == "native\nperennial")
idatC.non = filter(idatC, grass.group == "non-native\nannual")

# density table
idens <- tibble(nat.natdens.s = c(seq(min(idatT.nat$natdens.s), max(idatT.nat$natdens.s), length.out = 100), seq(min(idatC.nat$natdens.s), max(idatC.nat$natdens.s), length.out = 100)),
                nat.nondens.s = c(seq(min(idatT.nat$nondens.s), max(idatT.nat$nondens.s), length.out = 100), seq(min(idatC.nat$nondens.s), max(idatC.nat$nondens.s), length.out = 100)),
                nat.native.density = c(seq(min(idatT.nat$native.density), max(idatT.nat$native.density), length.out = 100), seq(min(idatC.nat$native.density), max(idatC.nat$native.density), length.out = 100)),
                nat.nonnative.density = c(seq(min(idatT.nat$nonnative.density), max(idatT.nat$nonnative.density), length.out = 100), seq(min(idatC.nat$nonnative.density), max(idatC.nat$nonnative.density), length.out = 100)),
                non.natdens.s = c(seq(min(idatT.non$natdens.s), max(idatT.non$natdens.s), length.out = 100), seq(min(idatC.non$natdens.s), max(idatC.non$natdens.s), length.out = 100)),
                non.nondens.s = c(seq(min(idatT.non$nondens.s), max(idatT.non$nondens.s), length.out = 100), seq(min(idatC.non$nondens.s), max(idatC.non$nondens.s), length.out = 100)),
                non.native.density = c(seq(min(idatT.non$native.density), max(idatT.non$native.density), length.out = 100), seq(min(idatC.non$native.density), max(idatC.non$native.density), length.out = 100)),
                non.nonnative.density = c(seq(min(idatT.non$nonnative.density), max(idatT.non$nonnative.density), length.out = 100), seq(min(idatC.non$nonnative.density), max(idatC.non$nonnative.density), length.out = 100)),
                exp.type = c(rep("Observational", 100), rep("Manipulated", 100)))

# extract coefficients
asim <- rbind(path_coef_fun(rpamodTa, "rpro", "Observational"), path_coef_fun(aiamodTa, "ainf", "Observational"), path_coef_fun(pcamodTa, "pcha", "Observational"), path_coef_fun(plamodTa, "plol", "Observational"), path_coef_fun(ptamodTa, "ptri", "Observational"), path_coef_fun(dramodTa, "dres", "Observational"), path_coef_fun(paamodTa, "pave", "Observational"), path_coef_fun(aiamodCa, "ainf", "Manipulated"), path_coef_fun(pcamodCa, "pcha", "Manipulated"), path_coef_fun(plamodCa, "plol", "Manipulated"), path_coef_fun(ptamodCa, "ptri", "Manipulated"), path_coef_fun(dramodCa, "dres", "Manipulated"), path_coef_fun(paamodCa, "pave", "Manipulated")) %>%
  spread(key = coef, value = value) %>%
  replace(is.na(.), 0) %>%
  full_join(idens) %>%
  mutate(nat_dens_nat_prev = logit2prob(`cond((Int))` + `cond(natdens.s)` * nat.natdens.s),
         nat_dens_non_prev = logit2prob(`cond((Int))`+ `cond(nonnative)` + `cond(natdens.s)` * non.natdens.s + `cond(natdens.s:nonnative)` * non.natdens.s),
         non_dens_nat_prev = logit2prob(`cond((Int))` + `cond(nondens.s)` * nat.nondens.s),
         non_dens_non_prev = logit2prob(`cond((Int))`+ `cond(nonnative)` + `cond(nondens.s)` * non.nondens.s + `cond(nondens.s:nonnative)` * non.nondens.s)) %>%
  select(exp.type, pathogen, nat.native.density, non.native.density, nat.nonnative.density, non.nonnative.density, nat_dens_nat_prev:non_dens_non_prev) %>%
  gather(key = dens_status, value = prev, -c(pathogen, exp.type, nat.native.density, non.native.density, nat.nonnative.density, non.nonnative.density)) %>%
  mutate(grass.group = substring(dens_status, 10, 12) %>% recode(nat = "native\nperennial", non = "non-native\nannual"),
         dens.type = substring(dens_status, 1, 3),
         exp.type = factor(exp.type, levels = c("Observational", "Manipulated"))) %>%
  left_join(otus) %>%
  filter(!(pathogen == "pcha" & grass.group == "native\nperennial") & !(pathogen == "ptri" & grass.group == "non-native\nannual"))

# split by density
asimnat <- filter(asim, dens.type == "nat") %>%
  mutate(native.density = case_when(grass.group == "native\nperennial" ~ nat.native.density,
                                    grass.group == "non-native\nannual" ~ non.native.density)) %>%
  rename(density = native.density)

asimnon <- filter(asim, dens.type == "non") %>%
  mutate(nonnative.density = case_when(grass.group == "native\nperennial" ~ nat.nonnative.density,
                                       grass.group == "non-native\nannual" ~ non.nonnative.density)) %>%
  rename(density = nonnative.density)

# manuscript figure function
dens_fig_fun <- function(rawdat, preddat, stat, exp){
  
  rawdat2 <- filter(rawdat, grass.group == stat & exp.type == exp)
  preddat2 <- filter(preddat, grass.group == stat & exp.type == exp)
  
  out_plot <- ggplot(rawdat2, aes(x = density, y = prev, color = path.abb, linetype = path.abb, shape = path.abb)) +
    geom_point(alpha = 0.5) +
    geom_line(data = preddat2) +
    scale_color_manual(values = col_pal, name = "Pathogen") +
    scale_linetype_manual(values = rep(c("solid", "dashed"), each = 4), name = "Pathogen") +
    scale_shape_manual(values = rep(c(19, 17), each = 4), name = "Pathogen") +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0, 1)) +
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
obs.nat.natdens.plot <- dens_fig_fun(idatFnat, asimnat, "native\nperennial", "Observational")
obs.non.natdens.plot <- dens_fig_fun(idatFnat, asimnat, "non-native\nannual", "Observational")
obs.nat.nondens.plot <- dens_fig_fun(idatFnon, asimnon, "native\nperennial", "Observational")
obs.non.nondens.plot <- dens_fig_fun(idatFnon, asimnon, "non-native\nannual", "Observational")
man.nat.natdens.plot <- dens_fig_fun(idatFnat, asimnat, "native\nperennial", "Manipulated")
man.non.natdens.plot <- dens_fig_fun(idatFnat, asimnat, "non-native\nannual", "Manipulated")
man.nat.nondens.plot <- dens_fig_fun(idatFnon, asimnon, "native\nperennial", "Manipulated")
man.non.nondens.plot <- dens_fig_fun(idatFnon, asimnon, "non-native\nannual", "Manipulated")

# max values
filter(idatFnat, grass.group == "native\nperennial" & exp.type == "Manipulated") %>% select(density) %>% max()
filter(idatFnat, grass.group == "non-native\nannual" & exp.type == "Manipulated") %>% select(density) %>% max()
filter(idatFnon, grass.group == "native\nperennial" & exp.type == "Manipulated") %>% select(density) %>% max()
filter(idatFnon, grass.group == "non-native\nannual" & exp.type == "Manipulated") %>% select(density) %>% max()

# plot settings
ymarg = 20
pathlist1 <- c("ainf", "dres", "pave")
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
    theme(legend.text = element_text(size = 9),
          legend.title = element_text(size = 10),
          legend.key.size = unit(3, "lines"),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.margin = margin(0, 0, 0, 0, unit="cm")) +
    guides(color = guide_legend(nrow = 1, title.position = "left"), linetype = guide_legend(nrow = 1, title.position = "left"), shape = guide_legend(nrow = 1, title.position = "left")))

leg.dens.plot2 <- get_legend(
  ggplot(filter(idatFnat, !(pathogen %in% pathlist1)), aes(x = density, y = prev, color = path.abb, linetype = path.abb, shape = path.abb)) +
    geom_point(size = 2) +
    geom_line(data = filter(asimnat, !(pathogen %in% pathlist1))) +
    scale_color_manual(values = col_pal[4:7], name = "Pathogen") +
    scale_linetype_manual(values = rep(c("solid", "dashed"), each = 4)[4:7], name = "Pathogen") +
    scale_shape_manual(values = rep(c(19, 17), each = 4)[4:7], name = "Pathogen") +
    ggtheme +
    theme(legend.text = element_text(size = 9),
          legend.title = element_blank(),
          legend.key.size = unit(3, "lines"),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.margin = margin(0, 0, 0, 0, unit="cm")) +
    guides(color = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)))


# paste plots
obs.natdens.plot <- cowplot::plot_grid(obs.nat.natdens.plot + 
                                         ggtitle("Observational") + 
                                         annotate(geom = "text", label = "A", x = 0, y = 1, color = "black", size = 4, fontface = "bold") +
                                         theme(plot.title = element_text(size = topLabel, hjust = 0.5)), 
                                  obs.non.natdens.plot + 
                                    annotate(geom = "text", label = "C", x = 0, y = 1, color = "black", size = 4, fontface = "bold"), 
                                  ncol = 1)

man.natdens.plot <- cowplot::plot_grid(man.nat.natdens.plot + 
                                    ggtitle("Manipulated") + 
                                      annotate(geom = "text", label = "B", x = 0, y = 1, color = "black", size = 4, fontface = "bold") + 
                                   annotate(geom = "text", label = "Native perennial hosts", x = 418.8, y = 0.5, angle = 270, color = "black", size = sideLabel) +
                                    coord_cartesian(xlim = c(0, 379), clip = 'off') + 
                                     theme(plot.title = element_text(size = topLabel, hjust = 0.5),
                                           plot.margin = unit(c(5.5, ymarg, 5.5, 5.5), "pt"),
                                           axis.text.y = element_blank()), 
                                  man.non.natdens.plot + 
                                    annotate(geom = "text", label = "D", x = 0, y = 1, color = "black", size = 4, fontface = "bold") + 
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
                                         annotate(geom = "text", label = "E", x = 0, y = 1, color = "black", size = 4, fontface = "bold"),
                                       obs.non.nondens.plot +
                                         annotate(geom = "text", label = "G", x = 0, y = 1, color = "black", size = 4, fontface = "bold"), 
                                       ncol = 1)

man.nondens.plot <- cowplot::plot_grid(man.nat.nondens.plot + 
                                         annotate(geom = "text", label = "F", x = 0, y = 1, color = "black", size = 4, fontface = "bold") + 
                                         annotate(geom = "text", label = "Native hosts", x = 250.9, y = 0.5, angle = 270, color = "black", size = sideLabel) +
                                         coord_cartesian(xlim = c(0, 227), clip = 'off') + 
                                         theme(plot.margin = unit(c(5.5, ymarg, 5.5, 5.5), "pt"),
                                               axis.text.y = element_blank()), 
                                       man.non.nondens.plot + 
                                         theme(plot.margin = unit(c(5.5, ymarg, 5.5, 5.5), "pt"),
                                               axis.text.y = element_blank()) + 
                                         annotate(geom = "text", label = "H", x = 0, y = 1, color = "black", size = 4, fontface = "bold") + 
                                         annotate(geom = "text", label = "Non-native annual hosts", x = 2848.3, y = 0.5, angle = 270, color = "black", size = sideLabel) +
                                         coord_cartesian(xlim = c(0, 2577), clip = 'off'), 
                                       ncol = 1)

nondens.plot <- cowplot::plot_grid(obs.nondens.plot, man.nondens.plot, nrow = 1) 

nondens.plot.fin <- ggdraw(add_sub(nondens.plot, expression(paste("Non-native annual grass density (", m^-1, ")", sep = "")), vpadding = grid::unit(0,"lines"), size = axisTitle))

dens.plot <- cowplot::plot_grid(natdens.plot.fin, nondens.plot.fin, nrow = 2)
leg.dens.plot <- cowplot::plot_grid(leg.dens.plot1, leg.dens.plot2, align = "v", nrow = 2) 
dens.plot.fin <- cowplot::plot_grid(dens.plot, leg.dens.plot, rel_heights = c(1, 0.08), rel_widths = c(1, 0.08), align = "v", nrow = 2)

pdf("./output/figure3_pathogen_relative_abundance_density.pdf", width = 7, height = 10.5)
grid.arrange(arrangeGrob(dens.plot.fin, left = textGrob("Pathogen relative abundance", gp=gpar(fontsize = 13), rot = 90)))
dev.off()


#### simulated density effects on damage ####

# ran the following code twice: once with March data and models and once with April, saved output separately

# change month
pdatT <- pdatTm
# pdatT <- pdatTa

# use plant-scale data because it's the same dataset as that used for mean.dam, but it includes the zeros
pdatT.nat <- filter(pdatT, grass.group == "native\nperennial")
pdatT.non <- filter(pdatT, grass.group == "non-native\nannual")
pdatC.nat <- filter(pdatC, grass.group == "native\nperennial")
pdatC.non <- filter(pdatC, grass.group == "non-native\nannual")

# density table
ddens <- tibble(nat.natdens.s = c(seq(min(pdatT.nat$natdens.s), max(pdatT.nat$natdens.s), length.out = 100), seq(min(pdatC.nat$natdens.s), max(pdatC.nat$natdens.s), length.out = 100)),
                nat.nondens.s = c(seq(min(pdatT.nat$nondens.s), max(pdatT.nat$nondens.s), length.out = 100), seq(min(pdatC.nat$nondens.s), max(pdatC.nat$nondens.s), length.out = 100)),
                nat.native.density = c(seq(min(pdatT.nat$native.density), max(pdatT.nat$native.density), length.out = 100), seq(min(pdatC.nat$native.density), max(pdatC.nat$native.density), length.out = 100)),
                nat.nonnative.density = c(seq(min(pdatT.nat$nonnative.density), max(pdatT.nat$nonnative.density), length.out = 100), seq(min(pdatC.nat$nonnative.density), max(pdatC.nat$nonnative.density), length.out = 100)),
                non.natdens.s = c(seq(min(pdatT.non$natdens.s), max(pdatT.non$natdens.s), length.out = 100), seq(min(pdatC.non$natdens.s), max(pdatC.non$natdens.s), length.out = 100)),
                non.nondens.s = c(seq(min(pdatT.non$nondens.s), max(pdatT.non$nondens.s), length.out = 100), seq(min(pdatC.non$nondens.s), max(pdatC.non$nondens.s), length.out = 100)),
                non.native.density = c(seq(min(pdatT.non$native.density), max(pdatT.non$native.density), length.out = 100), seq(min(pdatC.non$native.density), max(pdatC.non$native.density), length.out = 100)),
                non.nonnative.density = c(seq(min(pdatT.non$nonnative.density), max(pdatT.non$nonnative.density), length.out = 100), seq(min(pdatC.non$nonnative.density), max(pdatC.non$nonnative.density), length.out = 100)),
                exp.type = c(rep("Observational", 100), rep("Manipulated", 100)))

# change models depending on March or April
smodT <- smodTM
mod.avg <- F
imodT <- imodTMa

# smodT <- smodTAa
# mod.avg <- T
# imodT <- imodTAa

# extract coefficients
dsim <- rbind(dam_coef_fun(smodT, "surface", "Observational", mod.avg), dam_coef_fun(smodCa, "surface", "Manipulated", T), dam_coef_fun(imodT, "leaves", "Observational", T), dam_coef_fun(imodCa, "leaves", "Manipulated", T)) %>%
  spread(key = coef, value = value) %>%
  replace(is.na(.), 0) %>%
  full_join(ddens) %>%
  mutate(nat_dens_nat_dam = logit2prob(`cond((Int))` + `cond(natdens.s)` * nat.natdens.s),
         nat_dens_non_dam = logit2prob(`cond((Int))`+ `cond(nonnative)` + `cond(natdens.s)` * non.natdens.s + `cond(natdens.s:nonnative)` * non.natdens.s),
         non_dens_nat_dam = logit2prob(`cond((Int))` + `cond(nondens.s)` * nat.nondens.s),
         non_dens_non_dam = logit2prob(`cond((Int))`+ `cond(nonnative)` + `cond(nondens.s)` * non.nondens.s + `cond(nondens.s:nonnative)` * non.nondens.s)) %>%
  select(exp.type, dam.type, nat.native.density, nat.nonnative.density, non.native.density, non.nonnative.density, nat_dens_nat_dam:non_dens_non_dam) %>%
  gather(key = dens_status, value = dam, -c(dam.type, exp.type, nat.native.density:non.nonnative.density)) %>%
  mutate(grass.group = substring(dens_status, 10, 12) %>% recode(nat = "native\nperennial", non = "non-native\nannual"),
         dens.type = substring(dens_status, 1, 3),
         exp.type = factor(exp.type, levels = c("Observational", "Manipulated")))

# split by density
dsimnat <- filter(dsim, dens.type == "nat")  %>%
  mutate(native.density = case_when(grass.group == "native\nperennial" ~ nat.native.density,
                                    grass.group == "non-native\nannual" ~ non.native.density),
         host.group = recode(grass.group, "native\nperennial" = "native perennial", "non-native\nannual" = "non-native annual")) %>%
  rename(density = native.density)

dsimnon <- filter(dsim, dens.type == "non")  %>%
  mutate(nonnative.density = case_when(grass.group == "native\nperennial" ~ nat.nonnative.density,
                                       grass.group == "non-native\nannual" ~ non.nonnative.density),
         host.group = recode(grass.group, "native\nperennial" = "native perennial", "non-native\nannual" = "non-native annual")) %>%
  rename(density = nonnative.density)

# change month data
ldatT <- ldatTm
# ldatT <- ldatTa

# raw data: surface area
surfdat <- ldatT %>%
  filter(surface > 0) %>%
  group_by(year, experiment, subplot, host, grass.group, native.density, nonnative.density, plant) %>%
  summarise(mean.dam = mean(surface)) %>%
  ungroup() %>%
  select(experiment, grass.group, native.density, nonnative.density, mean.dam) %>%
  full_join(ldatC %>%
              filter(surface > 0) %>%
              group_by(experiment, subplot, host, grass.group, native.density, nonnative.density, plant) %>%
              summarise(mean.dam = mean(surface)) %>%
              ungroup() %>%
              select(experiment, grass.group, native.density, nonnative.density, mean.dam)) %>%
    rename(dam = mean.dam) %>%
    mutate(exp.type = case_when(experiment == "transect" ~ "Observational",
                                experiment == "competition" ~ "Manipulated") %>%
             factor(levels = c("Observational", "Manipulated")),
           dam.type = "surface")

# raw data: leaves
leafdat <- pdatT %>%
  select(experiment, grass.group, native.density, nonnative.density, prop.dam) %>%
  full_join(pdatC %>%
              select(experiment, grass.group, native.density, nonnative.density, prop.dam)) %>%
  rename(dam = prop.dam) %>%
  mutate(exp.type = case_when(experiment == "transect" ~ "Observational",
                              experiment == "competition" ~ "Manipulated") %>%
           factor(levels = c("Observational", "Manipulated")),
         dam.type = "leaves")

# split by density
rawdamdatnat <- full_join(surfdat, leafdat) %>%
  select(-nonnative.density) %>%
  rename(density = native.density) %>%
  mutate(host.group = recode(grass.group, "native\nperennial" = "native perennial", "non-native\nannual" = "non-native annual"))

rawdamdatnon <- full_join(surfdat, leafdat) %>%
  select(-native.density) %>%
  rename(density = nonnative.density) %>%
  mutate(host.group = recode(grass.group, "native\nperennial" = "native perennial", "non-native\nannual" = "non-native annual"))

# labels
damlabels <- tibble(dens.type = rep(c("native perennial", "non-native annual"), each = 2), exp.type = c("Observational", "Manipulated", "Observational", "Manipulated"), x = rep(0.5, 4), y = rep(1.1, 4), labels = c("A", "B", "C", "D")) %>%
  mutate(exp.type = factor(exp.type, levels = c("Observational", "Manipulated")))

# native density plot 
natdamplot <- ggplot(dsimnat, aes(x = density, y = dam)) +
  geom_line(size = 1, aes(linetype = host.group, color = dam.type)) +
  geom_point(data = rawdamdatnat, alpha = 0.5, aes(shape = host.group, color = dam.type)) +
  geom_text(data = filter(damlabels, dens.type == "native perennial"), aes(x = x, y = y, label = labels), size = 4, fontface = "bold") +
  facet_rep_wrap(~ exp.type, scales = "free_x") +
  scale_linetype_manual(values = c("solid", "dashed"), name = "Host group") +
  scale_shape_manual(values = c(19, 21), name = "Host group") +
  scale_color_manual(values = col_pal[1:2], name = "Metric") +
  xlab(expression(paste("Native perennial grass density (individuals ", m^-2, ")", sep = ""))) +
  ylab("Disease severity") +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  ggtheme +
  theme(legend.position = "none")

# non-native density plot  
nondamplot <- ggplot(dsimnon, aes(x = density, y = dam)) +
  geom_line(size = 1, aes(linetype = host.group, color = dam.type)) +
  geom_point(data = rawdamdatnon, alpha = 0.5, aes(shape = host.group, color = dam.type)) +
  geom_text(data = filter(damlabels, dens.type == "non-native annual"), aes(x = x, y = y, label = labels), size = 4, fontface = "bold") +
  facet_rep_wrap(~ exp.type, scales = "free_x") +
  scale_linetype_manual(values = c("solid", "dashed"), name = "Host group") +
  scale_shape_manual(values = c(19, 21), name = "Host group") +
  scale_color_manual(values = col_pal[1:2], name = "Metric") +
  xlab(expression(paste("Non-native annual grass density (individuals ", m^-2, ")", sep = ""))) +
  ylab("Disease severity") +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
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
dam.plot <- cowplot::plot_grid(natdamplot, nondamplot, nrow = 2, rel_heights = c(0.95, 1))

pdf("./output/figure5_pathogen_damage_density_march.pdf", width = 5, height = 6)
dam.plot
dev.off()

# pdf("./output/figureS5_pathogen_damage_density_april.pdf", width = 5, height = 6)
# dam.plot
# dev.off()


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
                               damage == 1 & is.na(infection) ~ "damage only",
                               is.na(damage) & infection == 1 ~ "infection only") %>%
           factor(levels = c("infection only", "damage only", "both")),
         native.rel = native.density / total.density,
         year.f = as.factor(year),
         exp.type = case_when(experiment == "transect" ~ "Observational",
                              experiment == "competition" ~ "Manipulated") %>% factor(levels = c("Observational", "Manipulated")))

# labels
denslabels <- tibble(exp.type = c("Observational", "Manipulated", "Observational", "Manipulated"), x = c(0, 0, -0.01, -0.01), y = c(1670, 2720, 1.05, 1.05), plot.type = c("density", "density", "relative", "relative"), labels = c("A", "B", "C", "D")) %>%
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
