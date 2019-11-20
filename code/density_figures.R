## Goal: manuscrpit figures

# update of figures.R


#### set up ####

# clear all existing data
rm(list=ls())

# load libraries
library(sjPlot)
library(tidyverse)
library(cowplot)

# figure settings
pal <- park_palette("Everglades")[1:5]

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
idatT <- read_csv("./output/infect_density_experiment_transect_data.csv")
idatC <- read_csv("./output/infect_density_experiment_competition_data.csv")

mdatT <- read_csv("./output/damage_density_experiment_meandam_transect16_data.csv")
mdatC <- read_csv("./output/damage_density_experiment_meandam_competition_data.csv")

ldatT <- read_csv("./output/damage_density_experiment_propdam_transect_data.csv")
ldatC <- read_csv("./output/damage_density_experiment_propdam_competition_data.csv")

pdatT <- read_csv("./output/damage_density_experiment_plant_transect_data.csv")
pdatC <- read_csv("./output/damage_density_experiment_plant_competition_data.csv")

fun <- read_csv("./output/taxonomy_species_origin_data.csv")

# import models
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

load("./output/damage_density_experiment_meandam_absolute_transect_avg_amodel.rda")
load("./output/damage_density_experiment_meandam_absolute_competition_avg_amodel.rda")
load("./output/damage_density_experiment_propdam_absolute_transect_avg_amodel.rda")
load("./output/damage_density_experiment_propdam_absolute_competition_avg_amodel.rda")

# functions
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

# function for extracting coefficients from pathogen models
coef_fun <- function(mod, pathogen, exp.type){
  out <- tibble(coef = names(mod$coefficients["full",]), value = mod$coefficients["full",]) %>% 
    mutate(pathogen = pathogen,
           exp.type = exp.type)
}

# function for extracting coefficients from damage models
dam_coef_fun <- function(mod, dam.type, exp.type){
  out <- tibble(coef = names(mod$coefficients["full",]), value = mod$coefficients["full",]) %>% 
    mutate(dam.type = dam.type,
           exp.type = exp.type)
}


#### table of pathogens ####

# combine possible host species across years
fun2 <- fun %>%
  group_by(pathogen, otu.id) %>%
  summarise(observed.host.species = toString(host.species)) %>%
  ungroup() %>%
  mutate(observed.host.species = case_when(otu.id == 1 ~ "Ab, Af, Bd, Bh, Eg, Fp, Pa, Sp",
                                           otu.id == 2 ~ "Ab, Bd, Bh, Eg, Sp",
                                           otu.id == 7 ~ "Bd, Bh, Eg, Sp",
                                           otu.id == 4 ~ "Ab, Af, Bd, Bh, Eg, Fp, Sp",
                                           otu.id == 5 ~ "Ab, Af, Bd, Bh, Eg, Fp, Sp",
                                           otu.id == 3 ~ "Ab, Bd, Bh, Eg, Fp, Pa, Sp",
                                           TRUE ~ observed.host.species)) %>% 
  rename(taxonomy = pathogen)

# otu's and code names
otus <- tibble(pathogen = c("ainf", "pcha", "plol", "ptri", "dres", "pave", "rpro"),
               otu.id = c(1, 4, 5, 8, 2, 7, 3),
               path.abb = c("A. inf.", "P. cha.", "P. lol.", "P. tri.", "Dres.", "P. ave.", "R. pro")) %>%
  left_join(select(fun2, otu.id, taxonomy, observed.host.species))


# edit data
prevdat <- idatT %>%
  select(experiment, origin, ainf, pcha, plol, ptri, dres, pave, rpro) %>%
  full_join(idatC %>% select(experiment, origin, ainf, pcha, plol, ptri, dres, pave, rpro)) %>%
  gather(key = pathogen, value = present, -c(experiment, origin)) %>%
  left_join(otus) %>%
  mutate(status = recode(origin, "non-native" = "exotic") %>% factor(levels = c("native", "exotic")),
         exp.type = case_when(experiment == "transect" ~ "Observational",
                              experiment == "competition" ~ "Manipulated") %>% factor(levels = c("Observational", "Manipulated"))) 

# make table
fun_tab <- prevdat %>%
  group_by(taxonomy, observed.host.species, status) %>%
  summarise(otus = sum(present)) %>%
  spread(key = status, value = otus) %>%
  rename(native.abundance = native, exotic.abundance = exotic)

tab_df(fun_tab)


#### scaled density plots ####

# Note that we can't merge the data using scaled values from each analysis because the same plot may have a different value in the damage and infection experiment. New scaled values are created.

# simplify datasets
iplotsT <- idatT %>%
  select(experiment, year, plot, subplot, nonnative.density, native.density, total.density, nonnative.rel) %>%
  unique() %>%
  mutate(infection = 1)

dplotsT <- pdatT %>%
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

# figure out plot uses
anti_join(dplotsC, iplotsC) # 52 only in damage
anti_join(iplotsC, dplotsC) # none only in infection

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

# figures
densplot <- ggplot(plots, aes(natdens.s, nondens.s, color = data.type, shape = year.f)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~exp.type) +
  scale_color_manual(values = pal[c(5, 2, 4)]) +
  scale_shape_manual(values = c(19, 17)) +
  xlab("Scaled native host density") +
  ylab("Scaled invasive\nhost density") +
  ggtheme +
  guides(color = "none", shape = "none")

relplot <- ggplot(plots, aes(native.rel, nonnative.rel, color = data.type, shape = year.f)) +
  geom_point(alpha = 0.6, position = position_jitter(0.02, 0.02, seed = 2)) +
  facet_wrap(~exp.type) +
  scale_color_manual(values = pal[c(5, 2, 4)], name = "Data type") +  
  scale_shape_manual(values = c(19, 17), name = "Year") +
  xlab("Native host relative abundance") +
  ylab("Invasive host\nrelative abundance") +
  ggtheme +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal")

leg <- get_legend(relplot + theme(legend.background = element_blank(), legend.box.background = element_rect(color = "black")))

pdf("./output/plant_community_gradients.pdf", width = 6, height = 6)
cowplot::plot_grid(densplot, relplot + theme(legend.position = "none"), leg,
                   nrow = 3,
                   rel_heights = c(1, 1, 0.15))
dev.off()


#### infection prevalence ####

# figure
pdf("./output/pathogen_species_prevalence_host_status.pdf", width = 6, height = 6)
ggplot(prevdat, aes(x = taxonomy, y = present, group = status, fill = status)) +
  scale_fill_manual(values = pal[c(4, 2)], name = "Host plant status") +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1, position = position_dodge(0.3), color = "black") +
  stat_summary(geom = "point", fun.y = "mean", size = 3, position = position_dodge(0.3), shape = 21) +
  facet_wrap(~ exp.type, nrow = 2, strip.position = "right") +
  ggtheme +
  xlab("Pathogen") +
  ylab("Pathogen community prevalence") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "italic"),
        legend.position = c(0.65, 0.93),
        legend.direction = "horizontal")
dev.off()
  

#### summarized density effects ####

# extract coefficients
acoef <- rbind(coef_fun(aiamodTa, "ainf", "Observational"), coef_fun(pcamodTa, "pcha", "Observational"), coef_fun(plamodTa, "plol", "Observational"), coef_fun(ptamodTa, "ptri", "Observational"), coef_fun(dramodTa, "dres", "Observational"), coef_fun(paamodTa, "pave", "Observational"), coef_fun(aiamodCa, "ainf", "Manipulated"), coef_fun(pcamodCa, "pcha", "Manipulated"), coef_fun(plamodCa, "plol", "Manipulated"), coef_fun(ptamodCa, "ptri", "Manipulated"), coef_fun(dramodCa, "dres", "Manipulated"), coef_fun(paamodCa, "pave", "Manipulated")) %>%
  spread(key = coef, value = value) %>%
  replace(is.na(.), 0) %>%
  mutate(nat_dens_nat_prev = logit2prob(`cond(natdens.s)`),
         nat_dens_inv_prev = logit2prob(`cond(natdens.s)` + `cond(natdens.s:nonnative)`),
         inv_dens_nat_prev = logit2prob(`cond(nondens.s)`),
         inv_dens_inv_prev = logit2prob(`cond(nondens.s)` + `cond(nondens.s:nonnative)`)) %>%
  select(exp.type, pathogen, nat_dens_nat_prev:inv_dens_inv_prev) %>%
  gather(key = dens_status, value = coef, -c(pathogen, exp.type)) %>%
  mutate(status = substring(dens_status, 10, 12) %>% recode(nat = "Native hosts", inv = "Invasive hosts") %>% factor(levels = c("Native hosts", "Invasive hosts")),
         density = substring(dens_status, 1, 3) %>% recode(nat = "Native grass density", inv = "Invasive grass density") %>% factor(levels = c("Native grass density", "Invasive grass density")),
         exp.type = factor(exp.type, levels = c("Observational", "Manipulated"))) %>%
  left_join(otus)

# remove specialists (origin is not in the model)
acoef2 <- filter(acoef, !(pathogen == "pcha" & status == "Native hosts") & !(pathogen == "ptri" & status == "Invasive hosts"))

# figure
ggplot(acoef2, aes(x = density, y = coef, group = status)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_point(aes(shape = path.abb, color = status), position = position_jitterdodge(dodge.width = 0.2, jitter.width = 0.1, seed = 3)) +
  facet_wrap(~ exp.type, strip.position = "right", nrow = 2) +
  scale_color_manual(values = pal[c(4, 2)], name = "Host plant status") +
  scale_shape_manual(values = c(17, 18, 15, 8, 10, 7), name = "Pathogen") +
  ggtheme +
  theme(axis.title.x = element_blank())


#### individual density effects ####

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

# number of individuals to increase by
ind <- 1000

# density table
inddens <- tibble(natdens.s = c(-natmuT/natsdT, (ind-natmuT)/natsdT, -natmuC/natsdC, (ind-natmuC)/natsdC),
                nondens.s = c(-nonmuT/nonsdT, (ind-nonmuT)/nonsdT, -nonmuC/nonsdC, (ind-nonmuC)/nonsdC),
                native.density = rep(c(0, ind), 2),
                nonnative.density = rep(c(0, ind), 2),
                exp.type = c(rep("Observational", 2), rep("Manipulated", 2)))

# check values
filter(idatC, native.density == min(native.density)) %>% select(native.density, natdens.s) %>% unique()
filter(idatT, native.density == min(native.density)) %>% select(native.density, natdens.s) %>% unique()
filter(idatC, nonnative.density == min(nonnative.density)) %>% select(nonnative.density, nondens.s) %>% unique()
filter(idatT, nonnative.density == min(nonnative.density)) %>% select(nonnative.density, nondens.s) %>% unique() # seem similar

# extract coefficients for absolute abundance
aind <- rbind(coef_fun(aiamodTa, "ainf", "Observational"), coef_fun(pcamodTa, "pcha", "Observational"), coef_fun(plamodTa, "plol", "Observational"), coef_fun(ptamodTa, "ptri", "Observational"), coef_fun(dramodTa, "dres", "Observational"), coef_fun(paamodTa, "pave", "Observational"), coef_fun(aiamodCa, "ainf", "Manipulated"), coef_fun(pcamodCa, "pcha", "Manipulated"), coef_fun(plamodCa, "plol", "Manipulated"), coef_fun(ptamodCa, "ptri", "Manipulated"), coef_fun(dramodCa, "dres", "Manipulated"), coef_fun(paamodCa, "pave", "Manipulated")) %>%
  spread(key = coef, value = value) %>%
  replace(is.na(.), 0) %>%
  full_join(inddens) %>%
  mutate(nat_dens_nat_prev = logit2prob(`cond((Int))` + `cond(natdens.s)` * natdens.s),
         nat_dens_inv_prev = logit2prob(`cond((Int))`+ `cond(nonnative)` + `cond(natdens.s)` * natdens.s + `cond(natdens.s:nonnative)` * natdens.s),
         inv_dens_nat_prev = logit2prob(`cond((Int))` + `cond(nondens.s)` * nondens.s),
         inv_dens_inv_prev = logit2prob(`cond((Int))`+ `cond(nonnative)` + `cond(nondens.s)` * nondens.s + `cond(nondens.s:nonnative)` * nondens.s)) %>%
  select(exp.type, pathogen, native.density, nonnative.density, nat_dens_nat_prev:inv_dens_inv_prev) %>%
  gather(key = dens_status, value = prev, -c(pathogen, exp.type, native.density, nonnative.density)) %>%
  mutate(status = substring(dens_status, 10, 12) %>% recode(nat = "Native hosts", inv = "Invasive hosts") %>% factor(levels = c("Native hosts", "Invasive hosts")),
         density = substring(dens_status, 1, 3),
         exp.type = factor(exp.type, levels = c("Observational", "Manipulated"))) %>%
  left_join(otus)

# split by density and calculate change in prevalence
aindnat <- filter(aind, density == "nat") %>% 
  select(-nonnative.density) %>%
  spread(key = native.density, value = prev) %>%
  mutate(prev_change = `1000` - `0`)

aindinv <- filter(aind, density == "inv") %>% 
  select(-native.density) %>%
  spread(key = nonnative.density, value = prev) %>%
  mutate(prev_change = `1000` - `0`)

# merge back together and remove specialists 
aind2 <- full_join(aindnat, aindinv) %>%
  mutate(density = recode(density, nat = "Native grass density", inv = "Invasive grass density")) %>%
  filter(!(pathogen == "pcha" & status == "Native hosts") & !(pathogen == "ptri" & status == "Invasive hosts"))

# figure
pdf("./output/prevalence_change_density_1000.pdf", width = 6, height = 4)
ggplot(aind2, aes(x = density, y = prev_change, fill = path.abb, shape = status, group = status)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_point(position = position_jitterdodge(dodge.width = 0.4, jitter.width = 0.1, seed = 8), size = 4) +
  facet_wrap(~ exp.type, strip.position = "right", nrow = 2) +
  ylab(expression(paste("Change in prevalence with 1000 individuals ", m^-2, sep = ""))) +
  scale_fill_manual(values = c(pal, "olivedrab4"), name = "Pathogen") +
  scale_shape_manual(values = c(22, 25), name = "Host plant status") +
  ggtheme +
  theme(axis.title.x = element_blank()) +
  guides(fill = guide_legend(override.aes = list(shape=21)))
dev.off()


#### simulated density effects ####

# density table
idens <- tibble(natdens.s = c(seq(min(idatT$natdens.s), max(idatT$natdens.s), length.out = 100), seq(min(idatC$natdens.s), max(idatC$natdens.s), length.out = 100)),
                nondens.s = c(seq(min(idatT$nondens.s), max(idatT$nondens.s), length.out = 100), seq(min(idatC$nondens.s), max(idatC$nondens.s), length.out = 100)),
                native.density = c(seq(min(idatT$native.density), max(idatT$native.density), length.out = 100), seq(min(idatC$native.density), max(idatC$native.density), length.out = 100)),
                nonnative.density = c(seq(min(idatT$nonnative.density), max(idatT$nonnative.density), length.out = 100), seq(min(idatC$nonnative.density), max(idatC$nonnative.density), length.out = 100)),
                exp.type = c(rep("Observational", 100), rep("Manipulated", 100)))

# extract coefficients for absolute abundance
asim <- rbind(coef_fun(aiamodTa, "ainf", "Observational"), coef_fun(pcamodTa, "pcha", "Observational"), coef_fun(plamodTa, "plol", "Observational"), coef_fun(ptamodTa, "ptri", "Observational"), coef_fun(dramodTa, "dres", "Observational"), coef_fun(paamodTa, "pave", "Observational"), coef_fun(aiamodCa, "ainf", "Manipulated"), coef_fun(pcamodCa, "pcha", "Manipulated"), coef_fun(plamodCa, "plol", "Manipulated"), coef_fun(ptamodCa, "ptri", "Manipulated"), coef_fun(dramodCa, "dres", "Manipulated"), coef_fun(paamodCa, "pave", "Manipulated")) %>%
  spread(key = coef, value = value) %>%
  replace(is.na(.), 0) %>%
  full_join(idens) %>%
  mutate(nat_prev = logit2prob(`cond((Int))`),
         inv_prev = logit2prob(`cond((Int))` + `cond(nonnative)`),
         nat_dens_nat_prev = logit2prob(`cond((Int))` + `cond(natdens.s)` * natdens.s),
         nat_dens_inv_prev = logit2prob(`cond((Int))`+ `cond(nonnative)` + `cond(natdens.s)` * natdens.s + `cond(natdens.s:nonnative)` * natdens.s),
         inv_dens_nat_prev = logit2prob(`cond((Int))` + `cond(nondens.s)` * nondens.s),
         inv_dens_inv_prev = logit2prob(`cond((Int))`+ `cond(nonnative)` + `cond(nondens.s)` * nondens.s + `cond(nondens.s:nonnative)` * nondens.s)) %>%
  select(exp.type, pathogen, native.density, nonnative.density, nat_dens_nat_prev:inv_dens_inv_prev) %>%
  gather(key = dens_status, value = prev, -c(pathogen, exp.type, native.density, nonnative.density)) %>%
  mutate(status = substring(dens_status, 10, 12) %>% recode(nat = "Native hosts", inv = "Invasive hosts") %>% factor(levels = c("Native hosts", "Invasive hosts")),
         density = substring(dens_status, 1, 3),
         exp.type = factor(exp.type, levels = c("Observational", "Manipulated"))) %>%
  left_join(otus) %>%
  filter(!(pathogen == "pcha" & status == "Native hosts") & !(pathogen == "ptri" & status == "Invasive hosts"))

# split by density
asimnat <- filter(asim, density == "nat")
asiminv <- filter(asim, density == "inv")

# figures
pdf("./output/native_density_prevalence.pdf", width = 6, height = 4)
ggplot(asimnat, aes(x = native.density, y = prev, color = path.abb)) +
  geom_line(size = 1.5) +
  facet_grid(status ~ exp.type, scales = "free_x") +
  scale_color_manual(values = c(pal, "olivedrab4"), name = "Pathogen") +
  xlab(expression(paste("Native grass density (individuals ", m^-2, ")", sep = ""))) +
  ylab("Pathogen community prevalence") +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  ggtheme
dev.off()

pdf("./output/invasive_density_prevalence.pdf", width = 6, height = 4)
ggplot(asiminv, aes(x = nonnative.density, y = prev, color = path.abb)) +
  geom_line(size = 1.5) +
  facet_grid(status ~ exp.type, scales = "free_x") +
  scale_color_manual(values = c(pal, "olivedrab4"), name = "Pathogen") +
  xlab(expression(paste("Invasive grass density (individuals ", m^-2, ")", sep = ""))) +
  ylab("Pathogen community prevalence") +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  ggtheme
dev.off()


#### simulated effects on damage ####

# extract coefficients for absolute abundance
dsim <- rbind(dam_coef_fun(mamodTa, "mean", "Observational"), dam_coef_fun(mamodCa, "mean", "Manipulated"), dam_coef_fun(pamodTa, "proportion", "Observational"), dam_coef_fun(pamodCa, "proportion", "Manipulated")) %>%
  spread(key = coef, value = value) %>%
  replace(is.na(.), 0) %>%
  full_join(idens) %>%
  mutate(nat_dens_nat_dam = logit2prob(`cond((Int))` + `cond(natdens.s)` * natdens.s),
         nat_dens_inv_dam = logit2prob(`cond((Int))`+ `cond(nonnative)` + `cond(natdens.s)` * natdens.s + `cond(natdens.s:nonnative)` * natdens.s),
         inv_dens_nat_dam = logit2prob(`cond((Int))` + `cond(nondens.s)` * nondens.s),
         inv_dens_inv_dam = logit2prob(`cond((Int))`+ `cond(nonnative)` + `cond(nondens.s)` * nondens.s + `cond(nondens.s:nonnative)` * nondens.s)) %>%
  select(exp.type, dam.type, native.density, nonnative.density, nat_dens_nat_dam:inv_dens_inv_dam) %>%
  gather(key = dens_status, value = dam, -c(dam.type, exp.type, native.density, nonnative.density)) %>%
  mutate(status = substring(dens_status, 10, 12) %>% recode(nat = "Native hosts", inv = "Invasive hosts") %>% factor(levels = c("Native hosts", "Invasive hosts")),
         density = substring(dens_status, 1, 3),
         exp.type = factor(exp.type, levels = c("Observational", "Manipulated")))

# split by density
dsimnat <- filter(dsim, density == "nat")
dsiminv <- filter(dsim, density == "inv")

# figures
pdf("./output/native_density_damage.pdf", width = 6, height = 4)
ggplot(dsimnat, aes(x = native.density, y = dam, linetype = dam.type)) +
  geom_line(size = 1.5) +
  facet_grid(status ~ exp.type, scales = "free_x") +
  scale_linetype_manual(values = c("dotted", "solid"), name = "Damage type") +
  xlab(expression(paste("Native grass density (individuals ", m^-2, ")", sep = ""))) +
  ylab("Damage") +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  ggtheme
dev.off()

pdf("./output/invasive_density_damage.pdf", width = 6, height = 4)
ggplot(dsiminv, aes(x = nonnative.density, y = dam, linetype = dam.type)) +
  geom_line(size = 1.5) +
  facet_grid(status ~ exp.type, scales = "free_x") +
  scale_linetype_manual(values = c("dotted", "solid"), name = "Damage type") +
  xlab(expression(paste("Invasive grass density (individuals ", m^-2, ")", sep = ""))) +
  ylab("Damage") +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  ggtheme
dev.off()