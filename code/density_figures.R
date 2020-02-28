## Goal: manuscrpit figures

# update of figures.R


#### set up ####

# clear all existing data
rm(list=ls())

# load libraries
library(sjPlot)
library(tidyverse)
library(cowplot)
library(nationalparkcolors)
library(grid)
library(gridExtra)
library(lemon)
library(MuMIn)

# figure settings
pal <- park_palette("Everglades")[1:5]
col_pal <- c("#000000", "#56B4E9", "#009E73", "#F0E442", "#000000", "#56B4E9", "#009E73")

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

bgsp <- read_csv("./output/background_species_table.csv")

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

# function for insets in plots with panels (https://stackoverflow.com/questions/32807665/removing-one-tablegrob-when-applied-to-a-box-plot-with-a-facet-wrap?answertab=votes#tab-top)
# didn't end up using, but saving for future - just make a ggplot of the inset and save it as a grob, add it to the main plot with this function
annotation_custom2 <- 
  function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data){
    layer(data = data, stat = StatIdentity, position = PositionIdentity, 
          geom = ggplot2:::GeomCustomAnn,
          params = list(grob = grob,
                        xmin = xmin, xmax = xmax,
                        ymin = ymin, ymax = ymax))
  }


#### table of background species ####
bgsp2 <- bgsp %>%
  mutate(grass_group = recode(grass_group, non.ann = "non-native annual", nat.per = "native perennial"),
         experiment = recode(experiment, competition = "manipulated", transect = "observational")) %>%
  rename(study = experiment) %>%
  filter(abundance > 0) %>%
  arrange(desc(study), year, species)

tab_df(bgsp2)


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
                                           otu.id == 3 ~ "Ab, Bd, Bh, Eg, Sp",
                                           TRUE ~ observed.host.species)) %>% 
  rename(taxonomy = pathogen)

# otu's and code names
otus <- tibble(pathogen = c("ainf", "pcha", "plol", "ptri", "dres", "pave", "rpro"),
               otu.id = c(1, 4, 5, 8, 2, 7, 3),
               path.abb = c("A. inf.", "P. cha.", "P. lol.", "P. tri.", "Drec.", "P. ave.", "R. pro.")) %>%
  left_join(select(fun2, otu.id, taxonomy, observed.host.species, host.num))


# edit data
prevdat <- idatT %>%
  select(experiment, grass.status, ainf, pcha, plol, ptri, dres, pave, rpro) %>%
  full_join(idatC %>% select(experiment, grass.status, ainf, pcha, plol, ptri, dres, pave, rpro)) %>%
  gather(key = pathogen, value = present, -c(experiment, grass.status)) %>%
  left_join(otus) %>%
  mutate(status = recode(grass.status, "non-native" = "non-native\nannual", "native" = "native\nperennial"),
         exp.type = case_when(experiment == "transect" ~ "Observational",
                              experiment == "competition" ~ "Manipulated") %>% factor(levels = c("Observational", "Manipulated"))) 

# make table
fun_tab <- prevdat %>%
  group_by(taxonomy, observed.host.species, host.num, status) %>%
  summarise(otus = sum(present),
            prop = round(otus/length(present), 2)) %>%
  mutate(out = paste(otus, prop, sep = " (")) %>%
  select(-c(otus, prop)) %>%
  spread(key = status, value = out) 

tab_df(fun_tab)


#### sample sizes ####

# transect isolates
idatT %>%
  group_by(year, grass_group) %>%
  summarise(n = length(unique(isolate.id)))

# competition isolates
idatC %>%
  group_by(grass_group) %>%
  summarise(n = length(unique(isolate.id)))

# competition plots
idatC %>%
  group_by(bg.species, competition.density, competition.type) %>%
  summarise(n = length(unique(subplot)))

# transect design
idatT %>%
  group_by(plot) %>%
  summarise(n = length(unique(subplot)))

idatT %>%
  group_by(year) %>%
  summarise(n = length(unique(subplot)))

idatT %>%
  group_by(year, grass.status) %>%
  summarise(n = length(unique(isolate.id)))

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


#### density plots ####

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

# figures for presentations
scaledplot <- ggplot(plots, aes(natdens.s, nondens.s, color = data.type, shape = year.f)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~exp.type) +
  scale_color_manual(values = pal[c(5, 2, 4)]) +
  scale_shape_manual(values = c(19, 17)) +
  xlab("Scaled native perennial grass density") +
  ylab("Scaled non-native annual\ngrass density") +
  ggtheme +
  guides(color = "none", shape = "none")

relplot <- ggplot(plots, aes(native.rel, nonnative.rel, color = data.type, shape = year.f)) +
  geom_point(alpha = 0.6, position = position_jitter(0.02, 0.02, seed = 2)) +
  facet_wrap(~exp.type) +
  scale_color_manual(values = pal[c(5, 2, 4)], name = "Data type") +  
  scale_shape_manual(values = c(19, 17), name = "Year") +
  xlab("Native perennial grass relative abundance") +
  ylab("Non-native annual\ngrass relative abundance") +
  ggtheme +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal")

leg <- get_legend(relplot + theme(legend.background = element_blank(), legend.box.background = element_rect(color = "black")))

pdf("./output/plant_community_gradients.pdf", width = 6, height = 6)
cowplot::plot_grid(scaledplot, relplot + theme(legend.position = "none"), leg,
                   nrow = 3,
                   rel_heights = c(1, 1, 0.15))
dev.off()

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

pdf("./output/figureS2_plant_community_density.pdf", width = 7, height = 7)
cowplot::plot_grid(densplot, relplot2 + theme(legend.position = "none"), leg2,
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
acoef <- rbind(coef_fun(rpamodTa, "rpro", "Observational"), coef_fun(aiamodTa, "ainf", "Observational"), coef_fun(pcamodTa, "pcha", "Observational"), coef_fun(plamodTa, "plol", "Observational"), coef_fun(ptamodTa, "ptri", "Observational"), coef_fun(dramodTa, "dres", "Observational"), coef_fun(paamodTa, "pave", "Observational"), coef_fun(aiamodCa, "ainf", "Manipulated"), coef_fun(pcamodCa, "pcha", "Manipulated"), coef_fun(plamodCa, "plol", "Manipulated"), coef_fun(ptamodCa, "ptri", "Manipulated"), coef_fun(dramodCa, "dres", "Manipulated"), coef_fun(paamodCa, "pave", "Manipulated")) %>%
  spread(key = coef, value = value) %>%
  replace(is.na(.), 0) %>%
  mutate(nat_dens_nat_prev = logit2prob(`cond(natdens.s)`),
         nat_dens_non_prev = logit2prob(`cond(natdens.s)` + `cond(natdens.s:nonnative)`),
         non_dens_nat_prev = logit2prob(`cond(nondens.s)`),
         non_dens_non_prev = logit2prob(`cond(nondens.s)` + `cond(nondens.s:nonnative)`)) %>%
  select(exp.type, pathogen, nat_dens_nat_prev:non_dens_non_prev) %>%
  gather(key = dens_status, value = coef, -c(pathogen, exp.type)) %>%
  mutate(status = substring(dens_status, 10, 12) %>% recode(nat = "Native perennial hosts", non = "Non-native annual hosts") %>% factor(levels = c("Native perennial hosts", "Non-native annual hosts")),
         density = substring(dens_status, 1, 3) %>% recode(nat = "Native perennial grass density", non = "Non-native annual grass density") %>% factor(levels = c("Native perennial grass density", "Non-native annual grass density")),
         exp.type = factor(exp.type, levels = c("Observational", "Manipulated"))) %>%
  left_join(otus)

# remove specialists (origin is not in the model)
acoef2 <- filter(acoef, !(pathogen == "pcha" & status == "Native perennial hosts") & !(pathogen == "ptri" & status == "Non-native annual hosts"))

# figure
ggplot(acoef2, aes(x = density, y = coef, group = status)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_point(aes(shape = path.abb, color = status), position = position_jitterdodge(dodge.width = 0.2, jitter.width = 0.1, seed = 3)) +
  facet_wrap(~ exp.type, strip.position = "right", nrow = 2) +
  scale_color_manual(values = pal[c(4, 2)], name = "Host plant status") +
  scale_shape_manual(values = c(17, 18, 15, 8, 10, 7, 5), name = "Pathogen") +
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

# number of perennial individuals to increase by
ind_p <- 50

# number of annual individuals
ind_a <- 5000

# density table
inddens <- tibble(natdens.s = c(-natmuT/natsdT, (ind_p-natmuT)/natsdT, -natmuC/natsdC, (ind_p-natmuC)/natsdC),
                nondens.s = c(-nonmuT/nonsdT, (ind_a-nonmuT)/nonsdT, -nonmuC/nonsdC, (ind_a-nonmuC)/nonsdC),
                native.density = rep(c(0, ind_p), 2),
                nonnative.density = rep(c(0, ind_a), 2),
                exp.type = c(rep("Observational", 2), rep("Manipulated", 2)))

# check values
filter(idatC, native.density == min(native.density)) %>% select(native.density, natdens.s) %>% unique()
filter(idatT, native.density == min(native.density)) %>% select(native.density, natdens.s) %>% unique()
filter(idatC, nonnative.density == min(nonnative.density)) %>% select(nonnative.density, nondens.s) %>% unique()
filter(idatT, nonnative.density == min(nonnative.density)) %>% select(nonnative.density, nondens.s) %>% unique() # seem similar

# extract coefficients for absolute abundance
aind <- rbind(coef_fun(rpamodTa, "rpro", "Observational"), coef_fun(aiamodTa, "ainf", "Observational"), coef_fun(pcamodTa, "pcha", "Observational"), coef_fun(plamodTa, "plol", "Observational"), coef_fun(ptamodTa, "ptri", "Observational"), coef_fun(dramodTa, "dres", "Observational"), coef_fun(paamodTa, "pave", "Observational"), coef_fun(aiamodCa, "ainf", "Manipulated"), coef_fun(pcamodCa, "pcha", "Manipulated"), coef_fun(plamodCa, "plol", "Manipulated"), coef_fun(ptamodCa, "ptri", "Manipulated"), coef_fun(dramodCa, "dres", "Manipulated"), coef_fun(paamodCa, "pave", "Manipulated")) %>%
  spread(key = coef, value = value) %>%
  replace(is.na(.), 0) %>%
  full_join(inddens) %>%
  mutate(nat_dens_nat_prev = logit2prob(`cond((Int))` + `cond(natdens.s)` * natdens.s),
         nat_dens_non_prev = logit2prob(`cond((Int))`+ `cond(nonnative)` + `cond(natdens.s)` * natdens.s + `cond(natdens.s:nonnative)` * natdens.s),
         non_dens_nat_prev = logit2prob(`cond((Int))` + `cond(nondens.s)` * nondens.s),
         non_dens_non_prev = logit2prob(`cond((Int))`+ `cond(nonnative)` + `cond(nondens.s)` * nondens.s + `cond(nondens.s:nonnative)` * nondens.s)) %>%
  select(exp.type, pathogen, native.density, nonnative.density, nat_dens_nat_prev:non_dens_non_prev) %>%
  gather(key = dens_status, value = prev, -c(pathogen, exp.type, native.density, nonnative.density)) %>%
  mutate(status = substring(dens_status, 10, 12) %>% recode(nat = "Native perennial\nhosts", non = "Non-native\nannual hosts") %>% factor(levels = c("Native perennial\nhosts", "Non-native\nannual hosts")),
         density = substring(dens_status, 1, 3),
         exp.type = factor(exp.type, levels = c("Observational", "Manipulated"))) %>%
  left_join(otus)

# split by density and calculate change in prevalence, remove specialists
aindnat <- filter(aind, density == "nat") %>% 
  select(-nonnative.density) %>%
  spread(key = native.density, value = prev) %>%
  mutate(prev_change = `50` - `0`) %>%
  filter(!(pathogen == "pcha" & status == "Native perennial\nhosts") & !(pathogen == "ptri" & status == "Non-native\nannual hosts"))

aindinv <- filter(aind, density == "non") %>% 
  select(-native.density) %>%
  spread(key = nonnative.density, value = prev) %>%
  mutate(prev_change = `5000` - `0`) %>%
  filter(!(pathogen == "pcha" & status == "Native perennial\nhosts") & !(pathogen == "ptri" & status == "Non-native\nannual hosts"))

# labels
indlabels <- tibble(dens.type = rep(c("native perennial", "non-native annual"), each = 2), exp.type = c("Observational", "Manipulated", "Observational", "Manipulated"), x = rep(0.5, 4), y = rep(0.93, 4), labels = c("A", "C", "B", "D")) %>%
  mutate(exp.type = factor(exp.type, levels = c("Observational", "Manipulated")))

# figure
changeplotnat <- ggplot(aindnat, aes(x = status, y = prev_change)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_point(aes(fill = path.abb, shape = path.abb), position = position_jitterdodge(dodge.width = 0.1, jitter.width = 0.1, seed = 8), size = 3) +
  geom_text(data = filter(indlabels, dens.type == "native perennial"), aes(x = x, y = y, label = labels), size = 4, fontface = "bold") +
  facet_wrap(~ exp.type, strip.position = "right", nrow = 2) +
  ylab(expression(paste("Response to 50 native perennial grasses ", m^-2, sep = ""))) +
  scale_fill_manual(values = col_pal, name = "Pathogen") +
  scale_shape_manual(values = rep(c(21, 24), each = 4), name = "Pathogen") +
  ggtheme +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        strip.text = element_blank()) +
  ylim(-0.55, 0.95)

changeplotnon <- ggplot(aindinv, aes(x = status, y = prev_change)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_point(aes(fill = path.abb, shape = path.abb), position = position_jitterdodge(dodge.width = 0.1, jitter.width = 0.1, seed = 8), size = 3) +
  geom_text(data = filter(indlabels, dens.type == "non-native annual"), aes(x = x, y = y, label = labels), size = 4, fontface = "bold") +
  facet_wrap(~ exp.type, strip.position = "right", nrow = 2) +
  ylab(expression(paste("Response to 5000 non-native annual grasses ", m^-2, sep = ""))) +
  scale_fill_manual(values = col_pal, name = "Pathogen") +
  scale_shape_manual(values = rep(c(21, 24), each = 4), name = "Pathogen") +
  ggtheme +
  theme(axis.title.x = element_blank()) +
  ylim(-0.55, 0.95)

pdf("./output/figure3_prevalence_change_density.pdf", width = 7, height = 4)
cowplot::plot_grid(changeplotnat, changeplotnon, rel_widths = c(0.7, 1))
dev.off()


#### simulated density effects ####

# raw data of focal pathogens
idatFnat <- idatT %>%
  select(isolate.id, experiment, grass.status, native.density, ainf, pcha, plol, ptri, dres, pave, rpro) %>%
  full_join(idatC %>% select(isolate.id, experiment, grass.status, native.density, ainf, pcha, plol, ptri, dres, pave, rpro)) %>%
  group_by(experiment, grass.status, native.density) %>%
  summarise(ainf = mean(ainf), pcha = mean(pcha), plol = mean(plol), ptri = mean(ptri), dres = mean(dres), pave = mean(pave), rpro = mean(rpro)) %>%
  ungroup() %>%
  gather(key = pathogen, value = prev, -c(experiment:native.density)) %>%
  left_join(otus) %>%
  mutate(status = recode(grass.status, "non-native" = "non-native annual", "native" = "native perennial"),
         exp.type = case_when(experiment == "transect" ~ "Observational",
                              experiment == "competition" ~ "Manipulated") %>% factor(levels = c("Observational", "Manipulated"))) %>%
  rename(density = native.density)

idatFnon <- idatT %>%
  select(isolate.id, experiment, grass.status, nonnative.density, ainf, pcha, plol, ptri, dres, pave, rpro) %>%
  full_join(idatC %>% select(isolate.id, experiment, grass.status, nonnative.density, ainf, pcha, plol, ptri, dres, pave, rpro)) %>%
  group_by(experiment, grass.status, nonnative.density) %>%
  summarise(ainf = mean(ainf), pcha = mean(pcha), plol = mean(plol), ptri = mean(ptri), dres = mean(dres), pave = mean(pave), rpro = mean(rpro)) %>%
  ungroup() %>%
  gather(key = pathogen, value = prev, -c(experiment:nonnative.density)) %>%
  left_join(otus) %>%
  mutate(status = recode(grass.status, "non-native" = "non-native annual", "native" = "native perennial"),
         exp.type = case_when(experiment == "transect" ~ "Observational",
                              experiment == "competition" ~ "Manipulated") %>% factor(levels = c("Observational", "Manipulated"))) %>%
  rename(density = nonnative.density)

# idat by grass.status
idatT.nat = filter(idatT, grass.status == "native")
idatT.non = filter(idatT, grass.status == "non-native")
idatC.nat = filter(idatC, grass.status == "native")
idatC.non = filter(idatC, grass.status == "non-native")

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

# extract coefficients for absolute abundance
asim <- rbind(coef_fun(rpamodTa, "rpro", "Observational"), coef_fun(aiamodTa, "ainf", "Observational"), coef_fun(pcamodTa, "pcha", "Observational"), coef_fun(plamodTa, "plol", "Observational"), coef_fun(ptamodTa, "ptri", "Observational"), coef_fun(dramodTa, "dres", "Observational"), coef_fun(paamodTa, "pave", "Observational"), coef_fun(aiamodCa, "ainf", "Manipulated"), coef_fun(pcamodCa, "pcha", "Manipulated"), coef_fun(plamodCa, "plol", "Manipulated"), coef_fun(ptamodCa, "ptri", "Manipulated"), coef_fun(dramodCa, "dres", "Manipulated"), coef_fun(paamodCa, "pave", "Manipulated")) %>%
  spread(key = coef, value = value) %>%
  replace(is.na(.), 0) %>%
  full_join(idens) %>%
  mutate(nat_dens_nat_prev = logit2prob(`cond((Int))` + `cond(natdens.s)` * nat.natdens.s),
         nat_dens_non_prev = logit2prob(`cond((Int))`+ `cond(nonnative)` + `cond(natdens.s)` * non.natdens.s + `cond(natdens.s:nonnative)` * non.natdens.s),
         non_dens_nat_prev = logit2prob(`cond((Int))` + `cond(nondens.s)` * nat.nondens.s),
         non_dens_non_prev = logit2prob(`cond((Int))`+ `cond(nonnative)` + `cond(nondens.s)` * non.nondens.s + `cond(nondens.s:nonnative)` * non.nondens.s)) %>%
  select(exp.type, pathogen, nat.native.density, non.native.density, nat.nonnative.density, non.nonnative.density, nat_dens_nat_prev:non_dens_non_prev) %>%
  gather(key = dens_status, value = prev, -c(pathogen, exp.type, nat.native.density, non.native.density, nat.nonnative.density, non.nonnative.density)) %>%
  mutate(status = substring(dens_status, 10, 12) %>% recode(nat = "native perennial", non = "non-native annual"),
         dens.type = substring(dens_status, 1, 3),
         exp.type = factor(exp.type, levels = c("Observational", "Manipulated"))) %>%
  left_join(otus) %>%
  filter(!(pathogen == "pcha" & status == "native perennial") & !(pathogen == "ptri" & status == "non-native annual"))

# split by density
asimnat <- filter(asim, dens.type == "nat") %>%
  mutate(native.density = case_when(status == "native perennial" ~ nat.native.density,
                                    status == "non-native annual" ~ non.native.density)) %>%
  rename(density = native.density)

asimnon <- filter(asim, dens.type == "non") %>%
  mutate(nonnative.density = case_when(status == "native perennial" ~ nat.nonnative.density,
                                       status == "non-native annual" ~ non.nonnative.density)) %>%
  rename(density = nonnative.density)

# manuscript figure function
dens_fig_fun <- function(rawdat, preddat, stat, exp){
  
  rawdat2 <- filter(rawdat, status == stat & exp.type == exp)
  preddat2 <- filter(preddat, status == stat & exp.type == exp)
  
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
obs.nat.natdens.plot <- dens_fig_fun(idatFnat, asimnat, "native perennial", "Observational")
obs.non.natdens.plot <- dens_fig_fun(idatFnat, asimnat, "non-native annual", "Observational")
obs.nat.nondens.plot <- dens_fig_fun(idatFnon, asimnon, "native perennial", "Observational")
obs.non.nondens.plot <- dens_fig_fun(idatFnon, asimnon, "non-native annual", "Observational")
man.nat.natdens.plot <- dens_fig_fun(idatFnat, asimnat, "native perennial", "Manipulated")
man.non.natdens.plot <- dens_fig_fun(idatFnat, asimnat, "non-native annual", "Manipulated")
man.nat.nondens.plot <- dens_fig_fun(idatFnon, asimnon, "native perennial", "Manipulated")
man.non.nondens.plot <- dens_fig_fun(idatFnon, asimnon, "non-native annual", "Manipulated")

# max values
filter(idatFnat, status == "native perennial" & exp.type == "Manipulated") %>% select(density) %>% max()
filter(idatFnat, status == "non-native annual" & exp.type == "Manipulated") %>% select(density) %>% max()
filter(idatFnon, status == "native perennial" & exp.type == "Manipulated") %>% select(density) %>% max()
filter(idatFnon, status == "non-native annual" & exp.type == "Manipulated") %>% select(density) %>% max()

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

pdf("./output/figure2_pathogen_relative_abundance_density.pdf", width = 7, height = 10.5)
grid.arrange(arrangeGrob(dens.plot.fin, left = textGrob("Pathogen relative abundance", gp=gpar(fontsize = 13), rot = 90)))
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

summary(mamodTa)
summary(pamodTa)
summary(mamodCa)
summary(pamodCa)

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

importance(mamodTa)
importance(pamodTa)
importance(mamodCa)
importance(pamodCa)



#### simulated effects on damage ####

# use plant-scale data because it's the same dataset as that used for mean.dam, but it includes the zeros
pdatT.nat <- filter(pdatT, origin == "native")
pdatT.non <- filter(pdatT, origin == "non-native")
pdatC.nat <- filter(pdatC, origin == "native")
pdatC.non <- filter(pdatC, origin == "non-native")

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

# extract coefficients for absolute abundance
dsim <- rbind(dam_coef_fun(mamodTa, "surface", "Observational"), dam_coef_fun(mamodCa, "surface", "Manipulated"), dam_coef_fun(pamodTa, "leaves", "Observational"), dam_coef_fun(pamodCa, "leaves", "Manipulated")) %>%
  spread(key = coef, value = value) %>%
  replace(is.na(.), 0) %>%
  full_join(ddens) %>%
  mutate(nat_dens_nat_dam = logit2prob(`cond((Int))` + `cond(natdens.s)` * nat.natdens.s),
         nat_dens_non_dam = logit2prob(`cond((Int))`+ `cond(nonnative)` + `cond(natdens.s)` * non.natdens.s + `cond(natdens.s:nonnative)` * non.natdens.s),
         non_dens_nat_dam = logit2prob(`cond((Int))` + `cond(nondens.s)` * nat.nondens.s),
         non_dens_non_dam = logit2prob(`cond((Int))`+ `cond(nonnative)` + `cond(nondens.s)` * non.nondens.s + `cond(nondens.s:nonnative)` * non.nondens.s)) %>%
  select(exp.type, dam.type, nat.native.density, nat.nonnative.density, non.native.density, non.nonnative.density, nat_dens_nat_dam:non_dens_non_dam) %>%
  gather(key = dens_status, value = dam, -c(dam.type, exp.type, nat.native.density:non.nonnative.density)) %>%
  mutate(status = substring(dens_status, 10, 12) %>% recode(nat = "native perennial", non = "non-native annual"),
         dens.type = substring(dens_status, 1, 3),
         exp.type = factor(exp.type, levels = c("Observational", "Manipulated")))

# split by density
dsimnat <- filter(dsim, dens.type == "nat")  %>%
  mutate(native.density = case_when(status == "native perennial" ~ nat.native.density,
                                    status == "non-native annual" ~ non.native.density)) %>%
  rename(density = native.density)

dsimnon <- filter(dsim, dens.type == "non")  %>%
  mutate(nonnative.density = case_when(status == "native perennial" ~ nat.nonnative.density,
                                       status == "non-native annual" ~ non.nonnative.density)) %>%
  rename(density = nonnative.density)

# raw data: surface area
surfdat <- mdatT %>%
  select(experiment, origin, native.density, nonnative.density, mean.dam) %>%
  full_join(mdatC %>%
              select(experiment, origin, native.density, nonnative.density, mean.dam)) %>%
    rename(dam = mean.dam) %>%
    mutate(exp.type = case_when(experiment == "transect" ~ "Observational",
                                experiment == "competition" ~ "Manipulated") %>%
             factor(levels = c("Observational", "Manipulated")),
           status = recode(origin, "non-native" = "non-native annual", "native" = "native perennial"),
           dam.type = "surface")

# raw data: leaves
leafdat <- pdatT %>%
  select(experiment, origin, native.density, nonnative.density, prop.dam) %>%
  full_join(pdatC %>%
              select(experiment, origin, native.density, nonnative.density, prop.dam)) %>%
  rename(dam = prop.dam) %>%
  mutate(exp.type = case_when(experiment == "transect" ~ "Observational",
                              experiment == "competition" ~ "Manipulated") %>%
           factor(levels = c("Observational", "Manipulated")),
         status = recode(origin, "non-native" = "non-native annual", "native" = "native perennial"),
         dam.type = "leaves")

# split by density
rawdamdatnat <- full_join(surfdat, leafdat) %>%
  select(-nonnative.density) %>%
  rename(density = native.density)

rawdamdatinv <- full_join(surfdat, leafdat) %>%
  select(-native.density) %>%
  rename(density = nonnative.density)

# labels
damlabels <- tibble(dens.type = rep(c("native perennial", "non-native annual"), each = 2), exp.type = c("Observational", "Manipulated", "Observational", "Manipulated"), x = rep(0.5, 4), y = rep(1.1, 4), labels = c("A", "B", "C", "D")) %>%
  mutate(exp.type = factor(exp.type, levels = c("Observational", "Manipulated")))

# native density plot 
natdamplot <- ggplot(dsimnat, aes(x = density, y = dam)) +
  geom_line(size = 1, aes(linetype = status,color = dam.type)) +
  geom_point(data = rawdamdatnat, alpha = 0.5, aes(shape = status, color = dam.type)) +
  geom_text(data = filter(damlabels, dens.type == "native"), aes(x = x, y = y, label = labels), size = 4, fontface = "bold") +
  facet_rep_wrap(~ exp.type, scales = "free_x") +
  scale_linetype_manual(values = c("solid", "dashed"), name = "Grass group") +
  scale_shape_manual(values = c(19, 21), name = "Grass group") +
  scale_color_manual(values = col_pal[1:2], name = "Metric") +
  xlab(expression(paste("Native perennial grass density (individuals ", m^-2, ")", sep = ""))) +
  ylab("Disease severity") +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  ggtheme +
  theme(legend.position = "none")

# exotic density plot  
nondamplot <- ggplot(dsimnon, aes(x = density, y = dam)) +
  geom_line(size = 1, aes(linetype = status, color = dam.type)) +
  geom_point(data = rawdamdatinv, alpha = 0.5, aes(shape = status, color = dam.type)) +
  geom_text(data = filter(damlabels, dens.type == "exotic"), aes(x = x, y = y, label = labels), size = 4, fontface = "bold") +
  facet_rep_wrap(~ exp.type, scales = "free_x") +
  scale_linetype_manual(values = c("solid", "dashed"), name = "Grass group") +
  scale_shape_manual(values = c(19, 21), name = "Grass group") +
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
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, -10, -10, -10),
        plot.margin = margin(3, 3, 3, 3)) +
  guides(shape = guide_legend(override.aes = list(alpha = 1)), 
         color = guide_legend(override.aes = list(alpha = 1)))

# combine plots
dam.plot <- cowplot::plot_grid(natdamplot, nondamplot, nrow = 2, rel_heights = c(0.96, 1))

pdf("./output/figure4_pathogen_damage_density.pdf", width = 7, height = 7)
dam.plot
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

# P lol with native dens, observational
plol_man_natdens_0 <- exp(plamodTa$coefficients["full",1]) / (1 + exp(plamodTa$coefficients["full",1]))
plol_man_natdens_1 <- exp(plamodTa$coefficients["full",1] + plamodTa$coefficients["full",2]) / (1 + exp(plamodTa$coefficients["full",1] + plamodTa$coefficients["full",2]))
(plol_man_natdens_1 - plol_man_natdens_0) / plol_man_natdens_0

# R pro with native dens, observational
rpro_man_natdens_0 <- exp(rpamodTa$coefficients["full",1]) / (1 + exp(rpamodTa$coefficients["full",1]))
rpro_man_natdens_1 <- exp(rpamodTa$coefficients["full",1] + rpamodTa$coefficients["full",2]) / (1 + exp(rpamodTa$coefficients["full",1] + rpamodTa$coefficients["full",2]))
(rpro_man_natdens_1 - rpro_man_natdens_0) / rpro_man_natdens_0

# proportion of leaves, observational
prop_obs_nat <- exp(pamodTa$coefficients["full",1]) / (1 + exp(pamodTa$coefficients["full",1]))
prop_obs_non <- exp(pamodTa$coefficients["full",1] + pamodTa$coefficients["full",2]) / (1 + exp(pamodTa$coefficients["full",1] + pamodTa$coefficients["full",2]))
(prop_obs_nat - prop_obs_non) / prop_obs_non

# proportion of leaves, manipulated
prop_man_nat <- exp(pamodCa$coefficients["full",1]) / (1 + exp(pamodCa$coefficients["full",1]))
prop_man_non <- exp(pamodCa$coefficients["full",1] + pamodCa$coefficients["full",2]) / (1 + exp(pamodCa$coefficients["full",1] + pamodCa$coefficients["full",2]))
(prop_man_nat - prop_man_non) / prop_man_non

# proportion with native dens, observational
prop_obs_natdens_0 <- exp(pamodTa$coefficients["full",1]) / (1 + exp(pamodTa$coefficients["full",1]))
prop_obs_natdens_1 <- exp(pamodTa$coefficients["full",1] + pamodTa$coefficients["full",3]) / (1 + exp(pamodTa$coefficients["full",1] + pamodTa$coefficients["full",3]))
(prop_obs_natdens_1 - prop_obs_natdens_0) / prop_obs_natdens_0

# proportion with native dens on non-native, observational
prop_obs_non_natdens_0 <- exp(pamodTa$coefficients["full",1] + pamodTa$coefficients["full",2]) / (1 + exp(pamodTa$coefficients["full",1] + pamodTa$coefficients["full",2]))
prop_obs_non_natdens_1 <- exp(pamodTa$coefficients["full",1] + pamodTa$coefficients["full",2] + pamodTa$coefficients["full",3] + pamodTa$coefficients["full",6]) / (1 + exp(pamodTa$coefficients["full",1] + pamodTa$coefficients["full",2] + pamodTa$coefficients["full",3] + pamodTa$coefficients["full",6]))
(prop_obs_non_natdens_1 - prop_obs_non_natdens_0) / prop_obs_non_natdens_0
