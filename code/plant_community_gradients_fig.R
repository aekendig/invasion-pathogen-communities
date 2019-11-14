## Goal: invasive and native density and relative abundance


#### set up ####

# clear everything except final data
rm(list = ls())

# load libraries
library(tidyverse)
library(nationalparkcolors)
library(cowplot)

# figure settings
pal <- park_palette("Everglades")[1:5]
axisText = 12
axisTitle = 14
legendText = 12
legendTitle = 12

# import data
dplots <- read_csv("./output/plot_scale_density_data.csv") # plots used in damage analysis
iplots15T <- read_csv("./output/infect_density_transect15_data.csv") %>%
  mutate(exp.year = "observational (2015)") # plots used in infection analysis
iplots16T <- read_csv("./output/infect_density_transect16_data.csv") %>%
  mutate(exp.year = "observational (2016)")
iplots16C <- read_csv("./output/infect_density_competition_data.csv") %>%
  mutate(exp.year = "manipulated (2016)")


#### edit data ####

# simplify infection data
iplots15Tb <- iplots15T %>% dplyr::select(year, experiment, exp.year, plot, subplot, native.density, nonnative.density, total.density)
iplots16Tb <- iplots16T %>% dplyr::select(year, experiment, exp.year, plot, subplot, native.density, nonnative.density, total.density)
iplots16Cb <- iplots16C %>% dplyr::select(year, experiment, exp.year, plot, subplot, native.density, nonnative.density, total.density)

# combine plots from infection data
iplots <- full_join(iplots15Tb, iplots16Tb) %>%
  full_join(iplots16Cb) %>%
  unique() %>%
  mutate(year.f = as.factor(year),
         infection = 1,
         exp.year = factor(exp.year, levels = c("observational (2015)", "observational (2016)", "manipulated (2016)")),
         exp.type = case_when(experiment == "transect" ~ "Observational",
                              experiment == "competition" ~ "Manipulated") %>% factor(levels = c("Observational", "Manipulated")))

# remove scaled densities from dplots
dplots2 <- select(dplots, -c(natdens.s, nondens.s, bg.species, competition.density)) %>%
  mutate(damage = 1,
         exp.year = case_when(experiment == "transect" & year == 2015 ~ "observational (2015)",
                              experiment == "transect" & year == 2016 ~ "observational (2016)",
                              experiment == "competition" & year == 2016 ~ "manipulated (2016)") %>%
           factor(levels = c("observational (2015)", "observational (2016)", "manipulated (2016)")),
         year.f = as.factor(year),
         exp.type = case_when(experiment == "transect" ~ "Observational",
                              experiment == "competition" ~ "Manipulated") %>% factor(levels = c("Observational", "Manipulated")))

# figure out plot uses
anti_join(dplots2, iplots) # 70 only in damage
anti_join(iplots, dplots2) # 11 only in infection

# combine the two plot datasets
plots <- full_join(dplots2, iplots) %>%
  mutate(data.type = case_when(damage == 1 & infection == 1 ~ "both",
                          damage == 1 & is.na(infection) ~ "damage",
                          is.na(damage) & infection == 1 ~ "infection") %>%
           factor(levels = c("damage", "infection", "both")),
         native.rel = native.density / total.density,
         nonnative.rel = nonnative.density / total.density)


#### visualizations ####

ggplot(plots, aes(log(native.density + 1), log(nonnative.density + 1), color = data.type)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~exp.year)

ggplot(plots, aes(native.rel, nonnative.rel, color = data.type)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~exp.year)

densplot <- ggplot(plots, aes(log(native.density + 1), log(nonnative.density + 1), color = data.type, shape = year.f)) +
  geom_point() +
  facet_wrap(~exp.type) +
  scale_color_manual(values = pal[c(4, 5, 2)]) +
  scale_shape_manual(values = c(19, 17)) +
  xlab("Log(native density + 1)") +
  ylab("Log(invasive density + 1)") +
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
  guides(color = "none", shape = "none")
  
relplot <- ggplot(plots, aes(native.rel, nonnative.rel, color = data.type, shape = year.f)) +
  geom_point() +
  facet_wrap(~exp.type) +
  scale_color_manual(values = pal[c(4, 5, 2)], name = "Data type") +
  scale_shape_manual(values = c(19, 17), name = "Year") +
  xlab("Native relative abundance") +
  ylab("Invasive relative abundance") +
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

leg <- get_legend(relplot + theme(legend.background = element_blank(), legend.box.background = element_rect(color = "black")))

pdf("./output/plant_community_gradients.pdf", width = 6, height = 6)
cowplot::plot_grid(densplot, relplot + theme(legend.position = "none"), leg,
          nrow = 3,
          rel_heights = c(1, 1, 0.15))
dev.off()


#### correlation ####

# raw density
cor.test(~ native.density + nonnative.density, data = plots) # slightly negative, not sig
cor.test(~ native.density + nonnative.density, data = filter(plots, exp.type == "Observational")) # slightly negative, not sig
cor.test(~ native.density + nonnative.density, data = filter(plots, exp.type == "Manipulated"))  # slightly negative, not sig
         
# log-transformed density
cor.test(~ log(native.density + 1) + log(nonnative.density + 1), data = plots) # slightly negative, not sig
cor.test(~ log(native.density + 1) + log(nonnative.density + 1), data = filter(plots, exp.type == "Observational")) # slightly positive, not sig
cor.test(~ log(native.density + 1) + log(nonnative.density + 1), data = filter(plots, exp.type == "Manipulated")) # slightly negative and sig at p = 0.04
