#### set up ####

# clear all existing data
rm(list=ls())

# load libraries
library(tidyverse)

# import data
dat <- read_csv("data/manipulated_experiment_plot_map.csv")


#### edit data ####

dat2 <- dat %>%
  mutate(X_adj = case_when(Species == "All" ~ X,
                           TRUE ~ X - 0.5),
         Y_adj = case_when(Species == "All" ~ Y,
                           TRUE ~ Y - 0.5) * -1,
         sp = case_when(Species == "All" ~ "all",
                        TRUE ~ paste(substr(Species, 1, 1), tolower(substr(Species, 2, 2)), sep = "")),
         sp_dens = case_when(sp == "all" ~ sp,
                             TRUE ~ paste(sp, "\n", Density, sep = "")),
         Age = case_when(Species %in% c("SP adult", "EG adults") ~ "adult",
                         TRUE ~ "seedling"),
         plot_size = case_when(sp == "all" ~ "4m2",
                               TRUE ~ "1m2"))


#### figure set-up ####

axisText = 10
axisTitle = 12
legendText = 10
legendTitle = 10

ggtheme <- theme_bw() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size=legendText),
        legend.title = element_text(size=legendTitle),
        legend.position = "bottom",
        legend.direction = "vertical",
        strip.background = element_blank(),
        strip.text = element_text(size = axisText)) 


#### figure ####

pdf("output/figureS2_manipulated_experiment_map.pdf", width = 6.5, height = 6.8)
ggplot(dat2, aes(X_adj, Y_adj)) +
  geom_point(aes(size = plot_size, fill = Treatment, color = Age), shape = 22) +
  geom_text(aes(label = sp_dens), size = 1.7) +
  scale_size_manual(values = c(6, 13), name = "Plot size", labels = c(expression(paste("1", m^2, sep = "")), expression(paste("4", m^2, sep = "")))) +
  scale_fill_manual(values = c("#56B4E9", "#009E73", "#E69F00")) +
  scale_color_manual(values = c("black", "gray")) +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  ggtheme
dev.off()