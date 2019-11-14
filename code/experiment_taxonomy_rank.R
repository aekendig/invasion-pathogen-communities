## Goal: assess taxonomy ranking within experiment


#### set up ####

# clear all existing data
rm(list=ls())

# load libraries
library(tidyverse)

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
fun <- read_csv("./output/taxonomy_species_origin_data.csv")


#### edit data ####

# transect data
tdatT <- idatT %>%
  group_by(otu.id, taxonomy) %>%
  summarise(abundance = length(isolate.id),
            host.species = paste(unique(host), collapse = ", "),
            host.origin = paste(unique(grass.status), collapse = ", ")) %>%
  ungroup() %>%
  mutate(rank = rank(-abundance, ties.method = "random"),
         exp.type = "Observational",
         taxonomy = case_when(taxonomy %in% fdat$taxonomy ~ taxonomy,
                              TRUE ~ ""))

# competition data
tdatC <- idatC %>%
  group_by(otu.id, taxonomy) %>%
  summarise(abundance = length(isolate.id),
            host.species = paste(unique(host), collapse = ", "),
            host.origin = paste(unique(grass.status), collapse = ", ")) %>%
  ungroup() %>%
  mutate(rank = rank(-abundance, ties.method = "random"),
         exp.type = "Manipulated",
         taxonomy = case_when(taxonomy %in% fdat$taxonomy ~ taxonomy,
                              TRUE ~ ""))

# abbreviate pathogen names
otus <- tibble(pathogen = c("ainf", "pcha", "plol", "ptri", "dres", "pave", "rpro"),
               otu.id = c(1, 4, 5, 8, 2, 7, 3),
               path.abb = c("A. inf.", "P. cha.", "P. lol.", "P. tri.", "Dres.", "P. ave.", "R. pro"))

# merge
tdat <- full_join(tdatT, tdatC) %>%
  left_join(otus) %>%
  mutate(exp.type = factor(exp.type, levels = c("Observational", "Manipulated")))


#### visualize ####
pdf("./output/experiment_taxonomy_rank.pdf", width = 5, height = 3)
tdat %>%
  ggplot(aes(x = rank, y = abundance)) +
  geom_point(size = 1.5) +
  geom_text(aes(label = path.abb, color = path.abb), hjust = -0.05, vjust = -0.05, size = 3) +
  facet_wrap(~exp.type) +
  scale_color_manual(values = c(pal, "olivedrab4", "black"), guide = F) +
  xlab("Rank") +
  ylab("Abundance") +
  ggtheme
dev.off()