## Goal: describe differences between pathogens associated with native and non-native plants

# notes: updated version is species_origin_community_differences.R


#### set up ####

# clear everything except final data
rm(list = ls())

# load libraries
library(tidyverse)
library(SpadeR)

# import data
dat <- read.csv("./data/fungal_pathogens_2015_2017.csv")


#### edit data ####

# remove JEF to make numbers between native and non-native more comparable
dat1 <- dat %>%
  filter(experiment != "JEF transect")


#### isolates ####

# by origin
dat %>% 
  group_by(grass.status) %>%
  summarise(n = n(),
            isolates = length(unique(isolate.id)))

# by origin without JEF (SP only)
dat1 %>% 
  group_by(grass.status) %>%
  summarise(isolates = length(unique(isolate.id)))

# by origin and year
dat1 %>%
  group_by(grass.status, year) %>%
  summarise(isolates = length(unique(isolate.id)))
# similar numbers


#### richness ####

# richness data
rdat <- dat1 %>%
  group_by(year, year.f, grass.status, otu.id) %>%
  summarise(abundance = length(isolate.id)) %>%
  ungroup() %>%
  group_by(year, year.f, grass.status) %>%
  summarise(isolates = sum(abundance),
            observed = length(otu.id),
            estimated = Diversity(abundance, datatype="abundance")$Species_richness[3,"Estimate"],
            lower = Diversity(abundance, datatype="abundance")$Species_richness[3,"95%Lower"],
            upper = Diversity(abundance, datatype="abundance")$Species_richness[3,"95%Upper"]) %>%
  gather(key = "type", value = "richness", -c(year, year.f, grass.status, isolates, lower, upper)) %>%
  mutate(lower = ifelse(type == "observed", NA, lower),
         upper = ifelse(type == "observed", NA, upper))

# effect of isolates
rdat %>%
  filter(type == "estimated") %>%
  ggplot(aes(x = isolates, y = richness, color = grass.status, shape = year.f)) +
  geom_point(size = 3) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 2) # not strong

# comparison by year
rdat %>%
  filter(type == "estimated") %>%
  ggplot(aes(x = year.f, y = richness, color = grass.status)) +
  geom_point(size = 3, position = position_dodge(0.5)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, position = position_dodge(0.5))

# cause of high errorbar?
dat1 %>%
  group_by(year, year.f, grass.status, otu.id) %>%
  summarise(abundance = length(isolate.id)) %>%
  ggplot(aes(x = abundance)) +
  geom_histogram() + 
  facet_grid(year.f ~ grass.status)
# two very high abundance pathogens for both native and non-native in year 2

# comparison by year without 2
rdat %>%
  filter(year != 2016 & type == "estimated") %>%
  ggplot(aes(x = year.f, y = richness, color = grass.status)) +
  geom_point(size = 3, position = position_dodge(0.5)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, position = position_dodge(0.5))
# similar numbers

# plot with estimated and observed
rdat %>%
  ggplot(aes(x = year.f, y = richness, group = grass.status, color = grass.status, shape = type)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, position = position_dodge(0.5), color = "black") +
  geom_point(aes(y = richness), size = 2, position = position_dodge(0.5)) +
  scale_shape_manual(values = c(19, 22))

# log-transform y-axis to make it easier to view
rdat %>%
  ggplot(aes(x = year.f, y = log(richness), group = grass.status, color = grass.status, shape = type)) +
  geom_errorbar(aes(ymin = log(lower), ymax = log(upper)), width = 0.1, position = position_dodge(0.5), color = "black") +
  geom_point(size = 3, position = position_dodge(0.5)) +
  scale_shape_manual(values = c(19, 22)) +
  theme(axis.title.x = element_blank())


#### diversity ####

# shannon diversity data
sdat <- dat1 %>%
  group_by(year, year.f, grass.status, otu.id) %>%
  summarise(abundance = length(isolate.id)) %>%
  ungroup() %>%
  group_by(year, year.f, grass.status) %>%
  summarise(isolates = sum(abundance),
            diversity = Diversity(abundance, datatype="abundance")$Shannon_diversity[4,"Estimate"],
            lower = Diversity(abundance, datatype="abundance")$Shannon_diversity[4,"95%Lower"],
            upper = Diversity(abundance, datatype="abundance")$Shannon_diversity[4,"95%Upper"])

# comparison by year
sdat %>%
  ggplot(aes(x = year.f, y = diversity, group = grass.status, color = grass.status)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, position = position_dodge(0.5), color = "black") +
  geom_point(size = 3, position = position_dodge(0.5), shape = 19) +
  theme(axis.title.x = element_blank()) + 
  geom_point(aes(y = isolates), size = 3, position = position_dodge(0.5), shape = 23)

# comparison by year without isolates
sdat %>%
  ggplot(aes(x = year.f, y = diversity, group = grass.status, color = grass.status)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, position = position_dodge(0.5), color = "black") +
  geom_point(size = 3, position = position_dodge(0.5), shape = 19) +
  theme(axis.title.x = element_blank())


#### fungal species ####

# make sure otu.id and taxonomy are 1:1
dat1 %>% 
  group_by(otu.id) %>% 
  summarise(taxa = length(unique(taxonomy))) %>%
  filter(taxa > 1)

dat1 %>% 
  group_by(taxonomy) %>% 
  summarise(otus = length(unique(otu.id))) %>%
  filter(otus > 1) # no - 14 taxa are assigned to multiple OTUs (usually genus level, not always)

sum(is.na(dat1$otu.id))
sum(is.na(dat1$taxanomy))

# taxonomy data
tdat <- dat1 %>%
  group_by(year, year.f, grass.status, otu.id, taxonomy) %>%
  summarise(abundance = length(isolate.id)) %>%
  ungroup() %>%
  group_by(year, year.f, grass.status) %>%
  mutate(rank = rank(-abundance, ties.method = "random")) %>%
  ungroup()

# visualize
tdat %>%
  ggplot(aes(x = rank, y = abundance, color = grass.status)) +
  geom_point() +
  geom_line() +
  facet_wrap(~year.f) +
  geom_text(aes(label = taxonomy), size = 3, hjust = 0)


#### outputs ####
write_csv(rdat, "./output/richness_origin_data.csv")
write_csv(sdat, "./output/diversity_origin_data.csv")
write_csv(tdat, "./output/taxonomy_origin_data.csv")
