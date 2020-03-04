## Goal: compare overall damage between the two groups

#### set up ####

# clear all existing data
rm(list=ls())

# load libraries
library(tidyverse)
library(cowplot)

# import data

# plant-scale data without zeros
mdatT <- read_csv("./output/damage_density_experiment_meandam_transect_data.csv")
mdatC <- read_csv("./output/damage_density_experiment_meandam_competition_data.csv")

# leaf-scale data without zeros
ldatT <- read_csv("./output/damage_density_experiment_propdam_transect_data.csv")
ldatC <- read_csv("./output/damage_density_experiment_propdam_competition_data.csv")

# plant-scale data with zeros
pdatT <- read_csv("./output/damage_density_experiment_plant_transect_data.csv")
pdatC <- read_csv("./output/damage_density_experiment_plant_competition_data.csv")

# models
load("./output/damage_density_experiment_meandam_absolute_transect_avg_amodel.rda")
load("./output/damage_density_experiment_meandam_absolute_competition_avg_amodel.rda")
load("./output/damage_density_experiment_propdam_absolute_transect_avg_amodel.rda")
load("./output/damage_density_experiment_propdam_absolute_competition_avg_amodel.rda")


#### edit raw data ####

# combine experiments
# add severity type
mdat <- mdatT %>%
  select(year:mean.dam, origin) %>%
  full_join(mdatC %>% 
              select(year:subplot, host, mean.dam, origin)) %>%
  rename(severity = mean.dam) %>%
  mutate(sev.type = "surface")

pdat <- pdatT %>%
  select(year:host, prop.dam, origin) %>%
  full_join(pdatC %>% 
              select(year:subplot, host, prop.dam, origin)) %>%
  rename(severity = prop.dam) %>%
  mutate(sev.type = "leaves")  
  
# combine data
# add life history and experiment type
dat <- full_join(mdat, pdat) %>%
  mutate(status = recode(origin, "non-native" = "non-native\nannual", "native" = "native\nperennial"),
         exp.type = case_when(experiment == "transect" ~ "Observational",
                              TRUE ~ "Manipulated") %>% 
           factor(levels = c("Observational", "Manipulated")),
         sev.type = fct_relevel(sev.type, "surface"))


#### extract model values ####

# prediction function
pred_fun <- function(dat, mod, exp.type, sev.type){
  
  # change the density values to average
  dat_new <- dat %>%
    mutate(natdens.s = 0,
           nondens.s = 0,
           othdens.s = 0)
  
  # predict mean
  dat_new$pred <- predict(mod, newdata = dat_new, re.form = NA, type = "response")
  
  # predict SE
  dat_new$pred.se <- predict(mod, newdata = dat_new, re.form = NA, se.fit = T, type = "response")$se.fit
  
  # simplify data
  dat_out <- dat_new %>%
    select(origin, pred, pred.se) %>%
    unique() %>%
    mutate(exp.type = exp.type,
           sev.type = sev.type)
  
  # export
  return(dat_out)
}

# evaluate models
mdatTM <- pred_fun(mdatT, mamodTa, "Observational", "surface")
mdatCM <- pred_fun(mdatC, mamodCa, "Manipulated", "surface")
ldatTM <- pred_fun(ldatT, pamodTa, "Observational", "leaves")
ldatCM <- pred_fun(ldatC, pamodCa, "Manipulated", "leaves")

# combine
mod_dat <- rbind(mdatTM, mdatCM, ldatTM, ldatCM) %>%
  mutate(status = recode(origin, "non-native" = "non-native\nannual", "native" = "native\nperennial"),
         exp.type = fct_relevel(exp.type, "Observational"),
         sev.type = fct_relevel(sev.type, "surface")) %>%
  rename(severity = pred)


#### figure ####

# labels
labs <- tibble(let = LETTERS[1:4], 
               sev.type = c("surface", "surface", "leaves", "leaves") %>% fct_relevel("surface"), 
               exp.type = c("Observational", "Manipulated", "Observational", "Manipulated") %>% fct_relevel("Observational"), 
               status = rep("native\nperennial", 4))

# figure settings
axisText = 10
axisTitle = 12
legendText = 10
legendTitle = 10

# figure
fig_raw <- ggplot(mod_dat, aes(x = status, y = severity, fill = status)) +
  geom_point(data = dat, alpha = 0.5, size = 0.5, aes(color = status)) +
  geom_errorbar(aes(ymin = severity - 2*pred.se, ymax = severity + 2*pred.se), width = 0.1) +
  geom_point(size = 3, shape = 21) +
  scale_fill_manual(values = c("black", "white"), guide = F) +
  scale_color_manual(values = c("gray80", "gray80"), guide = F) +
  facet_grid(sev.type ~ exp.type, scales = "free") + 
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
  ylab("Disease severity")

fig <- ggplot(mod_dat, aes(x = status, y = severity, fill = status)) +
  geom_errorbar(aes(ymin = severity - 2*pred.se, ymax = severity + 2*pred.se), width = 0.1) +
  geom_point(size = 3, shape = 21) +
  scale_fill_manual(values = c("black", "white"), guide = F) +
  facet_grid(sev.type ~ exp.type, scales = "free") + 
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
  ylab("Disease severity")

# save
pdf("./output/figure2_damage.pdf", width = 4, height = 4)
fig
dev.off()