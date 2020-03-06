## Goal: compare overall damage between the two groups

#### set up ####

# clear all existing data
rm(list=ls())

# load libraries
library(tidyverse)
library(cowplot)
library(glmmTMB)
library(DHARMa) # plot glmmTMB
library(MuMIn) # dredge

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


#### edit raw data ####

# combine experiments
# add severity type
mdat <- mdatT %>%
  select(year:mean.dam, nonnative, origin) %>%
  full_join(mdatC %>% 
              select(year:subplot, host, mean.dam, nonnative, origin)) %>%
  rename(severity = mean.dam) %>%
  mutate(sev.type = "surface")

pdat <- pdatT %>%
  select(year:host, prop.dam, nonnative, origin) %>%
  full_join(pdatC %>% 
              select(year:subplot, host, prop.dam, nonnative, origin)) %>%
  rename(severity = prop.dam) %>%
  mutate(sev.type = "leaves")  

ldat <- ldatT %>%
  select(year:infected, nonnative, origin) %>%
  full_join(ldatC %>% 
              select(year:subplot, host, plant, leaf, infected, nonnative, origin))
  
# combine data for figure
# add life history and experiment type
dat <- full_join(mdat, pdat) %>%
  mutate(status = recode(origin, "non-native" = "non-native\nannual", "native" = "native\nperennial"),
         exp.type = case_when(experiment == "transect" ~ "Observational",
                              TRUE ~ "Manipulated") %>% 
           factor(levels = c("Observational", "Manipulated")))


#### statistical models ####

lmod <- glmmTMB(infected ~ nonnative + (1|year/experiment/subplot/plant), data = ldat, family = binomial)
summary(lmod)
# plot(simulateResiduals(lmod))

mmod <- glmmTMB(severity ~ nonnative + (1|year/experiment/subplot), data = mdat, family = beta_family)
summary(mmod)
# plot(simulateResiduals(mmod))


#### extract model values ####

# prediction function
pred_fun <- function(dat, mod, sev.type){
  
  # predict mean
  dat$pred <- predict(mod, newdata = dat, re.form = NA, type = "response")
  
  # predict SE
  dat$pred.se <- predict(mod, newdata = dat, re.form = NA, se.fit = T, type = "response")$se.fit
  
  # simplify data
  dat_out <- dat %>%
    select(origin, pred, pred.se) %>%
    unique() %>%
    mutate(sev.type = sev.type)
  
  # export
  return(dat_out)
}

# evaluate models
ldatP <- pred_fun(ldat, lmod, "leaves")
mdatP <- pred_fun(mdat, mmod, "surface")

# combine
mod_dat <- rbind(mdatP, ldatP) %>%
  mutate(status = recode(origin, "non-native" = "non-native\nannual", "native" = "native\nperennial")) %>%
  rename(severity = pred)


#### figure ####

# labels
labs <- tibble(let = LETTERS[1:2], 
               sev.type = c("leaves", "surface"), 
               status = rep("native\nperennial", 2))

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
  facet_wrap(~ sev.type, ncol = 1, scales = "free") + 
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
  facet_wrap(~sev.type, ncol = 1, scales = "free") + 
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
pdf("./output/figure2_damage.pdf", width = 2.5, height = 4)
fig
dev.off()


#### percent change ####

# leaves
leaf_nat <- exp(fixef(lmod)$cond[1]) / (1 + exp(fixef(lmod)$cond[1]))
leaf_non <- exp(fixef(lmod)$cond[1] + fixef(lmod)$cond[2]) / (1 + exp(fixef(lmod)$cond[1] + fixef(lmod)$cond[2]))
(leaf_nat - leaf_non) / leaf_non

# surface
surf_nat <- exp(fixef(mmod)$cond[1]) / (1 + exp(fixef(mmod)$cond[1]))
surf_non <- exp(fixef(mmod)$cond[1] + fixef(mmod)$cond[2]) / (1 + exp(fixef(mmod)$cond[1] + fixef(mmod)$cond[2]))
(surf_nat - surf_non) / surf_non
