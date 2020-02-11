## Goal: determine whether plant species within the same origin group have more similar pathogen communities than between groups

# notes: this is an updated version of origin_community_differences.R


#### set up ####

# clear everything except final data
rm(list = ls())

# load libraries
library(tidyverse)
library(vegan)
library(nationalparkcolors)
library(sjPlot)
library(rusda)
library(mgsub)
library(cowplot)
library(ggrepel)

# figure settings
pal <- c(park_palette("Everglades")[1:5], "black", "mediumpurple3", "olivedrab4")
shapes <- c(0, 1, 2, 5, 13, 19, 17)
axisText = 12
axisTitle = 14
legendText = 12
legendTitle = 12

# import data
dat <- read.csv("./data/fungal_pathogens_2015_2017.csv")


#### edit data ####

# remove JEF to make numbers between native and non-native more comparable
dat1 <- dat %>%
  filter(experiment != "JEF transect" & 
           host %in% c("AB", "AF", "BRAD", "BD", "BH", "EG", "SP")) %>%
  mutate(grass_group = recode(grass.status, native = "native\nperennial", "non-native" = "non-native\nannual"),
         host_sp = recode(host, AB = "Avena barbata", AF = "Avena fatua", BRAD = "Brachypodium distachyonm", BD = "Bromus diandrus", BH = "Bromus hordeaceus", EG = "Elymus glaucus", SP = "Stipa pulchra"),
           host_sp = fct_relevel(host_sp, "Avena barbata", "Avena fatua", "Brachypodium distachyonm"))

# sample sizes
nrow(dat1)
unique(dat1$experiment)
unique(dat1$host)

# set minimum isolate number
min.iso <- 3

# create wide data
datw1 <- dat1 %>%
  group_by(year, year.f, host_sp, grass_group, otu.id) %>%
  summarise(abundance = length(isolate.id)) %>%
  group_by(year, year.f, host_sp, grass_group) %>%
  mutate(isolates = sum(abundance)) %>%
  filter(isolates > min.iso) %>%
  ungroup() %>%
  spread(otu.id, abundance)

# community matrix
cdat1 <- datw1 %>%
  select(-c(year:isolates)) %>%
  data.frame()
cdat1[is.na(cdat1)] = 0
cdat1 <- as.matrix(cdat1)

# environmental matrix
edat1 <- datw1 %>% select(c(year:grass_group))

# subset full data for those included in these matrices
dat2 <- edat1 %>%
  left_join(dat1)
# same number of observations - can use dat1

# sample sizes by group
dat2 %>%
  group_by(year.f, grass_group) %>%
  summarise(n = n())

# sample sizes by species 
(sp_table <- dat2 %>%
    group_by(year, grass_group, host_sp) %>%
    summarise(n = n()) %>%
    arrange(year, grass_group, host_sp))

write_csv(sp_table, "./output/host_species_sample_sizes.csv")


#### host specificity analysis ####

# pathogen data
pdat <- dat2 %>%
  group_by(grass_group, otu.id, taxonomy) %>%
  summarise(abundance = length(isolate.id)) %>%
  mutate(pathogen = mgsub(as.character(taxonomy), c(" cf.", " sp.", "unknown "), c("", "", "")),
         species = case_when(sapply(strsplit(pathogen, " "), length) == 2 ~ TRUE,
                             TRUE ~ FALSE)) %>%
  ungroup()

# look at non species
filter(pdat, species == F) %>%
  select(pathogen) %>%
  unique() %>%
  data.frame()
# genera

# number of unique pathogens
paths <- pdat %>%
  group_by(pathogen, species) %>%
  summarise(abundance = sum(abundance)) %>%
  unique()
# 58 unique

sum(paths$species)
# 29 species

paths %>%
  group_by(species) %>%
  summarise(abu = sum(abundance))

# pathogens with specificity data - couldn't get it to work within the case_when statment
# had to remove pathogens without records
sdat <- pdat %>%
  filter(species == T & !(pathogen %in% c("Parastagonospora cumpignensis", "Preussia lignicola", "Gilmaniella subornata", "Cystostereum heteromorphum", "Preussia minipascua", "Ochrocladosporium frigidarii"))) %>%
  select(pathogen) %>%
  unique() %>%
  mutate(host.num = NA)

# add specificity data to pathogens with loop (errors with mutate)
# for(i in 1:nrow(sdat)){
#  hosts <- length(unique(associations(sdat$pathogen[i], database = "FH", spec_type = "fungus")$associations$host))
#  sdat$host.num[i] <- hosts
#}

# save file
# write_csv(sdat3, "./output/pathogen_host_range_usda.csv")

# re-import file and replace host ranges for species not in the database as zeros
sdat2 <- read_csv("./output/pathogen_host_range_usda.csv") %>%
  full_join(paths %>%
              select(pathogen, species)) %>%
  mutate(host.num = case_when(is.na(host.num) & species == T ~ 0,
                              TRUE ~ host.num))

# check species with one host in database
filter(sdat2, host.num == 1)
# associations("Ramularia proteae", database = "FH", spec_type = "fungus") # not a host we have
# associations("Anthostomella proteae", database = "FH", spec_type = "fungus") # not a host we have
# associations("Cladosporium antarcticum", database = "FH", spec_type = "fungus") # not a host we have

# total abundance by group
gdat <- dat2 %>%
  group_by(grass_group) %>%
  summarise(n = n())

# merge specificity data back with pathogen data
pdat2 <- full_join(pdat, sdat2) %>%
  mutate(prop.abundance = case_when(grass_group == "native\nperennial" ~ abundance/gdat$n[1],
                                    grass_group == "non-native\nannual"  ~ abundance/gdat$n[2]),
         host.num.scaled = host.num * prop.abundance)
  
# figure
pdat2 %>%
  ggplot(aes(x = grass_group, y = host.num.scaled)) +
  stat_summary(geom = "point", fun.y = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se")

# merge specificity data with full dataset
sdat3 <- pdat2 %>%
  select(taxonomy, pathogen, host.num, species) %>%
  unique() %>%
  right_join(dat2) %>%
  mutate(species.known = as.numeric(species))

# summarise by using whole dataset
hostplot <- sdat3 %>%
  ggplot(aes(x = grass_group, y = host.num, fill = grass_group)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.1) +
  stat_summary(geom = "point", fun.y = "mean", size = 5, shape = 21) +
  scale_fill_manual(values = c("black", "white"), guide = F) +
  ylab("Host range of associated pathogens") +
  xlab("Host group") +
  theme_bw() +
  theme(axis.text = element_text(size = 10, color="black"),
        axis.title = element_text(size = 12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
# only takes pathogens identified to species into account

# proportion with species known
sdat3 %>%
  ggplot(aes(x = grass.status, y = species.known)) +
  stat_summary(geom = "point", fun.y = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot")

# total samples with species known
sdat3 %>%
  group_by(grass.status) %>%
  summarise(n = n(),
            database = sum(!is.na(host.num))) %>%
  mutate(prop = database / n)

# t-test
smod <- t.test(host.num ~ grass_group, data = filter(sdat3, !is.na(host.num)))
smod


#### output taxonomy table and figure ####

# abundance by grass status
adat <- sdat3 %>%
  group_by(year, otu.id, taxonomy, pathogen, grass_group) %>%
  summarise(abundance = length(isolate.id)) %>%
  spread(key = grass_group, value = abundance) %>%
  ungroup() %>%
  rename("native.abundance" = "native\nperennial", "nonnative.abundance" = "non-native\nannual")

# taxonomy data
tdat <- sdat3 %>%
  group_by(year, year.f, otu.id, taxonomy, pathogen, host.num) %>%
  summarise(abundance = length(isolate.id),
            host.species = paste(unique(host), collapse = ", ")) %>%
  ungroup() %>%
  group_by(year, year.f) %>%
  mutate(rank = rank(-abundance, ties.method = "random")) %>%
  ungroup() %>%
  left_join(adat) %>%
  mutate(native.abundance = replace_na(native.abundance, 0),
         nonnative.abundance = replace_na(nonnative.abundance, 0))

# look at all fungi
tdat %>%
  ggplot(aes(x = rank, y = abundance)) +
  geom_point(size = 2) +
  facet_wrap(~year)

# function to calulate break from next highest rank
break_fun <- function(dat){
  
  n = max(dat$rank)
  diff = c(0, rep(NA, n - 1))
  val = c(1:n)
  
  for(i in 2:n){
    diff[i] = filter(dat, rank == val[i-1])$abundance - filter(dat, rank == val[i])$abundance
  }
  
  out <- data.frame(rank = val, breaks = diff)
  return(out)
}

# calculate breaks (mutate didn't work)
tdat2015 <- filter(tdat, year == 2015)
tdat2015 <- break_fun(tdat2015) %>% full_join(tdat2015)
filter(tdat2015, breaks == max(breaks)) # rank 5 to 6

tdat2016 <- filter(tdat, year == 2016)
tdat2016 <- break_fun(tdat2016) %>% full_join(tdat2016)
# the max one is clearly too early in rank
sort(tdat2016$breaks)
filter(tdat2016, breaks == 15) # ranks 2 to 3
filter(tdat2016, breaks == 12) # ranks 5 to 6

tdat2017 <- filter(tdat, year == 2017)
tdat2017 <- break_fun(tdat2017) %>% full_join(tdat2017)
# the max one is clearly too early in rank
sort(tdat2017$breaks)
filter(tdat2017, breaks == 10) # ranks 4 to 5

# break points
breakdat <- tibble(year = c(2015, 2016, 2017), rank.break = c(5.5, 5.5, 4.5))

# extract top fungal species
fdat <- tdat %>%
  left_join(breakdat) %>%
  filter(rank < rank.break) %>%
  select(pathogen, year, host.species, abundance, native.abundance, nonnative.abundance, rank, host.num, otu.id) %>%
  mutate(pathogen = case_when(pathogen == "Drechslera" ~ "Drechslera sp.",
                              TRUE ~ pathogen),
         host.species = mgsub(host.species, c("AB", "BH", "SP", "BD", "EG", "AF"), c("Ab", "Bh", "Sp", "Bd", "Eg", "Af"))) %>%
  arrange(pathogen, year, rank)

# number of fungal species
length(unique(fdat$otu.id)) # 7

# save data
write_csv(fdat, "./output/taxonomy_species_origin_data.csv")

# fungal taxonomy figure for presentations and manuscripts
set.seed(4)
taxplot <- tdat %>%
  ggplot(aes(x = rank, y = abundance)) +
  geom_vline(data = breakdat, aes(xintercept = rank.break), linetype = "dashed", size = 0.2) +
  geom_text_repel(data = fdat, aes(x = rank, label = pathogen), fontface = "italic", size = 3, segment.color = "gray", min.segment.length = 0, hjust = 0, nudge_x = 1, direction = "y") +
  facet_wrap(~year, scales = "free")  +  
  geom_point(size = 0.5, color = "red") +
  xlab("Rank") +
  ylab("Abundance") +
  theme_bw() +
  theme(axis.text = element_text(size = 10, color="black"),
        axis.title = element_text(size = 12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 10, color="black"))

pdf("./output/species_origin_taxonomy_rank.pdf", width = 6, height = 3.5)
set.seed(4)
taxplot
dev.off()

pdf("./output/figureS1_fungal_taxonomy_rank.pdf", width = 7, height = 3.5)
set.seed(4)
taxplot
dev.off()


#### PERMANOVA ####
pmod1 <- adonis(cdat1 ~ grass_group + host_sp + year.f, data = edat1, method="chao")
pmod1$aov.tab
# all three are sig, status is the least important

tab_df(round(pmod1$aov.tab, 3))


#### NMDS ####

# models
nmds1 <- metaMDS(cdat1, distance = "chao", halfchange = F, expand = F, trymax = 100) # error: results may be meaningless with non-integer data in method “chao” (every value is an integer)
nmds1b <- metaMDS(cdat1, distance = "bray", halfchange = F, expand = F, trymax = 100)

# goodness of fit
gof1 <- goodness(nmds1)
gof1b <- goodness(nmds1b)

# extract data scores from nmds
ndat1 <- as.data.frame(scores(nmds1)) %>%
  cbind(edat1) %>%
  cbind(gof1)
head(ndat1) 

ndat1b <- as.data.frame(scores(nmds1b)) %>%
  cbind(edat1) %>%
  cbind(gof1b)
head(ndat1b) 

# function to extract ellipses
# edited version of : https://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplot
veganCovEllipse <- function (df){
  scale <- 1
  npoints <- 100
  cov <- with(df, cov.wt(cbind(NMDS1, NMDS2), wt=rep(1/length(NMDS1), length(NMDS1)))$cov)
  center <- with(df, c(mean(NMDS1), mean(NMDS2)))
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

# extract ellipses
# edited version of: https://stackoverflow.com/questions/56646820/functional-programming-use-broom-nest-tidy-unnest-and-map-within-a-function
ell1 <- ndat1 %>%
  group_by(grass_group) %>%
  nest() %>%
  mutate(out = map(data, veganCovEllipse),
         tidied = map(out, as_tibble)) %>%
  unnest(tidied)

# extract means for labels
# lab1 <- ndat1 %>%
#   group_by(grass_group) %>%
#   summarise(NMDS1 = mean(NMDS1),
#             NMDS2 = mean(NMDS2))
lab1 <- ndat1 %>%
  select(grass_group) %>%
  unique() %>%
  mutate(NMDS1 = c(0.37, 0.3),
         NMDS2 = c(-0.17, 0.15))

# abbreviate genus
ndat1 <- ndat1 %>%
  mutate(grass_sp = recode(host_sp, "Avena barbata" = "A. barbata", "Avena fatua" = "A. fatua", "Brachypodium distachyonm" = "B. distachyonm", "Bromus diandrus" = "B. diandrus", "Bromus hordeaceus" = "B. hordeaceus", "Elymus glaucus" = "E. glaucus", "Stipa pulchra" = "S. pulchra"))


#### visualize ####

# Chao nmds
# for presentation
pdf("./output/species_origin_community_differences.pdf", width = 5.5, height = 4)
ndat1 %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(fill = grass_sp, shape = grass_group, group = year.f), size = 3) +
  geom_path(data = ell1, aes(linetype = grass_group), size=1) +  
  geom_text(data = lab1, aes(label = grass_group)) +
  scale_fill_manual(values = pal, name = "Grass species") +
  scale_shape_manual(values = c(22, 25), name = "Grass group") +
  guides(fill=guide_legend(override.aes=list(shape=c(rep(25, 5), 22, 22)), label.theme = element_text(angle = 0, face = "italic")), linetype = "none") +
  theme_bw() +
  theme(axis.text = element_text(size = axisText, color="black"),
        axis.title = element_text(size = axisTitle),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size=legendText),
        legend.title = element_text(size=legendTitle))
dev.off()

# for publication
nmdsplot <- ndat1 %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(fill = grass_group, shape = grass_sp, group = year.f), size = 3) +
  geom_path(data = ell1, aes(linetype = grass_group), size=1) +  
  geom_text(data = lab1, aes(label = grass_group)) +
  scale_fill_manual(values = c("black", "white"), name = "Grass group") +
  scale_shape_manual(values = shapes, name = "Grass species") +
  guides(fill = guide_legend(override.aes = list(shape = 21)), linetype = "none", shape = guide_legend(label.theme = element_text(size = 10, angle = 0, face = "italic"))) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, color="black"),
        axis.title = element_text(size = 12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        legend.spacing.y = unit(-0.03, "cm"))

# bray-curtis nmds
ndat1b %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = host, shape = grass.status, size = year.f))


#### combine figures for manuscript
pdf("./output/figure1_pathogen_communities_host_status.pdf", width = 7.0, height = 3.5)
cowplot::plot_grid(nmdsplot, hostplot, nrow = 1, rel_widths = c(1, 0.4), labels = c("A", "B"), label_size = 12, hjust = c(-0.5, 1))
dev.off()