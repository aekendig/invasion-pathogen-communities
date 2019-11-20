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
shapes <- c(0, 1, 2, 5, 19, 13, 14, 15)
axisText = 12
axisTitle = 14
legendText = 12
legendTitle = 12

# import data
dat <- read.csv("./data/fungal_pathogens_2015_2017.csv")

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

# function to create wide data
datw_fun <- function(df, min.iso) {
  dout <- df %>%
    group_by(year, year.f, host, grass.status, otu.id) %>%
    summarise(abundance = length(isolate.id)) %>%
    group_by(year, year.f, host, grass.status) %>%
    mutate(isolates = sum(abundance)) %>%
    filter(isolates > min.iso) %>%
    ungroup() %>%
    spread(otu.id, abundance)
  
  return(dout)
}

# function to create community matrix
cdat_fun <- function(datw){
  cdat <- datw %>%
    select(-c(year:isolates)) %>%
    data.frame()
  
  cdat[is.na(cdat)] = 0
  cdat <- as.matrix(cdat)
  
  return(cdat)
}

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


#### edit data ####

# remove JEF to make numbers between native and non-native more comparable
dat1 <- dat %>%
  filter(experiment != "JEF transect") %>%
  mutate(grass.status = recode(grass.status, "non-native" = "exotic"))

# sample sizes
nrow(dat1)
unique(dat1$experiment)
unique(dat1$host)

# wide dataset
(datw <- datw_fun(dat1, 0))
# create one that excludes BRAD (much lower sampling than others)
(datw1 <- datw_fun(dat1, 4))

# community matrix
cdat1 <- cdat_fun(datw1)

# environmental matrix
edat1 <- datw1 %>% select(c(year:grass.status))


#### host specificity analysis ####

# pathogen data
pdat <- dat1 %>%
  group_by(grass.status, otu.id, taxonomy) %>%
  summarise(abundance = length(isolate.id)) %>%
  mutate(pathogen = mgsub(as.character(taxonomy), c(" cf.", " sp.", "unknown "), c("", "", "")),
         species = case_when(sapply(strsplit(pathogen, " "), length) == 2 ~ TRUE,
                             TRUE ~ FALSE)) %>%
  ungroup()

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
# write_csv(sdat, "./output/pathogen_host_range_usda.csv")
sdat <- read_csv("./output/pathogen_host_range_usda.csv")

# check species with one host in database
filter(sdat, host.num == 1)
# associations("Ramularia proteae", database = "FH", spec_type = "fungus") # not a host we have
# associations("Anthostomella proteae", database = "FH", spec_type = "fungus") # not a host we have
# associations("Cladosporium antarcticum", database = "FH", spec_type = "fungus") # not a host we have
# associations("Pyrenophora dactylidis", database = "FH", spec_type = "fungus") # not a host we have
# associations("Chaetomium pseudocochliodes", database = "FH", spec_type = "fungus") # not a host we have

# total abundance by group
gdat <- dat1 %>%
  group_by(grass.status) %>%
  summarise(n = n())

# merge specificity data back with pathogen data
pdat2 <- full_join(pdat, sdat) %>%
  mutate(prop.abundance = case_when(grass.status == "native" ~ abundance/gdat$n[1],
                                    grass.status == "exotic"  ~ abundance/gdat$n[2]),
         host.num.scaled = host.num * prop.abundance)
  
# figure
pdat2 %>%
  ggplot(aes(x = grass.status, y = host.num.scaled)) +
  stat_summary(geom = "point", fun.y = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se")

# merge specificity data with full dataset
sdat2 <- pdat2 %>%
  select(taxonomy, pathogen, host.num, species) %>%
  unique() %>%
  right_join(dat1) %>%
  mutate(species.known = as.numeric(species))

# summarise by using whole dataset
hostplot <- sdat2 %>%
  ggplot(aes(x = grass.status, y = host.num, fill = grass.status)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.1) +
  stat_summary(geom = "point", fun.y = "mean", size = 5, shape = 21) +
  scale_fill_manual(values = c("black", "white"), guide = F) +
  ylab("Host range of associated pathogens") +
  xlab("Host status") +
  theme_bw() +
  theme(axis.text = element_text(size = 10, color="black"),
        axis.title = element_text(size = 12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
# only takes USDA pathogens into account

# proportion with species known
sdat2 %>%
  ggplot(aes(x = grass.status, y = species.known)) +
  stat_summary(geom = "point", fun.y = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot")

# total samples with species known
sdat2 %>%
  group_by(grass.status) %>%
  summarise(n = n(),
            database = sum(!is.na(host.num))) %>%
  mutate(prop = database / n)

# t-test
smod <- t.test(host.num ~ grass.status, data = filter(sdat2, !is.na(host.num)))
smod


#### output taxonomy table and figure ####

# abundance by grass status
adat <- sdat2 %>%
  group_by(year, otu.id, taxonomy, pathogen, grass.status) %>%
  summarise(abundance = length(isolate.id)) %>%
  spread(key = grass.status, value = abundance) %>%
  ungroup() %>%
  rename(native.abundance = native, exotic.abundance = exotic)

# taxonomy data
tdat <- sdat2 %>%
  group_by(year, year.f, otu.id, taxonomy, pathogen, host.num) %>%
  summarise(abundance = length(isolate.id),
            host.species = paste(unique(host), collapse = ", ")) %>%
  ungroup() %>%
  group_by(year, year.f) %>%
  mutate(rank = rank(-abundance, ties.method = "random")) %>%
  ungroup() %>%
  left_join(adat) %>%
  mutate(native.abundance = replace_na(native.abundance, 0),
         exotic.abundance = replace_na(exotic.abundance, 0))

# look at all fungi
tdat %>%
  ggplot(aes(x = rank, y = abundance)) +
  geom_point(size = 2) +
  facet_wrap(~year)

# calculate breaks (mutate didn't work)
tdat2015 <- filter(tdat, year == 2015)
tdat2015 <- break_fun(tdat2015) %>% full_join(tdat2015)
filter(tdat2015, breaks == max(breaks)) # 5 - 6

tdat2016 <- filter(tdat, year == 2016)
tdat2016 <- break_fun(tdat2016) %>% full_join(tdat2016)
sort(tdat2016$breaks)
filter(tdat2016, breaks == 15) # 2 - 3
filter(tdat2016, breaks == 12) # 5 - 6

tdat2017 <- filter(tdat, year == 2017)
tdat2017 <- break_fun(tdat2017) %>% full_join(tdat2017)
sort(tdat2017$breaks)
filter(tdat2017, breaks == 10) # 4 - 5

# indicate breaks
breakdat <- tibble(year = c(2015, 2016, 2017), rank.break = c(5.5, 5.5, 4.5))

# extract top fungal species
fdat <- tdat %>%
  full_join(breakdat) %>%
  filter(rank < rank.break) %>%
  select(pathogen, year, host.species, abundance, native.abundance, exotic.abundance, rank, host.num, otu.id) %>%
  mutate(pathogen = case_when(pathogen == "Drechslera" ~ "Drechslera sp.",
                              TRUE ~ pathogen),
         host.species = mgsub(host.species, c("AB", "BH", "SP", "BD", "PA", "EG", "FP", "AF"), c("Ab", "Bh", "Sp", "Bd", "Pa", "Eg", "Fp", "Af"))) %>%
  arrange(pathogen, year, rank)

# number of fungal species
length(unique(fdat$otu.id)) # 7 (one only occurs on native species)

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
pmod1 <- adonis(cdat1 ~ grass.status + host + year.f, data = edat1, method="chao")
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

# extract ellipses
# edited version of: https://stackoverflow.com/questions/56646820/functional-programming-use-broom-nest-tidy-unnest-and-map-within-a-function
ell1 <- ndat1 %>%
  group_by(grass.status) %>%
  nest() %>%
  mutate(out = map(data, veganCovEllipse),
         tidied = map(out, as_tibble)) %>%
  unnest(tidied)

# extract means for labels
lab1 <- ndat1 %>%
  group_by(grass.status) %>%
  summarise(NMDS1 = mean(NMDS1),
            NMDS2 = mean(NMDS2))


#### visualize ####

# Chao nmds
# for presentation
pdf("./output/species_origin_community_differences.pdf", width = 5.5, height = 4)
ndat1 %>%
  mutate(host = recode(host, AB = "A. barbata", AF = "A. fatua", BD = "B. diandrus", BH = "B. hordeaceus", EG = "E. glaucus", FP = "F. perennis", PA = "P. aquatica", SP = "S. pulchra")) %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(fill = host, shape = grass.status, group = year.f), size = 3) +
  geom_path(data = ell1, aes(linetype = grass.status), size=1) +  
  geom_text(data = lab1, aes(label = grass.status)) +
  scale_fill_manual(values = pal, name = "Host species") +
  scale_shape_manual(values = c(22, 25), name = "Host status") +
  guides(fill=guide_legend(override.aes=list(shape=c(rep(25, 4), 22, 25, 25, 22)), label.theme = element_text(angle = 0, face = "italic")), linetype = "none") +
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
  mutate(host = recode(host, AB = "A. barbata", AF = "A. fatua", BD = "B. diandrus", BH = "B. hordeaceus", EG = "E. glaucus", FP = "F. perennis", PA = "P. aquatica", SP = "S. pulchra"),
         grass.status = recode(grass.status, "non-native" = "exotic")) %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(fill = grass.status, shape = host, group = year.f), size = 3) +
  geom_path(data = ell1, aes(linetype = grass.status), size=1) +  
  geom_text(data = lab1, aes(label = grass.status)) +
  scale_fill_manual(values = c("black", "white"), name = "Host status") +
  scale_shape_manual(values = shapes, name = "Host species") +
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