## Goal: How does foliar fungal pathogen community composition differ between native perennial and non-native annual grasses?


#### set up ####

# clear everything except final data
rm(list = ls())

# load packages
library(tidyverse)
library(vegan)
library(sjPlot)
library(rusda)
library(mgsub)
library(cowplot)
library(ggrepel)

# import data
dat <- read.csv("./data/fungal_pathogens_2015_2017.csv")

# figure settings
shapes <- c(0, 1, 2, 5, 19, 17)
axisText = 10
axisTitle = 12
legendText = 10
legendTitle = 10


#### edit data ####

# remove JEF to make numbers between native and non-native more comparable
dat1 <- dat %>%
  filter(experiment != "JEF transect" & 
           host %in% c("AB", "AF", "BD", "BH", "EG", "SP")) %>%
  mutate(host.sp = recode(host, AB = "Avena barbata", AF = "Avena fatua", BD = "Bromus diandrus", BH = "Bromus hordeaceus", EG = "Elymus glaucus", SP = "Stipa pulchra"),
         host.sp = fct_relevel(host.sp, "Avena barbata", "Avena fatua", "Bromus diandrus"))

# sample sizes
nrow(dat1)
unique(dat1$experiment)
unique(dat1$host)

# subplots for labelling map
unique(select(dat1, experiment, subplot)) %>% 
  arrange(experiment, subplot) %>%
  data.frame()

# set minimum isolate number
min.iso <- 4

# create wide data
datw1 <- dat1 %>%
  group_by(year, year.f, host.sp, grass.group, otu.id) %>%
  summarise(abundance = length(isolate.id)) %>%
  group_by(year, year.f, host.sp, grass.group) %>%
  mutate(isolates = sum(abundance)) %>%
  filter(isolates >= min.iso) %>%
  ungroup() %>%
  spread(otu.id, abundance)

# community matrix
cdat1 <- datw1 %>%
  select(-c(year:isolates)) %>%
  data.frame()
cdat1[is.na(cdat1)] = 0
cdat1 <- as.matrix(cdat1)

# environmental matrix
edat1 <- datw1 %>% select(c(year:grass.group))


#### host range ####

# pathogen data
pdat <- dat1 %>%
  group_by(grass.group, otu.id, taxonomy) %>%
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

# unique pathogens
paths <- pdat %>%
  group_by(otu.id, taxonomy, pathogen, species) %>%
  summarise(abundance = sum(abundance)) %>%
  ungroup() %>%
  unique()

# species
sum(paths$species)
# 36

# pathogen names
filter(paths, species == 1) %>%
  select(pathogen) %>%
  unique()
# 29

# see where the collapse comes from
filter(paths, species == 1) %>%
  group_by(pathogen) %>%
  summarise(taxa = paste(taxonomy, collapse  = ", "),
            otus = paste(otu.id, collapse = ", ")) %>%
  data.frame()
# combining different OTUs
# taxonomy are the same

# pathogens with host range data - couldn't get it to work within the case_when statement
# had to remove pathogens without records
sdat <- pdat %>%
  filter(species == T & !(pathogen %in% c("Parastagonospora cumpignensis", "Preussia lignicola", "Gilmaniella subornata", "Cystostereum heteromorphum", "Preussia minipascua", "Ochrocladosporium frigidarii"))) %>%
  select(pathogen) %>%
  unique() %>%
  mutate(host.num = NA)

# add host range data to pathogens with loop (errors with mutate)
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
gdat <- dat1 %>%
  group_by(grass.group) %>%
  summarise(n = n())

# merge host range data back with pathogen data
pdat2 <- full_join(pdat, sdat2) %>%
  mutate(prop.abundance = case_when(grass.group == "native\nperennial" ~ abundance/gdat$n[1],
                                    grass.group == "non-native\nannual"  ~ abundance/gdat$n[2]),
         host.num.scaled = host.num * prop.abundance)
  
# visualize
pdat2 %>%
  ggplot(aes(x = grass.group, y = host.num.scaled)) +
  stat_summary(geom = "point", fun.y = "mean") +
  stat_summary(geom = "errorbar", fun.data = "mean_se")

# merge host range data with full dataset
sdat3 <- pdat2 %>%
  select(taxonomy, pathogen, host.num, species) %>%
  unique() %>%
  right_join(dat1) %>%
  mutate(species.known = as.numeric(species))

# summarise by using whole dataset
hostplot <- sdat3 %>%
  ggplot(aes(x = grass.group, y = host.num, fill = grass.group)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.1) +
  stat_summary(geom = "point", fun.y = "mean", size = 5, shape = 21) +
  scale_fill_manual(values = c("black", "white"), guide = F) +
  ylab("Host range of associated pathogens") +
  xlab("Host group") +
  theme_bw() +
  theme(axis.text = element_text(size = axisText, color="black"),
        axis.title = element_text(size = axisTitle),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
hostplot

# t-test
smod <- t.test(host.num ~ grass.group, data = filter(sdat3, !is.na(host.num)))
smod


#### select focal pathogens ####

# abundance by grass status
adat <- sdat3 %>%
  group_by(year, otu.id, taxonomy, pathogen, grass.group) %>%
  summarise(abundance = length(isolate.id)) %>%
  spread(key = grass.group, value = abundance) %>%
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
  mutate(pathogen = case_when(pathogen == "Drechslera" ~ "Pyrenophora sp.",
                              TRUE ~ pathogen),
         host.species = mgsub(host.species, c("AB", "BH", "SP", "BD", "EG", "AF"), c("Ab", "Bh", "Sp", "Bd", "Eg", "Af"))) %>%
  arrange(pathogen, year, rank)

# number of fungal species
length(unique(fdat$otu.id)) # 7

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


#### PERMANOVA ####
pmod1 <- adonis(cdat1 ~ grass.group + host.sp + year.f, data = edat1, method="chao")
pmod1$aov.tab
# all three are sig, status is the least important




#### NMDS ####

# models
nmds1 <- metaMDS(cdat1, distance = "chao", halfchange = F, expand = F, trymax = 100) # error: results may be meaningless with non-integer data in method “chao” (every value is an integer)
nmds1b <- metaMDS(cdat1, distance = "bray", halfchange = F, expand = F, trymax = 100)

# goodness of fit
gof1 <- goodness(nmds1)
gof1b <- goodness(nmds1b)

# stress plot
stressplot(nmds1)
stressplot(nmds1b)

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
  group_by(grass.group) %>%
  nest() %>%
  mutate(out = map(data, veganCovEllipse),
         tidied = map(out, as_tibble)) %>%
  unnest(tidied)

# extract means for labels
# lab1 <- ndat1 %>%
#   group_by(grass.group) %>%
#   summarise(NMDS1 = mean(NMDS1),
#             NMDS2 = mean(NMDS2))
lab1 <- ndat1 %>%
  select(grass.group) %>%
  unique() %>%
  mutate(NMDS1 = c(0, -0.18),
         NMDS2 = c(-0.2, 0.25))

# abbreviate genus
ndat1 <- ndat1 %>%
  mutate(grass.sp = recode(host.sp, "Avena barbata" = "A. barbata", "Avena fatua" = "A. fatua", "Bromus diandrus" = "B. diandrus", "Bromus hordeaceus" = "B. hordeaceus", "Elymus glaucus" = "E. glaucus", "Stipa pulchra" = "S. pulchra"))

# visualize
nmdsplot <- ndat1 %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(fill = grass.group, shape = grass.sp, group = year.f), size = 3) +
  geom_path(data = ell1, aes(linetype = grass.group), size=1) +  
  geom_text(data = lab1, aes(label = grass.group)) +
  scale_fill_manual(values = c("black", "white"), name = "Host group") +
  scale_shape_manual(values = shapes, name = "Grass species") +
  guides(fill = guide_legend(override.aes = list(shape = 21)), linetype = "none", shape = guide_legend(label.theme = element_text(size = axisText, angle = 0, face = "italic"))) +
  theme_bw() +
  theme(axis.text = element_text(size = axisText, color="black"),
        axis.title = element_text(size = axisTitle),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = legendText),
        legend.title = element_text(size = legendTitle),
        legend.spacing.y = unit(0.05, "cm"))

# bray-curtis nmds
ndat1b %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = host.sp, shape = grass.group, size = year.f))


#### bar plot ####

# focal pathogens
fabb <- tibble(otu.id = c(1, 7, 4, 5, 8, 2, 3),
               path.abb = c("A. inf.", "P. ave.", "P. cha.", "P. lol.", "P. tri.", "Pyr. sp.", "R. pro."),
               color.pal = c("#000000", "#56B4E9", "#009E73", "#F0E442", "gray50", "#81ffff", "#00edad"))

# summarize by proportion
# add abbreviated genus
# add pathogen name
dat1sum <- dat1 %>%
  mutate(grass.sp = recode(host.sp, "Avena barbata" = "A. barbata", "Avena fatua" = "A. fatua", "Bromus diandrus" = "B. diandrus", "Bromus hordeaceus" = "B. hordeaceus", "Elymus glaucus" = "E. glaucus", "Stipa pulchra" = "S. pulchra")) %>%
  group_by(grass.status, grass.group, grass.sp, otu.id, taxonomy) %>%
  summarise(isolates = n()) %>%
  ungroup() %>%
  group_by(grass.status, grass.group, grass.sp) %>%
  mutate(tot = sum(isolates),
         prop = isolates / tot) %>%
  ungroup() %>%
  mutate(grass.sp = fct_rev(grass.sp)) %>%
  left_join(fabb) %>%
  mutate(path.abb = fct_rev(ifelse(is.na(path.abb), paste("z", otu.id, sep = ""), path.abb))) %>%
  arrange(path.abb)

# remaining colors
n = length(unique(dat1sum$otu.id)) - nrow(fabb)
reds <- colorRampPalette(c("#E69F00", "#6b2f00"))
color.red <- reds(n)

# bar plot
barplot <- dat1sum %>%
  ggplot(aes(x = grass.sp, y = prop, fill = path.abb)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = tot), check_overlap = T, y = 1.025, size = 3) +
  facet_grid(grass.group ~ ., scales = "free", space = "free") +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, color="black"),
        axis.text.y = element_text(size = 10, color="black", face = "italic"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 12),
        strip.background = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.text = element_text(size = 9, face = "italic"),
        legend.title = element_text(size = 12),
        legend.margin = margin(0, 0, 0, 0, unit="cm")) +
  ylim(0, 1.025) +
  ylab("Proportion of isolates") +
  scale_fill_manual(breaks = fabb$path.abb,
                    values = c(color.red, rev(fabb$color.pal)),
                    name = "Pathogen") +
  guides(fill = guide_legend(nrow = 1))
barplot

#### outputs ####

### sample sizes ###

# sample sizes by group
dat1 %>%
  group_by(year.f, grass.group) %>%
  summarise(n = n())

# sample sizes by species 
(sptable <- dat1 %>%
    group_by(year, grass.group, host.sp) %>%
    summarise(n = n()) %>%
    arrange(year, grass.group, host.sp))

write_csv(sptable, "./output/host_species_sample_sizes.csv")

### overlapping fungal species ###

# shared isolates
overlap <- dat1 %>%
  mutate(native = ifelse(grass.status == "native", 1, 0),
         non = ifelse(grass.status == "non-native", 1, 0)) %>%
  group_by(otu.id, taxonomy) %>%
  summarise(native = sum(native) > 0,
            non = sum(non) > 0) %>%
  mutate(native_only = ifelse(native == T & non == F, 1, 0),
         non_only = ifelse(native == F & non == T, 1, 0),
         both = ifelse(native == T & non == T, 1, 0))

length(unique(dat1$otu.id)) # 83 OTU's
sum(overlap$native_only) # 41
sum(overlap$non_only) # 18
sum(overlap$both) # 24

# merge overlap data back with isolate data
dat1 %>%
  left_join(select(overlap, otu.id, taxonomy, native_only, non_only, both)) %>%
  group_by(grass.group) %>%
  summarise(tot = n(),
            native_only = sum(native_only),
            non_only = sum(non_only),
            both = sum(both)) %>%
  mutate(nat_prop = native_only/tot,
         non_prop = non_only/tot,
         both_prop = both/tot)

# total isolates
534+427

### isolates identified to species ###

# samples with species known
sdat3 %>%
  group_by(grass.status) %>%
  summarise(n = n(),
            database = sum(!is.na(host.num))) %>%
  mutate(prop = database / n)

### fungal taxonomy rank ###

write_csv(fdat, "./output/fungal_taxonomy_rank.csv")

pdf("./output/figureS2_fungal_taxonomy_rank.pdf", width = 7, height = 3.5)
set.seed(4)
taxplot
dev.off()

### PERMANOVA table ###
tab_df(round(pmod1$aov.tab, 3))

### combined figures for manuscript ###

# lower plots
lowplots <- cowplot::plot_grid(nmdsplot, hostplot, nrow = 1, rel_widths = c(1, 0.4), labels = c("B", "C"), label_size = 12, hjust = c(-0.5, 1))

pdf("./output/figure1_pathogen_communities_host_status.pdf", width = 7.0, height = 7)
cowplot::plot_grid(barplot, lowplots, nrow = 2, labels = c("A", "", ""), label_size = 12, hjust = c(-0.5, 1))
dev.off()