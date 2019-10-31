### goal: create figures

# variables: density and relative abundance, frequency- and density-dependent transmission, intras- vs. interspecific transmission

# model: 
# Population dynamics of pathogens with multiple host species, Dobson 2004, The American Naturalist
# Edited model to make the populations at carrying capacity (deaths become births)
# Form of frequency-dependence (divide by transmitting population) from: Disease emergence in multi-host epidemic models, McCormack and Allen 2007, Mathematical Medicine and Biology
# Form of frequency-dependence (divide by receiving population) from: Mathematical Models for Microparasites of Wildlife, Heesterbeek and Roberts 1995, Ecology of Infectious Disease in Natural Populations
# closely related model: An immunization model for a heterogeneous population, Hethcote 1978, Theoretical Population Biology
# closely related model: Infectious Disease and Species Coexistence: A Model of Lotka-Volterra, Holt and Pickering 1985, The American Naturalist
# closely related model: Disease and community structure: The importance of host self-regulation in a host-host-pathogen model, Begon et al. 1992, The American Naturalist

# density/relative abundance scenarios: simulation_scenarios_102619.xlsx


#### set-up ####

# clear all existing data
rm(list=ls())

# load libraries
library(deSolve)
library(tidyverse)
library(nationalparkcolors)
library(cowplot)

# figure settings
pal <- park_palette("Everglades")[1:5]
axisText = 10
axisTitle = 12
legendText = 10
legendTitle = 12

gg_theme <- theme_bw() +
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


#### default parameters ####

# death rates
d1 <- 0.1
d2 <- 0.1

# infection-induced mortality  
alpha1 <- 0.01
alpha2 <- 0.01

# recover rate
theta1 <- 0.2
theta2 <- 0.2

# initial values for I and R
I10 <- 1
R10 <- 0

I20 <- 0
R20 <- 0

# simulation time
simtime <- 100

# simulation resolution
npoints <- 50


#### scenario descriptions ####
scen <- tibble(scenario = 1:8,
               description = c("1: no change", "2: sp. 1 density", "3: sp. 2 density", "4: sp. 1 rel. abundance", "5: sp. 2 rel. abundance", "6: all", "7: total density 1", "8: total density 2"))


#### density-dependent model ####
dd_mod <- function(d1, B11, B12, alpha1, theta1, d2, B22, B21, alpha2, theta2, S10, I10, R10, S20, I20, R20, simtime){
  
  # parameters
  parms<-list(d1 = d1, B11 = B11, B12 = B12, alpha1 = alpha1, theta1 = theta1, d2 = d2, B22 = B22, B21 = B21, alpha2 = alpha2, theta2 = theta2, S10 = S10, I10 = I10, R10 = R10, S20 = S20, I20 = I20, R20 = R20, simtime = simtime)
  
  # define model
  mymodel <- with(as.list(parms), function(t, x, parms){
    
    # State variables
    S1 <- x["S1"]
    I1 <- x["I1"]
    R1 <- x["R1"]
    
    S2 <- x["S2"]
    I2 <- x["I2"]
    R2 <- x["R2"]
    
    # species 1
    S1dot <- d1 * (S1 + R1) + (d1 + alpha1) * I1 - S1 * (B11 * I1 + B12 * I2) - d1 * S1
    I1dot <- S1 * (B11 * I1 + B12 * I2) - (d1 + alpha1 + theta1) * I1
    R1dot <- theta1 * I1 - d1 * R1
    
    # species 2
    S2dot <- d2 * (S2 + R2) + (d2 + alpha2) * I2 - S2 * (B22 * I2 + B21 * I1) - d2 * S2
    I2dot <- S2 * (B22 * I2 + B21 * I1) - (d2 + alpha2 + theta2) * I2
    R2dot <- theta2 * I2 - d2 * R2
    
    list(c(S1dot, I1dot, R1dot, S2dot, I2dot, R2dot))	
    
  })
  
  # start values
  xstart <- c(S1 = S10, I1 = I10, R1 = R10, S2 = S20, I2 = I20, R2 = R20)
  
  # vector of timesteps
  times  <- seq(0, simtime, length = simtime)
  
  # LSODA (default step size)
  out <- as.data.frame(lsoda(xstart, times, mymodel, parms, hmax=20))
  
  # output times series
  return(out)  
}


#### frequency-dependent model ####
fd_mod <- function(d1, B11, B12, alpha1, theta1, d2, B22, B21, alpha2, theta2, S10, I10, R10, S20, I20, R20, simtime){
  
  # parameters
  parms<-list(d1 = d1, B11 = B11, B12 = B12, alpha1 = alpha1, theta1 = theta1, d2 = d2, B22 = B22, B21 = B21, alpha2 = alpha2, theta2 = theta2, S10 = S10, I10 = I10, R10 = R10, S20 = S20, I20 = I20, R20 = R20, simtime = simtime)
  
  # define model
  mymodel <- with(as.list(parms), function(t, x, parms){
    
    # State variables
    S1 <- x["S1"]
    I1 <- x["I1"]
    R1 <- x["R1"]
    
    S2 <- x["S2"]
    I2 <- x["I2"]
    R2 <- x["R2"]
    
    # total population size
    N1 <- S1 + I1 + R1
    N2 <- S2 + I2 + R2
    
    # species 1
    S1dot <- d1 * (S1 + R1) + (d1 + alpha1) * I1 - S1 * (B11 * I1/N1 + B12 * I2/N2) - d1 * S1
    I1dot <- S1 * (B11 * I1/N1 + B12 * I2/N2) - (d1 + alpha1 + theta1) * I1
    R1dot <- theta1 * I1 - d1 * R1
    
    # species 2
    S2dot <- d2 * (S2 + R2) + (d2 + alpha2) * I2 - S2 * (B22 * I2/N2 + B21 * I1/N1) - d2 * S2
    I2dot <- S2 * (B22 * I2/N2 + B21 * I1/N1) - (d2 + alpha2 + theta2) * I2
    R2dot <- theta2 * I2 - d2 * R2
    
    list(c(S1dot, I1dot, R1dot, S2dot, I2dot, R2dot))	
    
  })
  
  # start values
  xstart <- c(S1 = S10, I1 = I10, R1 = R10, S2 = S20, I2 = I20, R2 = R20)
  
  # vector of timesteps
  times  <- seq(0, simtime, length = simtime)
  
  # LSODA (default step size)
  out <- as.data.frame(lsoda(xstart, times, mymodel, parms, hmax=20))
  
  # output times series
  return(out)  
}


#### simulation function ####

sim_fun <- function(S10_in, S20_in, points, filename, scenario, transmission, spillover){
  
  # empty dataframe
  out_df <- data.frame(S1_dens = S10_in, S2_dens = S20_in, time = NA, S1 = NA, I1 = NA, R1 = NA, S2 = NA, I2 = NA, R2 = NA)
  
  # start pdf
  pdf(file = filename)
  
  # cycle through densities
  for(i in 1:points){
    
    # set initial density
    S10 <- S10_in[i]
    S20 <- S20_in[i]
   
    # save model output
    ifelse(transmission == "dd", 
           temp_out <- dd_mod(d1, B11, B12, alpha1, theta1, d2, B22, B21, alpha2, theta2, S10, I10, R10, S20, I20, R20, simtime), 
           temp_out <- fd_mod(d1, B11, B12, alpha1, theta1, d2, B22, B21, alpha2, theta2, S10, I10, R10, S20, I20, R20, simtime))
    
    # make output long
    temp_out_l <- temp_out %>% 
      gather(key = "pop", value = "abundance", -time)
    
    # print simulations to pdf
    print(ggplot(temp_out_l, aes(time, abundance, color = pop)) + 
            geom_line())
    
    # add final density to output dataframe
    out_df[i, 3:9] <- filter(temp_out, time == simtime)
    
  }
  
  # close pdf
  dev.off()
  
  # calculate prevalence
  out_df_prev <- out_df %>%
    unique() %>%
    mutate(prev1 = I1/(S1 + I1 + R1),
           prev2 = I2/(S2 + I2 + R2),
           prev = (I1 + I2)/(S1 + I1 + R1 + S2 + I2 + R2)) %>%
    select(S1_dens, S2_dens, prev1, prev2, prev) %>%
    gather(key = "group", value = "prevalence", -c(S1_dens, S2_dens)) %>%
    mutate(scenario = scenario,
           transmission = transmission,
           spillover = spillover)
  
  # return dataframe
  return(out_df_prev)

}



#### dd spillover transmission rates ####

# species 1 is  the reservoir host
# species 2 is the nonreservoir host

# intraspecific transmission rates  
B11 <- 0.2
B22 <- 0   

# interspecific transmission rates  
B12 <- 0
B21 <- 0.1

# default simulation (change #)
# s#_dd_spill <- sim_fun(S10_s#, S20_s#, npoints, "./output/two_host_sim_scen#_dd_spill_102919.pdf", scenario = #, transmission = "dd", spillover = 1)


#### scenario 1 (no change),  dd, spillover ####

# initial values for S
S10_s1 <- rep(300, npoints)
S20_s1 <- rep(300, npoints)

# run simulation
s1_dd_spill <- sim_fun(S10_s1, S20_s1, npoints, "./output/two_host_sim_scen1_dd_spill_102919.pdf", scenario = 1, transmission = "dd", spillover = 1)


#### scenario 2 (spp 1 density),  dd, spillover ####

# initial values for S
S10_s2 <- seq(0, 300, length.out = npoints)
S20_s2 <- rep(10, npoints)

# run simulation
s2_dd_spill <- sim_fun(S10_s2, S20_s2, npoints, "./output/two_host_sim_scen2_dd_spill_102919.pdf", scenario = 2, transmission = "dd", spillover = 1)


#### scenario 3 (spp 2 density),  dd, spillover ####

# initial values for S
S10_s3 <- rep(10, npoints)
S20_s3 <- seq(0, 300, length.out = npoints)

# run simulation
s3_dd_spill <- sim_fun(S10_s3, S20_s3, npoints, "./output/two_host_sim_scen3_dd_spill_102919.pdf", scenario = 3, transmission = "dd", spillover = 1)


#### scenario 4 (spp 1 rel. abundance),  dd, spillover ####

# initial values for S
S10_s4 <- seq(0, 300, length.out = npoints)
S20_s4 <- rep(150, npoints)

# run simulation
s4_dd_spill <- sim_fun(S10_s4, S20_s4, npoints, "./output/two_host_sim_scen4_dd_spill_102919.pdf", scenario = 4, transmission = "dd", spillover = 1)


#### scenario 5 (spp 2 rel. abundance),  dd, spillover ####

# initial values for S
S10_s5 <- rep(150, npoints)
S20_s5 <- seq(0, 300, length.out = npoints)

# run simulation
s5_dd_spill <- sim_fun(S10_s5, S20_s5, npoints, "./output/two_host_sim_scen5_dd_spill_102919.pdf", scenario = 5, transmission = "dd", spillover = 1)


#### scenario 6 (all),  dd, spillover ####

# initial values for S
S10_s6 <- seq(0, 300, length.out = npoints)
S20_s6 <- seq(0, 300, length.out = npoints)

# run simulation
s6_dd_spill <- sim_fun(S10_s6, S20_s6, npoints, "./output/two_host_sim_scen6_dd_spill_102919.pdf", scenario = 6, transmission = "dd", spillover = 1)


#### scenario 7 (total density 1),  dd, spillover ####

# initial values for S
S10_s7 <- seq(0, 300, length.out = npoints)
S20_s7 <- seq(300, 0, length.out = npoints)

# run simulation
s7_dd_spill <- sim_fun(S10_s7, S20_s7, npoints, "./output/two_host_sim_scen7_dd_spill_102919.pdf", scenario = 7, transmission = "dd", spillover = 1)


#### scenario 8 (total density 2),  dd, spillover ####

# initial values for S
S10_s8 <- seq(300, 0, length.out = npoints)
S20_s8 <- seq(0, 300, length.out = npoints)

# run simulation
s8_dd_spill <- sim_fun(S10_s8, S20_s8, npoints, "./output/two_host_sim_scen8_dd_spill_102919.pdf", scenario = 8, transmission = "dd", spillover = 1)


#### density-dependent spillover figure ####

# combine data, add scenario descriptions
dd_spill <- rbind(s1_dd_spill, s2_dd_spill, s3_dd_spill, s4_dd_spill, s5_dd_spill, s6_dd_spill, s7_dd_spill, s8_dd_spill) %>% 
  as_tibble() %>%
  full_join(scen) %>%
  gather(key = "species", value = "density", -c(group, prevalence, scenario, transmission, spillover, description)) %>%
  mutate(species = recode(species, S1_dens = "species 1", S2_dens = "species 2"),
         group = recode(group, prev = "combined", prev1 = "species 1", prev2 = "species 2") %>% factor(levels = c("species 1", "species 2", "combined")))

# create figure
dd_spill_plot1 <- ggplot(data = filter(dd_spill , scenario %in% 1:4), aes(density, prevalence, color = group, linetype = group)) +
  geom_line(size = 1.5) +
  facet_grid(description ~ species)  +
  xlab("Density") +
  ylab("Prevalence") +
  scale_color_manual(values = pal[c(4, 5, 2)], name = "Group") +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), guide = F) +
  gg_theme

dd_spill_plot2 <- ggplot(data = filter(dd_spill , scenario %in% 5:8), aes(density, prevalence, color = group, linetype = group)) +
  geom_line(size = 1.5) +
  facet_grid(description ~ species)  +
  xlab("Density") +
  ylab("Prevalence") +
  scale_color_manual(values = pal[c(4, 5, 2)], name = "Group") +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), guide = F) +
  gg_theme

leg <- get_legend(dd_spill_plot1)

pdf("./output/two_host_sim_prev_dd_spill_102919.pdf")
plot_grid(dd_spill_plot1 + theme(legend.position = "none"), leg, nrow = 2, rel_heights = c(1, 0.05))

plot_grid(dd_spill_plot2 + theme(legend.position = "none"), leg, nrow = 2, rel_heights = c(1, 0.05))
dev.off()


#### fd spillover transmission rates ####

# max prevalence was 1.5e-6 when I used the same rates as dd
# species 1 is  the reservoir host
# species 2 is the nonreservoir host

# intraspecific transmission rates  
B11 <- 0.2 * 300
B22 <- 0  * 300 

# interspecific transmission rates  
B12 <- 0 * 300
B21 <- 0.1 * 300


#### scenario 1 (no change),  fd, spillover ####

# initial values for S
S10_s1 <- rep(300, npoints)
S20_s1 <- rep(300, npoints)

# run simulation
s1_fd_spill <- sim_fun(S10_s1, S20_s1, npoints, "./output/two_host_sim_scen1_fd_spill_102919.pdf", scenario = 1, transmission = "fd", spillover = 1)


#### scenario 2 (spp 1 density),  fd, spillover ####

# initial values for S
S10_s2 <- seq(0, 300, length.out = npoints)
S20_s2 <- rep(10, npoints)

# run simulation
s2_fd_spill <- sim_fun(S10_s2, S20_s2, npoints, "./output/two_host_sim_scen2_fd_spill_102919.pdf", scenario = 2, transmission = "fd", spillover = 1)


#### scenario 3 (spp 2 density),  fd, spillover ####

# initial values for S
S10_s3 <- rep(10, npoints)
S20_s3 <- seq(0, 300, length.out = npoints)

# run simulation
s3_fd_spill <- sim_fun(S10_s3, S20_s3, npoints, "./output/two_host_sim_scen3_fd_spill_102919.pdf", scenario = 3, transmission = "fd", spillover = 1)


#### scenario 4 (spp 1 rel. abundance),  fd, spillover ####

# initial values for S
S10_s4 <- seq(0, 300, length.out = npoints)
S20_s4 <- rep(150, npoints)

# run simulation
s4_fd_spill <- sim_fun(S10_s4, S20_s4, npoints, "./output/two_host_sim_scen4_fd_spill_102919.pdf", scenario = 4, transmission = "fd", spillover = 1)


#### scenario 5 (spp 2 rel. abundance),  fd, spillover ####

# initial values for S
S10_s5 <- rep(150, npoints)
S20_s5 <- seq(0, 300, length.out = npoints)

# run simulation
s5_fd_spill <- sim_fun(S10_s5, S20_s5, npoints, "./output/two_host_sim_scen5_fd_spill_102919.pdf", scenario = 5, transmission = "fd", spillover = 1)


#### scenario 6 (all),  fd, spillover ####

# initial values for S
S10_s6 <- seq(0, 300, length.out = npoints)
S20_s6 <- seq(0, 300, length.out = npoints)

# run simulation
s6_fd_spill <- sim_fun(S10_s6, S20_s6, npoints, "./output/two_host_sim_scen6_fd_spill_102919.pdf", scenario = 6, transmission = "fd", spillover = 1)


#### scenario 7 (total density 1),  fd, spillover ####

# initial values for S
S10_s7 <- seq(0, 300, length.out = npoints)
S20_s7 <- seq(300, 0, length.out = npoints)

# run simulation
s7_fd_spill <- sim_fun(S10_s7, S20_s7, npoints, "./output/two_host_sim_scen7_fd_spill_102919.pdf", scenario = 7, transmission = "fd", spillover = 1)


#### scenario 8 (total density 2),  fd, spillover ####

# initial values for S
S10_s8 <- seq(300, 0, length.out = npoints)
S20_s8 <- seq(0, 300, length.out = npoints)

# run simulation
s8_fd_spill <- sim_fun(S10_s8, S20_s8, npoints, "./output/two_host_sim_scen8_fd_spill_102919.pdf", scenario = 8, transmission = "fd", spillover = 1)


#### frequency-dependent spillover figure ####

# combine data, afd scenario descriptions
fd_spill <- rbind(s1_fd_spill, s2_fd_spill, s3_fd_spill, s4_fd_spill, s5_fd_spill, s6_fd_spill, s7_fd_spill, s8_fd_spill) %>% 
  as_tibble() %>%
  full_join(scen) %>%
  gather(key = "species", value = "density", -c(group, prevalence, scenario, transmission, spillover, description)) %>%
  mutate(species = recode(species, S1_dens = "species 1", S2_dens = "species 2"),
         group = recode(group, prev = "combined", prev1 = "species 1", prev2 = "species 2") %>% factor(levels = c("species 1", "species 2", "combined")))

# create figure
fd_spill_plot1 <- ggplot(data = filter(fd_spill , scenario %in% 1:4), aes(density, prevalence, color = group, linetype = group)) +
  geom_line(size = 1.5) +
  facet_grid(description ~ species)  +
  xlab("Density") +
  ylab("Prevalence") +
  scale_color_manual(values = pal[c(4, 5, 2)], name = "Group") +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), guide = F) +
  gg_theme

fd_spill_plot2 <- ggplot(data = filter(fd_spill , scenario %in% 5:8), aes(density, prevalence, color = group, linetype = group)) +
  geom_line(size = 1.5) +
  facet_grid(description ~ species)  +
  xlab("Density") +
  ylab("Prevalence") +
  scale_color_manual(values = pal[c(4, 5, 2)], name = "Group") +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), guide = F) +
  gg_theme

leg <- get_legend(fd_spill_plot1)

pdf("./output/two_host_sim_prev_fd_spill_102919.pdf")
plot_grid(fd_spill_plot1 + theme(legend.position = "none"), leg, nrow = 2, rel_heights = c(1, 0.05))

plot_grid(fd_spill_plot2 + theme(legend.position = "none"), leg, nrow = 2, rel_heights = c(1, 0.05))
dev.off()


#### dd specialist transmission rates ####

# species 1 is  the host
# species 2 is the non-host

# intraspecific transmission rates  
B11 <- 0.2
B22 <- 0   

# interspecific transmission rates  
B12 <- 0
B21 <- 0


#### scenario 1 (no change),  dd, specialist ####

# initial values for S
S10_s1 <- rep(300, npoints)
S20_s1 <- rep(300, npoints)

# run simulation
s1_dd_spec <- sim_fun(S10_s1, S20_s1, npoints, "./output/two_host_sim_scen1_dd_spec_102919.pdf", scenario = 1, transmission = "dd", spillover = 0)


#### scenario 2 (spp 1 density),  dd, specialist ####

# initial values for S
S10_s2 <- seq(0, 300, length.out = npoints)
S20_s2 <- rep(10, npoints)

# run simulation
s2_dd_spec <- sim_fun(S10_s2, S20_s2, npoints, "./output/two_host_sim_scen2_dd_spec_102919.pdf", scenario = 2, transmission = "dd", spillover = 0)


#### scenario 3 (spp 2 density),  dd, specialist ####

# initial values for S
S10_s3 <- rep(10, npoints)
S20_s3 <- seq(0, 300, length.out = npoints)

# run simulation
s3_dd_spec <- sim_fun(S10_s3, S20_s3, npoints, "./output/two_host_sim_scen3_dd_spec_102919.pdf", scenario = 3, transmission = "dd", spillover = 0)


#### scenario 4 (spp 1 rel. abundance),  dd, specialist ####

# initial values for S
S10_s4 <- seq(0, 300, length.out = npoints)
S20_s4 <- rep(150, npoints)

# run simulation
s4_dd_spec <- sim_fun(S10_s4, S20_s4, npoints, "./output/two_host_sim_scen4_dd_spec_102919.pdf", scenario = 4, transmission = "dd", spillover = 0)


#### scenario 5 (spp 2 rel. abundance),  dd, specialist ####

# initial values for S
S10_s5 <- rep(150, npoints)
S20_s5 <- seq(0, 300, length.out = npoints)

# run simulation
s5_dd_spec <- sim_fun(S10_s5, S20_s5, npoints, "./output/two_host_sim_scen5_dd_spec_102919.pdf", scenario = 5, transmission = "dd", spillover = 0)


#### scenario 6 (all),  dd, specialist ####

# initial values for S
S10_s6 <- seq(0, 300, length.out = npoints)
S20_s6 <- seq(0, 300, length.out = npoints)

# run simulation
s6_dd_spec <- sim_fun(S10_s6, S20_s6, npoints, "./output/two_host_sim_scen6_dd_spec_102919.pdf", scenario = 6, transmission = "dd", spillover = 0)


#### scenario 7 (total density 1),  dd, specialist ####

# initial values for S
S10_s7 <- seq(0, 300, length.out = npoints)
S20_s7 <- seq(300, 0, length.out = npoints)

# run simulation
s7_dd_spec <- sim_fun(S10_s7, S20_s7, npoints, "./output/two_host_sim_scen7_dd_spec_102919.pdf", scenario = 7, transmission = "dd", spillover = 0)


#### scenario 8 (total density 2),  dd, specialist ####

# initial values for S
S10_s8 <- seq(300, 0, length.out = npoints)
S20_s8 <- seq(0, 300, length.out = npoints)

# run simulation
s8_dd_spec <- sim_fun(S10_s8, S20_s8, npoints, "./output/two_host_sim_scen8_dd_spec_102919.pdf", scenario = 8, transmission = "dd", spillover = 0)


#### density-dependent specialist figure ####

# combine data, add scenario descriptions
dd_spec <- rbind(s1_dd_spec, s2_dd_spec, s3_dd_spec, s4_dd_spec, s5_dd_spec, s6_dd_spec, s7_dd_spec, s8_dd_spec) %>% 
  as_tibble() %>%
  full_join(scen) %>%
  gather(key = "species", value = "density", -c(group, prevalence, scenario, transmission, spillover, description)) %>%
  mutate(species = recode(species, S1_dens = "species 1", S2_dens = "species 2"),
         group = recode(group, prev = "combined", prev1 = "species 1", prev2 = "species 2") %>% factor(levels = c("species 1", "species 2", "combined")))

# create figure
dd_spec_plot1 <- ggplot(data = filter(dd_spec , scenario %in% 1:4), aes(density, prevalence, color = group, linetype = group)) +
  geom_line(size = 1.5) +
  facet_grid(description ~ species)  +
  xlab("Density") +
  ylab("Prevalence") +
  scale_color_manual(values = pal[c(4, 5, 2)], name = "Group") +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), guide = F) +
  gg_theme

dd_spec_plot2 <- ggplot(data = filter(dd_spec , scenario %in% 5:8), aes(density, prevalence, color = group, linetype = group)) +
  geom_line(size = 1.5) +
  facet_grid(description ~ species)  +
  xlab("Density") +
  ylab("Prevalence") +
  scale_color_manual(values = pal[c(4, 5, 2)], name = "Group") +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), guide = F) +
  gg_theme

leg <- get_legend(dd_spec_plot1)

pdf("./output/two_host_sim_prev_dd_spec_102919.pdf")
plot_grid(dd_spec_plot1 + theme(legend.position = "none"), leg, nrow = 2, rel_heights = c(1, 0.05))

plot_grid(dd_spec_plot2 + theme(legend.position = "none"), leg, nrow = 2, rel_heights = c(1, 0.05))
dev.off()


#### fd specialist transmission rates ####

# species 1 is  the host
# species 2 is the non-host

# intraspecific transmission rates  
B11 <- 0.2 * 300
B22 <- 0   

# interspecific transmission rates  
B12 <- 0
B21 <- 0


#### scenario 1 (no change),  fd, specialist ####

# initial values for S
S10_s1 <- rep(300, npoints)
S20_s1 <- rep(300, npoints)

# run simulation
s1_fd_spec <- sim_fun(S10_s1, S20_s1, npoints, "./output/two_host_sim_scen1_fd_spec_102919.pdf", scenario = 1, transmission = "fd", spillover = 0)


#### scenario 2 (spp 1 density),  fd, specialist ####

# initial values for S
S10_s2 <- seq(0, 300, length.out = npoints)
S20_s2 <- rep(10, npoints)

# run simulation
s2_fd_spec <- sim_fun(S10_s2, S20_s2, npoints, "./output/two_host_sim_scen2_fd_spec_102919.pdf", scenario = 2, transmission = "fd", spillover = 0)


#### scenario 3 (spp 2 density),  fd, specialist ####

# initial values for S
S10_s3 <- rep(10, npoints)
S20_s3 <- seq(0, 300, length.out = npoints)

# run simulation
s3_fd_spec <- sim_fun(S10_s3, S20_s3, npoints, "./output/two_host_sim_scen3_fd_spec_102919.pdf", scenario = 3, transmission = "fd", spillover = 0)


#### scenario 4 (spp 1 rel. abundance),  fd, specialist ####

# initial values for S
S10_s4 <- seq(0, 300, length.out = npoints)
S20_s4 <- rep(150, npoints)

# run simulation
s4_fd_spec <- sim_fun(S10_s4, S20_s4, npoints, "./output/two_host_sim_scen4_fd_spec_102919.pdf", scenario = 4, transmission = "fd", spillover = 0)


#### scenario 5 (spp 2 rel. abundance),  fd, specialist ####

# initial values for S
S10_s5 <- rep(150, npoints)
S20_s5 <- seq(0, 300, length.out = npoints)

# run simulation
s5_fd_spec <- sim_fun(S10_s5, S20_s5, npoints, "./output/two_host_sim_scen5_fd_spec_102919.pdf", scenario = 5, transmission = "fd", spillover = 0)


#### scenario 6 (all),  fd, specialist ####

# initial values for S
S10_s6 <- seq(0, 300, length.out = npoints)
S20_s6 <- seq(0, 300, length.out = npoints)

# run simulation
s6_fd_spec <- sim_fun(S10_s6, S20_s6, npoints, "./output/two_host_sim_scen6_fd_spec_102919.pdf", scenario = 6, transmission = "fd", spillover = 0)


#### scenario 7 (total density 1),  fd, specialist ####

# initial values for S
S10_s7 <- seq(0, 300, length.out = npoints)
S20_s7 <- seq(300, 0, length.out = npoints)

# run simulation
s7_fd_spec <- sim_fun(S10_s7, S20_s7, npoints, "./output/two_host_sim_scen7_fd_spec_102919.pdf", scenario = 7, transmission = "fd", spillover = 0)


#### scenario 8 (total density 2),  fd, specialist ####

# initial values for S
S10_s8 <- seq(300, 0, length.out = npoints)
S20_s8 <- seq(0, 300, length.out = npoints)

# run simulation
s8_fd_spec <- sim_fun(S10_s8, S20_s8, npoints, "./output/two_host_sim_scen8_fd_spec_102919.pdf", scenario = 8, transmission = "fd", spillover = 0)


#### frequency-dependent specialist figure ####

# combine data, afd scenario descriptions
fd_spec <- rbind(s1_fd_spec, s2_fd_spec, s3_fd_spec, s4_fd_spec, s5_fd_spec, s6_fd_spec, s7_fd_spec, s8_fd_spec) %>% 
  as_tibble() %>%
  full_join(scen) %>%
  gather(key = "species", value = "density", -c(group, prevalence, scenario, transmission, spillover, description)) %>%
  mutate(species = recode(species, S1_dens = "species 1", S2_dens = "species 2"),
         group = recode(group, prev = "combined", prev1 = "species 1", prev2 = "species 2") %>% factor(levels = c("species 1", "species 2", "combined")))

# create figure
fd_spec_plot1 <- ggplot(data = filter(fd_spec , scenario %in% 1:4), aes(density, prevalence, color = group, linetype = group)) +
  geom_line(size = 1.5) +
  facet_grid(description ~ species)  +
  xlab("Density") +
  ylab("Prevalence") +
  scale_color_manual(values = pal[c(4, 5, 2)], name = "Group") +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), guide = F) +
  gg_theme

fd_spec_plot2 <- ggplot(data = filter(fd_spec , scenario %in% 5:8), aes(density, prevalence, color = group, linetype = group)) +
  geom_line(size = 1.5) +
  facet_grid(description ~ species)  +
  xlab("Density") +
  ylab("Prevalence") +
  scale_color_manual(values = pal[c(4, 5, 2)], name = "Group") +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), guide = F) +
  gg_theme

leg <- get_legend(fd_spec_plot1)

pdf("./output/two_host_sim_prev_fd_spec_102919.pdf")
plot_grid(fd_spec_plot1 + theme(legend.position = "none"), leg, nrow = 2, rel_heights = c(1, 0.05))

plot_grid(fd_spec_plot2 + theme(legend.position = "none"), leg, nrow = 2, rel_heights = c(1, 0.05))
dev.off()


#### dd generalist transmission rates ####

# species 1 is  the host
# species 2 is the non-host

# intraspecific transmission rates  
B11 <- 0.2
B22 <- 0.2   

# interspecific transmission rates  
B12 <- 0.1
B21 <- 0.1


#### scenario 1 (no change),  dd, generalist ####

# initial values for S
S10_s1 <- rep(300, npoints)
S20_s1 <- rep(300, npoints)

# run simulation
s1_dd_gen <- sim_fun(S10_s1, S20_s1, npoints, "./output/two_host_sim_scen1_dd_gen_102919.pdf", scenario = 1, transmission = "dd", spillover = 2)


#### scenario 2 (spp 1 density),  dd, generalist ####

# initial values for S
S10_s2 <- seq(0, 300, length.out = npoints)
S20_s2 <- rep(10, npoints)

# run simulation
s2_dd_gen <- sim_fun(S10_s2, S20_s2, npoints, "./output/two_host_sim_scen2_dd_gen_102919.pdf", scenario = 2, transmission = "dd", spillover = 2)


#### scenario 3 (spp 2 density),  dd, generalist ####

# initial values for S
S10_s3 <- rep(10, npoints)
S20_s3 <- seq(0, 300, length.out = npoints)

# run simulation
s3_dd_gen <- sim_fun(S10_s3, S20_s3, npoints, "./output/two_host_sim_scen3_dd_gen_102919.pdf", scenario = 3, transmission = "dd", spillover = 2)


#### scenario 4 (spp 1 rel. abundance),  dd, generalist ####

# initial values for S
S10_s4 <- seq(0, 300, length.out = npoints)
S20_s4 <- rep(150, npoints)

# run simulation
s4_dd_gen <- sim_fun(S10_s4, S20_s4, npoints, "./output/two_host_sim_scen4_dd_gen_102919.pdf", scenario = 4, transmission = "dd", spillover = 2)


#### scenario 5 (spp 2 rel. abundance),  dd, generalist ####

# initial values for S
S10_s5 <- rep(150, npoints)
S20_s5 <- seq(0, 300, length.out = npoints)

# run simulation
s5_dd_gen <- sim_fun(S10_s5, S20_s5, npoints, "./output/two_host_sim_scen5_dd_gen_102919.pdf", scenario = 5, transmission = "dd", spillover = 2)


#### scenario 6 (all),  dd, generalist ####

# initial values for S
S10_s6 <- seq(0, 300, length.out = npoints)
S20_s6 <- seq(0, 300, length.out = npoints)

# run simulation
s6_dd_gen <- sim_fun(S10_s6, S20_s6, npoints, "./output/two_host_sim_scen6_dd_gen_102919.pdf", scenario = 6, transmission = "dd", spillover = 2)


#### scenario 7 (total density 1),  dd, generalist ####

# initial values for S
S10_s7 <- seq(0, 300, length.out = npoints)
S20_s7 <- seq(300, 0, length.out = npoints)

# run simulation
s7_dd_gen <- sim_fun(S10_s7, S20_s7, npoints, "./output/two_host_sim_scen7_dd_gen_102919.pdf", scenario = 7, transmission = "dd", spillover = 2)


#### scenario 8 (total density 2),  dd, generalist ####

# initial values for S
S10_s8 <- seq(300, 0, length.out = npoints)
S20_s8 <- seq(0, 300, length.out = npoints)

# run simulation
s8_dd_gen <- sim_fun(S10_s8, S20_s8, npoints, "./output/two_host_sim_scen8_dd_gen_102919.pdf", scenario = 8, transmission = "dd", spillover = 2)


#### density-dependent generalist figure ####

# combine data, add scenario descriptions
dd_gen <- rbind(s1_dd_gen, s2_dd_gen, s3_dd_gen, s4_dd_gen, s5_dd_gen, s6_dd_gen, s7_dd_gen, s8_dd_gen) %>% 
  as_tibble() %>%
  full_join(scen) %>%
  gather(key = "species", value = "density", -c(group, prevalence, scenario, transmission, spillover, description)) %>%
  mutate(species = recode(species, S1_dens = "species 1", S2_dens = "species 2"),
         group = recode(group, prev = "combined", prev1 = "species 1", prev2 = "species 2") %>% factor(levels = c("species 1", "species 2", "combined")))

# create figure
dd_gen_plot1 <- ggplot(data = filter(dd_gen , scenario %in% 1:4), aes(density, prevalence, color = group, linetype = group)) +
  geom_line(size = 1.5) +
  facet_grid(description ~ species)  +
  xlab("Density") +
  ylab("Prevalence") +
  scale_color_manual(values = pal[c(4, 5, 2)], name = "Group") +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), guide = F) +
  gg_theme

dd_gen_plot2 <- ggplot(data = filter(dd_gen , scenario %in% 5:8), aes(density, prevalence, color = group, linetype = group)) +
  geom_line(size = 1.5) +
  facet_grid(description ~ species)  +
  xlab("Density") +
  ylab("Prevalence") +
  scale_color_manual(values = pal[c(4, 5, 2)], name = "Group") +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), guide = F) +
  gg_theme

leg <- get_legend(dd_gen_plot1)

pdf("./output/two_host_sim_prev_dd_gen_102919.pdf")
plot_grid(dd_gen_plot1 + theme(legend.position = "none"), leg, nrow = 2, rel_heights = c(1, 0.05))

plot_grid(dd_gen_plot2 + theme(legend.position = "none"), leg, nrow = 2, rel_heights = c(1, 0.05))
dev.off()


#### fd generalist transmission rates ####

# species 1 is  the host
# species 2 is the non-host

# intraspecific transmission rates  
B11 <- 0.2 * 300
B22 <- 0.2 * 300   

# interspecific transmission rates  
B12 <- 0.1 * 300
B21 <- 0.1 * 300


#### scenario 1 (no change),  fd, generalist ####

# initial values for S
S10_s1 <- rep(300, npoints)
S20_s1 <- rep(300, npoints)

# run simulation
s1_fd_gen <- sim_fun(S10_s1, S20_s1, npoints, "./output/two_host_sim_scen1_fd_gen_102919.pdf", scenario = 1, transmission = "fd", spillover = 2)


#### scenario 2 (spp 1 density),  fd, generalist ####

# initial values for S
S10_s2 <- seq(0, 300, length.out = npoints)
S20_s2 <- rep(10, npoints)

# run simulation
s2_fd_gen <- sim_fun(S10_s2, S20_s2, npoints, "./output/two_host_sim_scen2_fd_gen_102919.pdf", scenario = 2, transmission = "fd", spillover = 2)


#### scenario 3 (spp 2 density),  fd, generalist ####

# initial values for S
S10_s3 <- rep(10, npoints)
S20_s3 <- seq(0, 300, length.out = npoints)

# run simulation
s3_fd_gen <- sim_fun(S10_s3, S20_s3, npoints, "./output/two_host_sim_scen3_fd_gen_102919.pdf", scenario = 3, transmission = "fd", spillover = 2)


#### scenario 4 (spp 1 rel. abundance),  fd, generalist ####

# initial values for S
S10_s4 <- seq(0, 300, length.out = npoints)
S20_s4 <- rep(150, npoints)

# run simulation
s4_fd_gen <- sim_fun(S10_s4, S20_s4, npoints, "./output/two_host_sim_scen4_fd_gen_102919.pdf", scenario = 4, transmission = "fd", spillover = 2)


#### scenario 5 (spp 2 rel. abundance),  fd, generalist ####

# initial values for S
S10_s5 <- rep(150, npoints)
S20_s5 <- seq(0, 300, length.out = npoints)

# run simulation
s5_fd_gen <- sim_fun(S10_s5, S20_s5, npoints, "./output/two_host_sim_scen5_fd_gen_102919.pdf", scenario = 5, transmission = "fd", spillover = 2)


#### scenario 6 (all),  fd, generalist ####

# initial values for S
S10_s6 <- seq(0, 300, length.out = npoints)
S20_s6 <- seq(0, 300, length.out = npoints)

# run simulation
s6_fd_gen <- sim_fun(S10_s6, S20_s6, npoints, "./output/two_host_sim_scen6_fd_gen_102919.pdf", scenario = 6, transmission = "fd", spillover = 2)


#### scenario 7 (total density 1),  fd, generalist ####

# initial values for S
S10_s7 <- seq(0, 300, length.out = npoints)
S20_s7 <- seq(300, 0, length.out = npoints)

# run simulation
s7_fd_gen <- sim_fun(S10_s7, S20_s7, npoints, "./output/two_host_sim_scen7_fd_gen_102919.pdf", scenario = 7, transmission = "fd", spillover = 2)


#### scenario 8 (total density 2),  fd, generalist ####

# initial values for S
S10_s8 <- seq(300, 0, length.out = npoints)
S20_s8 <- seq(0, 300, length.out = npoints)

# run simulation
s8_fd_gen <- sim_fun(S10_s8, S20_s8, npoints, "./output/two_host_sim_scen8_fd_gen_102919.pdf", scenario = 8, transmission = "fd", spillover = 2)


#### frequency-dependent generalist figure ####

# combine data, afd scenario descriptions
fd_gen <- rbind(s1_fd_gen, s2_fd_gen, s3_fd_gen, s4_fd_gen, s5_fd_gen, s6_fd_gen, s7_fd_gen, s8_fd_gen) %>% 
  as_tibble() %>%
  full_join(scen) %>%
  gather(key = "species", value = "density", -c(group, prevalence, scenario, transmission, spillover, description)) %>%
  mutate(species = recode(species, S1_dens = "species 1", S2_dens = "species 2"),
         group = recode(group, prev = "combined", prev1 = "species 1", prev2 = "species 2") %>% factor(levels = c("species 1", "species 2", "combined")))

# create figure
fd_gen_plot1 <- ggplot(data = filter(fd_gen , scenario %in% 1:4), aes(density, prevalence, color = group, linetype = group)) +
  geom_line(size = 1.5) +
  facet_grid(description ~ species)  +
  xlab("Density") +
  ylab("Prevalence") +
  scale_color_manual(values = pal[c(4, 5, 2)], name = "Group") +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), guide = F) +
  gg_theme

fd_gen_plot2 <- ggplot(data = filter(fd_gen , scenario %in% 5:8), aes(density, prevalence, color = group, linetype = group)) +
  geom_line(size = 1.5) +
  facet_grid(description ~ species)  +
  xlab("Density") +
  ylab("Prevalence") +
  scale_color_manual(values = pal[c(4, 5, 2)], name = "Group") +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), guide = F) +
  gg_theme

leg <- get_legend(fd_gen_plot1)

pdf("./output/two_host_sim_prev_fd_gen_102919.pdf")
plot_grid(fd_gen_plot1 + theme(legend.position = "none"), leg, nrow = 2, rel_heights = c(1, 0.05))

plot_grid(fd_gen_plot2 + theme(legend.position = "none"), leg, nrow = 2, rel_heights = c(1, 0.05))
dev.off()


#### save prevalence output ####

prev_save <- rbind(dd_spec, fd_spec, dd_spill, fd_spill, dd_gen, fd_gen)
write_csv(prev_save, "./output/two_host_dd_fd_prev_102919.csv")


#### examine output more closely ####

filter(dd_spill, scenario == 6 & density < 10 & species == "species 1")
