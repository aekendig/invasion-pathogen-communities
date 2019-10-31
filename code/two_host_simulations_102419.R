### goal: explore effects of varying density and relative abundance in a spillover model

# model: Population dynamics of pathogens with multiple host species, Dobson 2004, The American Naturalist

# notes: Edited model to make the populations at carrying capacity (deaths become births)


#### set-up ####

# clear all existing data
rm(list=ls())

# load libraries
library(deSolve)
library(tidyverse)


#### default parameters ####

# species 1 is  the reservoir host
# species 2 is the nonreservoir host

# death rates
d1 <- 0.1
d2 <- 0.1
  
# intraspecific transmission rates  
B11 <- 0.2
B22 <- 0   

# interspecific transmission rates  
B12 <- 0
B21 <- 0.1
  
# infection-induced mortality  
alpha1 <- 0.01
alpha2 <- 0.01

# recover rate
theta1 <- 0.2
theta2 <- 0.2

# initial values
S10 <- 100
I10 <- 1
R10 <- 0
  
S20 <- 100
I20 <- 0
R20 <- 0

# simulation time
simtime <- 100


#### model ####
mod <- function(d1, B11, B12, alpha1, theta1, d2, B22, B21, alpha2, theta2, S10, I10, R10, S20, I20, R20, simtime){
  
  # parameters
  parms<-list(d1 = d1, B11 = B11, B12 = B12, alpha1 = alpha1, theta1 = theta1, d2 = d2, B22 = B22, B12 = B12, alpha2 = alpha2, theta2 = theta2, S10 = S10, I10 = I10, R10 = R10, S20 = S20, I20 = I20, R20 = R20, simtime = simtime)
  
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
  

#### simulations ####

# default values (each S0 = 100)
out_def <- mod(d1, B11, B12, alpha1, theta1, d2, B22, B21, alpha2, theta2, S10, I10, R10, S20, I20, R20, simtime) 

out_def_l <- out_def %>% 
  gather(key = "pop", value = "abundance", -time)

ggplot(out_def_l, aes(time, abundance, color = pop)) + 
  geom_line()

# filter(out_def_l, time < 10) %>%
#   ggplot(aes(time, abundance, color = pop)) + 
#   geom_line()

tail(out_def)


# increase reservoir host (80% reservoir, 20% nonreservoir)
S10 <- 160
S20 <- 40

out_8020 <- mod(d1, B11, B12, alpha1, theta1, d2, B22, B21, alpha2, theta2, S10, I10, R10, S20, I20, R20, simtime)

out_8020_l <- out_8020 %>% 
  gather(key = "pop", value = "abundance", -time)

ggplot(out_8020_l, aes(time, abundance, color = pop)) + 
  geom_line()

tail(out_8020)


# increase total density (50% reservoir, 50% nonreservoir)
S10 <- 300
S20 <- 300

out_high <- mod(d1, B11, B12, alpha1, theta1, d2, B22, B21, alpha2, theta2, S10, I10, R10, S20, I20, R20, simtime)

out_high_l <- out_high %>% 
  gather(key = "pop", value = "abundance", -time)

ggplot(out_high_l, aes(time, abundance, color = pop)) + 
  geom_line()

tail(out_high)


# increase reservoir host relative abundance and density (80% reservoir, 20% nonreservoir)
S10 <- 480
S20 <- 120

out_high_8020 <- mod( d1, B11, B12, alpha1, theta1, d2, B22, B21, alpha2, theta2, S10, I10, R10, S20, I20, R20, simtime)

out_high_8020_l <- out_high_8020 %>% 
  gather(key = "pop", value = "abundance", -time)

ggplot(out_high_8020_l, aes(time, abundance, color = pop)) + 
  geom_line()

tail(out_high_8020)


# increase reservoir host relative abundance and density while keeping native the same (83% reservoir,17% nonreservoir)
S10 <- 500
S20 <- 100

out_8317 <- mod(d1, B11, B12, alpha1, theta1, d2, B22, B21, alpha2, theta2, S10, I10, R10, S20, I20, R20, simtime)

out_8317_l <- out_8317 %>% 
  gather(key = "pop", value = "abundance", -time)

ggplot(out_8317_l, aes(time, abundance, color = pop)) + 
  geom_line()

tail(out_8317)

# combine datasets
out_def$model = "5050"
out_8020$model = "8020"
out_high$model = "high 5050"
out_high_8020$model = "high 8020"
out_8317$model = "high 8317"

out_df = rbind(out_def, out_8020, out_high, out_high_8020, out_8317)

# calculate prevalence
out_df %>%
  filter(time == simtime) %>%
  group_by(model) %>%
  summarise(prev1 = I1/(S1 + I1 + R1),
            prev2 = I2/(S2 + I2 + R2),
            prev = (I1 + I2)/(S1 + I1 + R1 + S2 + I2 + R2))
# slight changes in prevalence


#### vary reservoir host density ####

# nonreservoir host density
S20 = 100

# reservoir host densities
S10_l = seq(1, 10000, length.out = 40)

# empty dataframe
out_res <- data.frame(S10 = S10_l, time = NA, S1 = NA, I1 = NA, R1 = NA, S2 = NA, I2 = NA, R2 = NA)

# start pdf
pdf(file = "./output/two_host_simulations_reservoir_density_102619.pdf")

# cycle through densities
for(i in 1:40){
  
  # set reservoir host density
  S10 <- S10_l[i]
  
  temp_out <- mod(d1, B11, B12, alpha1, theta1, d2, B22, B21, alpha2, theta2, S10, I10, R10, S20, I20, R20, simtime)
  
  temp_out_l <- temp_out %>% 
    gather(key = "pop", value = "abundance", -time)
  

  print(ggplot(temp_out_l, aes(time, abundance, color = pop)) + 
    geom_line())
  
  out_res[i, 2:8] <- filter(temp_out, time == simtime)
  
}

# close pdf
dev.off()

# calculate prevalence
out_res_prev <- out_res %>%
  group_by(S10) %>%
  summarise(prev1 = I1/(S1 + I1 + R1),
            prev2 = I2/(S2 + I2 + R2),
            prev = (I1 + I2)/(S1 + I1 + R1 + S2 + I2 + R2)) %>%
  ungroup() %>%
  gather(key = "group", value = "prevalence", -S10) %>%
  mutate(group = recode(group, prev1 = "reservoir", prev2 = "nonreservoir", prev = "combined"))

# figure of prevalence
ggplot(out_res_prev, aes(x = S10, y = prevalence, color = group)) +
  geom_line()


#### low reservoir host density ####

# nonreservoir host density
S20 = 100

# reservoir host densities
S10_low = seq(1, 200, length.out = 40)

# empty dataframe
out_res_low <- data.frame(S10 = S10_low, time = NA, S1 = NA, I1 = NA, R1 = NA, S2 = NA, I2 = NA, R2 = NA)

# start pdf
pdf(file = "./output/two_host_simulations_low_reservoir_density_102619.pdf")

# cycle through densities
for(i in 1:40){
  
  # set reservoir host density
  S10 <- S10_low[i]
  
  temp_out <- mod(d1, B11, B12, alpha1, theta1, d2, B22, B21, alpha2, theta2, S10, I10, R10, S20, I20, R20, simtime)
  
  temp_out_l <- temp_out %>% 
    gather(key = "pop", value = "abundance", -time)
  
  
  print(ggplot(temp_out_l, aes(time, abundance, color = pop)) + 
          geom_line())
  
  out_res_low[i, 2:8] <- filter(temp_out, time == simtime)
  
}

# close pdf
dev.off()

# calculate prevalence
out_res_prev_low <- out_res_low %>%
  group_by(S10) %>%
  summarise(prev1 = I1/(S1 + I1 + R1),
            prev2 = I2/(S2 + I2 + R2),
            prev = (I1 + I2)/(S1 + I1 + R1 + S2 + I2 + R2)) %>%
  ungroup() %>%
  gather(key = "group", value = "prevalence", -S10) %>%
  mutate(group = recode(group, prev1 = "reservoir", prev2 = "nonreservoir", prev = "combined"))

# figure of prevalence
ggplot(out_res_prev_low, aes(x = S10, y = prevalence, color = group)) +
  geom_line()


#### vary nonreservoir host density ####

# reservoir host density
S10 = 100

# reservoir host densities
S20_l = seq(1, 200, length.out = 40)

# empty dataframe
out_nonres <- data.frame(S20 = S20_l, time = NA, S1 = NA, I1 = NA, R1 = NA, S2 = NA, I2 = NA, R2 = NA)

# start pdf
pdf(file = "./output/two_host_simulations_nonreservoir_density_102619.pdf")

# cycle through densities
for(i in 1:40){
  
  # set reservoir host density
  S20 <- S20_l[i]
  
  temp_out <- mod(d1, B11, B12, alpha1, theta1, d2, B22, B21, alpha2, theta2, S10, I10, R10, S20, I20, R20, simtime)
  
  temp_out_l <- temp_out %>% 
    gather(key = "pop", value = "abundance", -time)
  
  
  print(ggplot(temp_out_l, aes(time, abundance, color = pop)) + 
          geom_line())
  
  out_nonres[i, 2:8] <- filter(temp_out, time == simtime)
  
}

# close pdf
dev.off()

# calculate prevalence
out_nonres_prev <- out_nonres %>%
  group_by(S20) %>%
  summarise(prev1 = I1/(S1 + I1 + R1),
            prev2 = I2/(S2 + I2 + R2),
            prev = (I1 + I2)/(S1 + I1 + R1 + S2 + I2 + R2)) %>%
  ungroup() %>%
  gather(key = "group", value = "prevalence", -S20) %>%
  mutate(group = recode(group, prev1 = "reservoir", prev2 = "nonreservoir", prev = "combined"))

# figure of prevalence
ggplot(out_nonres_prev, aes(x = S20, y = prevalence, color = group)) +
  geom_line()
