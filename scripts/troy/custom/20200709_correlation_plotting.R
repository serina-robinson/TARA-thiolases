# Install packages
pacman::p_load("tidyverse", "ggplot2", "GGally")

# Set working directory
setwd("C:/Users/tabie/OneDrive/Documents/GitHub/TARA-thiolases/")

# Load in data
dat <- read_csv("data/full50_allstats.csv")
dat$hydrophob_full

# RK ratio lowered in psychro
plot1 <- ggplot(dat) +
  geom_boxplot(aes(x = RKrat, col = temperature_range)) +
  theme_classic()
plot1


# NQ ratio lower in thermophile
plot2 <- ggplot(dat) +
  geom_point(aes(x = NQrat, y = temperature, col = temperature_range)) +
  theme_classic()
plot2


# disordered_res decreased in thermophile
plot3 <- ggplot(dat) +
  geom_point(aes(x = disordered_res, y = temperature, col = temperature_range)) +
  theme_classic()
plot3


# ordered_res hard to tell
plot4 <- ggplot(dat) +
  geom_point(aes(x = ordered_res, y = temperature, col = temperature_range)) +
  theme_classic()
plot4


# hydrophobicity hard to tell
plot5 <- ggplot(dat) +
  geom_point(aes(x = hydrophob_full, y = temperature, col = temperature_range)) +
  theme_classic()
plot5


# boman interactions hard to tell
plot6 <- ggplot(dat) +
  geom_point(aes(x = boman_interactions, y = temperature, col = temperature_range)) +
  theme_classic()
plot6


# hmoment hard to tell
plot7 <- ggplot(dat) +
  geom_point(aes(x = hmoment, y = temperature, col = temperature_range)) +
  theme_classic()
plot7

# instab hard to tell
plot8 <- ggplot(dat) +
  geom_point(aes(x = instab, y = temperature, col = temperature_range)) +
  theme_classic()
plot8


# pi_prot hard to tell
plot9 <- ggplot(dat) +
  geom_point(aes(x = pi_prot, y = temperature, col = temperature_range)) +
  theme_classic()
plot9


# Small hard to tell
plot10 <- ggplot(dat) +
  geom_point(aes(x = Small, y = temperature, col = temperature_range)) +
  theme_classic()
plot10


# Charged hard to tell
plot11 <- ggplot(dat) +
  geom_point(aes(x = Charged, y = temperature, col = temperature_range)) +
  theme_classic()
plot11


# Nonpolar seems to be higher in thermophiles
plot12 <- ggplot(dat) +
  geom_point(aes(x = Nonpolar, y = temperature, col = temperature_range)) +
  theme_classic()
plot12


# Polar seems to be lower in thermophiles
plot13 <- ggplot(dat) +
  geom_point(aes(x = Polar, y = temperature, col = temperature_range)) +
  theme_classic()
plot13


# 
plot14 <- ggplot(dat) +
  geom_point(aes(x = , y = temperature, col = temperature_range)) +
  theme_classic()
plot14



vars <- c("temperature", "temperature_range", "topt", "RKrat")
subdat <- dat[,colnames(dat) %in% vars]
ggpairs(subdat) +
  theme_bw()





