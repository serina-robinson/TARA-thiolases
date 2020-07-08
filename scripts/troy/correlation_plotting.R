# Install packages
pacman::p_load("tidyverse", "ggplot2")

# Set working directory
setwd("C:/Users/tabie/OneDrive/Documents/GitHub/TARA-thiolases/")

# Load in data
dat <- read_csv("data/full50_allstats.csv")


plot1 <- ggplot(dat) +
  geom_point(aes(x = depth_m, y = temperature, col = polar)) +
  theme_classic()
plot1
