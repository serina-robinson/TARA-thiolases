pacman::p_load("ggplot2", "tidyverse", "plotly")
setwd("C:/Users/tabie/OneDrive/Documents/GitHub/TARA-thiolases/")


dat <- read_csv("data/123_OleA_allstats.csv")


ggplot(dat) +
  geom_point(aes(RKrat, temperature))
