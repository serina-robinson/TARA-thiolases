# Load packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "ggseqlogo", "cowplot")

# Set working directory
setwd("C:/Users/tabie/OneDrive/Documents/GitHub/TARA-thiolases")

# Read in 123 OleA data
dat <- read_csv("data/84_OleA_temps_noNAs.csv")
colnames(dat)

# Histogram of temperatures
pdf("output/123_OleA_temperature_histogram.pdf", height = 4)
ggplot(dat) +
  geom_histogram(aes(dat$temperature)) +
  theme_classic() + 
  labs(x = "Temperature", y = "Count")
dev.off()
# Bar plot of temperature ranges
dat2 <- dat %>% 
  mutate(temp_category = case_when(temperature <= 15 ~ "< 15",
                                   temperature >= 30 ~ "> 30",
                                   TRUE ~ "15 - 30"))
dat2$temp_category <-  factor(dat2$temp_category, levels = c("< 15", "15 - 30", "> 30"))

pdf("output/123_OleA_temperature_barplot.pdf", height = 4)
ggplot(dat2) +
  geom_bar(aes(temp_category)) +
  theme_classic() + 
  labs(x = "Temperature", y = "Count")
dev.off()

