# Versions
# full50_raw: original combination of 30 tara sequences, 10 psychrophiles, and 10 thermophiles
# full50_1: designated tara sequences as psychrphile, thermophile, or mesophile
# full50_2: added jittered coords for 20 non-tara samples
# full50_3: designated temperature definition for tara samples with temp value as "sampled"
#           cleaned up column names
# full50_4: designated polar as |lat| > 60 and filled in NAs

library(tidyverse)
library(janitor)
setwd("C:/Users/tabie/OneDrive/Documents/GitHub/TARA-thiolases/")


# add sampling designation for temperature in tara samples
full50_2 <- read_csv("data/full50_2.csv")
view(full50_2)
full50_2[1:30, "temperature_definition"] <- "S"
full50_2 %>% 
  select("temperature", "temperature_definition", "temperature_range") %>% 
  view()

# remove sampling designations for samples with no recorded temperature

full50_3 <- full50_2 %>% 
 mutate(temperature_definition = case_when(is.na(temperature) ~ NA_character_,
           TRUE ~ temperature_definition))
view(select(full50_3, temperature_definition, temperature))

# clean up column names
full50_3 <- full50_3 %>% 
  janitor::clean_names()
view(full50_3)

# write new file
write_csv(full50_3, "data/full50_3.csv")

# designate polar or non polar
full50_3 <- read_csv("data/full50_3.csv")
full50_4 <- full50_3 %>% 
  mutate(polar = case_when(abs(lat) >= 60 ~ "polar",
                           abs(lat) < 60 ~ "non polar",
                           TRUE ~ polar))
view(select(full50_4, lat, polar))

# write new file
write_csv(full50_4, "data/full50_4.csv")
