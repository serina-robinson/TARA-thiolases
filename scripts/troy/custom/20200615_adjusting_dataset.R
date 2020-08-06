# Versions
# full50_raw: original combination of 30 tara sequences, 10 psychrophiles, and 10 thermophiles
# full50_1: designated tara sequences as psychrophile, thermophile, or mesophile
# full50_2: added jittered coords for 20 non-tara samples
# full50_3: designated temperature definition for tara samples with temp value as "sampled"
#           cleaned up column names
# full50_4: designated polar as |lat| > 60 and filled in NAs
# full50_5: set temperature ranges
# full50_6: added nams: combined labels for 2 datasets

library(tidyverse)
library(janitor)
library(Biostrings)
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

# set temperature ranges

full50_4 <- read_csv("data/full50_4.csv")
full50_5 <- full50_4 %>% 
  mutate(temperature_range = case_when(temperature >= 40 ~ "Thermophilic",
                           temperature < 15 ~ "Psychrophilic",
                           is.na(temperature) ~ NA_character_,
                           TRUE ~ "Mesophilic"))
view(select(full50_5, temperature, temperature_range))

# write new file
write_csv(full50_5, "data/full50_5.csv")

# add nams
full50_6 <- read_csv("data/full50_5.csv") %>% 
  mutate(nams = case_when(!is.na(label) ~ label,
                          is.na(label) ~ paste0(word(ole_a_patric, sep = 'fig\\|', -1),
                                                "_",
                                                genus,
                                                "_",
                                                case_when(temperature_range == "Psychrophilic" ~ "psychrophile",
                                                          temperature_range == "Mesophilic" ~ "mesophile",
                                                          temperature_range == "Thermophilic" ~ "thermophile",
                                                          TRUE ~ "thermophile"))))
write_csv(full50_6, "data/full50_6.csv")

# make file of 123 OleA that has only ones with temp data
pacman::p_load("sequinr")
all123 <- read_csv("data/123_OleA_temps.csv")
noNAs <- all123 %>% 
  filter(!is.na(temperature))
write_csv(noNAs, "data/84_OleA_temps_noNAs.csv")
noNAs2 <- noNAs %>% 
  select(nams, sequence)
writeXStringSet(noNAs2$sequence, noNAs2$nams, "data/123_OleA_no_NA_temps.fasta")
test <- AAStringSet(noNAs2$sequence)
names(test) <- noNAs2$nams
writeXStringSet(test, "data/123_OleA_no_NA_temps.fasta")
fasta123 <- readAAStringSet("data/123_OleA.fasta")


