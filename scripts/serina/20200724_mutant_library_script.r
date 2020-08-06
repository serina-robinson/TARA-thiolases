# Install packages
pacman::p_load("tidyverse", "readxl")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/TARA-thiolases/")

# Read in the mutant library
mut <- read_excel("data/Batch220D_Wackett_final_JGI_order.xlsx")

# Which are the Xanthomonas mutants
xan <- mut %>%
  dplyr::filter(grepl("NP_635607_Xanthomonas_campestris", user_id)) %>%
  dplyr::mutate(mutation = substr(user_id, 1, 5))
xan$mutation

# Challenges

# Make a table of the mutations