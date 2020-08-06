# Install packages
pacman::p_load("tidyverse", "readxl")

# Set working directory
setwd("C:/Users/tabie/OneDrive/Documents/GitHub/TARA-thiolases/")

# Read in the mutant library
mut <- read_excel("data/Batch220D_Wackett_final_JGI_order.xlsx")

# Which are the Xanthomonas mutants
xan <- mut %>%
  dplyr::filter(grepl("NP_635607_Xanthomonas_campestris", user_id)) %>%
  dplyr::mutate(mutation = substr(user_id, 1, 5))
xan$mutation

# Challenges

# retry categorical machine learning with 2/3 groups and channels
# try class imbalance things

# look at important variables for categorical and hope they line up with regression

# look at overlap between residues important in model and residues in mutant library

###### read goblirsch paper

# make seq logo of positions that are both important and 
# in mutant library for 84 seqs
# and contingency table

# look in pymol at any weird ones

# make slides
## intro: why machine learning
## global properties didn't work so well
## table for regression: channels worked best
## obs vs predicted temp plot too
## channels for categorical too
## maybe use mutant library for experiments
## show which we are interested in with sequence logo
## future: purify and look at thermostability/activity/any ideas?

# identify if any of the mutants in the library are worth testing
# or do we need to make new ones to test

# which ones do we want to test and why


# Make a table of the mutations