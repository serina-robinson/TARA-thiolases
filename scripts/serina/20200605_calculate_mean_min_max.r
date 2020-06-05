# Load packages
pacman::p_load("Biostrings", "tidyverse", "GGally")

# Set working directory
setwd("C:/Users/tabie/OneDrive/Documents/GitHub/TARA-thiolases/")

# Read in the TARA sequences only
tara <- readAAStringSet('data/TARA_marine_sequences_to_order.fasta') # these are the 30 sequences from TARA oceans
names(tara)

# Read in the sequence metadata
annot_df <- read_csv("data/TARA_tree_annotation_df.csv")
view(annot_df) # View the data frame in a different window
colnames(annot_df) # What sorts of variables are in this data frame
dim(annot_df) # Check the dimensions, 244 rows and 55 columns

# But wait! we only have 30 sequences we actually ordered 
# Subset the data frame...
twenty_nine <- annot_df[annot_df$label %in% names(tara),]
dim(twenty_nine) # 29 rows and 55 columns

### Challenge 1. Add the 30th TARA sequence (like yesterday) to make a data frame with 30 rows

# Now that we've got our data frame...let's calculate some statistics
# We'll start with some variables of interest
dat <- annot_df # you can change this to thirty, or twenty_nine if you just want to look at some seqs

# Subset the data frame only for variables you're interested in
vars <- c("temperature", "oxygen", "no3", "po4", "iron_5m") # subsetting for variables of interest
sub_df <- dat[,colnames(dat) %in% vars]

# Calculate the min, median, mean, max temperatures
summary(sub_df$temperature)

# If you want to just look at one of these
mean(sub_df$temperature, na.rm = T)

# This is the same as pulling out the 4th item from the summary 
summary(sub_df$temperature)[4]
mean(sub_df$temperature, na.rm = T) == summary(sub_df$temperature)[4]

# Calculate correlations between variables
ggpairs(sub_df) +
  theme_bw()

# Challenges..play around. 
# What does okubo_weiss mean? What about nitracline or brunt_vaisala?
# Hint: The literature reading for today is your friend. Also Google.
# Choose some interesting looking correlations or plots 
# Make some ggplots for lab meeting on Tuesday!
# Look up what some of hte 