# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "fastDummies", "caret", "protr")


# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/TARA-thiolases/")

# Read in 50 TARA seqs
sqs <- readAAStringSet("data/50_TARA_psychro_thermo_unaligned.fasta")
query_fils <- sapply(1:length(sqs), function(x) {tempfile(pattern = "", fileext = ".fasta")})
sapply(1:length(sqs), function(x) {writeXStringSet(sqs[x], query_fils[x])})

# Read in a separate function
source("src/extract_12angstrom_residues.R")

extract_84_list <- lapply(1:length(sqs), function(x) { extract_12angstrom(query_fils[x]) })
extract_84_df <- data.frame(matrix(unlist(extract_84_list), nrow = length(extract_84_list), byrow=T), 
                            stringsAsFactors=FALSE)

# Write un-encoded data frame to file
onehot <- fastDummies::dummy_cols(extract_84_df, remove_selected_columns = T, ignore_na = T)

# First, we need to learn about one-hot encoding!
# One-hot encoding is used to encode categorical variables as binary variables
# a new binary variable column is added for each unique categorical value.

onehot_df <- tibble(names(sqs)) %>%
  bind_cols(onehot)
colnames(onehot_df)[1] <- "nams"

write_csv(onehot_df, "data/12_angstrom_one_hot_aa_features_extracted.csv")

# Now let's try it with encoding each variable based on its physicochemical properties

source("src/convert_seq_15aap.R")x
extract_feat_list <- lapply(1:length(extract_84_list), function(x) { convert_seq_15aap(extract_84_list[[x]]) })
extract_feat_df <- data.frame(matrix(unlist(extract_feat_list), nrow = length(extract_feat_list), byrow=T), stringsAsFactors=FALSE)

feat_df <- tibble(names(sqs)) %>%
  bind_cols(extract_feat_df)
colnames(feat_df)[1] <- "nams"
head(feat_df)

write_csv(feat_df, "data/12_angstrom_physical_aa_features_extracted.csv")


## Challenges

## Challenge 1. 

# Read the NRPSPredictor2 paper

# Expand and calculate for all 80 sequences (or however many you have growth temp information for)

# Compare one-hot to physicochemical properties to categorical amino acids

# Rename amino acid features to correspond to what they are

# Test different feature sets with machine learning -> both classification and regression
