# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "fastDummies", "caret", "protr")


# Set working directory
setwd("C:/Users/tabie/OneDrive/Documents/GitHub/TARA-thiolases/")

# Read in 50 TARA seqs
# 123_OleA_no_NA_temps.fasta
sqs <- readAAStringSet("data/84_OleA_no_NA_temps.fasta")
query_fils <- sapply(1:length(sqs), function(x) {tempfile(pattern = "", fileext = ".fasta")})
sapply(1:length(sqs), function(x) {writeXStringSet(sqs[x], query_fils[x])})

# Read in a separate function
source("src/extract_10angstrom_residues.R")

extract_84_list <- lapply(1:length(sqs), function(x) { extract_10angstrom(query_fils[x]) })
extract_84_df <- data.frame(matrix(unlist(extract_84_list), nrow = length(extract_84_list), byrow=T), 
                            stringsAsFactors=FALSE)

# Write un-encoded data frame to file
separated_df <- tibble(names(sqs)) %>%
  bind_cols(extract_84_df)
colnames(separated_df)[1] <- "nams"

write_csv(separated_df, "data/residue_extraction/10_angstrom_84_OleA_aa_separated.csv")
onehot <- fastDummies::dummy_cols(extract_84_df, remove_selected_columns = T, ignore_na = T)

# First, we need to learn about one-hot encoding!
# One-hot encoding is used to encode categorical variables as binary variables
# a new binary variable column is added for each unique categorical value.

onehot_df <- tibble(names(sqs)) %>%
  bind_cols(onehot)
colnames(onehot_df)[1] <- "nams"

write_csv(onehot_df, "data/residue_extraction/10_angstrom_84_OleA_one_hot_aa_features_extracted.csv")

# Now let's try it with encoding each variable based on its physicochemical properties

source("src/convert_seq_5aap.R")
extract_feat_list <- lapply(1:length(extract_84_list), function(x) { convert_seq_5aap(extract_84_list[[x]]) })
extract_feat_df <- data.frame(matrix(unlist(extract_feat_list), nrow = length(extract_feat_list), byrow=T), stringsAsFactors=FALSE)

colnames(extract_feat_df) <- paste0("pos", "_", sort(rep(1:50, times = 5)), "_", rep(c("polarity", "secondary_structure", "size", "codon_diversity", "charge"), times = 50))

feat_df <- tibble(names(sqs)) %>%
  bind_cols(extract_feat_df)
colnames(feat_df)[1] <- "nams"
head(feat_df)

write_csv(feat_df, "data/residue_extraction/10_angstrom_84_OleA_physical_aa_features_extracted.csv")

## Challenges

# Read Atchley et al. 2005 PNAS (in lit/encoding_amino_acids). Inspect the file 5_aa_properties.csv
# Describe the encoding occuring in convert_seq_5aap.R

# How does convert_seq_5aap compare to convert_seq_15aap.R

# Rewrite code to calculate properties for all 80 sequences (or however many you have growth temp information for)

# Compare one-hot to physicochemical properties to categorical amino acids

# Rename amino acid features to correspond to what they are

colnames(extract_feat_df) <- paste0("pos", "_", sort(rep(1:50, times = 5)), "_", rep(c("polarity", "secondary_structure", "size", "codon_diversity", "charge"), times = 50))

?rep
# Test different feature sets with machine learning -> both classification and regression

