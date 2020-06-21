# Install packages
pacman::p_load("tidyverse", "Peptides", "seqinr")

# Set working directory
setwd("C:/....")

# Read in the protein sequences
sqs <-read.alignment("data/50_TARA_psychro_thermo_unaligned.fasta", format = "fasta")
sqs$seq <- toupper(sqs$seq)

## Calculate protein sequence properties using the Peptides package
# Look at the Peptides documentation for explanations
# Reference manual: https://cran.r-project.org/web/packages/Peptides/Peptides.pdf
# Github repo: https://github.com/dosorio/Peptides
# Also available in the docs/Peptides_package_protein_properties.xlsx

# Cruciani 
crciani <- sapply(sqs$seq, crucianiProperties)
cdf <- do.call(rbind.data.frame, crciani)
colnames(cdf) <- names(crciani[[1]])

# Fasgai
fsgai <- sapply(sqs$seq, fasgaiVectors)
fdf <- do.call(rbind.data.frame, fsgai)
colnames(fdf) <- names(fsgai[[1]])

# ST scale
st_scale <- sapply(sqs$seq, stScales)
stdf <- do.call(rbind.data.frame, st_scale)
colnames(stdf) <- names(st_scale[[1]])

# T-scale
t_scale <- sapply(sqs$seq, tScales)
tdf <- do.call(rbind.data.frame, t_scale)
colnames(tdf) <- names(t_scale[[1]])

# VHSE
vhse_scale <- sapply(sqs$seq, vhseScales)
vdf <- do.call(rbind.data.frame, vhse_scale)
colnames(vdf) <- names(vhse_scale[[1]])

# Z-scale
z_scale <- sapply(sqs$seq, zScales)
zdf <- do.call(rbind.data.frame, z_scale)
colnames(zdf) <- names(z_scale[[1]])

# Kidera
kdera <- sapply(sqs$seq, kideraFactors)
kdf <- do.call(rbind.data.frame, kdera)
colnames(kdf) <- names(kdera[[1]])

# Combine everything into one data frame
sqdf <- data.frame(cbind(sqs$nam, sqs$seq), stringsAsFactors = F)
colnames(sqdf) <- c("nams", "sqs")
sqs$nam

combdf <- sqdf %>%
  dplyr::mutate(acc = word(nams, sep = "\\.1", 1)) %>%
  dplyr::mutate(aliphat_full = aIndex(sqs)) %>%
  dplyr::mutate(hydrophob_full = hydrophobicity(sqs)) %>%
  dplyr::mutate(boman_interactions = boman(sqs)) %>%
  dplyr::mutate(hmoment = hmoment(sqs)) %>%
  dplyr::mutate(instab = instaIndex(sqs)) %>%
  dplyr::mutate(lngth_prot = lengthpep(sqs)) %>%
  dplyr::mutate(mw_prot = mw(sqs)) %>%
  dplyr::mutate(pi_prot = pI(sqs)) %>%
  bind_cols(., zdf, tdf, cdf, vdf, stdf, fdf, kdf)
combdf


# Read in the protein temperature predictions for your sequences
# Used the TOMER tool described in the following paper
# https://www.biorxiv.org/content/10.1101/2020.05.06.081737v1.full
tomerdf <- read_csv("data/tomer/50_TARA_psychro_thermo_results.csv") %>%
  janitor::clean_names() %>%
  dplyr::select(-x1, -std_err)
# topt is the predicted optimal temperature for the proteins (in degrees Celsius)

# Challenge 1. Read about the Peptides package and do your best to understand some of the
# indices that are being calculated.

# Challenge 2. Join the tomerdf with combdf. Take a look at the distribution of topt. 
# What do you notice?

# Challenge 3. Read in full50 (or some version of it) and join with df from Challenge 2.

# Challenge 4. Plot relationships between protein predicted optimal temperature and 
# different protein properties like hydrophobicity (Kyte-Doolittle) and instability indices

# Challenge 5. Use Ggally function from 20200605 to investigate associations between all protein properties
# and temperature properties (both optimal growth temp of organism from full50 and protein optimum from tomerdf)

# Super challenge 6. Make a new page in your web app for interactive plotting of protein properties
# relative to temperature (both optimal growth temp of organism from full50 and protein optimum from tomerdf)

