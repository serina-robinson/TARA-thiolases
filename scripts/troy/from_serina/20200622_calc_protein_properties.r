# Install packages
pacman::p_load("tidyverse", "Peptides", "seqinr", "ggplot2", "GGally")

# Set working directory
setwd("C:/Users/tabie/OneDrive/Documents/GitHub/TARA-thiolases/")

# Read in the protein sequences
sqs <-read.alignment("data/50_TARA_psychro_thermo_unaligned.fasta", format = "fasta")
sqs$seq <- toupper(sqs$seq)
sqs$nam
sqs$seq
## Calculate protein sequence properties using the Peptides package
# Look at the Peptides documentation for explanations
# Reference manual: https://cran.r-project.org/web/packages/Peptides/Peptides.pdf
# Github repo: https://github.com/dosorio/Peptides
# Also available in the docs/Peptides_package_protein_properties.xlsx

# Cruciani 
crciani <- sapply(sqs$seq, crucianiProperties)
cdf <- do.call(rbind.data.frame, crciani)
colnames(cdf) <- names(crciani[[1]])
crucianiProperties(sqs$seq)
cdf
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
view(combdf)


# Read in the protein temperature predictions for your sequences
# Used the TOMER tool described in the following paper
# https://www.biorxiv.org/content/10.1101/2020.05.06.081737v1.full
tomerdf <- read_csv("data/tomer/50_TARA_psychro_thermo_results.csv") %>%
  janitor::clean_names() %>%
  dplyr::select(-x1, -std_err)
# topt is the predicted optimal temperature for the proteins (in degrees Celsius)
view(tomerdf)
# Challenge 1. Read about the Peptides package and do your best to understand some of the
# indices that are being calculated.

# Challenge 2. Join the tomerdf with combdf. Take a look at the distribution of topt. 
# What do you notice?
combined <- combdf %>% 
  left_join(tomerdf, by = c("nams" = "sequence"))
view(combined)
combined[grep("Psychrobacter", combined$nams), "nams"] <- "282669.3.peg.967_Kocuria_mesophile"
combined[grep("2026735.166.peg.848_Deltaproteobacteria", combined$nams), "nams"] <- "2026735.166.peg.848_Deltaproteobacteria_psychrophile"
combined[grep("2026735.148.peg.1837_Deltaproteobacteria", combined$nams), "nams"] <- "2026735.148.peg.1837_Deltaproteobacteria_psychrophile"
combined[grep("2026809.16.peg.2167_Epsilonproteobacteri", combined$nams), "nams"] <- "2026809.16.peg.2167_Epsilonproteobacteria_psychrophile"
combined[grep("2026799.183.peg.457_Verrucomicrobia", combined$nams), "nams"] <- "2026799.183.peg.457_Verrucomicrobia_psychrophile"
combined[grep("2483033.3.peg.4057_Sedimentitalea", combined$nams), "nams"] <- "2483033.3.peg.4057_Sedimentitalea_mesophile"

write_csv(combined, "data/50_protein_props.csv")
# Challenge 3. Read in full50 (or some version of it) and join with df from Challenge 2.
full50 <- read_csv("data/full50_5.csv")
# full50$label == full50$newnams
# full50$label[30]
# full50$newnams[30]
# combined[50,1]

combined1 <- combined %>% 
  mutate(genome = word(nams, sep = "_", 1))
full50join <- full50 %>% 
  mutate(genome = case_when(!is.na(label) ~ word(label, sep = "_", 1),
                          is.na(label) ~ paste0(word(ole_a_patric, sep = "\\|", 2))))
  

  
fullcomb <- combined1 %>% 
  left_join(full50join, by = "genome")
view(fullcomb)
fullcomb$topt
cbind(fullcomb$topt, fullcomb$genome)
# Challenge 4. Plot relationships between protein predicted optimal temperature and 
# different protein properties like hydrophobicity (Kyte-Doolittle) and instability indices
ggplot(fullcomb) + 
  geom_point(aes(x = topt, y = hydrophob_full)) +
  theme_classic() 
ggplot(fullcomb) + 
  geom_point(aes(x = topt, y = instab)) +
  theme_classic() 
# Challenge 5. Use Ggally function from 20200605 to investigate associations between all protein properties
# and temperature properties (both optimal growth temp of organism from full50 and protein optimum from tomerdf)
ggpairs(fullcomb, columns = c(4:6, 57), cardinality_threshold = 50) +
  theme_classic()
grep("topt", colnames(fullcomb))
colnames(fullcomb)
# Super challenge 6. Make a new page in your web app for interactive plotting of protein properties
# relative to temperature (both optimal growth temp of organism from full50 and protein optimum from tomerdf)


write_csv(fullcomb, "data/fullcomb.csv")
###### work on this -> make inputs for columns
