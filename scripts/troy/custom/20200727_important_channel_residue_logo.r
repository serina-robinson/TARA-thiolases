# Load packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "ggseqlogo", "cowplot")

# Set working directory
setwd("C:/Users/tabie/OneDrive/Documents/GitHub/TARA-thiolases")

# Read channel residue important residue data
dat <- read_csv("output/residue_extraction/important_channel_residues.csv")
colnames(dat)
dat2 <- dat %>% 
  filter(important_residue == TRUE,
         in_mutant_library == TRUE)
dat2$'4ku5_residue'

# Read in 123 OleA sequences
oleas <- readAAStringSet("data/123_OleA.fasta")
ref_4ku5 <- readAAStringSet("data/4KU5.fasta")
names(ref_4ku5) <- '4ku5'

# Get function to get residues from query sequence
source("src/get_residues.r")

# Get residues for 1 sequence
getresidues(oleas[16], ref_4ku5, dat2$'4ku5_residue')


# Get data frames of residues for multiple sequences
result_df <- data.frame(names(oleas))
colnames(result_df) <- "nams"
for(i in 1:length(oleas)) {
  query <- oleas[i]
  result_df$imp_residues[i] <- getresidues(query = oleas[i],
                                       ref = ref_4ku5,
                                       aa_inds = dat2$'4ku5_residue')
}
view(result_df)

# Make logo of all sequences together
pdf("output/residue_extraction/channel_logo_impres_together.pdf")
p <- ggplot() + 
  geom_logo(result_df[,2], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p
dev.off()
