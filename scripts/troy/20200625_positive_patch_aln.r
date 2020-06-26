# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "ggseqlogo")

# Set working directory
setwd("C:/Users/tabie/OneDrive/Documents/GitHub/TARA-thiolases/")

# Read in the query and reference sequences
query_50 <- readAAStringSet("data/50_TARA_psychro_thermo_unaligned.fasta")
query <- query_50[2]
names(query) # 1075768.5.peg.2737_Janthinobacterium_psychrophile
ref <- readAAStringSet("data/4KU5.fasta")
names(ref) <- "4KU5"
  
# Align the query and the reference
alned <- DECIPHER::AlignSeqs(c(ref, query), verbose = TRUE)
BrowseSeqs(alned, "output/alnedseq.html")
ref_aln <- alned["4KU5"]
query_aln <- alned[length(alned)]
length(alned)
width(alned)
# Read in the indices (from 4KU5 from Pymol)
# pos_patch <- c() # this could be where you put your indices you discover for the positive patch!
channel_b <- sort(c(176, 173, 172, 242, 243, 112, 111, 171, 117, 316, 203, 246))
channel_b
channel_a <- sort(c(253, 258, 261, 284, 291, 292, 295, 343, 345, 349, 351, 353))
channel_a

aa_inds <- c(channel_a, channel_b)
aa_inds

# Exract the 34 amino acid positions
poslist <- list()
position = 1

for(i in 1:width(ref_aln)) {
  if (substr(ref_aln, i, i) != "-") {
    if (position %in% aa_inds) {
      poslist[[i]] <- i
    }
  position = position + 1
  }
}

# Get the new indices
new_inds <- unlist(poslist)
new_inds
  
# Get complete code
query_pos <- as.character(unlist(lapply(1:length(new_inds), function(x) {
  substr(query_aln, new_inds[x], new_inds[x]) })))
olea_chans <- paste0(query_pos, collapse = "")
olea_chans

# Visualize 
pdf("output/channel_logo.pdf", width = 10, height = 2)
p <- ggplot() + 
  geom_logo(olea_chans, method = "p", col_scheme = 'hydrophobicity') + 
  theme_logo() 
  #theme(legend.position = 'none', axis.text.x = element_blank()) 
p
dev.off()

## Challenge 1. Work through the ggseqlogo tutorial:
# https://omarwagih.github.io/ggseqlogo/#colour_schemes
# Play around with changing the color. 
# What does the chemistry coloring look like if you add a postively charged residue?
# What does hydrophobicity coloring look like?

## Challenge 2. Look at the DECIPHER documentation:
# http://www2.decipher.codes/Documentation/ArtOfAlignmentInR.pdf
# How does DECIPHER compare to other tools like Muscle?
# Pay particular attention to the flowchart in Fig. 6.
# What function would you run after AlignSeqs if you were doing this for phylogenetics?
# What would you do if you had nucleotide sequences?

## Challenge 3. (hard) Write a function to extract the channel residues from all 50 TARA sequences
# This writing functions tutorial might be helpful: https://swcarpentry.github.io/r-novice-inflammation/02-func-R/

## Challenge 4. (hard) Modify your function so that it has an input allowing you to change which residues you search for.
# For example, allow it to only search for channel_a residues? Or channel_b residues? (and next week, positive patch residues)

## Challenge 5. (hard) Make ggseqlogos for each of your TARA sequences. 
# Use cowplot to write them to the same pdf. 
# If possible, make one PDF for Channel A residues in TARA sequences
# And one PDF for Channel B residues.
# You might want to use a for loop here.
# Are they the same, or different? Where are they different?

## Challenge 6. (even harder) Try to follow the tutorial to first align all 50 sequences 
# and then visualize in one ggseqlogo what the alignment looks like. What residues are conserved?
# which positions are variable?

## Super challenge 7. Make a sequence logo, sequence alignment and bar plot like the last part of the tutorial.
# What other types of plots could you make? Could you some include temperature information?
# What about PUFA presence/absence? Do the PUFA-containing organisms have different channel residues?

## For next week: Identify 'positive patch' residues and run this same sort of code to visualize and compare
# positive patches among your 50 sequences.

## read about apply family of functions (lapply in line 49ish)
# maybe replace with for loop

## Mutate the serine back to a C in 4ku5.fasta (write to new file)

## read papers about positive patch
  