# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "ggseqlogo", "cowplot")

# Set working directory
setwd("C:/Users/tabie/OneDrive/Documents/GitHub/TARA-thiolases/")

# Read in the query and reference sequences
seqs50 <- readAAStringSet("data/50_TARA_psychro_thermo_unaligned.fasta")
seqs123 <- readAAStringSet("data/123_OleA.fasta")
ref_1ebl <- readAAStringSet("data/1EBL_A.fasta")
ref_4ku5 <- readAAStringSet("data/4KU5.fasta")
names(ref_1ebl) <- '1ebl'
names(ref_4ku5) <- '4ku5'
# BrowseSeqs(AlignSeqs(c(ref_1ebl, ref_4ku5)), "output/alnedseq2.html")
# Make function to get residues from query sequence
getresidues <- function(query, ref, aa_inds) {
  
  alned <- DECIPHER::AlignSeqs(c(ref, query), verbose = FALSE)
  ref_aln <- alned[1]
  query_aln <- alned[length(alned)]
  
  poslist <- list()
  position = 1
  
  for(j in 1:width(ref_aln)) {
    if (substr(ref_aln, j, j) != "-") {
      if (position %in% aa_inds) {
        poslist[[j]] <- j}
      position = position + 1}}
  
  new_inds <- unlist(poslist)
  
  query_pos <- as.character(unlist(lapply(1:length(new_inds), function(x) {
    substr(query_aln, new_inds[x], new_inds[x]) })))
  
  olea_chans <- paste0(query_pos, collapse = "")
  return(olea_chans)
}

# Define parameters 
pospatch_1ebl <- c(214, 249, 253, 256, 257)
channelA <- sort(c(253, 258, 261, 284, 291, 292, 295, 343, 345, 349, 351, 353))
channelB <- sort(c(176, 173, 172, 242, 243, 112, 111, 171, 117, 316, 203, 246))
rad8 <- read_csv("data/8_angstrom_radius.csv")
rad8vec <-  rad8$ang8.resi
rad10 <- read_csv("data/10_angstrom_radius.csv")
rad10vec <-  rad10$ang10.resi
rad12 <- read_csv("data/12_angstrom_radius.csv")
rad12vec <-  rad12$ang12.resi
rad14 <- read_csv("data/14_angstrom_radius.csv")
rad14vec <-  rad14$ang14.resi

seqs <- seqs123
# seqs <- seqs50
# Get residues for 1 sequence
getresidues(seqs[35], ref_4ku5, channelA)

# Get list of residues for multiple sequences
result_df <- data.frame(names(seqs))
colnames(result_df) <- "nams"
for(i in 1:length(seqs)) {
  query <- seqs[i]
  result_df$channelA[i] <- getresidues(query = query,
                                  ref = ref_4ku5,
                                  aa_inds = channelA)
}
view(result_df)

# write_csv(result_df, "data/123_OleA_important_residues.csv")
# BrowseSeqs(AlignSeqs(c(ref_1ebl, seqs[1])), "output/alnedseq3.html")

# make a function to separate the sequences into their own column

AAseq <- result_df[1,2]
nchar(AAseq)
sepresidues <- function(AAseq, ){

}
result_df2 <- data.frame(names(seqs))
for(i in 1:nchar(AAseq)) {
  residue[i] <- substr(AAseq, i, i)
}
cbind(result_df2, residue)
