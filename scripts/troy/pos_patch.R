# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "ggseqlogo", "cowplot")

# Set working directory
setwd("C:/Users/tabie/OneDrive/Documents/GitHub/TARA-thiolases/")

# Read in the query and reference sequences
seqs50 <- readAAStringSet("data/50_TARA_psychro_thermo_unaligned.fasta")
ref_4ku5 <- readAAStringSet("data/1EBL_A.fasta")
names(ref_4ku5) <- '4ku5'
ref_1ebl <- readAAStringSet("data/4KU5.fasta")
names(ref_1ebl) <- '1ebl'

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

# Get residues for 1 sequence
getresidues(ref_4ku5, ref_1ebl, pospatch_1ebl)

# Get list of residues for multiple sequences
resultlist <- list()
for(i in 1:length(seqs50)) {
  query <- seqs50[i]
  resultlist[[i]] <- getresidues(query = query,
                                  ref = ref_1ebl,
                                  aa_inds = pospatch_1ebl)
}
view(resultlist)
