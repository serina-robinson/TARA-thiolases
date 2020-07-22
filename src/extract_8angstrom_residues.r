extract_8angstrom <- function(query_fil) {

  # Read in the reference sequence
  query <- readAAStringSet(query_fil)
  ref <- readAAStringSet("data/4KU5.fasta")
  names(ref) <- "4KU5"
  
  # Align the query and the reference
  alned <- AlignSeqs(c(ref, query), verbose = TRUE)
  # alned <- AAStringSet(muscle(c(ref, query), in1 = "data/A_domains_muscle.fasta", in2 = query_fil, profile = T)) # If you want to align with a profile
  query_aln <- alned[length(alned)]
  ref_aln <- alned["4KU5"]
  
  # Read in the 8 angstrom residues
  aa33_inds <- read_csv("data/residue_extraction/8_angstrom_radius.csv") %>%
    pull()
  aa33_inds_adj <- aa33_inds # Depends on ref

  # Exract the 33 amino acid positions
  poslist <- list()
  position = 1
 
  for(i in 1:width(ref_aln)) {
    if (substr(ref_aln, i, i) != "-") {
      if (position %in% aa33_inds_adj) {
        poslist[[i]] <- i
      }
    position = position + 1
    }
  }
  
  # Get the new indices
  new_33inds <- unlist(poslist)
  new_33inds
  
  # Get 33 aa code
  query_pos <- as.character(unlist(lapply(1:length(new_33inds), function(x) {
    substr(query_aln, new_33inds[x], new_33inds[x]) })))
  return(query_pos)
}



  