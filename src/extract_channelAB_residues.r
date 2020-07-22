extract_channelAB <- function(query_fil) {

  # Read in the reference sequence
  query <- readAAStringSet(query_fil)
  ref <- readAAStringSet("data/4KU5.fasta")
  names(ref) <- "4KU5"
  
  # Align the query and the reference
  alned <- AlignSeqs(c(ref, query), verbose = TRUE)
  # alned <- AAStringSet(muscle(c(ref, query), in1 = "data/A_domains_muscle.fasta", in2 = query_fil, profile = T)) # If you want to align with a profile
  query_aln <- alned[length(alned)]
  ref_aln <- alned["4KU5"]
  
  # Read in the 12 angstrom residues
  aachannel_inds <- sort(c(176, 173, 172, 242, 243, 112, 111, 171, 117, 316, 203, 246, 253, 258, 261, 284, 291, 292, 295, 343, 345, 349, 351, 353))
  aachannel_inds_adj <- aachannel_inds # Depends on ref

  # Exract the 34 amino acid positions
  poslist <- list()
  position = 1
 
  for(i in 1:width(ref_aln)) {
    if (substr(ref_aln, i, i) != "-") {
      if (position %in% aachannel_inds_adj) {
        poslist[[i]] <- i
      }
    position = position + 1
    }
  }
  
  # Get the new indices
  new_channel_inds <- unlist(poslist)
  new_channel_inds
  
  # Get 84 aa code
  query_pos <- as.character(unlist(lapply(1:length(new_channel_inds), function(x) {
    substr(query_aln, new_channel_inds[x], new_channel_inds[x]) })))
  return(query_pos)
}



  