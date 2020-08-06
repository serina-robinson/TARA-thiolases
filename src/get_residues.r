# Gets residues on query that match aa_inds from a reference

getresidues <- function(query, ref, aa_inds) {
# Align sequences
  alned <- DECIPHER::AlignSeqs(c(ref, query), verbose = FALSE)
  ref_aln <- alned[1]
  query_aln <- alned[length(alned)]
  
  poslist <- list()
  position = 1
# Count through alignment to account for gaps  
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