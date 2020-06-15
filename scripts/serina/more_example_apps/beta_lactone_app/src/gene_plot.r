gene_plot<-function(cys){

# Function that reads in a list of DNA sequences
# Returns a list of data frames for plotting using the packages genoPlotR
  
  #Calls the get_dna_segs function
  source("src/get_dna_segs.R")
  
  ll <- list(cys)
  dna <- get_dna_segs(ll)
  #dnal <- list(cys)
  #print(head(dnal[[1]]))
  dnal <- list(dna)
  
  dnal[[1]]$strand <- 1
  colored <- colorRampPalette(brewer.pal(12, "Paired"))(length(cys))
  
  dnal[[1]]$col <- colored
  dnal[[1]]$fill <- colored
  dnal[[1]]$x <- ((dnal[[1]]$start + dnal[[1]]$end)/2) / (dnal[[1]]$end[length(dnal[[1]]$end)])
  dnal[[1]]$y <- 0.4

  proteins <- dnal[[1]]$name # Names of each gene
  print(head(data.frame(dnal[[1]])))
  write.csv(x = data.frame(dnal[[1]]), file = paste0("data/gene_clusters/", dnal[[1]]$name[1],".csv"), quote = F, row.names = F)
  return(list(dnal, proteins))

}
