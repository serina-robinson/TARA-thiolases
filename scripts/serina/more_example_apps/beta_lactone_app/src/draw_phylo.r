draw_phylo <- function(tab, trefil){
  ### Function that reads in a table and tree file, 
  # Converts to ggtree format and annotates with pubchem ID (in future event of matching to the images of compuonds)
  
  # Pull the genera and clean up the taxonomy names
  gen <- unlist(lapply(1:length(tab$Organism), function(x) { strsplit(tab$Organism, " ")[[x]][1]}))
  spec <- unlist(lapply(1:length(tab$Organism), function(x) { strsplit(tab$Organism, " ")[[x]][2]}))
  spec[is.na(spec)] <- ""
  genspec <- paste0(gen, " ", spec)
  
  # Read in the tree
  p <- ggtree(trefil)

  # Match to the PubChem compound IDs
  cv <- read.csv("data/cid_compound_matchup_key.csv", stringsAsFactors = F)  # Match up compound with its pubchem ID
  genspec <- trimws(genspec)
  genspec <- gsub(" ", "_", genspec)
  genspec <- gsub("\\,", "", genspec)
  dtf <- data.frame(cbind(genspec, cv), stringsAsFactors= F)
  dtf$url <- paste0(dtf$cid, ".png")
  imgurl <- c("data/np_dwnloads/")
  lenn <- length(dtf$url[dtf$url %in% list.files(imgurl)])

  ll <- list()
  for(i in 1:lenn){
     ll[[i]] <- paste0(imgurl,"/",dtf$url[dtf$url %in% list.files(imgurl)][i])
     names(ll[[i]]) <- as.character(dtf$genspec)[dtf$url %in% list.files(imgurl)][i]
   }
  dtf$url[grep("NA.png", dtf$url)] <- "white.png"

  dtfcut <- dtf[dtf$genspec %in% p$data$label[p$data$isTip],]
  dtfcut <- dtfcut[!duplicated(as.character(dtfcut$genspec)),]
  p$data$url[p$data$isTip] <- p$data$label[p$data$isTip]
  p$data$prod[p$data$isTip] <- p$data$label[p$data$isTip]
  
  for(i in 1:length(dtfcut$url)){
    if(length(grep(p$data$url[i], dtfcut$genspec)>0)){
       p$data$url[i] <- as.character(dtfcut$url)[grep(p$data$label[i], dtfcut$genspec)]
       p$data$prod[i] <- as.character(dtfcut$nam)[grep(p$data$label[i], dtfcut$genspec)]
     }
   }

  p$data$label <- gsub("_", " ", p$data$label)

  return(p)
}