# Install packages
pacman::p_load("tidyverse", "Biostrings", "seqinr", "DECIPHER", "Biostrings", "ggseqlogo", "cowplot", "gggenes", "ggtree")

# Set working directory
setwd("C:/Users/tabie/OneDrive/Documents/GitHub/TARA-thiolases/")

# Read in the sequences
seqs <- read.fasta("data/cah_project/11_CAH_sequences_no_his_tag.fasta", seqtype = "AA", strip.desc = TRUE)
seqs
# Create a data frame
genes <- attr(seqs, "name") #generates a gene list
seq.df <- data.frame(nams = attr(seqs, "name")) %>%
  dplyr::mutate(seqs = seqs)
##Calculating AA statistics for each sequence using AAstat() in the "seqinr" package

aastat.df <- data.frame(sapply(seq.df$seqs, AAstat)) # check out the plot
# aastat.df[[1]]
# class(aastat.df)
# view(aastat.df)
# ?AAstat
# AAstat(seq.df$seqs[[11]])

new.df <- data.frame(nams = attr(seqs, "name"))
for(i in 1:length(genes)){
  new.df$len[i] <- length(seq.df$seqs[[i]])-1 #aa sequence length
  new.df$Ala[i] <- aastat.df[[i]][1][[1]][which(names(aastat.df[[i]][1][[1]])=="A")] / new.df$len[i] 
  new.df$Cys[i] <- aastat.df[[i]][1][[1]][which(names(aastat.df[[i]][1][[1]])=="C")] / new.df$len[i]
  new.df$Asp[i] <- aastat.df[[i]][1][[1]][which(names(aastat.df[[i]][1][[1]])=="D")] / new.df$len[i]
  new.df$Glu[i] <- aastat.df[[i]][1][[1]][which(names(aastat.df[[i]][1][[1]])=="E")] / new.df$len[i]
  new.df$Phe[i] <- aastat.df[[i]][1][[1]][which(names(aastat.df[[i]][1][[1]])=="F")] / new.df$len[i]
  new.df$Gly[i] <- aastat.df[[i]][1][[1]][which(names(aastat.df[[i]][1][[1]])=="G")] / new.df$len[i]
  new.df$His[i] <- aastat.df[[i]][1][[1]][which(names(aastat.df[[i]][1][[1]])=="H")] / new.df$len[i]
  new.df$Ile[i] <- aastat.df[[i]][1][[1]][which(names(aastat.df[[i]][1][[1]])=="I")] / new.df$len[i]
  new.df$Lys[i] <- aastat.df[[i]][1][[1]][which(names(aastat.df[[i]][1][[1]])=="K")] / new.df$len[i]
  new.df$Leu[i] <- aastat.df[[i]][1][[1]][which(names(aastat.df[[i]][1][[1]])=="L")] / new.df$len[i] 
  new.df$Met[i] <- aastat.df[[i]][1][[1]][which(names(aastat.df[[i]][1][[1]])=="M")] / new.df$len[i]
  new.df$Asn[i] <- aastat.df[[i]][1][[1]][which(names(aastat.df[[i]][1][[1]])=="N")] / new.df$len[i]
  new.df$Pro[i] <- aastat.df[[i]][1][[1]][which(names(aastat.df[[i]][1][[1]])=="P")] / new.df$len[i]
  new.df$Gln[i] <- aastat.df[[i]][1][[1]][which(names(aastat.df[[i]][1][[1]])=="Q")] / new.df$len[i]
  new.df$Arg[i] <- aastat.df[[i]][1][[1]][which(names(aastat.df[[i]][1][[1]])=="R")] / new.df$len[i]
  new.df$Ser[i] <- aastat.df[[i]][1][[1]][which(names(aastat.df[[i]][1][[1]])=="S")] / new.df$len[i]
  new.df$Thr[i] <- aastat.df[[i]][1][[1]][which(names(aastat.df[[i]][1][[1]])=="T")] / new.df$len[i]
  new.df$Val[i] <- aastat.df[[i]][1][[1]][which(names(aastat.df[[i]][1][[1]])=="V")] / new.df$len[i]
  new.df$Trp[i] <- aastat.df[[i]][1][[1]][which(names(aastat.df[[i]][1][[1]])=="W")] / new.df$len[i] 
  new.df$Tyr[i] <- aastat.df[[i]][1][[1]][which(names(aastat.df[[i]][1][[1]])=="Y")] / new.df$len[i] 
  new.df$RKrat[i] <- (new.df$Arg[i]*new.df$len[i]) /(new.df$Arg[i]*new.df$len[i] + new.df$Lys[i]*new.df$len[i]) 
  new.df$NQrat[i] <- (new.df$Asn[i]*new.df$len[i] + new.df$Gln[i]*new.df$len[i]) / (new.df$Asn[i]*new.df$len[i] + new.df$Glu[i]*new.df$len[i] + new.df$Gln[i]*new.df$len[i] + new.df$Arg[i]*new.df$len[i])
  new.df$disordered_res[i] <- (new.df$Gln[i]*new.df$len[i] + new.df$Ser[i]*new.df$len[i] + new.df$Pro[i]*new.df$len[i] + new.df$Glu[i]*new.df$len[i] + new.df$Lys[i]*new.df$len[i]) / new.df$len[i]
  new.df$ordered_res[i] <- (new.df$Val[i]*new.df$len[i] + new.df$Leu[i]*new.df$len[i] + new.df$Ile[i]*new.df$len[i] + new.df$Met[i]*new.df$len[i] + new.df$Phe[i]*new.df$len[i] + new.df$Trp[i]*new.df$len[i] + new.df$Tyr[i]*new.df$len[i]) / new.df$len[i]
  new.df$Tiny[i] <- aastat.df[[i]][2][[1]][which(names(aastat.df[[i]][2][[1]])=="Tiny")]
  new.df$Small[i] <- aastat.df[[i]][2][[1]][which(names(aastat.df[[i]][2][[1]])=="Small")]
  new.df$Aliphatic[i] <- aastat.df[[i]][2][[1]][which(names(aastat.df[[i]][2][[1]])=="Aliphatic")]
  new.df$Aromatic[i] <- aastat.df[[i]][2][[1]][which(names(aastat.df[[i]][2][[1]])=="Aromatic")]
  new.df$Nonpolar[i] <- aastat.df[[i]][2][[1]][which(names(aastat.df[[i]][2][[1]])=="Non.polar")]
  new.df$Polar[i] <- aastat.df[[i]][2][[1]][which(names(aastat.df[[i]][2][[1]])=="Polar")]
  new.df$Charged[i] <- aastat.df[[i]][2][[1]][which(names(aastat.df[[i]][2][[1]])=="Charged")]
  new.df$Basic[i] <- aastat.df[[i]][2][[1]][which(names(aastat.df[[i]][2][[1]])=="Basic")]
  new.df$Acidic[i] <- aastat.df[[i]][2][[1]][which(names(aastat.df[[i]][2][[1]])=="Acidic")]
}





## Align sequences
ref <- readAAStringSet("data/cah_project/pseudolabrys.fasta")
query1 <- readAAStringSet("data/cah_project/Azorhizobium caulinodans ORS 571 (4NQ3_A) no X.fasta")
query2 <- readAAStringSet("data/cah_project/Enterobacter cloacae (5T13_A) no H tag.fasta")
query3 <- readAAStringSet("data/cah_project/Frankia sp. Eul1b (5HY0_A) no H tag.fasta")
query4 <- readAAStringSet("data/cah_project/Moorella thermoacetica ATCC 39073 (6BUM_A).fasta")
query5 <- readAAStringSet("data/cah_project/Moorella thermoacetica ATCC 39073 (6DHJ_A) no H tag.fasta")
query6 <- readAAStringSet("data/cah_project/Pseudomonas sp. ADP (4BVQ_A) no H tag.fasta")
query7 <- readAAStringSet("data/cah_project/Rhodococcus erythropolis (5HWE_A) no H tag.fasta")
alned <- DECIPHER::AlignSeqs(c(ref, query1, query2, query3, query4, query5, query6, query7
                               ), verbose = FALSE)
BrowseSeqs(alned, "output/alnedseq.html")


## Make logo

seqs1 <- readAAStringSet("data/cah_project/10_CAH_sequences_no_his_tag.fasta")

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
resultlist <- list()
for(i in 1:length(seqs1)) {
  query <- seqs1[i]
  resultlist[[i]] <- getresidues(query = query,
                                  ref = ref,
                                  aa_inds = c())
}



## to do
# remove his tags #
# long format #
# Pseudolabrys on top #
# make heatmap #
# fix his tags and use new file #
# blast frankia one #
# xtalpred without his tags #
# clean up table from xtalpred #
##### later: Ca binding loop below (moorella as ref)


# Heat map and such
for(i in 2:35){
  new.df[i] <- as.numeric(unlist(new.df[i]))
}
class(new.df$Polar)
new.df$nams2 <- c("Pseudolabrys", "Hydrogenophaga", "Bradyrhizobium", "Variovorax", "Azorhizobium", "Enterobacter", "Frankia", "Moorella 1", "Moorella 2", "Pseudomonas", "Rhodococcus")
new.df$nams2
colnames(new.df)
# wide to long
new.df_long <- new.df %>% 
  gather(key = stat, value = value, -nams, -len, -nams2)
new.df_long_large <- new.df %>% 
  gather(key = stat, value = value, RKrat, Tiny, Small, Nonpolar, Polar)
new.df_long_med <- new.df %>% 
  gather(key = stat, value = value, NQrat, ordered_res, disordered_res, Aliphatic, Charged)
new.df_long_small <- new.df %>% 
  gather(key = stat, value = value, -nams, -len, -nams2, -RKrat, -NQrat, -ordered_res, -disordered_res, -Tiny, -Small, -Aliphatic, -Nonpolar, -Polar, -Charged)

# Large values
pdf("output/CAH_heatmap_large.pdf", height = 4, width = 8)
ggplot(new.df_long_large, aes(stat, nams2)) +
  geom_tile(aes(fill = value)) + 
 # geom_text(aes(label = round(value, 4))) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 50, size = 8,
                                   vjust = 1.05, hjust = 1.2)) +
  scale_fill_gradient(low = "white", high = "darkgreen")

dev.off()


# Small values
pdf("output/CAH_heatmap_small.pdf", height = 4, width = 8)
ggplot(new.df_long_small, aes(stat, nams2)) +
  geom_tile(aes(fill = value)) + 
  # geom_text(aes(label = round(value, 4))) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 50, size = 8,
                                   vjust = 1.05, hjust = 1.2)) +
  scale_fill_gradient(low = "white", high = "darkgreen")

dev.off()

# Med values
pdf("output/CAH_heatmap_med.pdf", height = 4, width = 8)
ggplot(new.df_long_med, aes(stat, nams2)) +
  geom_tile(aes(fill = value)) + 
  # geom_text(aes(label = round(value, 4))) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 50, size = 8,
                                   vjust = 1.05, hjust = 1.2)) +
  scale_fill_gradient(low = "white", high = "darkgreen")

dev.off()


# All values
pdf("output/CAH_heatmap.pdf", height = 4, width = 8)
ggplot(new.df_long, aes(stat, nams2)) +
  geom_tile(aes(fill = value)) + 
  # geom_text(aes(label = round(value, 4))) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 50, size = 8,
                                   vjust = 1.05, hjust = 1.2)) +
  scale_fill_gradient(low = "white", high = "darkgreen")

dev.off()



# Tree stuff

ml <- read.tree("data/CAH_project/10_CAH_tree_100_bootstrap.nwk") # maximum-likelihood method
mlt <- ggtree(ml)
view(mlt$data)
mlt$data <- mlt$data %>% 
  mutate(label2 = c("Pseudolabrys sp. Root1462 (WP_056911810.1)", "Azorhizobium caulinodans ORS 571 (4NQ3_A)", "Hydrogenophaga (WP_009515690.1)", "Variovorax paradoxus (WP_018906567.1)", "Enterobacter cloacae (5T13_A)", "Moorella thermoacetica ATCC 39073 (6DHJ_A)", "Bradyrhizobium sp. WSM1253 (WP_007596559.1)", "Pseudomonas sp. ADP (4BVQ_A)", "Frankia sp. Eul1b (5HY0_A)", "Rhodococcus erythropolis (5HWE_A)", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_),
         color = case_when(label2 == "Pseudolabrys sp. Root1462 (WP_056911810.1)" ~ "blue",
                           TRUE ~ "black"),
         label = round(as.numeric(label), 2))
mlt_boot <- mlt + 
  geom_tiplab(aes(label = label2), color = mlt$data$color[1:10]) +
  xlim(0, 2) +
  geom_nodelab(aes(label = label), hjust = 1.2, vjust = -0.25)
pdf("data/CAH_project/CAH_tree.pdf")
mlt_boot
dev.off()

