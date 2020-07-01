# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "ggseqlogo", "cowplot")

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
olea_chans2
olea_chans2 <- paste0(olea_chans, "HHRKRH")

# Visualize 
pdf("output/channel_logo2.pdf", width = 10, height = 2)
p <- ggplot() + 
  geom_logo(olea_chans2, method = "p", col_scheme = 'chemistry') + 
  theme_logo() 
  #theme(legend.position = 'none', axis.text.x = element_blank()) 
p
dev.off()
?geom_logo
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

# if seqs is a text file of all of the unaligned fasta olea sequences
seqs <- readAAStringSet("data/50_TARA_psychro_thermo_unaligned.fasta")
query <- seqs[1]
length(seqs)
ref <- readAAStringSet("data/4KU5.fasta")
names(ref) <- "4KU5"
channelB <- sort(c(176, 173, 172, 242, 243, 112, 111, 171, 117, 316, 203, 246))
channelA <- sort(c(253, 258, 261, 284, 291, 292, 295, 343, 345, 349, 351, 353))
getresidues(query, ref, c(channelA,channelB))

result34 <- getresidues(seqs[34], ref, c(channelA,channelB))

resultlistA <- c(1:50)
resultlistB <- c(1:50)
for(i in 1:length(seqs)) {
  query <- seqs[i]
  resultlistA[[i]] <- getresidues(query = query,
                                 ref = ref,
                                 aa_inds = channelA)
}
for(i in 1:length(seqs)) {
  query <- seqs[i]
  resultlistB[[i]] <- getresidues(query = query,
                                  ref = ref,
                                  aa_inds = channelB)
}
view(resultlist)
resultlist[34] == result34
resultdat <- data.frame(resultlist)
colnames(resultdat) <- "sequence"
resultdatA <- data.frame(resultlistA)
colnames(resultdatA) <- "sequence"
resultdatB <- data.frame(resultlistB)
colnames(resultdatB) <- "sequence"

p = list()

for(i in 1:length(resultdatA)){
 p[[i]] <- ggplot() + 
  geom_logo(resultdatA[i,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
 return(p)}

pdf("output/channel_logo_togetherB.pdf")
p <- ggplot() + 
  geom_logo(resultdatB[,i], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p
dev.off()
#theme(legend.position = 'none', axis.text.x = element_blank()) 
pdf("output/channel_logo_allB.pdf", height = 100)
plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p23,p24,p25,p26,p27,p28,p29,p30,p31,p32,p33,p34,p35,p36,p37,p38,p39,p40,p41,p42,p43,p44,p45,p46,p47,p48,p49,p50,ncol = 1)
dev.off()
p
p1 <- ggplot() + 
  geom_logo(resultdatB[1,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p2 <- ggplot() + 
  geom_logo(resultdatB[2,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p3 <- ggplot() + 
  geom_logo(resultdatB[3,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p4 <- ggplot() + 
  geom_logo(resultdatB[4,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p5 <- ggplot() + 
  geom_logo(resultdatB[5,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p6 <- ggplot() + 
  geom_logo(resultdatB[6,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p7 <- ggplot() + 
  geom_logo(resultdatB[7,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p8 <- ggplot() + 
  geom_logo(resultdatB[8,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p9 <- ggplot() + 
  geom_logo(resultdatB[9,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p10 <- ggplot() + 
  geom_logo(resultdatB[10,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p11 <- ggplot() + 
  geom_logo(resultdatB[11,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p12 <- ggplot() + 
  geom_logo(resultdatB[12,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p13 <- ggplot() + 
  geom_logo(resultdatB[13,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p14 <- ggplot() + 
  geom_logo(resultdatB[14,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p15 <- ggplot() + 
  geom_logo(resultdatB[15,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p16 <- ggplot() + 
  geom_logo(resultdatB[16,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p17 <- ggplot() + 
  geom_logo(resultdatB[17,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p18 <- ggplot() + 
  geom_logo(resultdatB[18,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p19 <- ggplot() + 
  geom_logo(resultdatB[19,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p20 <- ggplot() + 
  geom_logo(resultdatB[20,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p21 <- ggplot() + 
  geom_logo(resultdatB[21,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p22 <- ggplot() + 
  geom_logo(resultdatB[22,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p23 <- ggplot() + 
  geom_logo(resultdatB[23,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p24 <- ggplot() + 
  geom_logo(resultdatB[24,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p25 <- ggplot() + 
  geom_logo(resultdatB[25,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p26 <- ggplot() + 
  geom_logo(resultdatB[26,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p27 <- ggplot() + 
  geom_logo(resultdatB[27,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p28 <- ggplot() + 
  geom_logo(resultdatB[28,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p29 <- ggplot() + 
  geom_logo(resultdatB[29,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p30 <- ggplot() + 
  geom_logo(resultdatB[30,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p31 <- ggplot() + 
  geom_logo(resultdatB[31,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p32 <- ggplot() + 
  geom_logo(resultdatB[32,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p33 <- ggplot() + 
  geom_logo(resultdatB[33,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p34 <- ggplot() + 
  geom_logo(resultdatB[34,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p35 <- ggplot() + 
  geom_logo(resultdatB[35,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p36 <- ggplot() + 
  geom_logo(resultdatB[36,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p37 <- ggplot() + 
  geom_logo(resultdatB[37,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p38 <- ggplot() + 
  geom_logo(resultdatB[38,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p39 <- ggplot() + 
  geom_logo(resultdatB[39,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p40 <- ggplot() + 
  geom_logo(resultdatB[40,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p41 <- ggplot() + 
  geom_logo(resultdatB[41,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p42 <- ggplot() + 
  geom_logo(resultdatB[42,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p43 <- ggplot() + 
  geom_logo(resultdatB[43,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p44 <- ggplot() + 
  geom_logo(resultdatB[44,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p45 <- ggplot() + 
  geom_logo(resultdatB[45,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p46 <- ggplot() + 
  geom_logo(resultdatB[46,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p47 <- ggplot() + 
  geom_logo(resultdatB[47,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p48 <- ggplot() + 
  geom_logo(resultdatB[48,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p49 <- ggplot() + 
  geom_logo(resultdatB[49,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()
p50 <- ggplot() + 
  geom_logo(resultdatB[50,1], method = "p", col_scheme = 'chemistry') + 
  theme_logo()

plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, ncol = 1)







p
getresidues <- function(query, ref, aa_inds) {
  

  # channel_b <- sort(c(176, 173, 172, 242, 243, 112, 111, 171, 117, 316, 203, 246))
  # channel_a <- sort(c(253, 258, 261, 284, 291, 292, 295, 343, 345, 349, 351, 353))
  # aa_inds <- c(channel_a, channel_b)

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
olea_chans
olea_chans2
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
###### maybe replace with for loop

## Mutate the serine back to a C in 4ku5.fasta (write to new file)

## read papers about positive patch

## make aa_inds variable
  