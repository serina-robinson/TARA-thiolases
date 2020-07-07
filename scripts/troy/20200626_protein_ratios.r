# Install packages
pacman::p_load("tidyverse", "Biostrings", "seqinr")

# Set working directory
setwd("C:/Users/tabie/OneDrive/Documents/GitHub/TARA-thiolases/")

# Read in the sequences
seqs <- read.fasta("data/50_TARA_psychro_thermo_unaligned.fasta", seqtype = "AA", strip.desc = TRUE)
# seqs
# Create a data frame
genes <- attr(seqs, "name") #generates a gene list
seq.df <- data.frame(nams = attr(seqs, "name")) %>%
  dplyr::mutate(seqs = seqs)
##Calculating AA statistics for each sequence using AAstat() in the "seqinr" package
aastat.df <- data.frame(sapply(seq.df$seqs, AAstat)) # check out the plot
aastat.df[[1]]
class(aastat.df)
# view(aastat.df)
# ?AAstat
# seq.df$seqs[[2]]
# class(seq.df$seqs)
# Wow this is ugly...I'm sure there is a better way to do this haha
new.df_og <- data.frame(nams = attr(seqs, "name"))
# rm(new.df)
new.df <- new.df_og
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
# aastat.df[[i]][2][[1]]
aastat.df2 <- data.frame(t(aastat.df))
aastat.df3 <- aastat.df2 %>% 
  mutate(nams = rownames(aastat.df2))
aastat.df4 <- aastat.df3 %>% 
  mutate(nams = case_when(substr(nams, 1, 1) == "X" ~ substr(nams, 2, 100),
                          TRUE ~ nams))
all.df <- full_join(aastat.df4, seq.df, by = "nams")
# which(names(aastat.df[[2]][1][[1]])=="R")
# view(all.df)
## to do list:
# did the join work?
# keep doing this mutate thing
# read some of the articles
# try the next challenges
# ?t

# allstat.df <- all.df %>% 
  # mutate(aa_length = lengths(seqs)
         # )

## These are challenges for the WEEK not just for Friday!
# cold.df = new.df
# Challenge 0. Create a new column in cold.df calculating the arginine-to-lysine ratio 
# This ratio was shown to be decreased in some psychrophilic enzymes
# Why is that? To find out:
# See this paper (also in your lit folder as Siddiqui_2006) 
# https://onlinelibrary-wiley-com.ezp1.lib.umn.edu/doi/epdf/10.1002/prot.20989
# This ratio is calculated using the following formula: Arg /(Lys + Arg) e.g. R/(R + K)
# You can also find a description of this, and other properties in this section of a Google book:
# https://books.google.com/books?id=LANEAgAAQBAJ&pg=PA82&lpg=PA82&dq=arginine-to-lysine+ratio+cold-adapted+enzymes+temperature&source=bl&ots=rC2L68FUye&sig=ACfU3U10pk5Pdz-B7B43ZfUNGFcW8Vf_ww&hl=en&sa=X&ved=2ahUKEwic8I_RwJ7qAhUMbs0KHTpeARcQ6AEwAnoECAgQAQ#v=onepage&q=arginine-to-lysine%20ratio%20cold-adapted%20enzymes%20temperature&f=false

# Challenge 1. Create a new column in cold.df calculating the ratio of asparagine and glutamine against all amino acids
# and specifically against Asn, Glu, Gln, and Asp. You can use the formula: Asn + Gln / (Asn + Glu + Gln + Asp)
# Why are we testing this ratio? Do you expect it to be higher or lower in thermophilic enzymes?
# See Bhanuramanand et al. 2014 for some clues

# Challenge 2. Add new column in cold.df for properties that were excluded including
# Basic, Non.polar, Small etc....incidentally, what is the difference between small and tiny?
# Which residues are included in the small and tiny groups?
# To figure that out you will also need to add columns to cold.df for the other amino acids not currently 
# included e.g. cysteine...check that in the end you include columns for all 20 amino acids
# double check that adding up ratios = 1

# Challenge 3 (long-term). From perusing the papers in the lit folder and reading other papers they cite,
# calculate other amino acid ratios that or properties may be interesting or involved in temperature adaptation.
# For example, order-promoting residues (V, L, I, M, F, W, and Y)
# and disorder-promoting residues (Q, S, P, E, K)

# Challenge 4. Use ggplot boxplots to visualize the ratios of amino acids in your proteins
# See Figure 1 in Yang_2015_GenBiolEvol for an example
# Also here: https://www.r-graph-gallery.com/265-grouped-boxplot-with-ggplot2.html
# Then use the geom_boxplot(fill = ) to color the boxplots by their temperature status
# Try making a different type of plot to visualize this. For example, trying using geom_bar(position = "dodge")
# What else can you come up with? 
full50 <- read.csv("data/full50_6.csv")

new.df[grep("Psychrobacter", new.df$nams), "nams"] <- "282669.3.peg.967_Kocuria_mesophile"
new.df[grep("2026735.166.peg.848_Deltaproteobacteria", new.df$nams), "nams"] <- "2026735.166.peg.848_Deltaproteobacteria_psychrophile"
new.df[grep("2026735.148.peg.1837_Deltaproteobacteria", new.df$nams), "nams"] <- "2026735.148.peg.1837_Deltaproteobacteria_psychrophile"
new.df[grep("2026809.16.peg.2167_Epsilonproteobacteri", new.df$nams), "nams"] <- "2026809.16.peg.2167_Epsilonproteobacteria_psychrophile"
new.df[grep("2026799.183.peg.457_Verrucomicrobia", new.df$nams), "nams"] <- "2026799.183.peg.457_Verrucomicrobia_psychrophile"
new.df[grep("2483033.3.peg.4057_Sedimentitalea", new.df$nams), "nams"] <- "2483033.3.peg.4057_Sedimentitalea_mesophile"
dat <- full_join(full50, new.df, by = "nams")
# ggplot(dat, aes(x=Arg, y=temperature_range)) + 
#   geom_boxplot()
# colnames(dat)
# table(dat$temperature_range)
 ## find topt and add to dat and then write to file 
## Challenge 5. Now make scatter plots comparing the amino acid ratios to optimal/isolation temperature
# and protein topt. You can either use Ggally or other make them individually

## Challenge 6. You guessed it...add a page to your growing web app that allows you to 
# interactively explore these amino acid ratios relative to temperature. 
# You can also combine this with the overall protein properties like hydrophobicity etc.
# that you calculated on Tuesday.

## Challenge 7. Now that you've done all this work and looked at all these plots..can you find any
# correlations in your data? You can use ggally to start, but then make individual plots for
# The most interesting correlations to show in lab meeting. You can use geom_cor to calculate
# correlations and p-values. https://rdrr.io/bioc/DEGreport/man/geom_cor.html
# What is the strongest correlation with isolation temperature? With topt?
# Are any of the predictors themselves correlated?
# (e.g. hydrophobicity and order-promoting residues for example)

## Challenge 8. If you're up for it, you could present one or more of these plots above at lab meeting on Tuesday. 

## Challenge 9. Based on the results from challenge 7, design some experiments using site-directed mutagenesis
# in which you would mutate residues in your OleAs to make them more thermostable/thermophilic, for example
# Write out which residues you would mutate in particular...perhaps surface residues so as not to disrupt
# the active site? When you disrupt the active site it often kills proteins. Use Pymol to visualize
# for example, all arginines, or all cysteines, and choose candidates not near the active site.
# You can look up some papers that describe protein engineering or rational mutagenesis to improve thermostability.
# The review in lit by Modares_2016 is a start.
# Propose few different experiments (actually write it in a Word doc). 
# We might even be able to try one of them this summer or fall!

## Challenge 10. We've been focusing on temperature up to this point. But what about looking at PUFA vs.
# non-PUFA sequences with respect to these ratios and properties that you have calculated?
# Use full50 has_pufa column to find which sequences have pfa synthases associated.

## Challenge 11. Now that you've looked at these properties in the full length protein sequence
# you could also look at them in just different regions of the protein such as the positive patch 
# or the channel A and B residues. Could you make this interactive in flexdashboard e.g. with a drop-down box for
# if you wanted to specify certain regions of the protein or the full length protein?

## Challenge 12. (hard) See if there is any way to clean up the monstrosity of the for loop code above. Can
# you use the %>% and mutate instead somehow? # Maybe something in the 'apply' family? 
# I don't know the answer to this one, but give it a shot :)

## Challenge 13. If you get super bored, it's time to bring in the 73 extra OleAs we already have.
# 73_OleA_seqs_for_lookup.csv are ALWAYS fair game, and we already have them in the freezer
# so they're ready to go!
# Could you look up if any of them are thermophiles or psychrophiles? 
# Can you find isolation temperatures and locations?
# Do the properties of their sequences support trends you're observing with the 50 you've been working with?
# stuff <- read_csv("data/73_OleA_seqs_for_lookup.csv")
## Challenge 14. Have a fun week, don't work to hard, and take at least Friday off!

## adding topt to dat
props <- read_csv("data/50_protein_props.csv")
full50_allstats <- full_join(props, dat, by = "nams")
view(dat)
view(full50_allstats)
### to do
# make plots
# add temp stat for 73 sequences 
# merge datasets
# write new dat file
# write down expected 4ku5 sequences
# 