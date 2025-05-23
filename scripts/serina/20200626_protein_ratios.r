# Install packages
pacman::p_load("tidyverse", "Biostrings", "seqinr")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/TARA-thiolases/")

# Read in the sequences
seqs <- read.fasta("data/50_TARA_psychro_thermo_unaligned.fasta", seqtype = "AA", strip.desc = TRUE)

# Create a data frame
genes <- attr(seqs, "name") #generates a gene list
seq.df <- data.frame(nams = attr(seqs, "name")) %>%
  dplyr::mutate(seqs = seqs)

##Calculating AA statistics for each sequence using AAstat() in the "seqinr" package
aastat.df <- data.frame(sapply(seq.df$seqs, AAstat)) # check out the plot
aastat.df[[1]]

# Wow this is ugly...I'm sure there is a better way to do this haha
new.df <- data.frame(nams = attr(seqs, "name"))
for(i in 1:length(genes)){
  new.df$len[i] <- length(seq.df$seqs[[i]])-1 #aa sequence length
  new.df$Arg[i] <- aastat.df[[i]][1][[1]][which(names(aastat.df[[i]][1][[1]])=="R")] / new.df$len[i] 
  new.df$Lys[i] <- aastat.df[[i]][1][[1]][which(names(aastat.df[[i]][1][[1]])=="K")] / new.df$len[i]
  new.df$Asn[i] <- aastat.df[[i]][1][[1]][which(names(aastat.df[[i]][1][[1]])=="N")] / new.df$len[i]
  new.df$Gln[i] <- aastat.df[[i]][1][[1]][which(names(aastat.df[[i]][1][[1]])=="Q")] / new.df$len[i]
  new.df$Glu[i] <- aastat.df[[i]][1][[1]][which(names(aastat.df[[i]][1][[1]])=="E")] / new.df$len[i]
  new.df$Asp[i] <- aastat.df[[i]][1][[1]][which(names(aastat.df[[i]][1][[1]])=="N")] / new.df$len[i]
  new.df$Pro[i] <- aastat.df[[i]][1][[1]][which(names(aastat.df[[i]][1][[1]])=="D")] / new.df$len[i]
  new.df$Thr[i] <- aastat.df[[i]][1][[1]][which(names(aastat.df[[i]][1][[1]])=="T")] / new.df$len[i]
  new.df$Leu[i] <- aastat.df[[i]][1][[1]][which(names(aastat.df[[i]][1][[1]])=="L")] / new.df$len[i]
  new.df$Aliphatic[i] <- aastat.df[[i]][2][[1]][which(names(aastat.df[[i]][2][[1]])=="Aliphatic")]
  new.df$Acidic[i] <- aastat.df[[i]][2][[1]][which(names(aastat.df[[i]][2][[1]])=="Acidic")]
  new.df$Polar[i] <- aastat.df[[i]][2][[1]][which(names(aastat.df[[i]][2][[1]])=="Polar")]
}

## These are challenges for the WEEK not just for Friday!

# Challenge 0. Create a new column in cold.df calculating the arginine-to-lysine ratio 
# This ratio was shown to be decreased in some psychrophilic enzymes
# Why is that? To find out:
# See this paper (also in your lit folder as Siddiqui_2006) 
# https://onlinelibrary-wiley-com.ezp1.lib.umn.edu/doi/epdf/10.1002/prot.20989
# This ratio is calculated using the following formula: Arg /(Lys + Arg) e.g. R/(R + K)
# You can also find a description of this, and other properties in this section of a Google book:
# https://books.google.com/books?id=LANEAgAAQBAJ&pg=PA82&lpg=PA82&dq=arginine-to-lysine+ratio+cold-adapted+enzymes+temperature&source=bl&ots=rC2L68FUye&sig=ACfU3U10pk5Pdz-B7B43ZfUNGFcW8Vf_ww&hl=en&sa=X&ved=2ahUKEwic8I_RwJ7qAhUMbs0KHTpeARcQ6AEwAnoECAgQAQ#v=onepage&q=arginine-to-lysine%20ratio%20cold-adapted%20enzymes%20temperature&f=false

# Challenge 1. Create a new column in cold.df calculating the ratio of asparagine and glutamine against all amino acids
# and specifically against Asn, Glu, Gln, and Asp. You can use the formula: Asn + Glu / (Asn + Glu + Gln + Asp)
# Why are we testing this ratio? Do you expect it to be higher or lower in thermophilic enzymes?
# See Bhanuramanand et al. 2014 for some clues

# Challenge 2. Add new column in cold.df for properties that were excluded including
# Basic, Non.polar, Small etc....incidentally, what is the difference between small and tiny?
# Which residues are included in the small and tiny groups?
# To figure that out you will also need to add columns to cold.df for the other amino acids not currently 
# included e.g. cysteine...check that in the end you include columns for all 20 amino acids

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

## Challenge 14. Have a fun week, don't work to hard, and take at least Friday off!
