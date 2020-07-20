# Install packages
pacman::p_load("tidyverse", "Biostrings", "seqinr", "DECIPHER", "Biostrings", "ggseqlogo", "cowplot", "gggenes", "ggtree")

# Set working directory
setwd("C:/Users/tabie/OneDrive/Documents/GitHub/TARA-thiolases/")

# Read in the sequences
# seqs0 <- read_csv("data/123_OleA_temps.csv")
# seqs0
# fastaseq <- AAStringSet(seqs0$sequence)
# names(fastaseq) <- seqs0$nams
# 
# writeXStringSet(fastaseq, "data/123_OleA.fasta")

seqs <- read.fasta("data/123_OleA.fasta", seqtype = "AA", strip.desc = TRUE)
seqs
# Create a data frame
genes <- attr(seqs, "name") #generates a gene list
seq.df <- data.frame(nams = attr(seqs, "name")) %>%
  dplyr::mutate(seqs = seqs)

##Calculating AA statistics for each sequence using AAstat() in the "seqinr" package

aastat.df <- data.frame(sapply(seq.df$seqs, AAstat)) # check out the plot

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

for(i in 2:35){
  new.df[i] <- as.numeric(unlist(new.df[i]))
}
class(new.df$Polar)
# new.df$nams2 <- 
# new.df$nams2
colnames(new.df)
# wide to long
new.df_long <- new.df %>% 
  gather(key = stat, value = value, -nams, -len)
new.df_long_large <- new.df %>%
  gather(key = stat, value = value, Tiny, ordered_res, NQrat, RKrat, Small, Nonpolar, Polar)
new.df_long_med <- new.df %>%
  gather(key = stat, value = value, disordered_res, Aliphatic, Charged)
new.df_long_small <- new.df %>%
  gather(key = stat, value = value, -nams, -len, -RKrat, -NQrat, -ordered_res, -disordered_res, -Tiny, -Small, -Aliphatic, -Nonpolar, -Polar, -Charged)

# All values
pdf("output/123_OleA__protprops_heatmap.pdf", height = 20, width = 12)
ggplot(new.df_long, aes(stat, nams)) +
  geom_tile(aes(fill = value)) + 
  # geom_text(aes(label = round(value, 4))) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 50, size = 8,
                                   vjust = 1.05, hjust = 1.2)) +
  scale_fill_gradient(low = "white", high = "darkgreen")

dev.off()

# Large values
pdf("output/123_OleA__protprops_heatmap_large.pdf", height = 20, width = 8)
ggplot(new.df_long_large, aes(stat, nams)) +
  geom_tile(aes(fill = value)) + 
  # geom_text(aes(label = round(value, 4))) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 50, size = 8,
                                   vjust = 1.05, hjust = 1.2)) +
  scale_fill_gradient(low = "white", high = "darkgreen")

dev.off()


# Small values
pdf("output/123_OleA__protprops_heatmap_small.pdf", height = 20, width = 8)
ggplot(new.df_long_small, aes(stat, nams)) +
  geom_tile(aes(fill = value)) + 
  # geom_text(aes(label = round(value, 4))) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 50, size = 8,
                                   vjust = 1.05, hjust = 1.2)) +
  scale_fill_gradient(low = "white", high = "darkgreen")

dev.off()

# Med values
pdf("output/123_OleA__protprops_heatmap_med.pdf", height = 20, width = 8)
ggplot(new.df_long_med, aes(stat, nams)) +
  geom_tile(aes(fill = value)) + 
  # geom_text(aes(label = round(value, 4))) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 50, size = 8,
                                   vjust = 1.05, hjust = 1.2)) +
  scale_fill_gradient(low = "white", high = "darkgreen")

dev.off()

temps123 <- read_csv("data/123_OleA_temps.csv")

allstats123 <- temps123 %>% 
  full_join(new.df, by = "nams")

plot1 <- ggplot(allstats123) +
  geom_boxplot(aes(x = temperature_range, y = NQrat))
plot1


write.csv(allstats123, "data/123_OleA_allratios.csv")

# Install packages
pacman::p_load("tidyverse", "Peptides", "seqinr", "ggplot2", "GGally")

# Set working directory
setwd("C:/Users/tabie/OneDrive/Documents/GitHub/TARA-thiolases/")

# Read in the protein sequences
sqs <-read.alignment("data/123_OleA.fasta", format = "fasta")
sqs$seq <- toupper(sqs$seq)
sqs$nam
sqs$seq
## Calculate protein sequence properties using the Peptides package
# Look at the Peptides documentation for explanations
# Reference manual: https://cran.r-project.org/web/packages/Peptides/Peptides.pdf
# Github repo: https://github.com/dosorio/Peptides
# Also available in the docs/Peptides_package_protein_properties.xlsx

# Cruciani 
crciani <- sapply(sqs$seq, crucianiProperties)
cdf <- do.call(rbind.data.frame, crciani)
colnames(cdf) <- names(crciani[[1]])
crucianiProperties(sqs$seq)
cdf
# Fasgai
fsgai <- sapply(sqs$seq, fasgaiVectors)
fdf <- do.call(rbind.data.frame, fsgai)
colnames(fdf) <- names(fsgai[[1]])

# ST scale
st_scale <- sapply(sqs$seq, stScales)
stdf <- do.call(rbind.data.frame, st_scale)
colnames(stdf) <- names(st_scale[[1]])

# T-scale
t_scale <- sapply(sqs$seq, tScales)
tdf <- do.call(rbind.data.frame, t_scale)
colnames(tdf) <- names(t_scale[[1]])

# VHSE
vhse_scale <- sapply(sqs$seq, vhseScales)
vdf <- do.call(rbind.data.frame, vhse_scale)
colnames(vdf) <- names(vhse_scale[[1]])

# Z-scale
z_scale <- sapply(sqs$seq, zScales)
zdf <- do.call(rbind.data.frame, z_scale)
colnames(zdf) <- names(z_scale[[1]])

# Kidera
kdera <- sapply(sqs$seq, kideraFactors)
kdf <- do.call(rbind.data.frame, kdera)
colnames(kdf) <- names(kdera[[1]])

# Combine everything into one data frame
sqdf <- data.frame(cbind(sqs$nam, sqs$seq), stringsAsFactors = F)
colnames(sqdf) <- c("nams", "sqs")
sqs$nam

combdf <- sqdf %>%
  # dplyr::mutate(acc = word(nams, sep = "\\.1", 1)) %>%
  dplyr::mutate(aliphat_full = aIndex(sqs)) %>%
  dplyr::mutate(hydrophob_full = hydrophobicity(sqs)) %>%
  dplyr::mutate(boman_interactions = boman(sqs)) %>%
  dplyr::mutate(hmoment = hmoment(sqs)) %>%
  dplyr::mutate(instab = instaIndex(sqs)) %>%
  dplyr::mutate(lngth_prot = lengthpep(sqs)) %>%
  dplyr::mutate(mw_prot = mw(sqs)) %>%
  dplyr::mutate(pi_prot = pI(sqs)) %>%
  bind_cols(., zdf, tdf, cdf, vdf, stdf, fdf, kdf)
view(combdf)
ratios <- read_csv("data/123_OleA_allratios.csv") %>% 
  select(-X1)
allstatsfinal <- ratios %>% 
  full_join(combdf, by = "nams")
duplicated(allstatsfinal)
table(allstatsfinal$temperature_range)
write_csv(allstatsfinal, "data/123_OleA_moststats.csv")

# oops add topt
allstatsfinal <- read_csv("data/123_OleA_moststats.csv")
full50_topt_only <- read_csv("data/full50_allstats.csv") %>% 
  select(nams, topt)
full73_topt_only <- read_csv("data/73_OleA_tomer_results.csv") %>% 
  mutate(topt = Topt,
         nams = word(Sequence, sep = "_", start = 3, end = -1)) %>%
  select(nams, topt)

full73_topt_only$nams == allstatsfinal$nams[1:73]
full73_topt_only[1,1] <- "Xanthomonas_campestris_pv._campestris_str._ATCC_33913"
full73_topt_only[46,1] <- "Streptococcus_pneumoniae"
full73_topt_only[51,1] <- "Pelobacter_propionicus_DSM_2379"
full73_topt_only[63,1] <- "Chloroflexi_bacterium_RBG_13_51_36"
full73_topt_only[20,1] <- "Lysobacter_tolerans"
allstatsfinal[49,1] <- "Kocuria_varians_strain_G6"
full73_topt_only[52,1] <- "Pseudoxanthomonas_sp._NML171200"
allstatsfinal[52,1] <- "Pseudoxanthomonas_sp._NML171200"
allstatsfinal[54,1] <- "Leifsonia_sp._Leaf325"
full73_topt_only[54,1] <- "Leifsonia_sp._Leaf325"
allstatsfinal[57,1] <- "Enhygromyxa_salina_strain_SWB007"
allstatsfinal[59,1] <- "Halobacteriovorax_marinus_strain_BE01"
full73_topt_only[59,1] <- "Halobacteriovorax_marinus_strain_BE01"
allstatsfinal[69,1] <- "Xanthomonas_translucens_pv._graminis"
full73_topt_only[69,1] <- "Xanthomonas_translucens_pv._graminis"
all_topt <- full_join(full50_topt_only, full73_topt_only)
full123_allstats <- full_join(all_topt, allstatsfinal)
write_csv(full123_allstats, "data/123_OleA_allstats.csv")
