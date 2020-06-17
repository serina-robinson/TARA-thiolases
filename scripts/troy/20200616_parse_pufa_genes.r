# Install packages
pacman::p_load("tidyverse", 
               "Biostrings",
               "janitor",
               "viridis")

# Set working directory
 setwd("C:/Users/tabie/OneDrive/Documents/GitHub/TARA-thiolases/") # fill in your working directory
 
# Today we're going to be focused on
# PUFA genes and genome neighborhoods 

# In order to look for PUFA genes, 
# You start by BLASTing the Pfa synthase gene
# A type I polyketide synthase called PfaA (below)
pfaA <- readAAStringSet("data/query_seqs/PfaA.fasta")
head(pfaA)
width(pfaA)

# Go to your most favorite site: PATRIC!
# Create a genome group of just your psychro and thermophiles
# of interest! 
# This is described in the Challenges below (Challenge 2)

# Once you've created your genome group
# You can now BLAST your query gene (pfaA) within PATRIC
# https://www.patricbrc.org/app/BLAST
# Download the results as a CSV

# Read in the BLAST results
blast <- read_csv("data/blast_output/20200615_TARA_psychrophiles_PUFA_blast.csv") %>%
  janitor::clean_names() # one of my favorite functions of all times! Can you tell what it does?
view(blast)
spec(blast)

# Read in your mapping dataset

full50 <- read_csv("data/full50_4.csv")

full50$genome_x # this is the column we are interested in!
view(full50)
# Challenge 1.
# Combine your BLAST results with genome.x in your full50 dataset... how would you do this?
# Hint: the argument 'by' in the 'join' family of functions is helpful.
# If you do it correctly you should get 1 TRUE and 49 FALSE
merged <- full50 %>% 
  left_join(blast, by = c("genome_x" = "genome"))
merged$start
# Challenge 2. Create your own psychrophile/thermophile genome
# group using in PATRIC and search for your 20 genome IDs to create the groups. 
# Follow the steps here: 
# https://docs.patricbrc.org/tutorial/genome_groups/creating_genome_groups.html

# Challenge 3. BLAST the PfaA gene against the genome group you've created and download the results.

# Challenge 4. Also BLAST your PATRIC genomes with the OleB gene in 
# the query_seqs/oleABCD.fasta file
# You may want to separate the OleB into a separate FASTA file.
# You can use a text editor like Notepad to create new FASTA files.
# At some point you will want to choose and download a text editor you like. 
# Here are some options: https://kinsta.com/blog/best-text-editors/

# Challenge 5. Read in the new BLAST results and join with full50 again. Compare the results to Challenge 1.
# Are there any 'thermophiles' that appear to have PUFA genes? What about OleBs?
OleB <- read_csv("data/BLAST_OleB_genome_group.csv") %>% 
  janitor::clean_names() %>% 
  filter(query_cover > 75, identity >= 25) %>% 
  group_by(genome) %>% 
  arrange(identity) %>% 
  slice(1)
OleB
duplicated(OleB$genome)

merged_OleB <- full50 %>% 
  left_join(OleB, by = c("genome_x" = "genome"))
merged_OleB1 <- merged_OleB %>% 
  mutate(ole_b_patric = merged_OleB$patric_id)

  
pfa <- read_csv("data/BLAST_pufa_genome_group.csv") %>% 
  janitor::clean_names() %>% 
 # filter() %>% 
  group_by(genome) %>% 
  arrange(identity) %>% 
  slice(1)
view(pfa)


merged_pfa <- full50 %>% 
  left_join(pfa, by = c("genome_x" = "genome"))
merged_pfa1 <- merged_OleB %>% 
  mutate(pfa_patric = merged_pfa$patric_id)
# mutate new column "pufa"
fullmerged <- left_join(merged_OleB1, merged_pfa1)
view(fullmerged)
# Challenge 6. Make a plot of whether each genome has a hit for PfaA, OleA, B, C, and D
# You can use ggplot with geom_tile() See example code below. Note this is an abbreviated
# version using only the OleACD hits from the full50 dataset - you will want to expand this
# to also include a column if there is a PfaA or OleB hit from your BLAST search.

# Pull out the genome positions for OleA, C, and D using the
# function in the stringr package called 'word'
# You can use ?word to look it up
?word
dat_wide <- fullmerged %>%
  dplyr::filter(!is.na(ole_a_patric)) %>%
  dplyr::mutate(OleA_pos = stringr::word(ole_a_patric, sep = "\\.", -1),
                OleB_pos = stringr::word(ole_b_patric, sep = "\\.", -1),
                OleC_pos = stringr::word(ole_c_patric, sep = "\\.", -1),
                OleD_pos = stringr::word(ole_d_patric, sep = "\\.", -1),
                Pfa_pos = stringr::word(pfa_patric, sep = "\\.", -1)) %>%
  dplyr::select(genome_x, OleA_pos, OleB_pos, OleC_pos, OleD_pos, Pfa_pos)

# Convert to long format for ggplot  
dat_long <- dat_wide %>%
  gather(key = gene.id, value = gene.pos, -genome_x) %>%
  mutate(gene.pos = as.numeric(gene.pos))


# Plot data
pdf("output/OleACD_heatmap.pdf", height = 20, width = 10)
ggplot(dat_long, aes(gene.id, genome_x)) +
  geom_tile(aes(fill = gene.pos)) + 
  geom_text(aes(label = round(gene.pos, 1))) + 
  scale_fill_gradient(low = "dodgerblue", high = "red") +
  theme_classic()
dev.off()

abc <- read_csv("data/blast_output/PfaA_all_sequences.csv") %>% 
  janitor::clean_names()

view(abc)

stringr::word(merged_OleB1$ole_a_patric, sep = "|", -1)



# Try messing around with the plot. You can change the color gradient for example.
# You can also make other types of plots using the geom_* functions