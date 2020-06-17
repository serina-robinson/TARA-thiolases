# Load packages
pacman::p_load("gggenes", "ggplot2", "tidyverse")

# Set working directory
setwd("C:/Users/tabie/OneDrive/Documents/GitHub/TARA-thiolases/") # fill in your working directory

# gggenes is a package built on ggplot2
# For drawing gene arrow maps
# First try an example with the built-in dataset called example_genes
exgenes <- example_genes
ggplot(exgenes, 
       aes(xmin = start, xmax = end, y = molecule, fill = gene)) +
  geom_gene_arrow() +
  facet_wrap(~ molecule, scales = "free", ncol = 1) + # what does this do?
  scale_fill_brewer(palette = "Set3") +
  theme_genes()

# Let's say we want to align by genE
# We'll use make_alignment_dummies
dummies <- make_alignment_dummies(
  exgenes,
  aes(xmin = start, xmax = end, y = molecule, id = gene),
  on = "genF"
)
dummies
exgenes
class(dummies)
ggplot(exgenes, aes(xmin = start, xmax = end, y = molecule, fill = gene)) +
  geom_gene_arrow() +
  geom_blank(data = dummies) +  #what does geom_blank do?
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3") +
  theme_genes()
?geom_blank
# You can find and try out the rest of the gggenes tutorial here:
# https://cran.r-project.org/web/packages/gggenes/vignettes/introduction-to-gggenes.html

# Now let's try it with our data!
# First we need to get our data into a format that matches example_genes...

# To do this go to PATRIC
# Run a BLAST, for example for OleA
# against your custom genome group
# Then 'SELECT ALL'
# Then in the right-hand panel click
# Create a 'Feature list'
# Click the black plus sign to expand features
# Add start, end, and strand
# These are the features we need for plotting!
# Then click download as CSV

# Challenge 1. BLAST OleA, B, C, and D sequences
# against your custom gene groups and download the feature lists
# Then run the code below with your new feature lists
# The example_features.csv files below are incomplete

# Here are some example feature lists
olea <- read_csv("data/feature_lists/OleA_example_features.csv") %>%
  mutate(gene = "oleA")
oleb <- read_csv("data/feature_lists/OleB_example_features.csv") %>%
  mutate(gene = "oleB")
olec <- read_csv("data/feature_lists/OleC_example_features.csv") %>%
  mutate(gene = "oleC")
oled <- read_csv("data/feature_lists/OleD_example_features.csv") %>%
  mutate(gene = "oleD")
view(olea)




# Combine everything using the bind_rows function
combined <- olea %>%
  bind_rows(oleb, olec, oled) %>%
  janitor::clean_names() %>%
  dplyr::rename(molecule = genome) %>% 
  group_by(gene, molecule) %>% 
  dplyr::slice(1)

# Filter for only the sequences that are in your dataset
full50 <- read_csv("data/full50_4.csv")
gens_to_keep <- full50$genome_x
gens_to_keep
plotdat <- combined %>%
  dplyr::filter(molecule %in% gens_to_keep) %>%
  dplyr::select(molecule, gene, start, end, strand) %>%
  dplyr::mutate(direction = case_when(strand == "+" ~ 1,
                                       strand == "-" ~ -1)) %>%
  dplyr::mutate(strand = case_when(direction == 1 ~ "forward",
                                   direction == -1 ~ "reverse")) %>%
  arrange(molecule)
plotdat
exgenes
# Make a gene plot
pdf("output/gene_diagram.pdf")
ggplot(plotdat, 
       aes(xmin = start, xmax = end, fill = gene, y = molecule)) +
  geom_gene_arrow() +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3") +
  theme_genes()
dev.off()
# You can see that in some of the genomes the OleA is missing, in others
# OleC only...at least Moritella sp. JT01 looks good! Everything else is 
# looking a little rough.

# Challenge 2. Recreate the plot above with the complete set
# of 20 gene neighborhoods for your thermophiles/psychrophiles

# Challenge 3. Use the make_alignment_dummies() function to
# align the gene neighborhoods by oleA

# Challenge 4. Use geom_gene_label() function to label 
# the gene clusters. See gggenes tutorial for more info. 
# https://cran.r-project.org/web/packages/gggenes/vignettes/introduction-to-gggenes.html

# Challenge 5 (Required). Download and install the phylogenetic tree software
# MEGA X. You'll want the GUI version for Windows. We'll be using it on Thursday June 18.
# If you have tons of extra time you can look over the documentation.

# Super/optional Challenge 6. Try adding the PfaA gene to
# the gene neighborhood diagrams. You will first need to pull the
# feature list from BLASTing your genome group with PfaA in data/query_seqs.
# The diagram is probably going to look pretty rough because in most cases
# The pfaA genes won't be adjacent to the oleABCD genes. That's ok!

