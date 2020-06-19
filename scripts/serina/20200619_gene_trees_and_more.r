# Install packages
pacman::p_load("gggenes", "ggtree", "ggplot2", "tidyverse", "plotly")

### All of this code below is calculate a phylogenetic tree for example_genes
# The geom_motif is defined in ggtree and it is a wrapper layer of gggenes::geom_gene_arrow. 
# The geom_motif can automatically adjust genomic alignment by selective gene (via the on parameter) 
# and can label genes via the label parameter.

# As the example_genes dataset only provide genomic coordinations of a set of genes, 
# a phylogeny for the genomes need to be firstly constructed. 
# We calculate jaccard similarity based on the ratio of overlapping genes
# among genomes and correspondingly determine genome distance. 


##### All of this is just to generate a phylogenetic tree ###
get_genes <- function(data, genome) {
  filter(data, molecule == genome) %>% pull(gene)
}

g <- unique(example_genes[,1])
n <- length(g)
d <- matrix(nrow = n, ncol = n)
rownames(d) <- colnames(d) <- g
genes <- lapply(g, get_genes, data = example_genes)

for (i in 1:n) {
  for (j in 1:i) {
    jaccard_sim <- length(intersect(genes[[i]], genes[[j]])) / 
      length(union(genes[[i]], genes[[j]]))
    d[j, i] <- d[i, j] <- 1 - jaccard_sim
  }
}
d
tree <- ape::bionj(d) 
#########

tree # This looks just like any newick tree that you read in from MEGAX

# Example gene plot
p <- ggtree(tree, branch.length='none') + 
  geom_tiplab() + xlim_tree(5.5) + 
  geom_facet(mapping = aes(xmin = start, xmax = end, fill = gene),
             data = example_genes, geom = geom_motif, panel = 'Alignment',
             on = 'genE', label = 'gene', align = 'left') +
  scale_fill_brewer(palette = "Set3") + 
  scale_x_continuous(expand=c(0,0)) +
  theme(strip.text=element_blank(),
        panel.spacing=unit(0, 'cm'))

## in case the facet panels were not ordered properly
p <- p + facet_grid(cols = vars(factor(.panel, levels = c("Tree", "Alignment"))),
                    scales = 'free_x')
p


#### Now let's do it with our genes
# This code is copied from the 20200617 script

# Here are some example feature lists
olea <- read_csv("data/feature_lists/OleA_complete_features.csv") %>%
  mutate(gene = "oleA")
blastolea <- read_csv("data/blast_output/BLAST_OleA_genome_group.csv")
combolea <- olea %>% 
  left_join(blastolea) %>% 
  janitor::clean_names() %>% 
  filter(query_cover > 75, identity > 25) %>% 
  arrange(desc(identity)) %>% 
  group_by(genome, gene) %>% 
  slice(1)
combolea


oleb <- read_csv("data/feature_lists/OleB_complete_features.csv") %>%
  mutate(gene = "oleB")
blastoleb <- read_csv("data/blast_output/BLAST_OleB_genome_group.csv")
# view(blastoleb)
# view(oleb)
# dim(oleb)
# setdiff(oleb$Genome, blastoleb$Genome)
comboleb <- oleb %>% 
  left_join(blastoleb) %>% 
  janitor::clean_names() %>% 
  filter(query_cover > 75, identity > 25) %>% 
  arrange(desc(identity)) %>% 
  group_by(genome, gene) %>% 
  slice(1)
comboleb
table(comboleb$genome, comboleb$gene)


olec <- read_csv("data/feature_lists/OleC_complete_features.csv") %>%
  mutate(gene = "oleC")
blastolec <- read_csv("data/blast_output/BLAST_OleC_genome_group.csv")
combolec <- olec %>% 
  left_join(blastolec) %>% 
  janitor::clean_names() %>% 
  filter(query_cover > 75, identity > 25) %>% 
  arrange(desc(identity)) %>% 
  group_by(genome, gene) %>% 
  slice(1)
combolec


oled <- read_csv("data/feature_lists/OleD_complete_features.csv") %>%
  mutate(gene = "oleD")
blastoled <- read_csv("data/blast_output/BLAST_OleD_genome_group.csv")
comboled <- oled %>% 
  left_join(blastoled) %>% 
  janitor::clean_names() %>% 
  filter(query_cover > 75, identity > 25) %>% 
  arrange(desc(identity)) %>% 
  group_by(genome, gene) %>% 
  slice(1)
comboled



view(olea)
# Combine everything using the bind_rows function
combined <- combolea %>%
  bind_rows(comboleb, combolec, comboled) %>%
  janitor::clean_names() %>%
  dplyr::rename(molecule = genome) %>% 
  group_by(molecule, gene) %>% 
  slice(1) %>% 
  ungroup()
view(combined)
# Filter for only the sequences that are in your dataset
full50 <- read_csv("data/full50_4.csv")
gens_to_keep <- full50$genome_x
# plotdat <- combined
plotdat <- combined %>%
 # dplyr::filter(molecule %in% gens_to_keep) %>%
  dplyr::select(molecule, gene, start, end, strand) %>%
  dplyr::mutate(direction = case_when(strand == "+" ~ 1,
                                      strand == "-" ~ -1)) %>%
  dplyr::mutate(strand = case_when(direction == 1 ~ "forward",
                                   direction == -1 ~ "reverse")) %>%
  arrange(molecule)

# Make a gene plot - check it is still working
gpl <- ggplot(plotdat, 
       aes(xmin = start, xmax = end, fill = gene, y = molecule)) +
  geom_gene_arrow() +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3") +
  theme_genes()
pdf("output/tree3.pdf")
gpl
dev.off()
# Read in the phylogenetic tree
ml <- read.tree("data/trees/20200617_20_psychro_thermo_ML_500boot.nwk") # maximum-likelihood method
mlt <- ggtree(ml)
view(plotdat)
# Uh oh...we have a problem...the tip labels don't match the molecule names!
mlt$data$label %in% plotdat$molecule # all FALSE
mlt$data$label
# How to fix???
# Do a join between the tree tip labels and the gene diagram plotting coordinates...by "genus"
tree_df <- data.frame(label = mlt$data$label) %>%
  dplyr::mutate(genus = stringr::word(label, sep = "_", 2)) %>%
  dplyr::filter(!is.na(genus))
duplicated(tree_df) # are any duplicated??
tree_df$genus[grep("Deltaproteobacteria", tree_df$genus)][1] <- "Deltaproteobacteria1"
tree_df$genus[grep("Deltaproteobacteria", tree_df$genus)][2] <- "Deltaproteobacteria2"
plotdat_df <- plotdat %>%
  dplyr::mutate(genus = stringr::word(molecule, sep = " ", 1))
view(plotdat_df)
table(plotdat_df$genus)
view(tree_df)
# Join them by = "genus"
# What sort of gymnastics is happening here?
plotdat_fixed <- plotdat_df %>%
  left_join(., tree_df, by = "genus") %>%
  dplyr::select(-molecule, -genus) %>%
  dplyr::rename(molecule = label) %>%
  dplyr::select(molecule, gene, start, end, strand, direction) %>%
  arrange(molecule, gene)  %>%
  group_by(molecule) %>%
  dplyr::mutate(difference = start - first(start)) %>% # what the heck am I doing here???
  dplyr::filter(abs(difference) < 10000) # what about here
view(plotdat_fixed)

mltr2 <- mlt + 
  geom_tiplab() + 
  xlim_tree(15) + 
  geom_facet(mapping = aes(xmin = start, xmax = end, fill = gene),
             data = plotdat_fixed, geom = geom_motif, panel = 'Alignment',
             on = "oleA", label = 'gene', align = 'left') +
  scale_fill_brewer(palette = "Set3") + 
  scale_x_continuous(expand=c(0,0)) +
  theme(strip.text=element_blank(),
        panel.spacing=unit(0, 'cm'))
# well that looks bad

## in case the facet panels were not ordered properly
mltr3 <- mltr2 + facet_grid(cols = vars(factor(.panel, 
                    levels = c("Tree", "Alignment"))),
                    scales = 'free_x')
pdf("output/tree2.pdf", width = 10, height = 7)
mltr3 
dev.off()

# Challenge 1. Add the gene diagrams for all 20 leaves (not just the 4 plotted currently...)

# Challenge 2. You will probably get pretty terrible looking gene diagrams...
# To fix it will be tricky...but I believe you can do it!
# What we want to do is keep ALL oleAs and then only plot the oleBs, oleCs, and oleDs...
# if they're within say an approximately 10,000 base pair window

# One possible way to do this is first to group by the molecule
# Then create an additional column called 'difference'
# That calculates the difference between oleB, oleC, and oleD and oleA, which is the 'first'
# entry per molecule (see code above, also here):
plotdat2 <- plotdat %>%
  arrange(molecule, gene)  %>%
  group_by(molecule) %>%
  dplyr::mutate(difference = start - first(start)) %>% # what the heck am I doing here???
  dplyr::filter(abs(difference) < 10000)
plotdat2

# Bonus challenge; is there a way to select the oleA rather than just the first 'start' entry?

# Challenge 3. 
# Fix the tree like you did yesterday, with improved names, fixed psychrophile/thermophile status
# Use geom_point to add a point for the temperature (challenge from yesterday)
# If you run into overlapping labels try adjusting xlim_tree()
# Add the bootstrap values
# Clean up the tree however you see fit! 
# Write the tree to file and prepare to show in lab meeting on Tuesday

# Super Challenge 4.
# Add a new tab to your web app for the phylogenetic tree
# Unfortunately the tree won't be interactive but you can still display it 
# Remember you'll need to add library(ggtree) to the r setup section
# And you need to navigate to the right folder to read in the tree (as shown below)
# Can you read in and display your final version of your tree in your app?

# ```{r }
# ml <- read.tree("../../data/trees/20200617_20_psychro_thermo_ML_500boot.nwk") # maximum-likelihood method
# mlt <- ggtree(ml)
# mlt
# ```

# Challenge 5. Read the article for journal club!

# Preparation for next week!

# We're going to need to download some more open-source software :)
# Unless you already have PyMOL downloaded (that's great)

# PyMOL
# https://pymol.org/2/

#  Also create a user account on Phyre2:
# http://www.sbg.bio.ic.ac.uk/phyre2/html/login.html





