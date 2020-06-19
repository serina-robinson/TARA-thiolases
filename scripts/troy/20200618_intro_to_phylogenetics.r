# Install packages
pacman::p_load("ggtree", "patchwork", "tidyverse")

# Set working directory
setwd("C:/Users/tabie/OneDrive/Documents/GitHub/TARA-thiolases/") # fill in your working directory

# First let's practice reading in a tree in Newick format 
nj <- read.tree("data/trees/20_psychro_thermo_seqs_NJ_500boot.nwk") # neighbor-joining method
ml <- read.tree("data/trees/20_olea_ML_50_psychro_thermo_tree.nwk") # maximum-liklihood method

# Convert the trees into a ggplot format using the amazing package ggtree
# It has it's own book! https://yulab-smu.github.io/treedata-book/index.html
njt <- ggtree(nj) # neighbor-joining
mlt <- ggtree(ml) # maximum-likelihood

# Using the patchwork package you can plot them side by side using a '+'
# https://www.datanovia.com/en/blog/ggplot-multiple-plots-made-ridiculuous-simple-using-patchwork-r-package/
njt + mlt
# But where are the words?

# Add tip labels 
njt2 <- njt +
  geom_tiplab() 
njt2

# Resize tree so it doesn't go off the page
njt3 <- njt2 +
  xlim(0, 10) # 4 is an arbitrary unit from the right side of the page
# try playing around with this number
njt3
mlt$data$label
# Putting it all together
# There is so much you can do! For example, let's display the bootstrap values
mlt_boot <- mlt + 
  geom_tiplab() +
  xlim(0, 4) +
  geom_nodelab()
mlt_boot
# Compare to the bootstrapped tree to check it works as expected

# What if we want to change anything/append new variables
# For example if we wanted to fix the names
# First create a data frame with the first column called label
# Where label corresponding to the phylogenetic tree labels
# That are stored in tree$data$label (see below)
# Then add whatever additional columns you want
mlt_df <- data.frame(label = mlt_boot$data$label) %>%
  dplyr::mutate(fixnam = paste0(word(label, sep = "_", 2)))  # WHOA what is happening here??
?paste0

mlt_df2 <- mlt_df %>% 
  mutate(temp_status = word(label, sep = "_", -1))

mlt_df3 <- mlt_df2 %>% 
  mutate(fixnam2 = as.numeric(label)*100)
mlt_df3[21,"fixnam2"] <- " "

mlt_df4 <- mlt_df3 %>% 
  dplyr::mutate(temp_status2 = case_when(grepl("Deltaproteobacteria", label) ~ "psychrophile",
                                         grepl("Epsilonproteobacteria", label) ~ "psychrophile",
                                         grepl("Verrucomicrobia", label) ~ "psychrophile",
                                         grepl("Psychrobacter", label) ~ "mesophile",
                                         grepl("Sedimentitalea", label) ~ "mesophile",
                                         TRUE ~ temp_status))


full50 <- read_csv("data/full50_4.csv")
full50[43, "genus"] <- "Psychrobacter"
full50[39, "genus"] <- "Deltaproteobacteria1"
full50[40, "genus"] <- "Deltaproteobacteria2"
mlt_df4[17, "fixnam"] <- "Deltaproteobacteria1"
mlt_df4[13, "fixnam"] <- "Deltaproteobacteria2"
view(full50)
mlt_df5 <- mlt_df4 %>% 
  left_join(full50, by = c("fixnam" = "genus"))

mlt_df6 <- mlt_df5 %>% 
  mutate(temp_color = case_when(temp_status2 == "thermophile" ~ "red",
                                temp_status2 == "mesophile" ~ "purple",
                                temp_status2 == "psychrophile" ~ "blue",
                                TRUE ~ "gray"))
view(mlt_df6)
#### fix rhodonellum, sedimentitalea, deltaproteobacteria replicates
# Then merge the data frame with the original tree using %<+%
mlt_append <- mlt %<+% mlt_df6 #oh my goodness what a ridiculous operator %<+% is
table(mlt_append$data$temp_status2)
table(mlt_append$data$temp_color)
# Now plot with updated labels
pdf("output/tree.pdf", width = 6, height = 5)
mlt_append +
  xlim(0, 4) +
  scale_color_gradient(low = "yellow", high = "red") +
  geom_tippoint(x = 3.5, aes(color = temperature)) +
  geom_tiplab(aes(label = fixnam), color = mlt_df6$temp_color[1:20]) +
  geom_nodelab(aes(label = fixnam2), hjust = 1.2, vjust = -0.25) 
  
  
dev.off()
as.factor(mlt_df5$temp_status2[1:20])
mlt_df6$temp_color[mlt_append$data$isTip]
?geom_point
?geom_nodelab
view(mlt_df6)
view(full50)
?geom_nodelab
view(mlt_append$data$fixnam2)
# Challenge 1. Color the tree leaves by whether or not the sequence is a psychrophile 
# or a thermophile. Hint: you can create a new column in mlt_df that is called temp_status
# or something like that. Then you use the 'color' option within the aes() of geom_tiplab
# to color by the status

  
# Challenge 1.5 If you wanted to change the colors of the tip labels
# to be say red and blue..how would you do that? 
# Scale_color_manual is helpful here
# for example try adding scale_color_manual(values = c("blue", "red)) 
# or any other colors of choice

# Challenge 2. Some of the so called "thermophiles" in this dataset are not actually
# thermophiles after you looked up their sampling temperature in the hydrothermal vent.
# Let's say you want to replace those incorrectly labeled thermophiles with psychrophiles
# How would you do this? Hint: case_when() could be useful here. Try adding a new column
# to mlt_append and using the 'grepl' function that searches for partial word matches and 
# evaluates to TRUE/FALSE to search for the incorrect 
# thermophiles and convert them to psychrophiles

view(mlt_df4)
# Challenge 3. The node labels (bootstrap values) are really ugly...way too many decimal places..
# and could you make them range between 0 and 100 instead of 0 and 1?
# Hint: make a new column where you transform the node labels
# Then plot these within geom_nodelab(label = ...)

           # case_when(is.numeric(label) ~ mlt_df2$label * 100,
           #                  is.character(label) ~ fixnam,
           #                  TRUE ~ fixnam))
# Challenge 4. What happens if you add geom_point(x = 4) to your tree ?
# Could you color the points to correspond to the optimal temperatures (if known)


view(mlt_df5)
view(full50)
# Challenge 5. Write your finished tree to a pdf with dimensions that look good.
# You might need to change the width and height as follows:
# pdf("output/tree.pdf", width = 5, height = 10)

# Extra fun: just try skimming the ggtree book and playing around with all the cool things
# you can do. Try adding a picture. Or zooming in on one clade. Or changing 
# the angle or layout of the tree. The options are endless!

# Super challenge 6. If you're feeling ambitious, you could try messing around 
# with matching up the gggenes with the tree. We'll get to that Friday if not
# Here is the tutorial
# https://yulab-smu.github.io/treedata-book/chapter11.html#genome-locus










