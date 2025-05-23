# Install packages
pacman::p_load("ggtree", "patchwork", "tidyverse")

# Set working directory
setwd("C:/...") # fill in your working directory

# First let's practice reading in a tree in Newick format 
nj <- read.tree("data/trees/20_psychro_thermo_seqs_NJ_500boot.nwk") # neighbor-joining method
ml <- read.tree("data/trees/20200617_20_psychro_thermo_ML_500boot.nwk") # maximum-liklihood method

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
  xlim(0, 4) # 4 is an arbitrary unit from the right side of the page
# try playing around with this number
njt3

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
  dplyr::mutate(fixnam = paste0(word(label, sep = "_", 2), " ",
                                word(label, sep = "_", 3)))  # WHOA what is happening here??

# Then merge the data frame with the original tree using %<+%
mlt_append <- mlt %<+% mlt_df #oh my goodness what a ridiculous operator %<+% is

# Now plot with updated labels
mlt_append +
  xlim(0, 4) +
  geom_tiplab(aes(label = fixnam, color = temp_status)) +
  geom_nodelab() 

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

# Challenge 3. The node labels (bootstrap values) are really ugly...way too many decimal places..
# and could you make them range between 0 and 100 instead of 0 and 1?
# Hint: make a new column where you transform the node labels
# Then plot these within geom_nodelab(label = ...)

# Challenge 4. What happens if you add geom_point(x = 4) to your tree ?
# Could you color the points to correspond to the optimal temperatures (if known)

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










