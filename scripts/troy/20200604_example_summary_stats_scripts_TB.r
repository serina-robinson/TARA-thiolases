# Install packages
# First time run:
 if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
 BiocManager::install("Biostrings")

# Load packages
pacman::p_load("Biostrings", "tidyverse")

# Set working directory
# On a a Mac:
# setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/TARA-thiolases/")
# On a PC:
 setwd("C:/Users/tabie/OneDrive/Documents/GitHub/TARA-thiolases")

# Read in the TARA sequences only
tara <- readAAStringSet('data/TARA_marine_sequences_to_order.fasta') # these are the 30 sequences from TARA oceans
names(tara)

# Read in the sequence metadata
annot_df <- read_csv("data/TARA_tree_annotation_df.csv")
view(annot_df) # View the data frame in a different window
colnames(annot_df) # What sorts of variables are in this data frame
dim(annot_df) # Check the dimensions, 244 rows and 55 columns

# But wait! we only have 30 sequences we actually ordered 
# Subset the data frame...
twenty_nine <- annot_df[annot_df$label %in% names(tara),]
dim(twenty_nine) # 29 rows and 55 columns

### Challenge 1. Can you find which sequence in the TARA sequences (tara) wasn't found with the above command? 
## There are only 29 rows but there are 30 TARA sequences...
## The following command will help you:
## names(tara)[!names(tara) %in% annot_df$label]
## Why wasn't the missing row found?
## Bonus: Can you make the data frame complete with thirty rows
## Hint: you can use the function bind_rows to bind a row to an existing data frame
## thirty <- bind_rows(twenty_nine, your_missing_row)

# Now that we've got our data frame...let's make some plots
dat <- twenty_nine # if you find the missing row you can change this to:
# dat <- thirty 
# If you wanted to look at all the data (not just the thirty), you could do:
# dat <- annot_df
colnames(dat) # these variables are available to play with

# First we can look at histograms of different continuous variables
# For example:
ggplot(dat) +
  geom_histogram(aes(temperature)) +
  theme_classic()

# Challenge 2. Can you adust the histogram bin size?


# Now we can look at relationships between continuous variables
# For example, a scatterplot of oxygen vs. temperature
ggplot(dat) +
  geom_point(aes(x = oxygen, y = temperature)) +
  theme_classic()

# Let's say you want to write it to file. Here is one way to do it:
pdf("output/temp_vs_oxygen.pdf")
ggplot(dat) +
  geom_point(aes(x = oxygen, y = temperature)) +
  theme_classic()
dev.off()

# Challenge 3. Can you color the points by depth?

# We can also look at categorical variables using boxplots
ggplot(dat) +
  geom_boxplot(aes(x = polar, y = temperature)) +
  geom_point(aes(x = polar, y = temperature)) +
  theme_classic() 


# Challenge 4. Can you remove/fix the NA? 


# BONUS challenges: What other variables can you look at? What other types of plots  can you make? 
