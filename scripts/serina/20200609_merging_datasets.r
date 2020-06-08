# Load packages
pacman::p_load("tidyverse")

# Set working directory
setwd("C:/Users/tabie/OneDrive/Documents/GitHub/TARA-thiolases/")

# Code from the following tutorial
# https://psu-psychology.github.io/r-bootcamp/tutorials/joins_tutorial.html

# Create the superheroes dataset
superheroes <-"
name, alignment, gender,         publisher
Magneto,       bad,   male,            Marvel
Storm,      good, female,            Marvel
Mystique,       bad, female,            Marvel
Batman,      good,   male,                DC
Joker,       bad,   male,                DC
Catwoman,       bad, female,                DC
Hellboy,      good,   male, Dark Horse Comics
"
superheroes <- read_csv(superheroes, trim_ws = TRUE, skip = 1)

# Create the publishers dataset
publishers <- "
publisher, yr_founded
DC,       1934
Marvel,       1939
Image,       1992
"
publishers <- read_csv(publishers, trim_ws = TRUE, skip = 1)

########
# Try different types of 'joins'
# Function: inner_join
# Require smatch in both datasets (non-matching rows are dropped)
ijsp <- inner_join(superheroes, publishers, by="publisher")

########
# Function: left_join
# Keep all rows in left-hand ‘x’ dataset (i.e., superheroes). 
# Add columns from publishers where there is a match. Fill in NA for non-matching observations.
ljsp <- left_join(superheroes, publishers)

########
# Function: right_join
# Keep all rows in right-hand ‘y’ dataset (i.e., publishers). 
# Add columns from superheroes where there is a match. 
rjsp <- superheroes %>%  # WHOA what is THIS '%>%' ??? see below
  right_join(publishers)

# The '%>%' is called a 'pipe' 
# "Pipes let you take the output of one function and send it directly to the next"
# This is useful when you need to many things to the same data set
#  You can find a massive tutorial on pipes and their history here:
# https://www.datacamp.com/community/tutorials/pipe-r-tutorial

# So using the pipe was another way to write this:
rjsp2 <- right_join(superheroes, publishers)
rjsp == rjsp2 # remember this '=='? You can check if the two objects match
# Do they?


########
# Function: full_join
# Keep all rows in left-hand ‘x’ (superheroes) and right-hand ‘y’ (publishers) datasets. 
# Resulting dataset will have all columns of both datasets, but filling in NA for any non-matches on either side.
fjsp <- full_join(superheroes, publishers)


# Challenge 1. Complete the rest of the tutorial after full_join here:
# https://psu-psychology.github.io/r-bootcamp/tutorials/joins_tutorial.html


# Challenge 2. Joining is just part of the art of 'data wrangling.' How do you subset, filter, create new columns?
# The packages dplyr and tidyr in the tidyverse are great for this.
# I really love the tutorial here: https://datacarpentry.org/R-ecology-lesson/03-dplyr.html

# You'll be using the following dataset to complete the tutorial
surveys <- read_csv("data/tutorials/surveys.csv")
view(surveys)
# Go through the tutorial and complete the challenges
# https://datacarpentry.org/R-ecology-lesson/03-dplyr.html
