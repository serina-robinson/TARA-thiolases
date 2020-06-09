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
superheroes[3,3]
pull(superheroes, alignment)
superheroes$alignment
# Create the publishers dataset
publishers <- "
publisher, yr_founded
DC,       1934
Marvel,       1939
Image,       1992
"
publishers <- read_csv(publishers, trim_ws = TRUE, skip = 1)

publishers$yr_founded <- as.integer(publishers$yr_founded)
publishers
########
# Try different types of 'joins'
# Function: inner_join
# Requires match in both datasets (non-matching rows are dropped)
ijsp <- inner_join(superheroes, publishers, by="publisher")
superheroes
ijsp
########
# Function: left_join
# Keep all rows in left-hand dataset (i.e., superheroes). 
# Add columns from publishers where there is a match. Fill in NA for non-matching observations.
ljsp <- left_join(superheroes, publishers)
ljsp
########
thing <- 10
# Function: right_join
# Keep all rows in right-hand ‘y’ dataset (i.e., publishers). 
# Add columns from superheroes where there is a match. 
rjsp <- superheroes %>%  # WHOA what is THIS '%>%' ??? see below
  right_join(., publishers) %>% 
  mutate(thing = rep(c(10,11), 3))
rjsp
?rep
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
# Keep all rows in left-hand (superheroes) and right-hand ‘y’ (publishers) datasets. 
# Resulting dataset will have all columns of both datasets, but filling in NA for any non-matches on either side.
fjsp <- full_join(superheroes, publishers)
fjsp

# Challenge 1. Complete the rest of the tutorial after full_join here:
# https://psu-psychology.github.io/r-bootcamp/tutorials/joins_tutorial.html
## one to many join
df1 <- tibble(x = c(1, 1, 2), y = 1:3)
df2 <- tibble(x = c(1, 1, 2), z = c("a", "b", "a"))

df1 %>% left_join(df2)

## semi_join: keep things that are the same in y as x
semi_join(x = publishers, y = superheroes)
semi_join(x = superheroes, y = publishers)

## anti_join: things that are not the same in y as x
anti_join(x = publishers, y = superheroes)
anti_join(x = superheroes, y = publishers)

## Joining multiple datasets
df1 <- data.frame(id=1:10, x=rnorm(10), y=runif(10))
df2 <- data.frame(id=1:11, z=rnorm(11), a=runif(11))
df3 <- data.frame(id=2:10, b=rnorm(9), c=runif(9))

mergedf = df1 %>% full_join(df2) %>% full_join(df3)
mergedf
# Challenge 2. Joining is just part of the art of 'data wrangling.' How do you subset, filter, create new columns?
# The packages dplyr and tidyr in the tidyverse are great for this.
# I really love the tutorial here: https://datacarpentry.org/R-ecology-lesson/03-dplyr.html

# You'll be using the following dataset to complete the tutorial
surveys <- read_csv("data/tutorials/surveys.csv")
view(surveys)
dim(surveys)
# Go through the tutorial and complete the challenges
# https://datacarpentry.org/R-ecology-lesson/03-dplyr.html
