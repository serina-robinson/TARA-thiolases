# Load packages
pacman::p_load("tidyverse", "leaflet")

# NOTE: Before you start it will be helpful to read the leaflet documentation
# https://rstudio.github.io/leaflet/
  
# Set working directory
setwd("C:/Users/tabie/OneDrive/Documents/GitHub/TARA-thiolases/")

# Read in environmental dataset for metagenome-assembled genomes with oleACD clusters in TARA
dat <- read_csv("data/TARA_tree_annotation_df.csv")

# Add the 'content' for pop-up text
content <- paste("Genus:", dat$genus, "<br>",
                 "Family:", dat$family, "<br>",
                 "Temperature:", round(dat$temperature, 1), "<br>",
                 "Depth:", dat$depth_m, "meters")

# Set the custom color palette
pal <- colorNumeric(
  palette = "RdYlBu", # name of a color palette 
  # full list of color palettes available here
  # https://www.nceas.ucsb.edu/sites/default/files/2020-04/colorPaletteCheatsheet.pdf
  na.color = "gray70",
  reverse = TRUE,
  domain = dat$temperature)

# Make an interactive map!
leaflet(data = dat) %>%
  addProviderTiles(providers$CartoDB.PositronNoLabels) %>% 
  #  addProviderTiles(providers$Esri.WorldPhysical) %>% # example of changing the map style
  # you can find the names of different map tiles here: 
  # http://leaflet-extras.github.io/leaflet-providers/preview/
  addCircleMarkers(lng = ~lon, 
                   lat = ~lat, 
                   popup = content,
                   radius = 4,      
                   fillOpacity = 0.5, 
                   stroke = FALSE,
                   color = ~pal(temperature))

# Challenge 1. Change the map style to a physical world map rather than just country outlines

# Challenge 2. Change the point coloring to be a different variable than temperature 
# e.g. depth, oxygen, iron etc.

# Challenge 3. Change the radius of the points to visualize a different variable
# Note, you may run into issues just using the variable directly. 
# Log transformations will be your friend.

# Challenge 4. Change the popup to display not only the genus and family but also the order, class and phylum.
 
# Challenge 5. You will notice some of the points overlap so you can't see them all.
# What happens if you use lat_jitter and lon_jitter (existing columns) for plotting instead?
# I've created these columns using the following function: 
# https://rdrr.io/github/lmullen/mullenMisc/man/jitter_latlong.html

# Challenge 6. The current map is displaying all 244 samples. 
# Try subsetting the data to only display the 30 points for the sequences we are actually getting!

# Challenge 7. Open to showing your map in lab meeting?

# Super challenge 8: Can you also plot the locations of your thermophiles and psychrophiles? 
# (This is part of the long-term goal for the week so no need to accomplish today)
