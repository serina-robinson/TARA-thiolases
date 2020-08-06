# Load packages
pacman::p_load("tidyverse", "leaflet", "Biostrings")

# NOTE: Before you start it will be helpful to read the leaflet documentation
# https://rstudio.github.io/leaflet/
  
# Set working directory
setwd("C:/Users/tabie/OneDrive/Documents/GitHub/TARA-thiolases/")

# Read in environmental dataset for metagenome-assembled genomes with oleACD clusters in TARA
dat <- read_csv("data/TARA_tree_annotation_df.csv")
dat
tara <- readAAStringSet('data/TARA_marine_sequences_to_order.fasta') # these are the 30 sequences from TARA oceans
names(tara)
view(dat)
dat[218,1] <- gsub("-","_", dat[218,1])
dat2 <- dat[dat$label %in% names(tara),]
dim(dat2)

PTM <- read_csv("data/Psychrophiles_Thermophiles_Mesophiles.csv")
dim(PTM)
view(PTM)
full50 <- full_join(dat2, PTM)
dim(full50)
view(full50)

write_csv(full50, "data/full50_raw")

# select(dat, -lat, -lon)
# Add the 'content' for pop-up text
content <- paste("Genus:", full50$genus, "<br>",
                 "Strain:", full50$genome.x, "<br>",
                 "Temperature:", round(full50$temperature, 1), "<br>",
                 "Depth:", full50$depth_m, "meters")
class(content)
# Set the custom color palette
pal <- colorNumeric(
  palette = "RdYlBu", # name of a color palette 
  # full list of color palettes available here
  # https://www.nceas.ucsb.edu/sites/default/files/2020-04/colorPaletteCheatsheet.pdf
  na.color = "gray70",
  reverse = TRUE,
  domain = dat$temperature)
 # rad <- ifelse(is.na(full50[,"depth_m"]), 2, log(full50$depth_m))
rad
new <- full50 %>% 
  mutate(rad = case_when((log(depth_m) < 3) ~ 3,
                         !is.na(depth_m) ~ log(depth_m),
                         TRUE ~ 3))
new$rad
view(new)
?ifelse
?case_when
# Make an interactive map!
leaflet(data = new) %>%
  addProviderTiles(providers$Stamen.Terrain) %>% 
  #  addProviderTiles(providers$OpenStreetMap.Mapnik) %>% # example of changing the map style
  # you can find the names of different map tiles here: 
  # http://leaflet-extras.github.io/leaflet-providers/preview/
  addCircleMarkers(lng = ~lon, 
                   lat = ~lat, 
                   popup = content,
                   radius = new$rad,      
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

jitter_latlong(coord, type = c("lat", "long"), latitude, km = 50)
?jitter_latlong
