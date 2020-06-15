# Install packages
pacman::p_load("tidyverse")

# Set working directory
setwd("C:/Users/tabie/OneDrive/Documents/GitHub/TARA-thiolases")

# Read in your full50 dataset
full50 <- read_csv("data/full50_raw.csv")
# Credit for these functions goes to:
# https://github.com/lmullen/mullenMisc

#' Jitter latitude and longitude coordinates
#'
#' Given a latitude or longitude coordinate, return that coordinate jittered
#' plus or minus a certain maximum amount. This function assumes a spherical
#' earth when calculating the maximum amount to jitter.
#'
#' @param coord The coordinate in degrees
#' @param type Whether the coordinate is latitude or longitude
#' @param latitude If the coordinate is longitude, then the latitude to use for
#'   calculating the maximum amount to jitter.
#' @param km The maximum number of kilometers to jitter a point plus or minus.
#' @return The jittered coordinate in degrees.
#' @examples
#' jitter_latlong(-73, type = "long", lat = 43, km = 1)
#' jitter_latlong(42, type = "lat", km = 1)
#' @export
#' @seealso \code{\link{length_of_degree}}
jitter_latlong <- function(coord, type = c("lat", "long"), latitude, km = 1) { # making a function of coord, type (default c("lat", "long")), latitude, and km (default 1)
  type = match.arg(type) # allowing inputs for lat and long to be sorted to the correct one even if they don't input exactly "lat" or "long"?
  
  if(missing(latitude) & type == "lat") { 
    latitude <- coord } # if there isn't a value given for latitude, but coord is a latitude, then use that

  km_per_degree <- length_of_degree(latitude, type = type) # finding the length of a degree from the function explained below
  degree_per_km <- 1 / km_per_degree # blah
  coord + (runif(1, min = -1, max = 1) * degree_per_km * km)
} #^ taking the coord, adding or subtracting a randomized amount, and then scaling it by the length a degree should be and how many km we want to allow it to move

# What is going on here???
jitter_latlong <- Vectorize(jitter_latlong,
                           vectorize.args = c("coord", "latitude"))
#^ turning the coord and latitude into a vector? I don't know--if I comment this out nothing seems different


#' Length of a degree of latitude or longitude
#'
#' Calculates the length of a degree of latitude or longitude in kilometers,
#' assuming an spherical earth.
#'
#' @param degree The degree to calculate the length for
#' @param type Whether to return the length of a degree for latitude or longitude
#' @return Length of the degree in kilometers
#'
#' @export
# Helper function
length_of_degree <- function(degree, type = c("lat", "long")) { # making a function that depends on degree and type 
  type <- match.arg(type) #making sure your input gets correctly sorted at lat or long again?
  length_at_equator <- 110.5742727 # in kilometers /// specifying how long 1 degree at the equator is
  if (type == "long") { # if degree is a longitude then we need to convert it
    cos(degree * (2 * pi) / 360) * length_at_equator # convert the degrees to radians and take the cos to do something
  } else if (type == "lat") { # if degree is a latitude then leave it be
    length_at_equator
  }
}

# This is a little different challenge. Instead of writing your own code you're 
# going to first practice reading and understanding someone else's code!

# Challenge 1. Walk through the functions and look up any commands you don't know.
# For example, what is match.arg? runif? Vectorize? 
?match.arg
  # https://alistaire.rbind.io/blog/match.arg/
  # checks the input against a set of choices, can partial match
?runif
  # seems like it will generate random values in a distribution specified
?Vectorize
  # vectorizes the action of its argument
  # maybe takes all loose values of coord and lat and combines them into a respective vector?
# Challenge 2. Write a comment on each line to describe what is happening. 
# Basically pretend you want to explain this to someone who has never seen it before
# (You can actually explain this to me on Monday).
# You might find the following and other resources on the spherical law of cosines
# to be helpful:
# http://www.math.ucdenver.edu/~hartkes/teaching/2011m896/SphericalLawOfCosines.pdf
# https://www.theoremoftheday.org/GeometryAndTrigonometry/SphericalCos/TotDSphericalCos.pdf
 
# Challenge 3. Run the functions. You should see them appear in your environment. 
# Now use jitter_latlong to 'jitter' both the lattitude and longitude 
# For your psychrophiles and thermophiles by 50 km.
# You will need to use mutate to create two new columns: lat_jitter and lon_jitter
# This documentation will be helpful:
# https://rdrr.io/github/lmullen/mullenMisc/man/jitter_latlong.html

full50_1 <- full50 %>% 
  mutate(lon_jitter = jitter_latlong(full50$lon, type = "long", latitude = full50$lat, km = 50),
         lat_jitter = jitter_latlong(full50$lat, type = "lat", km = 50))
# Why does it jitter them all the same?

select(full50_1, lon_jitter, lat_jitter) %>% 
  view()

write_csv(full50_1, "data/full50_1.csv")

### to do:
# challenges 4 and 5
# fix the mutate thing ^
# read article and take notes on figures
# finish comments on rows for helper function
# remove NAs dynamically

# Challenge 4. Now plot your jittered lat and long coordinates in leaflet! 
# Add it your flexdashboard

# Super challenge 5. This is really tricky...ok if you don't do it. Prioritize journal club.
# But what if you wanted to have the option to either show either jittered 
# coordinates or not jittered coordinates on the map. How would you implement this?

