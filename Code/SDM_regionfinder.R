## shapefile extraction code
## used when global projections cannot be done for a species
## identifies and saves a shapefile for the target countries for a focal species

# clean ----
rm(list=ls())

# packages ----
library(rnaturalearth)
library(sf)
library(dplyr)

# directory ----
directory <- "C:/Users/mstew/Documents/Modelling"

# Shapefile for Neogale vison ----
## Countries chosen from non-native abundances dataset
## N.vison <- Spain, UK, Poland, Finland, Argentina, Norway, Russia
n_vison_countries <- c("spain", "united kingdom", "poland", "finland", "argentina",
                       "norway", "russia")
world <- ne_countries(scale = "large", returnclass = "sf") ## 1km grain size
target <- tolower(n_vison_countries)
region <- world %>%
  dplyr::filter(
    tolower(name) %in% target |
      tolower(dplyr::coalesce(name_long, "")) %in% target |
      tolower(dplyr::coalesce(formal_en, "")) %in% target
  )
setwd(file.path(directory, "Vectors"))
st_write(region, "n_vison_region.gpkg", delete_dsn=TRUE)

# shapefile for R rattus ----
r_rattus_countries <- c("japan", "madagascar", "ecuador", "australia", "portugal", 
                        "chile", "new zealand", "new caledonia", "united states", 
                        "puerto rico", "sierra leone", "guinea", "france", "benin", 
                        "gabon", "niger", "iran")

us_states_1 <- c("Alabama","Alaska","Arizona","Arkansas","California","Colorado",
                 "Connecticut","Delaware","Florida","Georgia","Hawaii","Idaho",
                 "Illinois","Indiana","Iowa","Kansas","Kentucky","Louisiana")
## 18
us_states_2 <- c("Maine", "Maryland", "Massachusetts","Michigan","Minnesota","Mississippi",
                 "Missouri", "Montana","Nebraska","Nevada","New Hampshire",
                 "New Jersey","New Mexico","New York","North Carolina","North Dakota")
## 16
us_states_3 <- c("Ohio","Oklahoma","Oregon","Pennsylvania","Rhode Island",
                 "South Carolina","South Dakota","Tennessee","Texas",
                 "Utah","Vermont","Virginia","Washington","West Virginia",
                 "Wisconsin","Wyoming")
## 16
## us states
us_states <- ne_states("United States of America", returnclass = "sf")
states_1 <- us_states[us_states$name %in% us_states_1,]
states_2 <- us_states[us_states$name %in% us_states_2,]
states_3 <- us_states[us_states$name %in% us_states_3,]

## save
setwd(file.path(directory, "Vectors"))
st_write(states_1, "us_states_1.gpkg", delete_dsn=TRUE)
st_write(states_2, "us_states_2.gpkg", delete_dsn=TRUE)
st_write(states_3, "us_states_3.gpkg", delete_dsn=TRUE)

# countries from the world map using target_countries ----
world <- ne_countries(scale = "large", returnclass = "sf") ## 1km grain size
target <- tolower(target_countries)
region <- world %>%
  dplyr::filter(
    tolower(name) %in% target |
      tolower(dplyr::coalesce(name_long, "")) %in% target |
      tolower(dplyr::coalesce(formal_en, "")) %in% target
  )
setwd(file.path(directory, "Vectors"))
st_write(region, ".gpkg", delete_dsn=TRUE)

# europe, remove large/unconnected landmasses ----
## continents
europe <- ne_countries(scale="large", continent = "europe", returnclass = "sf")
s_am <- ne_countries(scale="large", continent = "south america", returnclass = "sf")
n_am <- ne_countries(scale="large", continent = "north america", returnclass = "sf")
oce <- ne_countries(scale="large", continent = "oceania", returnclass = "sf")
asia  <- ne_countries(scale="large", continent = "asia", returnclass = "sf")
africa <- ne_countries(scale="large", continent = "africa", returnclass = "sf")
  
europe <- europe %>%
  filter(!name %in% c("Russia", "Iceland", "Greece"))

?st_write
setwd(file.path(directory, "Vectors"))
st_write(europe, "europe_red.gpkg", delete_dsn=TRUE)
st_write(s_am, "s_am.gpkg", delete_dsn=TRUE)
st_write(n_am, "n_am.gpkg", delete_dsn=TRUE)
st_write(asia, "asia.gpkg", delete_dsn=TRUE)
st_write(africa, "africa.gpkg", delete_dsn=TRUE)

## leave shapefile directory
setwd(directory)
