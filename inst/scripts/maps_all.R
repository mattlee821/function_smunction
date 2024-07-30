# http://joshuamccrain.com/tutorials/maps/streets_tutorial.html

library(remotes)
remotes::install_github("ropensci/osmdata")
library(tidyverse)
library(osmdata) # package for working with streets
library(showtext) # for custom fonts
library(ggmap)
library(rvest)
library(sf)

#
location <- "lyon france"
available_tags("highway")
getbb(location)

# grab map data ====
## big streets
big_streets <- getbb(location)%>%
  opq()%>%
  add_osm_feature(key = "highway", 
                  value = c("motorway", "primary", "motorway_link", "primary_link")) %>%
  osmdata_sf()

## medium streets
med_streets <- getbb(location)%>%
  opq()%>%
  add_osm_feature(key = "highway", 
                  value = c("secondary", "tertiary", "secondary_link", "tertiary_link")) %>%
  osmdata_sf()

# small streets
small_streets <- getbb(location)%>%
  opq()%>%
  add_osm_feature(key = "highway", 
                  value = c("residential", "living_street",
                            "unclassified",
                            "service", "footway"
                  )) %>%
  osmdata_sf()

## river
river <- getbb(location)%>%
  opq()%>%
  add_osm_feature(key = "waterway", value = "river") %>%
  osmdata_sf()

## railway
railway <- getbb(location)%>%
  opq()%>%
  add_osm_feature(key = "railway", value="rail") %>%
  osmdata_sf()

# plot ====
## big streets
ggplot() +
  geom_sf(data = big_streets$osm_lines,
          inherit.aes = FALSE,
          color = "black")

## river
ggplot() +
  geom_sf(data = river$osm_lines,
          inherit.aes = FALSE,
          color = "black")

## big streets and river zoomed ====
ggplot() +
  geom_sf(data = big_streets$osm_lines,
          inherit.aes = FALSE,
          color = "black") +
  geom_sf(data = river$osm_lines,
          inherit.aes = FALSE,
          color = "black") +
  coord_sf(ylim = c(getbb(location)[2,1], getbb(location)[2,2]), 
           xlim = c(getbb(location)[1,1], getbb(location)[1,2]),
           expand = FALSE)

## everything ====
ggplot() +
  geom_sf(data = river$osm_lines,
          inherit.aes = FALSE,
          color = "steelblue",
          size = .8,
          alpha = .3) +
  geom_sf(data = railway$osm_lines,
          inherit.aes = FALSE,
          color = "black",
          size = .2,
          linetype="dotdash",
          alpha = .5) +
  geom_sf(data = med_streets$osm_lines,
          inherit.aes = FALSE,
          color = "black",
          size = .3,
          alpha = .5) +
  geom_sf(data = small_streets$osm_lines,
          inherit.aes = FALSE,
          color = "#666666",
          size = .2,
          alpha = .3) +
  geom_sf(data = big_streets$osm_lines,
          inherit.aes = FALSE,
          color = "black",
          size = .5,
          alpha = .6) +
  coord_sf(ylim = c(getbb(location)[2,1], getbb(location)[2,2]), 
           xlim = c(getbb(location)[1,1], getbb(location)[1,2]),
           expand = FALSE)
 

