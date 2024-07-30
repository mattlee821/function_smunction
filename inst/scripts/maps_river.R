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
## river
river <- getbb(location)%>%
  opq()%>%
  add_osm_feature(key = "waterway", value = "river") %>%
  osmdata_sf()


# plot ====
## river zoomed ====
a <- ggplot() +
  geom_sf(data = river$osm_lines,
          color = "black",
          size = 1,
          alpha = 1) +
  coord_sf(ylim = c(getbb(location)[2,1], getbb(location)[2,2]), 
           xlim = c(getbb(location)[1,1], getbb(location)[1,2]),
           expand = FALSE) +
  theme_void()

ggsave(a, file="~/Desktop/lyon_river.eps", device="eps")


