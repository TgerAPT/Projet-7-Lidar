library(lidR)
library(happign)
library(sf)
library(tmap); tmap_mode("view") # Set map to interactive
library(dplyr)
library(ggplot2);sf_use_s2(FALSE) # Avoid problem with spherical geometry
library(purrr)
library(stars)
library(terra)
library(jsonlite)

draw.area = function(){
  zone = st_coordinates(st_transform(mapedit::drawFeatures(), crs = 2154))
  y1 = min(zone[, 2])
  x1 = min(zone[, 1])
  y2 = max(zone[, 2])
  x2 = max(zone[, 1])
  return(c(x1,y1,x2,y2))
}
