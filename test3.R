library(lidR)
library(happign)
library(sf)
library(tmap); tmap_mode("view") # Set map to interactive
library(dplyr)
library(ggplot2);sf_use_s2(FALSE) # Avoid problem with spherical geometry
library(purrr)
library(stars)
library(terra)

borders <- read_sf("C:/Users/mathi/Desktop/projet troyes/cadastre_foret.gpkg")
layers <- get_layers_metadata("wms-r", "altimetrie")
mnt_layer <- layers[3,1] # "ELEVATION.ELEVATIONGRIDCOVERAGE.HIGHRES"
mns_layer <- layers[4,1] # "ELEVATION.ELEVATIONGRIDCOVERAGE.HIGHRES.MNS"

mnt <- get_wms_raster(borders, mnt_layer, res = 5, crs = 2154, rgb = FALSE)
mns <- get_wms_raster(borders, mns_layer, res = 5, crs = 2154, rgb = FALSE)

level_curve <- get_wfs(borders, "ELEVATION.CONTOUR.LINE:courbe",
                       spatial_filter = "intersects") |> 
  st_intersection(borders)

# Calculate digital height model i.e. tree height
mnh <- mns - mnt
mnh[mnh < 0] <- NA  # Remove negative value 
mnh[mnh > 50] <- 40 # Remove height more than 50m
plot(mnh)
plot(mnt)
plot(mns)

las <- readLAS("C:/Users/mathi/Downloads/LHD_FXX_1003_6798_PTS_C_LAMB93_IGN69.copc.laz")
plot(las)
clip_bb = st_bbox(las)-500
clip_las = clip_roi(las, clip_bb)
plot(clip_las, color="Classification")
# Khosravipour et al. pitfree algorithm
thr <- c(0,2,5,10,2000)
edg <- c(0, 1.5)
chm <- rasterize_canopy(clip_las, 1, pitfree(thr, edg))
sol <- rasterize_terrain(clip_las, 1, tin())
plot(chm)
plot(sol)
plot(chm-sol)
