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

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
dir <- getwd()

lien <- fromJSON(txt = "https://data.geopf.fr/private/wfs/?service=WFS&version=2.0.0&apikey=interface_catalogue&request=GetFeature&typeNames=IGNF_LIDAR-HD_TA:nuage-dalle&outputFormat=application/json&bbox=766276,6353692.630280115,769016.9951275728,6355112.417896504")
lien <- lien[["features"]][["properties"]][["url"]][1]

draw.area = function(){
  zone = st_coordinates(st_transform(mapedit::drawFeatures(), crs = 2154))
  y1 = min(zone[, 2])
  x1 = min(zone[, 1])
  y2 = max(zone[, 2])
  x2 = max(zone[, 1])
  return(c(x1,y1,x2,y2))
}


download.lidar = function(x1, y1, x2, y2) {
  # Créer une séquence de coordonnées avec un pas de 10000
  x = seq(x1, x2, 10000)
  y = seq(y1, y2, 10000)
  
  # Boucle à travers toutes les combinaisons de x et y
  for (i in x) {
    for (j in y) {
      # Construire l'URL pour obtenir les données JSON
      json_url <- paste0("https://data.geopf.fr/private/wfs/?service=WFS&version=2.0.0&apikey=interface_catalogue&request=GetFeature&typeNames=IGNF_LIDAR-HD_TA:nuage-dalle&outputFormat=application/json&bbox=", 
                         i, ",", j, ",", i, ",", j)
      print(paste("Fetching JSON from:", json_url))
      
      # Essayer de récupérer les données JSON
      tryCatch({
        json = fromJSON(txt = json_url)
        
        # Récupérer le lien du fichier .laz à partir des propriétés du JSON
        lien = json[["features"]][["properties"]][["url"]][1]
        
        if (!is.null(lien)) {  # Vérifier que le lien n'est pas nul
          print(paste("Downloading .laz file from:", lien))
          
          # Téléchargement du fichier .laz avec mode binaire
          download.file(lien, 
                        destfile = paste0("dir", i, "_", j, ".laz"),
                        mode = "wb")  # "wb" pour mode binaire
          
          print(paste("Saved file:", paste0(i, "_", j, ".laz")))
        } else {
          print("No valid link found in the JSON response.")
        }
      }, error = function(e) {
        print(paste("Error fetching or downloading data for bbox:", i, j))
        print(e)  # Affiche l'erreur rencontrée
      })
    }
  }
}

cut.area = function(las, liste){
  clip = clip_rectangle(las, liste[1], liste[2], liste[3], liste[4])
  return(clip)
}


coord <- draw.area()
download.lidar(coord[1],coord[2],coord[3],coord[4])

laz_dir <- list.files(dir, 
                  full.names = T, 
                  pattern= '.laz')

laz <- readLAS(laz_dir)
laz <- cut.area(laz, coord)
plot(laz)


las <- laz

laz_soil <- lidR::filter_ground(laz)
mycsf <- csf(TRUE, 1, 1, time_step = 1)
laz_soil <- classify_ground(laz, mycsf)
laz_soil <-lidR::filter_poi(laz_soil, Classification == 2L & ReturnNumber > 5L)

plot(laz_soil)

test = segment_shapes(LAS_soil, shp_hline(th1 = 100, th2 = 2, k = 3), attribute = "Shape")

plot(test, color = "Shape")

ground <- classify_ground(laz, csf())
plot(ground, color = "Classification")
chm <- grid_canopy(laz, res = 1, pitfree())
plot(chm)
density <- grid_density(laz, res = 1)
plot (density)
################################################################################
#comparaison d'algo : aucune différence 

mnt_tin <- rasterize_terrain(las, 1, tin())
plot(mnt_tin)
writeRaster(mnt_tin,"mnt_tin.tif")

mnt_knnidw <- rasterize_terrain(las, 1, knnidw())
plot(mnt_knnidw)
writeRaster(mnt_knnidw,"mnt_knnidw.tif")

mnt_kriging <- rasterize_terrain(las, 1, kriging())
plot(mnt_kriging)
writeRaster(mnt_kriging,"mnt_kriging.tif")

test_1 <- mnt_tin-mnt_knnidw
plot (test_1)

test_2 <- mnt_tin-mnt_kriging
plot (test_2)
##################################################################################
las_norm <- normalize_height(las, mnt)
las_filtered <- filter_poi(las_norm, Z >= 0 & Z <= 0.2)
mnt_filtered <- rasterize_terrain(las_filtered, 1, tin())
plot(las_filtered)
writeRaster(mnt_filtered, "las_filtered.tif")
mnt_filtered2 <- rasterize_canopy(las, 1, pitfree(thr, edg))
writeRaster(mnt_filtered2,"las_filtered2.tif")


las_file <- laz
las <- laz
detect_cloisonnements <- function(las, resolution = 1, threshold = 0.1, output_file = "cloisonnements.gpkg") {
  
  if (is.null(las) || npoints(las) == 0) {
    stop("Le fichier LAS est vide ou n'a pas été chargé correctement.")
  }
  
  # 1. Calculer la densité des points
  density <- grid_density(las, res = resolution)
  
  # 2. Calculer le nombre de retours
  returns <- grid_metrics(las, ~length(Z), res = resolution)
  names(returns) <- "num_returns"
  
  #3. MNH
  chm <- grid_canopy(las, res = resolution, algorithm = pitfree())
  names(chm) <- "chm_height"
  
  # 4. Empiler les rasters pour les combiner
  combined_stack <- raster::stack(density, returns, chm)
  
  # Convertir en data frame pour traitement
  combined_df <- as.data.frame(combined_stack, xy = TRUE)
  
  # Vérifier la présence de valeurs NA et appliquer na.rm = TRUE
  combined_df <- na.omit(combined_df)  # Suppression des lignes avec des NA
  
  # 5. Détecter les zones à faible densité, faible nombre de retours et faible CHM
  print(combined_df)
  combined_df$low_density <- ifelse(combined_df$density < quantile(combined_df$density, threshold, na.rm = TRUE), 1, 0)
  combined_df$low_returns <- ifelse(combined_df$num_returns < quantile(combined_df$num_returns, threshold, na.rm = TRUE)-15, 1, 0)
  combined_df$low_chm <- ifelse(combined_df$chm_height < quantile(combined_df$chm_height, threshold, na.rm = TRUE)-6, 1, 0)
  
  # 6. Détecter les cloisonnements
  combined_df$cloisonnement <- ifelse(combined_df$low_density == 1 & combined_df$low_returns == 1 & combined_df$low_chm == 1, 1, 0)
  
  # 7. Filtrer les zones de cloisonnements
  cloisonnements_points <- combined_df[combined_df$cloisonnement == 1, ]
  print("i")
  # 8. Convertir en objet spatial
  if (nrow(cloisonnements_points) > 0) {
    cloisonnements_sf <- st_as_sf(cloisonnements_points, coords = c("x", "y"), crs = st_crs(las))
    
    # 9. Exporter au format GPKG
    st_write(cloisonnements_sf, output_file, driver = "GPKG")
    message("Cloisonnements exportés vers ", output_file)
  } else {
    message("Aucun cloisonnement détecté.")
  }
  
  return(cloisonnements_sf)
}


# Exemple d'utilisation
cloisonnements <- detect_cloisonnements(las)





############################# Gab ###########################
dtm_tin <-rasterize_terrain(laz, res = 1, algorithm = tin())
plot_dtm3d(dtm_tin, bg="black")

###################### tri points entre 0 20 cm ################

writeRaster(dtm_tin,"blackandwhite_sol.tif")
mnh <- rasterize_terrain(las, 1, tin())


carte <- as.numeric(laz)
triangulation <- deldir(laz)
 x <- laz$X
 y <- laz$Y
 z <- laz$Z
 
tritest <-deldir(laz$X, laz$Y)
plot(tritest)





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
las <- laz

las <- readLAS("C:/Users/mathi/Downloads/LHD_FXX_0638_6971_PTS_C_LAMB93_IGN69.copc.laz")
#plot(las)
clip_bb = st_bbox(las)
clip_las = clip_roi(las, clip_bb)
#plot(clip_las, color="Classification")
# Khosravipour et al. pitfree algorithm
thr <- c(0,2,5,10,2000)
edg <- c(0, 1.5)
chm <- rasterize_canopy(las, 1, pitfree(thr, edg))
sol <- rasterize_terrain(las, 1, tin())
plot(chm)
plot_dtm3d(sol)
test_mnh <- chm-sol
writeRaster(test_mnh, "mnh.tif")
