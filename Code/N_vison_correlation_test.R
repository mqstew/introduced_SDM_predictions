## stats extraction script for Neogale vison
## With thanks to Felipe Espinoza who's work provides a basis for this code
## and for his assistance with errors and providing advice.

## this code takes the ensemble projection(s) and extracts habitat suitability data
## it then produces a comparison table with the suitability values for each study
## alongside the relative local abundances
## this comparison table is then passed to both Spearman's Rank and Pearson 
## correlation tests, which are also saved

# clean ----
rm(list=ls())
gc()

# packages ----
library(biomod2)
library(terra)
library(tidyterra)
library(exactextractr)
library(sf)

# directory ----
directory <- "C:/Users/mstew/Documents/Modelling"
setwd(directory)

# variables ----
n_vison_threshold <- 718
n_vison_points <- 21175

resp_name <- "neogalevison"
species_name <- "neogale_vison"
n_vis_target_countries <- c("United Kingdom", "Spain", "Poland",
                            "Finland", "Norway", "Argentina")

## mollweide crs
moll_crs <- "+proj=moll +lon_0=0 +units=m +datum=WGS84 +no_defs"
n_vis_abund <- read.csv(file.path("Data", "neogale_vison_dataset.csv"))
## clim only or footprint
scenario_name <- "Footprint"
# scenario_name <- "ClimateOnly"

# loop through country folders to extract projections ----
raster_list <- list()
for(i in seq_along(n_vis_target_countries)){
  # name (and safe version)
  ctry_name <- n_vis_target_countries[i]
  safe_name <- gsub("[^A-Za-z0-9]", "_", ctry_name)
  # locate folder
  ensemble_dir <- file.path(getwd(), resp_name, 
                            paste0("proj_", safe_name, "_", 
                                   scenario_name, "_ensemble"))
  emmean_files <- list.files(ensemble_dir, pattern = "\\.tif$",
                            recursive = TRUE, full.names = TRUE)
  if (length(emmean_files) == 0) { 
    message("No EMmean for ", ctry_name, " in ", ensemble_dir); next
  }
  r <- terra::rast(emmean_files[1])
  raster_list[[ctry_name]] <- r
}

# merge us to rest of projections for global projection for r rattus ----
# n_vis_list <- raster_list
n_vis_list <- sprc(raster_list)
n_vis_projection <- merge(n_vis_list)
terra::writeRaster(n_vis_projection,
                   file.path("Projections",
                             paste0(species_name, "_", scenario_name,"_projection.tif")),
                   overwrite=TRUE)

# merge Europe to Argentina for n_vis ----
europe_dir <- file.path(getwd(), resp_name, 
                        paste0("proj_Europe_", scenario_name, "_ensemble"))
europe_files <- list.files(europe_dir, pattern = "\\.tif$", recursive = TRUE,
                           full.names = TRUE)
europe <- terra::rast(europe_files[1])

## check crs, res and extents, then merge into a sinle raster layer ----
# terra::crs(europe, describe = TRUE)
# terra::crs(raster_list[["Argentina"]], describe = TRUE)
# 
# terra::res(europe);   terra::res(raster_list[["Argentina"]])      # resolution
# terra::ext(europe);   terra::ext(raster_list[["Argentina"]])      # extent

# merge and save ----
n_vis_projection <- terra::merge(europe, raster_list[["Argentina"]])
terra::writeRaster(n_vis_projection,
            file.path("Projections",
                      paste0(species_name, "_", scenario_name,"_projection.tif")),
            overwrite=TRUE)

# threshold reproject to extract ----
## Mollweide CRS already defined
## reproject
n_vis_region <- vect(file.path("Vectors", "n_vison_region.gpkg"))
## remove Russia from shapefile before reprojecting
n_vis_region <- n_vis_region[n_vis_region$name_en %in% n_vis_target_countries, ]
n_vis_moll <- project(n_vis_region, moll_crs)
n_vis_sf_moll <- st_as_sf(n_vis_moll)

## convert the projection to metres too
n_vis_projection <- project(n_vis_projection, moll_crs)
## Apply threshold to raster
n_vis_thresh <- terra::ifel(n_vis_projection > n_vison_threshold, 
                            n_vis_projection, NA)

# Extract abundance stats and range
stats <- exact_extract(n_vis_thresh, n_vis_sf_moll, 
                       fun = function(values, coverage_fraction) {
                         ## country area
                         cell_area_m2 <- res(n_vis_thresh)[1] * res(n_vis_thresh)[2]
                         total_area_km2 <- sum(coverage_fraction) * cell_area_m2 / 1e6
                         
                         ## Only keep cells where species is present 
                         present <- which(values > 0 & !is.na(values))
                         if (length(present) == 0) {
                           return(data.frame(mean = NA, max = NA, range_km2 = 0,
                                             total_area_km2 = total_area_km2,
                                             prop_suitable = 0))
                         }
                         cell_area_m2 <- 
                           res(n_vis_thresh)[1] * res(n_vis_thresh)[2]
                         range_km2 <- 
                           sum(coverage_fraction[present]) * cell_area_m2 / 1e6
                         
                         data.frame(
                           mean = mean(values[present]),
                           max = max(values[present]),
                           suitable_area_km2 = range_km2,
                           total_area_km2 = total_area_km2,
                           prop_suitable = range_km2 / total_area_km2
                         )
                       })
stats$region <- n_vis_sf_moll$name_en
stats

# tables (to then merge) ----
## country data
country_table <- data.frame(
  country = stats$region,
  country_suitability = stats$mean / 1000,
  coutry_area = stats$total_area_km2,
  suitable_area_km2 = stats$suitable_area_km2,
  suitable_prop = stats$prop_suitable
)

## remove russian points from dataset and any points without N/km2 measures
n_vis_abund <- n_vis_abund[n_vis_abund$sampling_country %in% 
                                       n_vis_target_countries, ]
n_vis_abund <- n_vis_abund[!is.na(n_vis_abund$ind_per_km2), ]

## abundance data, reproject points to moll_crs
abundance_points <- vect(n_vis_abund, geom = 
                           c("sampling_longitude", "sampling_latitude"), 
                         crs = "EPSG:4326")
abundance_points <- project(abundance_points, moll_crs)

## buffers
buff_10 <- buffer(abundance_points, width = 10000)
buff_50 <- buffer(abundance_points, width = 50000)
buff_100 <- buffer(abundance_points, width = 100000)
## convert study area into a buffer radius
study_radius <- sqrt(n_vis_abund$area_km2 * 1e6 / pi) ## 1e6 converts km to meters
study_area <- buffer(abundance_points, width = study_radius)

## suitability for study location
study_suitability <- exact_extract(n_vis_thresh, st_as_sf(study_area), "mean")
## suitability at exact point
point_suitability <- terra::extract(n_vis_thresh, abundance_points)[, 2]
## suitability in buffer zones
buff_10_suit <- terra::extract(n_vis_thresh, buff_10, fun = mean, 
                               na.rm = TRUE)[,2]
buff_50_suit <- terra::extract(n_vis_thresh, buff_50, fun = mean, 
                               na.rm = TRUE)[,2]
buff_100_suit <- terra::extract(n_vis_thresh, buff_100, fun = mean, 
                                na.rm = TRUE)[,2]

## location and buffer suitability
points_table <- data.frame(
  country          = n_vis_abund$sampling_country,
  actual_abundance  = n_vis_abund$ind_per_km2,
  point_suitability = point_suitability / 1000,
  study_area_suitability  = study_suitability / 1000,
  
  buffer_10k = buff_10_suit / 1000,
  buffer_50k = buff_50_suit / 1000,
  buffer_100k = buff_100_suit / 1000
)
# --- add country average to each point ---
comparison_table <- merge(points_table, country_table, by = "country", all.x = TRUE)

write.csv(comparison_table,
          file.path(directory, "Outputs",
                    paste0(species_name, "_", scenario_name, "_comparison_table.csv")))

# running correlation across each column ----
## spearman
scale_columns <- c("point_suitability", "study_area_suitability", 
                   "country_suitability", "buffer_10k", "buffer_50k", "buffer_100k")
spearman_results <- lapply(scale_columns, function(column){
  corr_test <- cor.test(comparison_table[[column]], 
                        comparison_table$actual_abundance,
                        method = "spearman")
  data.frame(scale = column, rho = corr_test$estimate, p = corr_test$p.value)
})
do.call(rbind, spearman_results)
spear_res <- do.call(rbind, spearman_results)

## pearson
pearson_results <- lapply(scale_columns, function(column){
  corr_test <- cor.test(comparison_table[[column]], 
                        comparison_table$actual_abundance,
                        method = "pearson")
  data.frame(scale = column, rho = corr_test$estimate, p = corr_test$p.value)
})
do.call(rbind, pearson_results)
pears_res <- do.call(rbind, pearson_results)
## merge and save
res_table <- merge(spear_res, pears_res, by = "scale", all.x = TRUE)
write.csv(res_table,
          file.path("Outputs", paste0(species_name, "_",
                                      scenario_name, "_correlation_results.csv")))
