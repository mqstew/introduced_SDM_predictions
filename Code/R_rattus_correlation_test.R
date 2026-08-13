## stats extraction script for Rattus rattus
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
r_rattus_threshold <- 484
r_rattus_points <- 16917

resp_name <- "rattusrattus"
species_name <- "rattus_rattus"

r_rat_target_countries <- c("Japan", "Madagascar", "Ecuador", 
                            "United States of America", "Australia",
                            "Portugal", "Chile", "New Zealand", "New Caledonia",
                            "Puerto Rico", "France", "Gabon", "Niger")
## mollweide crs
moll_crs <- "+proj=moll +lon_0=0 +units=m +datum=WGS84 +no_defs"
r_rat_abund <- read.csv(file.path("Data", "rattus_rattus_dataset.csv"))

## clim only or footprint
scenario_name <- "Footprint"
# scenario_name <- "ClimateOnly"

# us states loop for R. rattus ----
us_states <- c("Alabama","Alaska","Arizona","Arkansas","California","Colorado",
               "Connecticut","Delaware","Florida","Georgia","Hawaii","Idaho",
               "Illinois","Indiana","Iowa","Kansas","Kentucky","Louisiana",
               "Maine", "Maryland", "Massachusetts","Michigan","Minnesota",
               "Mississippi","Missouri", "Montana","Nebraska","Nevada",
               "New Hampshire","New Jersey","New Mexico","New York",
               "North Carolina","North Dakota","Ohio","Oklahoma","Oregon",
               "Pennsylvania","Rhode Island","South Carolina","South Dakota",
               "Tennessee","Texas","Utah","Vermont","Virginia","Washington",
               "West Virginia","Wisconsin","Wyoming")

states_rasters <- list()
for(s in seq_along(us_states)){
  # name (and safe version)
  state_name <- us_states[s]
  safe_name <- gsub("[^A-Za-z0-9]", "_", state_name)
  # locate folder
  ensemble_dir <- file.path(getwd(), resp_name, 
                            paste0("proj_", safe_name, "_", 
                                   scenario_name, "_ensemble"))
  emmean_files <- list.files(ensemble_dir, pattern = "\\.tif$",
                             recursive = TRUE, full.names = TRUE)
  if (length(emmean_files) == 0) { 
    message("No EMmean for ", state_name, " in ", ensemble_dir); next
  }
  r <- terra::rast(emmean_files[1])
  states_rasters[[state_name]] <- r
}
## merge us states
rlist <- states_rasters
rsc <- sprc(rlist)
r_rat_us <- merge(rsc)
## save
terra::writeRaster(r_rat_us, 
                   file.path(resp_name, paste0("proj_United_States_of_America_", 
                                               scenario_name, "_ensemble"),
                             paste0("proj_United_States_of_America_", scenario_name, "_ensemble_", 
                                    resp_name, "_ensemble.tif")))

# loop through country folders to extract projections ----
raster_list <- list()
for(i in seq_along(r_rat_target_countries)){
  # name (and safe version)
  ctry_name <- r_rat_target_countries[i]
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
# r_rat_list <- raster_list
r_rat_list <- sprc(raster_list)
r_rat_projection <- merge(r_rat_list)
terra::writeRaster(r_rat_projection,
                   file.path("Projections",
                             paste0(species_name, "_", scenario_name,"_projection.tif")),
                   overwrite=TRUE)

# threshold reproject to extract ----
## Mollweide CRS already defined
## reproject
r_rat_region <- vect(file.path("Vectors", "r_rattus_region.gpkg"))
## remove any countries from shapefile that are not in the target countries
r_rat_region <- r_rat_region[r_rat_region$name_en %in% r_rat_target_countries, ]
# plot(r_rat_region)
r_rat_moll <- project(r_rat_region, moll_crs)
r_rat_sf_moll <- st_as_sf(r_rat_moll)
## convert projection to m too
r_rat_projection <- project(r_rat_projection, moll_crs)
## appply threshold
r_rat_thresh <- terra::ifel(r_rat_projection > r_rattus_threshold, 
                            r_rat_projection, NA)

# Extract abundance stats and range
stats <- exact_extract(r_rat_thresh, r_rat_sf_moll, 
                       fun = function(values, coverage_fraction) {
                         ## country area
                         cell_area_m2 <- res(r_rat_thresh)[1] * res(r_rat_thresh)[2]
                         total_area_km2 <- sum(coverage_fraction) * cell_area_m2 / 1e6
                         
                         ## Only keep cells where species is present 
                         present <- which(values > 0 & !is.na(values))
                         if (length(present) == 0) {
                           return(data.frame(mean = NA, max = NA, range_km2 = 0,
                                             total_area_km2 = total_area_km2,
                                             prop_suitable = 0))
                         }
                         cell_area_m2 <- 
                           res(r_rat_thresh)[1] * res(r_rat_thresh)[2]
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
stats$region <- r_rat_sf_moll$name_en
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
## or area NAs
r_rat_abund <- r_rat_abund[r_rat_abund$sampling_country %in% 
                             r_rat_target_countries, ]
r_rat_abund <- r_rat_abund[!is.na(r_rat_abund$ind_per_km2), ]
r_rat_abund <- r_rat_abund[!is.na(r_rat_abund$area_km2), ]

## abundance data, reproject points to moll_crs
abundance_points <- vect(r_rat_abund, geom = 
                           c("sampling_longitude", "sampling_latitude"), 
                         crs = "EPSG:4326")
abundance_points <- project(abundance_points, moll_crs)

# buffers ----
buff_10 <- buffer(abundance_points, width = 10000)
buff_50 <- buffer(abundance_points, width = 50000)
buff_100 <- buffer(abundance_points, width = 100000)
## convert study area into a buffer radius
study_radius <- sqrt(r_rat_abund$area_km2 * 1e6 / pi) ## 1e6 converts km2 to m2
study_area <- buffer(abundance_points, width = study_radius)

## suitability for study location
study_suitability <- exact_extract(r_rat_thresh, st_as_sf(study_area), "mean")
## suitability at exact point
point_suitability <- terra::extract(r_rat_thresh, abundance_points)[, 2]
## suitability in buffer zones
buff_10_suit <- terra::extract(r_rat_thresh, buff_10, fun = mean, 
                               na.rm = TRUE)[,2]
buff_50_suit <- terra::extract(r_rat_thresh, buff_50, fun = mean, 
                               na.rm = TRUE)[,2]
buff_100_suit <- terra::extract(r_rat_thresh, buff_100, fun = mean, 
                                na.rm = TRUE)[,2]

## location and buffer suitability
points_table <- data.frame(
  country          = r_rat_abund$sampling_country,
  actual_abundance  = r_rat_abund$ind_per_km2,
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
