## sdm by country or region
## individual country projections

## With thanks to Felipe Espinoza who's work provides a basis for this code
## and for his assistance with errors and providing advice.

## project to target regions separately as my system cannot handle a global projection

# clean ----
rm(list=ls())
gc()

# packages ----
library(terra)
library(usdm)
library(biomod2)
library(readr)
library(dplyr)
library(stringr)
library(rnaturalearth)
library(sf)

# ---- Resolution control ----
# 30s native ~900M cells. Pick your test resolution:
#   agg_factor = 5 -> 2.5 arc-min (~37M)  [safe/fast]
#   agg_factor = 4 -> 2   arc-min (~56M)
#   agg_factor = 3 -> 1.5 arc-min (~100M) [finer, heavier]
#   agg_factor = 2 -> 1   arc-min (~225M) [heavy for a laptop]
#   agg_factor = 1 -> native 30s (~900M)  [fine for a few countries]
# NOTE: now that we clip to a few countries, agg_factor = 1 is usually fine.
agg_factor <- 1

# directories ----
directory <- "C:/Users/mstew/Documents/introduced_SDM_predictions"
data_dir <- "C:/Users/mstew/Documents/introduced_SDM_predictions/Data"
out_dir <- "C:/Users/mstew/Documents/introduced_SDM_predictions/Outputs"
eval_dir <- "C:/Users/mstew/Documents/introduced_SDM_predictions/Evaluations"
setwd(directory)


# dir.create("C:/Users/mstew/Documents/Modelling/terra_tmp", 
#            showWarnings = FALSE)
terraOptions(tempdir = "C:/Users/mstew/Documents/introduced_SDM_predictions/terra_tmp", 
             memfrac = 0.5)

# variables ----
## selected bioclimatic variables
n_selected_vars <- 5
## boyce threshold for choosing models to pass to Ensemble step
boyce_threshold <- 0.75
#######################################################################################
### SCENARIO NAME -- MAKE SURE THIS IS AS DESIRED BEFORE CONTINUING
#######################################################################################
## "ClimateOnly" or "Footprint"
# scenario_name <- "ClimateOnly"
scenario_name <- "Footprint"

# identify focal species from occ_files ----
occurrence_files <-
  list.files("Data", pattern = "_dataset_clean\\.csv$", full.names = TRUE)
species_names <- basename(occurrence_files) %>% str_remove("_dataset_clean\\.csv$")
# species name variables for the first species in the occurrence files
genus <- str_split(species_names[1], "_")[[1]][1]
epithet <- str_split(species_names[1], "_")[[1]][2]
## name with _ between genus and species
species_name <- paste0(genus, "_", epithet)
## no break between genus and species
resp_name <- paste0(genus, epithet)

occ <- read.csv(occurrence_files[1])
occ$occurrenceStatus <- 1
resp <- as.numeric(occ$occurrenceStatus)

# coords <- occ[, c("selected_variables.Longitude", "selected_variables.Latitude")]
coords <- as.matrix(occ[, c("selected_variables.Longitude", 
                            "selected_variables.Latitude")])


# environmental stacks ----
environmental_stack <- rast(list.files("Layers", 
                                       pattern="\\.tif$", full.names=TRUE))
## footprint ----
footprint_name <- names(environmental_stack)[
  grepl("human-footprint|human.footprint", names(environmental_stack))]
## climate ----
clim_names <- names(environmental_stack)[
  grepl("^wc2.1_30s_bio_", names(environmental_stack))]
clim_stack <- environmental_stack[[clim_names]]

## recall selected variables
selected_clim <- read.csv(file.path(directory, "Outputs",
                                    paste0(resp_name, "_selected_variables.csv")))

# full climates ----
## based on scenario_name
if(scenario_name=="ClimateOnly"){
  env_full <- environmental_stack[[selected_clim]]
} else if(scenario_name=="Footprint"){
  env_full <- environmental_stack[[c(selected_clim, footprint_name)]]
} else{
  message("\n> Scenario name not established")
}

# ---- Target countries (.gpkg) ----
countries_file   <- vect(file.path("Vectors", "r_rattus_region.gpkg"))
# europe <- vect(file.path("Vectors", "europe_red.gpkg"))
# <-- your .gpkg path
target_countries <- c("Benin", "Chile", "Ecuador", "France", "Gabon", "Guinea", 
                      "Iran", "Niger", "Portugal", "Sierra Leone", "Japan", 
                      "Madagascar", "Australia", "New Zealand", "New Caledonia", 
                      "Puerto Rico")
## "United States of America" being done by state
# target_countries <- c("Ohio","Oklahoma","Oregon","Pennsylvania","Rhode Island",
#                       "South Carolina","South Dakota","Tennessee","Texas",
#                       "Utah","Vermont","Virginia","Washington","West Virginia",
#                       "Wisconsin","Wyoming")

# <-- or NULL to use ALL features
name_col         <- "name_en"
## for US ones
# name_col <- "woe_name"
# <-- attribute column with country names

# ---- Load layers ----
# all_layers <- rast(list.files("Layers", pattern = "\\.tif$", full.names = TRUE))

# ---- Load + filter countries, then CROP extent (before aggregation = big saving) ----
countries <- (countries_file)
if (!is.null(target_countries)) {
  keep <- values(countries)[[name_col]] %in% target_countries
  stopifnot("No matching countries — check name_col / target_countries" = any(keep))
  countries <- countries[keep, ]
}
if (!same.crs(countries, env_full)) countries <-
  project(countries, crs(env_full))
env_full <- crop(env_full, countries)
message("Cropped to ", nrow(countries), " country feature(s)")

# env_full <- crop(env_full, europe)

# ---- Aggregate (if requested), then materialise to ONE compressed file ----
if (agg_factor > 1) {
  env_full <- aggregate(env_full, fact = agg_factor, fun = "mean")
  env_full <- writeRaster(env_full, "layers_agg_tmp.tif", overwrite = TRUE,
                               gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"))
  terra::tmpFiles(orphan = TRUE, old = TRUE, remove = TRUE)  # free aggregate temp
}
message("Cells per layer: ", ncell(env_full))
available_layers <- names(env_full)

# ---- Load model + ensemble ----
model_file    <- file.path(eval_dir, paste0(species_name, "_", 
                                             scenario_name, "_model_out.rds"))
ensemble_file <- file.path(eval_dir, paste0(species_name, "_", 
                                             scenario_name, "_ensemble_out.rds"))
stopifnot(file.exists(model_file), file.exists(ensemble_file))
model_out <- (readRDS(model_file))
ensemble  <- (readRDS(ensemble_file))

# ---- Which single models does the ENSEMBLE actually need? ----
needed <- tryCatch(unique(unlist(get_kept_models(ensemble))),
                   error = function(e) "all")
if (length(needed) == 0 || all(is.na(needed))) needed <- "all"
message("Ensemble needs ",
        if (identical(needed, "all")) "ALL models" else paste(length(needed), "models"),
        ": ", paste(needed, collapse = ", "))

# ---- Variable/layer name matching ----
expl_vars <- model_out@expl.var.names
norm <- function(x) gsub("[^a-z0-9]", "", tolower(x))
layer_match <- sapply(expl_vars, function(v) {
  if (v %in% available_layers) return(v)
  idx  <- which(make.names(available_layers) == v);          if (length(idx)  == 1) return(available_layers[idx])
  idx2 <- which(endsWith(norm(available_layers), norm(v)));  if (length(idx2) == 1) return(available_layers[idx2])
  NA_character_
})
stopifnot(!any(is.na(layer_match)))

env_stack <- env_full[[layer_match]]
names(env_stack) <- expl_vars

# loop for country here, using mask to limit projections
env_by_country <- list()
for(i in seq_len(nrow(countries))){
  ctry <- countries[i, ]
  ctry_name <- values(ctry)[[name_col]]
  safe_name <- gsub("[^A-Za-z0-9]", "_", ctry_name)
  message("\n> Processing: ", ctry_name)
  ## crop then mask by country
  env_ctry <- crop(env_stack, ctry)
  env_ctry <- mask(env_ctry, ctry)
  message("\n> Cropped to ", ctry_name)
  env_by_country[[safe_name]] <- env_ctry
}

## mask to europe
# env_stack <- mask(env_stack, europe)

## project to europe
# message("\n> Projecting to Europe")
# proj <- BIOMOD_Projection(
#   bm.mod              = model_out,
#   proj.name           = paste0("Europe_", scenario_name),
#   new.env             = env_stack,
#   models.chosen       = needed,          # <<< the big saving
#   build.clamping.mask = FALSE,
#   do.stack            = TRUE,
#   on_0_1000           = TRUE,
#   nb.cpu              = 1)
# 
# ## ensemble
# message("\n> Ensemble for Europe")
# BIOMOD_EnsembleForecasting(
#   bm.em = ensemble,
#   bm.proj = proj,
#   proj.name = paste0("Europe_", scenario_name, "_ensemble"),
#   nb.cpu = 1)

## try to loop now that countries work
for(y in seq_len(nrow(countries))){
  ctry <- countries[y, ]
  ctry_name <- values(ctry)[[name_col]]
  safe_name <- gsub("[^A-Za-z0-9]", "_", ctry_name)
  message("\n> Projecting to: ", ctry_name)
  proj <- BIOMOD_Projection(
    bm.mod              = model_out,
    proj.name           = paste0(safe_name, "_", scenario_name), 
    new.env             = env_by_country[[safe_name]],
    models.chosen       = needed,          # <<< the big saving
    build.clamping.mask = FALSE,
    do.stack            = TRUE,
    on_0_1000           = TRUE,
    nb.cpu              = 1)
  message("\n> Ensemble projection for: ", ctry_name)
  BIOMOD_EnsembleForecasting(
    bm.em = ensemble,
    bm.proj = proj,
    proj.name = paste0(safe_name, "_", scenario_name, "_ensemble"),
    nb.cpu=1)
}
