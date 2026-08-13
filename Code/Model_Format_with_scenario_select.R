## SDM Protocol from model format to output
## With thanks to Felipe Espinoza who's work provides a basis for this code
## and for his assistance with errors and providing advice.

## VIF select 5 bioclim variables for a focal species
## format the predictor layers, pseudo-absences and occurrences for modelling
## save model outputs

# clean and libraries ----
rm(list=ls())
gc()

library(terra)
library(usdm)
library(biomod2)
library(readr)
library(dplyr)
library(stringr)

# file structure ----
## file structure should be as follows:
# .../Modelling
#       /Data         -> global occurrence data and relative local abundances 
#       /Layers       -> predictor data as .tif
#       /Outputs      -> location for selected variables, comparison and correlation results
#       /Projections  -> Ensemble projections destination
#       /Vectors      -> location for target region shapefiles
#       /terra_tmp    -> terra package's temporary directory (needs frequent cleaning)

# wd ----
## working directory as "Modelling" folder
directory <- "C:/Users/mstew/Documents/introduced_SDM_predictions"
data_dir <- "C:/Users/mstew/Documents/introduced_SDM_predictions/Data"
out_dir <- "C:/Users/mstew/Documents/introduced_SDM_predictions/Outputs"
eval_dir <- "C:/Users/mstew/Documents/introduced_SDM_predictions/Evaluations"
setwd(directory)


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

## data (occs and bioclim)
# occurrences ----
## clean version of the global occurrences
occurrence_files <-
  list.files("Data", pattern = "_dataset_clean\\.csv$", full.names = TRUE)
## get species names directly from the file names
species_names <- basename(occurrence_files) %>% str_remove("_dataset_clean\\.csv$")

## target a single species from the occurrence files
occ <- read.csv(occurrence_files[1])
species_name <- species_names[1]
genus <- str_split(species_name, "_")[[1]][1]
epithet <- str_split(species_name, "_")[[1]][2]
## response name as genus/species with no name to prevent pathway issues
resp_name <- paste0(genus, epithet)
## assign presences
occ$occurrenceStatus <- 1
## set as response
resp <- as.numeric(occ$occurrenceStatus)
## matrix coordinates for model formatting
coords <- as.matrix(occ[, c("selected_variables.Longitude", 
                            "selected_variables.Latitude")])

# environmental stacks ----
environmental_stack <- rast(list.files("Layers", 
                                       pattern="\\.tif$", full.names=TRUE))
## footprint
names(environmental_stack)[names(environmental_stack) == 
                             "layer_wildareas-2009-human-footprint"] <-
  "layer_wildareas_2009_human_footprint"
footprint_name <- names(environmental_stack)[
  grepl("human-footprint|human.footprint", names(environmental_stack))]
## bioclim layers
clim_names <- names(environmental_stack)[
  grepl("^wc2.1_30s_bio_", names(environmental_stack))]
## stack
clim_stack <- environmental_stack[[clim_names]]

# VIF selection ----
env_vals <- terra::extract(clim_stack, coords)
env_vals <- env_vals[,-1]
env_vals <- env_vals[complete.cases(env_vals),]

vif_res <- vifstep(env_vals, th=5)
selected_clim <- head(vif_res@results$Variables, n_selected_vars)
## save selected variables, these can be recalled later instead of rerunning selection
write_csv(data.frame(selected_variables = selected_clim),
          file.path("Outputs", paste0(resp_name, "_selected_variables.csv")))
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

# model format ----
## full global environmental layers
biomod_data <- BIOMOD_FormatingData(
  resp.var = resp,
  expl.var = env_full,
  resp.xy = coords,
  resp.name = resp_name,
  PA.nb.rep = 1,
  PA.nb.absences = 10000,
  PA.strategy = "random",
  filter.raster = TRUE)

# build with cross-validation ----
message("\n> Model with cross validation")
model_out <- BIOMOD_Modeling(
  bm.format = biomod_data,
  modeling.id = scenario_name,
  models = c("RF", "GLM", "MAXNET"),
  CV.strategy = "random",    ## fixed, proper cross-validation
  CV.nb.rep = 5,            ## 5-fold cross-valid; 
  CV.perc = 0.7,            ## 70% for training
  OPT.strategy = "default",
  metric.eval = c("AUCroc", "TSS", "BOYCE"),
  var.import = 0,           ## 0 to save time
  seed.val = 42)

# ensemble model ----
message("\n> Ensemble")
ensemble <- BIOMOD_EnsembleModeling(
  bm.mod = model_out,
  models.chosen = "all",
  em.by = "all",
  em.algo = "EMmean",
  metric.select = "BOYCE",
  metric.select.thresh = boyce_threshold,
  metric.eval = c("AUCroc", "TSS", "BOYCE"),
  var.import = 0,
  seed.val = 42
)

# save model and ensemble outputs ----
message("\n> Saving model_out for later use")
saveRDS(
  model_out,
  file.path(eval_dir, paste0(resp_name, "_", scenario_name, "_model_out.rds")))
message("\n> Saving ensemble for later use")
saveRDS(
  ensemble,
  file.path(eval_dir, paste0(resp_name, "_", scenario_name, "_ensemble_out.rds")))

# save model and ensemble evaluations ----
message("\n> Saving model evaluations")
model_evals<- get_evaluations(model_out)
saveRDS(model_evals, 
        file.path(eval_dir, paste0(resp_name, "_", scenario_name,
                                    "_model_evaluations.rds")))
message("\n> Saving ensemble evaluations")
ensemble_evals <- get_evaluations(ensemble)
saveRDS(ensemble_evals, 
        file.path(eval_dir, paste0(resp_name, "_", scenario_name, 
                                    "_ensemble_evaluations.rds")))
# cross-validation metrics ----
## as a saved .csv
cv_metrics_list <- list()
for(model_name in names(model_evals)){
  for(metric_name in names(model_evals[[model_name]])){
    cv_values <- model_evals[[model_name]][[metric_name]]
    cv_values <- cv_values[!is.na(cv_values)]
    if(length(cv_values) > 0){
      cv_metrics_list[[length(cv_metrics_list) + 1]] <- data.frame(
        Model = model_name,
        Metric = metric_name,
        Mean_Value = mean(cv_values),
        SD_Value = sd(cv_values),
        Min_Value = min(cv_values),
        Max_Value = max(cv_values)
      )
    }
  }
}
cv_metrics_df <- do.call(rbind, cv_metrics_list)
write.csv(cv_metrics_df,
          file.path(directory, "Outputs", paste0(species_name, "_", scenario_name,
                                                 "_crossvalidated_metrics.csv")),
          row.names = FALSE)