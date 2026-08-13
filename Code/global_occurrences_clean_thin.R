## Coordinate cleaner and thinning
## With thanks to Felipe Espinoza who's work provides a basis for this code
## and for his assistance with errors and providing advice.

## This code takes the occurrence data (downloaded from GBIF) and cleans and thins
## saves as species_dataset_clean.csv to the Data folder

# clean ----
rm(list=ls())

# packages ----
library(CoordinateCleaner)
library(spThin)
library(dplyr)
library(stringr)

# set directory ----
directory <- "C:/Users/mstew/Documents/introduced_SDM_predictions"
data_dir <- "C:/Users/mstew/Documents/introduced_SDM_predictions/Data"
out_dir <- "C:/Users/mstew/Documents/introduced_SDM_predictions/Outputs"
eval_dir <- "C:/Users/mstew/Documents/introduced_SDM_predictions/Evaluations"
setwd(directory)

# occ data clean and save ----
## list files in the folder to confirm number of species
occurrence_files <-
  list.files("Data", pattern = "_occ\\.csv$", full.names = TRUE)
## get species names directly from the file names
species_names <- basename(occurrence_files) %>% str_remove("_occ\\.csv$")

## clean function ----
clean_thin <- function(file_path, species_name, max_n=5000){
  ## extract model name from the file
  genus <- str_split(species_name, "_")[[1]][1]
  epithet <- str_split(species_name, "_")[[1]][2]
  resp_name <- paste0(genus, "_", epithet)
  ## pull occ data in
  ## this is the download from gbif
  gbif_data <- read.delim(file_path)
  gbif_data$occurrenceStatus <- 1
  ## ensure correctness of coords
  gbif_data <- subset(gbif_data,
                      !is.na(decimalLatitude) & !is.na(decimalLongitude) &
                        decimalLatitude >= -90 & decimalLatitude <= 90 &
                        decimalLongitude >= -180 & decimalLongitude <= 180 &
                        !(decimalLatitude == 0 & decimalLongitude == 0))
  ## clean coords
  cleaned_data <- CoordinateCleaner::clean_coordinates(gbif_data, 
                                                       "decimalLongitude", 
                                                       "decimalLatitude",
                                                       species="species")
  ## cleaned data
  valid_data <- cleaned_data %>%
    select(decimalLongitude, decimalLatitude, occurrenceStatus, species) %>%
    mutate(occurrenceStatus = as.factor(occurrenceStatus))
  ## occ points for thinning steps
  occ_points <- valid_data[, c("decimalLongitude", "decimalLatitude")]
  occ_df <- data.frame(
    Longitude = occ_points[,1],
    Latitude  = occ_points[,2],
    Species   = valid_data$occurrenceStatus)
  
  ## thinning prep
  n <- nrow(occ_df); k <- max(1, ceiling(n / max_n))
  if (k == 1) return(list(occ_df))
  cl <- stats::kmeans(occ_df[, c("Longitude","Latitude")], 
                      centers = k, iter.max = 50)$cluster
  parts <- split(occ_df, cl)
  
  ## thin tiles, loop
  occ_thinned <- NULL
  for(tile in 1:length(parts)){
    thin_results <- thin(
      loc.data = parts[[tile]],
      lat.col = "Latitude",
      long.col = "Longitude",
      spec.col = "Species",
      thin.par = 1,          # distance in km
      reps = 1,              # number of replicates
      locs.thinned.list.return = TRUE,
      write.files = FALSE,
      write.log.file = FALSE)
    thin_results[[1]]
    occ_thinned <- rbind(occ_thinned, thin_results[[1]])
  }
  
  # file_name <- paste0(gsub(" ", "_", species_name), "_dataset_clean.csv")
  # write.csv(occ_thinned, file_name, row.names = FALSE)
  
  write_csv(data.frame(selected_variables = occ_thinned),
            file.path("Data", paste0(resp_name, "_dataset_clean.csv")))
}

# run cleaner ----
# clean_thin(occurrence_files, species_names)
walk2(occurrence_files, species_names, ~ clean_thin(.x, .y))