## plot layers
## plot predictive layers to a visual form
## includes legend

rm(list=ls())

library(terra)
library(dplyr)
library(raster)

# folders/working directory ----
directory <- "C:/Users/mstew/Documents/Modelling"
setwd(directory)

# environment layers ----
bio_1 <- rast("Layers/wc2.1_30s_bio_1.tif")
bio_2 <- rast("Layers/wc2.1_30s_bio_2.tif")
bio_3 <- rast("Layers/wc2.1_30s_bio_3.tif")
bio_4 <- rast("Layers/wc2.1_30s_bio_4.tif")
bio_5 <- rast("Layers/wc2.1_30s_bio_5.tif")
bio_6 <- rast("Layers/wc2.1_30s_bio_6.tif")
bio_7 <- rast("Layers/wc2.1_30s_bio_7.tif")
bio_8 <- rast("Layers/wc2.1_30s_bio_8.tif")
bio_9 <- rast("Layers/wc2.1_30s_bio_9.tif")
bio_10 <- rast("Layers/wc2.1_30s_bio_10.tif")
bio_11 <- rast("Layers/wc2.1_30s_bio_11.tif")
bio_12 <- rast("Layers/wc2.1_30s_bio_12.tif")
bio_13 <- rast("Layers/wc2.1_30s_bio_13.tif")
bio_14 <- rast("Layers/wc2.1_30s_bio_14.tif")
bio_15 <- rast("Layers/wc2.1_30s_bio_15.tif")
bio_16 <- rast("Layers/wc2.1_30s_bio_16.tif")
bio_17 <- rast("Layers/wc2.1_30s_bio_17.tif")
bio_18 <- rast("Layers/wc2.1_30s_bio_18.tif")
bio_19 <- rast("Layers/wc2.1_30s_bio_19.tif")
elev <- rast("Layers/wc2.1_elev.tif")

footprint <- rast("Layers/layer_wildareas_2009_human_footprint.tif")
# global(footprint, range, na.rm = TRUE)
footprint[footprint>201] <- NA
plot(footprint, main = "Human Footprint Layer", col=map.pal("viridis"))

# plot ----
par(mfrow=c(3,2))
## r rattus bios: 2, 14, 15, 18, 19
# par(mfrow=c(1,3))
plot(bio_2, main = "A. Mean diurnal range", col=map.pal("viridis"))
plot(bio_14, main = "B. Precipitation of driest month", col=map.pal("viridis"))
plot(bio_15, main = "C. Precipitation seasonality", col=map.pal("viridis"))
# par(mfrow=c(1,2))
plot(bio_18, main = "D. Precipitation of warmest quarter", col=map.pal("viridis"))
plot(bio_19, main = "E. Precipitation of coldest quarter", col=map.pal("viridis"))
## n vison bios: 2, 10, 13, 15, 18
plot(bio_2, main = "A. Mean diurnal range", col=map.pal("viridis"))
plot(bio_10, main = "B. Mean temperature of warmest quarter", col=map.pal("viridis"))
plot(bio_13, main = "C. Precipitation of wettest month", col=map.pal("viridis"))
plot(bio_15, main = "D. Precipitation seasonality", col=map.pal("viridis"))
plot(bio_18, main = "E. Precipitation of warmest quarter", col=map.pal("viridis"))

