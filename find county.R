library(raster)
library(sf)
library(here)
library(ggplot2)
library(tmap)
library(patchwork)
library(measurements)
library(scam)
library(parallel)
library(future.apply)
library(dplyr)
library(tidyverse)
library(modelsummary)
library(jsonlite)
library(data.table)
library(exactextractr)
library(future.apply)
library(sf)
library(here)
library(grf)
library(dplyr)
library(ggplot2)
library(data.table)
library(mgcv)
library(gratia)
library(spdep)
library(spatialreg)
library(FedData)


source(
  "https://github.com/tmieno2/OnFarmExperiments/blob/master/Functions/prepare.R?raw=TRUE",
  local = TRUE
)

#--- github ---#
source("https://raw.githubusercontent.com/brittanikedge/DIFM/main/Functions.R")
source(here("./Code/functions.R"))

#--- define field-year ---#
ffy <- "Wendte_LaueLib80_2017"

get_county(ffy)
