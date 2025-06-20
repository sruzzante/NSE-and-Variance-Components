

library(ncdf4)
library(hydroGOF)
library(dplyr)
library(ggplot2)
library(hydroGOF)
library(sf)
library(tmap)
library(tictoc)
library(lubridate)
library(stringr)
setwd("/home/ruzzante/projects/def-tgleeson/ruzzante/climatological_benchmarks/lstm-camelsbr/")

fls<-list.files("../../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-br/preprocessed/")%>%
  str_remove(".csv")

write.table(fls,
            "allbasins.txt",
            col.names = FALSE,
            append = FALSE,quote = FALSE,
            row.names = FALSE)

forcings<-read.csv("../../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-br/preprocessed/10500000.csv")

static<-read.csv("../../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-br/01_CAMELS_BR_attributes/CAMELS_BR_attributes_description.xlsx")