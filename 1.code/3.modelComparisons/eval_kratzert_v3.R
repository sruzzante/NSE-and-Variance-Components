# climatological benchmarks and variance components
# Author: Sacha Ruzzante
# sachawruzzante@gmail.com
# Last Update: 2025-06-18

# Evaluate Kratzert LSTM
# Kratzert, F., Gauch, M., Klotz, D., & Nearing, G. (2024). HESS Opinions: Never train a Long Short-Term Memory (LSTM) network on a single basin. Hydrology and Earth System Sciences, 28(17), 4187–4201. https://doi.org/10.5194/hess-28-4187-2024

# This is the regional model, with a 10-member ensemble, which was found to have the best performance

library(ncdf4)
library(hydroGOF)
library(dplyr)
library(ggplot2)
library(hydroGOF)
library(sf)
library(tmap)
library(tictoc)
library(lubridate)
library(tidyr)
library(reticulate)
setwd("/home/ruzzante/projects/def-tgleeson/ruzzante/climatological_benchmarks/")
use_python("/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/python/3.11.5/bin/python3", required = TRUE)
source("1.code/5.utils/utils.R")


Sys.setenv(UV_OFFLINE=1)
reticulate::py_require("pandas")
reticulate::py_require("xarray")
pd <- import("pandas")



# get IDs used in Kratzert model
pickle_data <- pd$read_pickle("../DATA/1.Spatial_data/regional/USA/sw_surfacewater_streamflow_runoff_river_network_waterstress/Kratzert2024/single-basin-vs-regional-model/run_dirs/regional_model/531_basins_multi_forcings_temporal_split_ensemble_member_0_0510_195856/test/model_epoch029/test_results.p")
IDs<-names(pickle_data)

# Initialize dataframe with stations
stns<-read.delim("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-us/camels_topo.txt",sep = ";",colClasses = c("character","numeric","numeric","numeric","numeric","numeric"))%>%
  st_as_sf(coords = c("gauge_lon","gauge_lat"),crs = "EPSG:4326")%>%
  left_join(read.delim("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-us/camels_name.txt",sep = ";",colClasses = c("character","character","character")))%>%
  filter(gauge_id %in% IDs)


# read in pickle files with simulated streamflow for the test period (1989-1999)
pickle_data<-list()

pickle_data[[1]]  <- pd$read_pickle("../DATA/1.Spatial_data/regional/USA/sw_surfacewater_streamflow_runoff_river_network_waterstress/Kratzert2024/single-basin-vs-regional-model/run_dirs/regional_model/531_basins_multi_forcings_temporal_split_ensemble_member_0_0510_195856/test/model_epoch029/test_results.p")
pickle_data[[2]]  <- pd$read_pickle("../DATA/1.Spatial_data/regional/USA/sw_surfacewater_streamflow_runoff_river_network_waterstress/Kratzert2024/single-basin-vs-regional-model/run_dirs/regional_model/531_basins_multi_forcings_temporal_split_ensemble_member_1_0510_195906/test/model_epoch026/test_results.p")
pickle_data[[3]]  <- pd$read_pickle("../DATA/1.Spatial_data/regional/USA/sw_surfacewater_streamflow_runoff_river_network_waterstress/Kratzert2024/single-basin-vs-regional-model/run_dirs/regional_model/531_basins_multi_forcings_temporal_split_ensemble_member_2_0510_195904/test/model_epoch028/test_results.p")
pickle_data[[4]]  <- pd$read_pickle("../DATA/1.Spatial_data/regional/USA/sw_surfacewater_streamflow_runoff_river_network_waterstress/Kratzert2024/single-basin-vs-regional-model/run_dirs/regional_model/531_basins_multi_forcings_temporal_split_ensemble_member_3_0510_195855/test/model_epoch022/test_results.p")
pickle_data[[5]]  <- pd$read_pickle("../DATA/1.Spatial_data/regional/USA/sw_surfacewater_streamflow_runoff_river_network_waterstress/Kratzert2024/single-basin-vs-regional-model/run_dirs/regional_model/531_basins_multi_forcings_temporal_split_ensemble_member_4_0510_195910/test/model_epoch022/test_results.p")
pickle_data[[6]]  <- pd$read_pickle("../DATA/1.Spatial_data/regional/USA/sw_surfacewater_streamflow_runoff_river_network_waterstress/Kratzert2024/single-basin-vs-regional-model/run_dirs/regional_model/531_basins_multi_forcings_temporal_split_ensemble_member_5_0510_195859/test/model_epoch030/test_results.p")
pickle_data[[7]]  <- pd$read_pickle("../DATA/1.Spatial_data/regional/USA/sw_surfacewater_streamflow_runoff_river_network_waterstress/Kratzert2024/single-basin-vs-regional-model/run_dirs/regional_model/531_basins_multi_forcings_temporal_split_ensemble_member_6_0510_195906/test/model_epoch026/test_results.p")
pickle_data[[8]]  <- pd$read_pickle("../DATA/1.Spatial_data/regional/USA/sw_surfacewater_streamflow_runoff_river_network_waterstress/Kratzert2024/single-basin-vs-regional-model/run_dirs/regional_model/531_basins_multi_forcings_temporal_split_ensemble_member_7_0510_195908/test/model_epoch025/test_results.p")
pickle_data[[9]]  <- pd$read_pickle("../DATA/1.Spatial_data/regional/USA/sw_surfacewater_streamflow_runoff_river_network_waterstress/Kratzert2024/single-basin-vs-regional-model/run_dirs/regional_model/531_basins_multi_forcings_temporal_split_ensemble_member_8_0510_195855/test/model_epoch025/test_results.p")
pickle_data[[10]] <- pd$read_pickle("../DATA/1.Spatial_data/regional/USA/sw_surfacewater_streamflow_runoff_river_network_waterstress/Kratzert2024/single-basin-vs-regional-model/run_dirs/regional_model/531_basins_multi_forcings_temporal_split_ensemble_member_9_0510_195901/test/model_epoch027/test_results.p")


# initialize list to store performance metrics
stnPerf<-vector(mode = "list", length = nrow(stns))

# initialize list to store seasonality metrics
stnSeas<-vector(mode = "list", length = nrow(stns))



# loop through stations
for(it in 1:nrow(stns)){
  tic(sprintf("station %d",it))
  
  # Read in data from each member of the 10-member ensemble model
  dat_ls<-list()
  for(it_ens in 1:10){
    xarray_dat<-pickle_data[[it_ens]][[it]]$`1D`$xr
    dat_ls[[it_ens]]<-
      data.frame(
        QObs = xarray_dat$`QObs(mm/d)_obs`$data,
        QSim = xarray_dat$`QObs(mm/d)_sim`$data,
        dt = as.Date(xarray_dat$date$data)
        
      )
    
  }
  # take mean of ensemble members
  dat<-dat_ls%>%
    bind_rows()%>%
    group_by(dt)%>%
    summarize(across(QObs:QSim,~mean(.x)))
  
  if(any(is.na(dat$QSim))){
    print(sprintf("THERE ARE NA VALUES IN STATION %s (number %d), SKIPPING",stns$ID[it],it))
    {next}
  }
  
  # calculate goodness of fit statistics
  perfMetrics<-calcAllGOF(dat,waterYearStart = 10,minYears = 8)
  perfMetrics$gauge_id<-stns$gauge_id[it]
  stnPerf[[it]]<-data.frame(perfMetrics)
  
  # calculate seasonality statistics
  seas<-calcSeasonality(dat)
  stnSeas[[it]] = data.frame(gauge_id = stns$gauge_id[it],
                             CoV = seas["CoV"],
                             QCI = seas["QCI"])
  
  
  toc()
}


# save data

stns$mdl<-"lstm-kratzert2024"


x<-bind_rows(stnPerf)%>%
  left_join(bind_rows(stnSeas))


stns<-left_join(stns,x)
write.csv(stns%>%st_drop_geometry(),"2.data/highLowBenchmarkGOF/GOF_kratzert_lstm.csv")


