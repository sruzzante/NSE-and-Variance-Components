# climatological benchmarks and variance components
# Author: Sacha Ruzzante
# sachawruzzante@gmail.com
# Last Update: 2025-06-18

# Evaluate LSTM for Camels-BR dataset
# The model runs are stored in 2.data/4.lstm-camelsbr


library(ncdf4)
library(hydroGOF)
library(dplyr)
library(ggplot2)
library(hydroGOF)
library(sf)
# library(tmap)
library(tictoc)
library(lubridate)
library(stringr)
library(reticulate)

# reticulate is a bit finnicky, this worked on compute canada - Narval
Sys.setenv(UV_OFFLINE=1)
use_python("/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/python/3.11.5/bin/python3", required = TRUE)

reticulate::py_require("pandas")
reticulate::py_require("xarray")
pd <- import("pandas")

setwd("/home/ruzzante/projects/def-tgleeson/ruzzante/climatological_benchmarks/")
source('1.code/utils.R')


# read in pickle files with simulated streamflow for the test period (1990-2009)
pickle_data<-list()
pickle_data[[1]] <- pd$read_pickle("lstm-camelsbr/runs/test_run_2605_194845/test/model_epoch004/test_results.p")
pickle_data[[2]] <- pd$read_pickle("lstm-camelsbr/runs/ens1_2705_130421/test/model_epoch002/test_results.p")
pickle_data[[3]] <- pd$read_pickle("lstm-camelsbr/runs/ens2_2705_164646/test/model_epoch006/test_results.p")
pickle_data[[4]] <- pd$read_pickle("lstm-camelsbr/runs/ens3_2705_154103/test/model_epoch003/test_results.p")
pickle_data[[5]] <- pd$read_pickle("lstm-camelsbr/runs/ens4_2705_154103/test/model_epoch003/test_results.p")

IDs<-names(pickle_data[[1]])

stns<-data.frame(
  ID =  IDs
)


# initialize list to store performance metrics
stnPerf<-vector(mode = "list", length = nrow(stns))

# initialize list to store seasonality metrics
stnSeas<-vector(mode = "list", length = nrow(stns))


# loop through stations

for(it in 1:nrow(stns)){
  tic(sprintf("station %d",it))
  dat_obs<-read.csv(sprintf("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-br/preprocessed/%s.csv",stns$ID[it]))
  # Read in data from each member of the 10-member ensemble model
  dat_ls<-list()
  for(it_ens in 1:5){
    
    # Use reticulate to extract the datetime64[ns] array from Python as a list
    py_dates <- xarray_dat$date$data
    
    # Manually convert each datetime to R Date
    dates <- sapply(py_dates, function(d) as.Date(py_to_r(d)))
    
    
    xarray_dat<-pickle_data[[it_ens]][[it]]$`1D`$xr
    
    dates <- sapply(xarray_dat$date         $data, function(d) as.Date(py_to_r(d)))
    
    dat_ls[[it_ens]]<-
      data.frame(
        QObs = xarray_dat$streamflow_mm_obs  $data,
        QSim = xarray_dat$streamflow_mm_sim  $data,
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
  perfMetrics<-calcAllGOF(dat,waterYearStart = 10,minYears = 7)
  perfMetrics$ID<-stns$ID[it]
  stnPerf[[it]]<-data.frame(perfMetrics)
  
  # calculate seasonality statistics
  seas<-calcSeasonality(dat)
  stnSeas[[it]] = data.frame(ID = stns$ID[it],
                             CoV = seas["CoV"],
                             QCI = seas["QCI"])
  
  
  
  toc()
}

stns<-data.frame(
  ID =  IDs
)
# is the benchmark efficiency worse in high-NSE benchmark stations
stns$mdl<-"lstm-camels-br"

x<-bind_rows(stnPerf)%>%
  left_join(bind_rows(stnSeas))%>%
  left_join(stn_var)


stns<-left_join(stns,x)

write.csv(stns,"2.data/highLowBenchmarkGOF/GOF_camels-br.csv")








