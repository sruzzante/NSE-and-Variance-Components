# climatological benchmarks and variance components
# Author: Sacha Ruzzante
# sachawruzzante@gmail.com
# Last Update: 2025-08-06

# Evaluate PREVAH model from Kraft (2025)
# Kraft, B., Schirmer, M., Aeberhard, W. H., Zappa, M., Seneviratne, S. I., & Gudmundsson, L. (2025). CH-RUN: A deep-learning-based spatially contiguous runoff reconstruction for Switzerland. Hydrology and Earth System Sciences, 29(4), 1061–1082. https://doi.org/10.5194/hess-29-1061-2025
# Gauge-based simulations were share privately via email



library(ncdf4)
library(hydroGOF)
library(dplyr)
library(ggplot2)
library(hydroGOF)
library(sf)
library(tmap)
library(tictoc)
library(lubridate)
setwd("/home/ruzzante/projects/def-tgleeson/ruzzante/climatological_benchmarks/")
source("1.code/5.utils/utils.R")

# load ncdf data
stns<-st_read("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-ch/camels_ch/catchment_delineations/CAMELS_CH_gauging_stations.shp")%>%
  filter(type == "stream")

# initialize list to store performance metrics
stnPerf<-vector(mode = "list", length = nrow(stns))

# initialize list to store seasonality metrics
stnSeas<-vector(mode = "list", length = nrow(stns))

# initialize dataframe to store variance components
stn_var<-data.frame(ID = stns$gauge_id,
                    varSeas_stl = NA,
                    varInterannual_stl = NA,
                    varRem_stl = NA,
                    varSeas_clas = NA,
                    varInterannual_clas = NA,
                    varRem_clas = NA,
                    varSeas_fourier = NA,
                    varInterannual_fourier = NA,
                    varRem_fourier = NA
)


stnCompNSE<- vector(mode = "list", length = nrow(stns))
# loop through stations

for(it in 1:nrow(stns)){
  tic(sprintf("station %d",it))
  
  # read data in
  
  dat_obs<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-ch/camels_ch/timeseries/observation_based/CAMELS_CH_obs_based_",
                         stns$gauge_id[it],".csv"),
                  # skip  =7,
                  sep  =",")%>%
    mutate(dt = ymd(date),
           QObs = discharge_vol.m3.s.,
           yday = pmin(yday(dt),365),
           year = year(dt))
  
  
  
  dat_sim<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-ch/camels_ch/timeseries/simulation_based/CAMELS_CH_sim_based_",
                         stns$gauge_id[it],".csv"),
                  # skip  =7,
                  sep  =",")%>%
    mutate(dt = ymd(date),
           QSim = discharge_vol_sim.m3.s.)
  
  dat<-left_join(dat_obs%>%select(dt,QObs,yday,year,QObs),
                 dat_sim%>% select(dt,QSim))
  

  
  if(any(is.na(dat$QSim))){
    print(sprintf("THERE ARE NA VALUES IN STATION %s (number %d), SKIPPING",stns$gauge_id[it],it))
    {next}
  }
  
  # calculate goodness of fit statistics
  perfMetrics<-calcAllGOF(dat,waterYearStart = 10)
  perfMetrics$ID<-stns$gauge_id[it]
  stnPerf[[it]]<-data.frame(perfMetrics)
  
  # calculate seasonality statistics
  seas<-calcSeasonality(dat)
  stnSeas[[it]] = data.frame(ID = stns$gauge_id[it],
                             CoV = seas["CoV"],
                             QCI = seas["QCI"])
  
  # Goodness of fit of variance components
  NSE_comps<- gof_components(dat)
  NSE_comps$ID<-stns$gauge_id[it]
  stnCompNSE[[it]] = data.frame(NSE_comps)
  
  
  dat<-mutate(dat,
              q = QObs,
              year  =year(dt),
              yday = pmin(yday(dt),365))
  
  # calculate variance components
  
  # STL decomposition
  stl_var_x<-stl_var(dat,s.window = 7,t.window = 365) 
  stn_var$varSeas_stl[it] <-stl_var_x["varSeas"] 
  stn_var$varInterannual_stl[it] <-stl_var_x["varInterannual"] 
  stn_var$varRem_stl[it] <-stl_var_x["varRem"]    
  
  #classical decomposition
  clas_var_x<-clas_var(dat)    
  stn_var$varSeas_clas[it] <-clas_var_x["varSeas"]  
  stn_var$varInterannual_clas[it] <-clas_var_x["varInterannual"]  
  stn_var$varRem_clas[it] <-clas_var_x["varRem"]
  
  # Fourier decomposition
  fourier_var_x<-fourier_var(dat)    
  stn_var$varSeas_fourier[it] <-fourier_var_x["varSeas"]  
  stn_var$varInterannual_fourier[it] <-fourier_var_x["varInterannual"]  
  stn_var$varRem_fourier[it] <-fourier_var_x["varRem"]
  
  toc()
}

# is the benchmark efficiency worse in high-NSE benchmark stations
stns$mdl<-"prevah-camels_ch"


x<-bind_rows(stnPerf)%>%
  left_join(bind_rows(stnSeas),by = c("ID" = "ID"))%>%
  left_join(bind_rows(stnCompNSE),by = c("ID"),suffix = c("",".maxRun"))%>%
  left_join(stn_var)

stns<-left_join(stns,x,by = c("gauge_id" = "ID"))

stns<-st_drop_geometry(stns)

write.csv(stns,"2.data/highLowBenchmarkGOF/GOF_camels_ch_Prevah.csv")


