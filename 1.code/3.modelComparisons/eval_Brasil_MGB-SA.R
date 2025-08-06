# climatological benchmarks and variance components
# Author: Sacha Ruzzante
# sachawruzzante@gmail.com
# Last Update: 2025-08-06

# Evaluate MGB-SA for Camels-BR dataset
# Chagas, V. B. P., Chaffe, P. L. B., Addor, N., Fan, F. M., Fleischmann, A. S., Paiva, R. C. D., & Siqueira, V. A. (2020). CAMELS-BR: Hydrometeorological time series and landscape attributes for 897 catchments in Brazil. Earth System Science Data, 12(3), 2075–2096. https://doi.org/10.5194/essd-12-2075-2020



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


setwd("/home/ruzzante/projects/def-tgleeson/ruzzante/climatological_benchmarks/")
source("1.code/5.utils/utils.R")


# read in pickle files with simulated streamflow for the test period (1990-2009)

stns<-st_read("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-br/13_CAMELS_BR_gauge_location/location_gauges_streamflow.gpkg")

sim_ids<-list.files("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-br/04_CAMELS_BR_streamflow_simulated/")%>%
  str_remove("_simulated_streamflow.txt")
stns<-filter(stns,gauge_id %in% sim_ids)

# initialize list to store performance metrics
stnPerf<-vector(mode = "list", length = nrow(stns))

# initialize list to store seasonality metrics
stnSeas<-vector(mode = "list", length = nrow(stns))

# initialize list to store goodness of fit for variance components
stnCompNSE<- vector(mode = "list", length = nrow(stns))

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


# loop through stations
for(it in 1:nrow(stns)){
  tic(sprintf("station %d",it))
  dat_obs<-read.delim(sprintf("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-br/02_CAMELS_BR_streamflow_all_catchments/%s_streamflow.txt",stns$gauge_id[it]),
                      sep = " ")%>%
    filter(qual_flag == 1)
  
  dat_sim<-read.delim(sprintf("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-br/04_CAMELS_BR_streamflow_simulated/%s_simulated_streamflow.txt",stns$gauge_id[it]),
                      sep = " ")
  
  dat<-left_join(dat_obs,dat_sim)%>%
    mutate(dt = ymd(paste(year,month,day,sep = "/")),
           
           yday = pmin(yday(dt),365),
           QObs =streamflow_m3s,
           QSim = simulated_streamflow_m3s)%>%
    filter(year%in% (1980:2014))
  
 if(nrow(dat)<(365)){next}
  
  
  if(any(is.na(dat$QSim))){
    print(sprintf("THERE ARE NA VALUES IN STATION %s (number %d), SKIPPING",stns$gauge_id[it],it))
    {next}
  }
  
  # calculate goodness of fit statistics
  perfMetrics<-calcAllGOF(dat,waterYearStart = 10,minYears = 7)
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
  
  # calculate variance components
  dat<-mutate(dat,
              q = QObs,
              year  =year(dt),
              yday = pmin(yday(dt),365))
  
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

# stns<-data.frame(
#   ID =  IDs
# )
# is the benchmark efficiency worse in high-NSE benchmark stations
stns$mdl<-"MGB-SA-camels-br"

x<-bind_rows(stnPerf)%>%
  left_join(bind_rows(stnSeas))%>%
  left_join(bind_rows(stnCompNSE),by = c("ID"),suffix = c("",".maxRun"))%>%
  left_join(stn_var)


stns<-left_join(stns,x,by =c("gauge_id" = "ID"))
stns<-st_drop_geometry(stns)
write.csv(stns,"2.data/highLowBenchmarkGOF/GOF_camels-br-MGB-SA.csv")








