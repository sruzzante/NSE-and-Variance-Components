# climatological benchmarks and variance components
# Author: Sacha Ruzzante
# sachawruzzante@gmail.com
# Last Update: 2025-08-06

# Evaluate VIC-Gl model from Schnorbus (2018)
# Schnorbus, M. (2018). VIC Glacier (VIC-GL)—Description of VIC model changes and upgrades, VIC Generation 2 Deployment Report (p. 40). Pacific Climate Impacts Consortium, University of Victoria. https://dspace.library.uvic.ca/server/api/core/bitstreams/850d5363-e326-4b13-bddc-f5cd6227c928/content




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

stns<-read.csv("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/VIC-Gl/station_metadata/stations.csv")

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
  
  
  
  
  dat<-readRDS(paste0(
      "../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/VIC-Gl/streamflow/",
    stns$gauge_id[it],".rds"
  ))
  
  

  
  if(any(is.na(dat$QSim))){
    print(sprintf("THERE ARE NA VALUES IN STATION %s (number %d), SKIPPING",stns$ID[it],it))
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
stns$mdl<-"vic-gl"


x<-bind_rows(stnPerf)%>%
  left_join(bind_rows(stnSeas),by = c("ID" ))%>%
  left_join(bind_rows(stnCompNSE),by = c("ID"),suffix = c("",".maxRun"))%>%
  left_join(stn_var)

stns<-left_join(stns,x,by = c("gauge_id" = "ID"))

# stns<-st_drop_geometry(stns)

write.csv(stns,"2.data/highLowBenchmarkGOF/GOF_vic_gl.csv")


