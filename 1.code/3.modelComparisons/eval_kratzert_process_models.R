# climatological benchmarks and variance components
# Author: Sacha Ruzzante
# sachawruzzante@gmail.com
# Last Update:2025-08-06


# evaluate the camels benchmarks compiled by Kratzert (2019)
# Kratzert, F. (2019). CAMELS benchmark models [Dataset]. HydroShare. https://doi.org/10.4211/hs.474ecc37e7db45baa425cdb4fc1b61e1

library(ncdf4)
library(hydroGOF)
library(dplyr)
library(tidyr)
library(ggplot2)
library(hydroGOF)
library(sf)
# library(tmap)
library(tictoc)
library(lubridate)
library(stringr)
setwd("/home/ruzzante/projects/def-tgleeson/ruzzante/climatological_benchmarks/")
source("1.code/5.utils/utils.R")

# path to dataset save directory
pth = "../DATA/1.Spatial_data/regional/USA/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels_benchmark/netcdf/"

IDs<-list.files(pth)%>%
  str_remove("_benchmark_models.nc")

# all models included in the dataset
mdls<-c("SAC_SMA","VIC_basin","HBV_lb","HBV_ub","q_sim_fuse_900","q_sim_fuse_902","q_sim_fuse_904","VIC_conus","mHm_basin","mHm_conus")

# select only one of each type (the best)
mdls<-c("SAC_SMA","VIC_basin","HBV_ub","q_sim_fuse_900","mHm_basin")

# initialize datafrrame of station attributes from CAMELS dataset
stns<-read.delim("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-us/camels_topo.txt",sep = ";",colClasses = c("character","numeric","numeric","numeric","numeric","numeric"))%>%
  st_as_sf(coords = c("gauge_lon","gauge_lat"),crs = "EPSG:4326")%>%
  left_join(read.delim("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-us/camels_name.txt",sep = ";",colClasses = c("character","character","character")))%>%
  filter(gauge_id %in% IDs)

# store blank dataframe
stns_na<-stns

# stns_ls<-list()
for(it_mdl in 1:length(mdls)){
  
  stns<-stns_na
  
  # initialize list to store performance metrics
  stnPerf<-vector(mode = "list", length = nrow(stns))
  
  # initialize list to store seasonality metrics
  stnSeas<-vector(mode = "list", length = nrow(stns))
  
  # initialize dataframe to store variance components
  stn_var<-data.frame(gauge_id = stns$gauge_id,
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
  
  # initialize list to store goodness of fit for variance components
  stnCompNSE<- vector(mode = "list", length = nrow(stns))
  
  
  for(it in 1:nrow(stns)){
    tic(sprintf("station %d",it))
    
    # load streamflow data
    dat_nc<-nc_open(sprintf("%s%s_benchmark_models.nc",pth,stns$gauge_id[it]))
    
    if(!mdls[it_mdl] %in% names(dat_nc$var)){next} # some models are missing some gauges
    
    dat<-
      data.frame(
        QObs =ncvar_get(dat_nc,"QObs"),
        QSim = ncvar_get(dat_nc,mdls[it_mdl]),
        ind = ncvar_get(dat_nc,"index")
      )%>%
      mutate(dt = ymd("1989-10-01")+ind)
    
    if(any(is.na(dat$QSim))){
      print(sprintf("THERE ARE NA VALUES IN STATION %s (number %d), SKIPPING",stns$gauge_id[it],it))
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
    
    # Goodness of fit of variance components
    NSE_comps<- gof_components(dat)
    NSE_comps$gauge_id<-stns$gauge_id[it]
    stnCompNSE[[it]] = data.frame(NSE_comps)
    
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
  
  stns$mdl<-mdls[it_mdl]
  
  x<-bind_rows(stnPerf)%>%
    left_join(bind_rows(stnSeas))%>%
    left_join(bind_rows(stnCompNSE),by = c("gauge_id"),suffix = c("",".maxRun"))%>%
    left_join(stn_var)
  



stns<-left_join(stns,x)

write.csv(stns%>%st_drop_geometry(),sprintf("2.data/highLowBenchmarkGOF/GOF_kratzert_%s.csv",mdls[it_mdl]))
stn_result<-compareMetricsHighLow(stns,splitVar = "NSEB",thresh = 0.5)

write.csv(stn_result,sprintf("2.data/highLowBenchmarkGOF/highLowBenchmarkGOF_kratzert_%s.csv",mdls[it_mdl]))

}



