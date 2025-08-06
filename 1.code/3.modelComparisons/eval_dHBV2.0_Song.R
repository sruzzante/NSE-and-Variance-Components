# climatological benchmarks and variance components
# Author: Sacha Ruzzante
# sachawruzzante@gmail.com
# Last Update: 2025-08-06

# Evaluate Song et al model (dHBV)
# Song, Y., Bindas, T., Shen, C., Ji, H., Knoben, W. J. M., Lonzarich, L., Clark, M. P., Liu, J., van Werkhoven, K., Lamont, S., Denno, M., Pan, M., Yang, Y., Rapp, J., Kumar, M., Rahmani, F., Thébault, C., Adkins, R., Halgren, J., … Lawson, K. (2025). High-Resolution National-Scale Water Modeling Is Enhanced by Multiscale Differentiable Physics-Informed Machine Learning. Water Resources Research, 61(4), e2024WR038928. https://doi.org/10.1029/2024WR038928


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



dat_obs_all<-read.delim("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/dHBV2.0_Song2025/1980-2020_daily_streamflow_gages/observations.txt",
                        header = FALSE,
                        sep = ",")



dat_sim_all<-read.delim("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/dHBV2.0_Song2025/1980-2020_daily_streamflow_gages/dhbv_uh_sim.txt",
                        header = FALSE,
                        sep = ",")

# Initialize dataframe with stations
stns<-read.delim("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/dHBV2.0_Song2025/1980-2020_daily_streamflow_gages/station_list.txt",sep = "s",
                 header = FALSE)
names(stns)<-"gauge_id"

tm<-read.delim("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/dHBV2.0_Song2025/1980-2020_daily_streamflow_gages/gagesii_time.txt",sep = "s",
               header = FALSE)%>%
  mutate(dt = ymd(substr(V1,1,10)))%>%
  pull(dt)


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

stnCompNSE<- vector(mode = "list", length = nrow(stns))


# loop through stations
for(it in 1:nrow(stns)){
  tic(sprintf("station %d",it))
  
  # Read in data 
  dat<-data.frame(dt = tm,
                  QObs = dat_obs_all[it,]%>%t()%>%as.numeric(),
                  QSim = dat_sim_all[it,]%>%t()%>%as.numeric())
  
  
  if(any(is.na(dat$QSim))){
    print(sprintf("THERE ARE NA VALUES IN STATION %s (number %d), SKIPPING",stns$ID[it],it))
    {next}
  }
  
  # calculate goodness of fit statistics
  perfMetrics<-calcAllGOF(dat,waterYearStart = 10)
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


# save data

stns$mdl<-"dHBV2.0_Song"


x<-bind_rows(stnPerf)%>%
  left_join(bind_rows(stnSeas),by = c( "gauge_id"))%>%
  left_join(bind_rows(stnCompNSE),by = c("gauge_id"),suffix = c("",".maxRun"))%>%
  left_join(stn_var)

stns<-left_join(stns,x)
write.csv(stns%>%st_drop_geometry(),"2.data/highLowBenchmarkGOF/GOF_Song_dHBV.csv")


