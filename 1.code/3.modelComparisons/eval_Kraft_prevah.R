# climatological benchmarks and variance components
# Author: Sacha Ruzzante
# sachawruzzante@gmail.com
# Last Update: 2025-06-18

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
source('1.code/utils.R')

# load ncdf data
dat_nc<-ncdf4::nc_open("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/Kraft2025/chrun_obs_grid.nc")

dates<-ymd("1961-01-01")+
  ncvar_get(dat_nc,
            "time")
stns<-data.frame(
  ID =  ncvar_get(dat_nc,"station")
)




# initialize list to store performance metrics
stnPerf<-vector(mode = "list", length = nrow(stns))

# initialize list to store seasonality metrics
stnSeas<-vector(mode = "list", length = nrow(stns))

# initialize dataframe to store variance components
stn_var<-data.frame(ID = stns$ID,
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


for(it in 1:nrow(stns)){
  tic(sprintf("station %d",it))
  
  # read data in
  dat<-data.frame(
    dt = dates,
    QObs = ncvar_get(dat_nc,"Qmm",start = c(1,it),count = c(22645,1) ),
    QSim =  ncvar_get(dat_nc,"Qmm_prevah",start = c(1,it),count = c(22645,1) )
  )
  # Prevah simulations are only from 1981-2022
  dat<-dat%>%filter(year(dt)>=1981&
                      year(dt)<=2022)
  
  if(any(is.na(dat$QSim))){
    print(sprintf("THERE ARE NA VALUES IN STATION %s (number %d), SKIPPING",stns$ID[it],it))
    {next}
  }
  
  # calculate goodness of fit statistics
  perfMetrics<-calcAllGOF(dat,waterYearStart = 10)
  perfMetrics$ID<-stns$ID[it]
  stnPerf[[it]]<-data.frame(perfMetrics)
  
  # calculate seasonality statistics
  seas<-calcSeasonality(dat)
  stnSeas[[it]] = data.frame(gauge_id = stns$ID[it],
                             CoV = seas["CoV"],
                             QCI = seas["QCI"])
  
  
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
stns$mdl<-"prevah-Kraft2025"

x<-bind_rows(stnPerf)%>%
  left_join(bind_rows(stnSeas),by = c("ID" = "gauge_id"))%>%
  left_join(stn_var)



stns<-left_join(stns,x)
write.csv(stns,"2.data/highLowBenchmarkGOF/GOF_Kraft_Prevah.csv")


