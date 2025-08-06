# climatological benchmarks and variance components
# Author: Sacha Ruzzante
# sachawruzzante@gmail.com
# Last Update: 2025-08-06

# Evaluate COSERO model from Klinger et al (2021)
# Klingler, C., Schulz, K., and Herrnegger, M.: LamaH-CE: LArge-SaMple DAta for
# Hydrology and Environmental Sciences for Central Europe, Earth Syst. Sci. Data, 
# 13, 4529–4565, https://doi.org/10.5194/essd-13-4529-2021, 2021.



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

stns<-st_read("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/lamah-ce/D_gauges/3_shapefiles/Gauges.shp")

# initialize list to store performance metrics
stnPerf<-vector(mode = "list", length = nrow(stns))

# initialize list to store seasonality metrics
stnSeas<-vector(mode = "list", length = nrow(stns))

# initialize dataframe to store variance components
stn_var<-data.frame(ID = stns$govnr,
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
  
  
  
  if(!file.exists(paste0(
    "../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/lamah-ce/F_hydrol_model/2_timeseries/ID_",
    stns$ID[it],".csv"
  ))){next}
  
  dat<-read.delim(paste0(
    "../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/lamah-ce/F_hydrol_model/2_timeseries/ID_",
    stns$ID[it],".csv"
  ),
  sep = ";")%>%
    
    mutate(QSim = Qsim_A,
           QObs = Qobs,
           dt = ymd(paste(YYYY,MM,DD)),
           year = year(dt),
           yday = pmin(yday(dt),365))
  
  
  

  
  if(any(is.na(dat$QSim))){
    print(sprintf("THERE ARE NA VALUES IN STATION %s (number %d), SKIPPING",stns$ID[it],it))
    {next}
  }
  
  # calculate goodness of fit statistics
  perfMetrics<-calcAllGOF(dat,waterYearStart = 10)
  perfMetrics$ID<-stns$govnr[it]
  stnPerf[[it]]<-data.frame(perfMetrics)
  
  # calculate seasonality statistics
  seas<-calcSeasonality(dat)
  stnSeas[[it]] = data.frame(ID = stns$govnr[it],
                             CoV = seas["CoV"],
                             QCI = seas["QCI"])
  
  # Goodness of fit of variance components
  NSE_comps<- gof_components(dat)
  NSE_comps$ID<-stns$govnr[it]
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
stns$mdl<-"cosero-lamah-ce"


x<-bind_rows(stnPerf)%>%
  left_join(bind_rows(stnSeas),by = c("ID" ))%>%
  left_join(bind_rows(stnCompNSE),by = c("ID"),suffix = c("",".maxRun"))%>%
  left_join(stn_var)

stns<-left_join(stns,x,by = c("govnr" = "ID"))

stns<-st_drop_geometry(stns)

write.csv(stns,"2.data/highLowBenchmarkGOF/GOF_lamah_ce_cosero.csv")


