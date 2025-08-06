# climatological benchmarks and variance components
# Author: Sacha Ruzzante
# sachawruzzante@gmail.com
# Last Update: 2025-08-06

# Evaluate Arsenault LSTM
# Arsenault, R., Martel, J.-L., Brunet, F., Brissette, F., & Mai, J. (2023). Continuous streamflow prediction in ungauged basins: Long short-term memory neural networks clearly outperform traditional hydrological models. Hydrology and Earth System Sciences, 27(1), 139–157. https://doi.org/10.5194/hess-27-139-2023

# Note that you must download and rerun the codes at https://osf.io/3s2pq/ to get the simulated time series data (HYSETS_NENA_Q_basin_*.csv)

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


# observed discharge data
dat_nc<-ncdf4::nc_open("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/Arsenault2023/Data/NortheastNA_regionalization_data.nc")
x<-st_read("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/Arsenault2023/Shapefiles/shapefiles_148_hysets_regionalization_paper.shp")
dates<-ymd("0000-01-01")+
  ncvar_get(dat_nc,
            "time")-1
ncvar_get(dat_nc,
          "watersheds")

stns<-data.frame(ID = 0:147,
                 lat = ncvar_get(dat_nc,"lat"),
                 lon = ncvar_get(dat_nc,"lon"),
                 area = ncvar_get(dat_nc,"area"),
                 elevation = ncvar_get(dat_nc,"elevation")
                 
)%>%
  st_as_sf(coords = c("lon","lat"),remove = FALSE,crs = "EPSG:4326")

plot(stns$area,x$Area)
stns$ID<-x$OfficialID

# stns<-left_join(stns%>%mutate(area = round(area)),
#                 x%>%select(Area,OfficialID)%>%st_drop_geometry()%>%mutate(Area = round(Area)),
#                 by = c("area" = "Area"))


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
# initialize list to store goodness of fit for variance components
stnCompNSE<- vector(mode = "list", length = nrow(stns))


# loop through stations

for(it in 1:nrow(stns)){
  tic(sprintf("station %d",it))
  
  # read simulated discharge for station
  dat<-read.csv(sprintf("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/Arsenault2023/model_predictions/HYSETS_NENA_Q_basin_%s.csv",it-1))%>%
    mutate(dt = seq.Date(from = ymd("1981-12-31"),
                         to = ymd("2018-12-29"),
                         by = "1 day")
    )
  
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
  stnSeas[[it]] = data.frame(ID = stns$ID[it],
                             CoV = seas["CoV"],
                             QCI = seas["QCI"])
  
  
  # Goodness of fit of variance components
  NSE_comps<- gof_components(dat)
  NSE_comps$ID<-stns$ID[it]
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
stns$mdl<-"lstm-arsenault2022"

x<-bind_rows(stnPerf)%>%
  left_join(bind_rows(stnSeas))%>%
  left_join(bind_rows(stnCompNSE),by = c("ID"),suffix = c("",".maxRun"))%>%



  left_join(stn_var)

stns<-left_join(stns,x)
write.csv(stns%>%st_drop_geometry(),"2.data/highLowBenchmarkGOF/GOF_Arsenault2022.csv",row.names = FALSE)
