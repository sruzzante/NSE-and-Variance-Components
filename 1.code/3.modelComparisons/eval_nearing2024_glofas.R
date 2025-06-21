# climatological benchmarks and variance components
# Author: Sacha Ruzzante
# sachawruzzante@gmail.com
# Last Update: 2025-06-18

# Evaluate GLOFAS from Nearing et al (2024)
# Nearing, G., Cohen, D., Dube, V., Gauch, M., Gilon, O., Harrigan, S., Hassidim, A., Klotz, D., Kratzert, F., Metzger, A., Nevo, S., Pappenberger, F., Prudhomme, C., Shalev, G., Shenzis, S., Tekalign, T. Y., Weitzner, D., & Matias, Y. (2024). Global prediction of extreme floods in ungauged watersheds. Nature, 627(8004), 559–563. https://doi.org/10.1038/s41586-024-07145-1

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

pth = "../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/Nearing2024/"

stn_meta<-read.csv(paste0(pth,"metadata/basin_attributes.csv"))%>%
  mutate(gauge_ID = as.numeric(str_remove(X,"GRDC_")))

gauge_IDs<-list.files("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/Nearing2024/model_data/google/dual_lstm/hydrologically_separated//")%>%
  str_remove(".nc")%>%
  str_remove("GRDC_")
dat_nc<-nc_open("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/Nearing2024/model_data/GRDCstattions_GloFASv40/dis24h_GLOFAS4.0_3arcmin_197901-202212_statsgoogle20230918.nc")

timeUnits_glofas<-dat_nc$dim$time$units

stns_grdc<-readxl::read_xlsx("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/GRDC_Stations.xlsx")

stns<-data.frame(gauge_ID = ncvar_get(dat_nc,"statid"),
                 latitude = ncvar_get(dat_nc,"statlat"),
                 longitude = ncvar_get(dat_nc,"statlon"))%>%
  left_join(stn_meta%>%select(gauge_ID,calculated_drain_area))%>%
  left_join(stns_grdc%>%mutate(gauge_ID = (grdc_no))%>%select(gauge_ID,country),
            by = c("gauge_ID"))%>%
  # st_as_sf(coords = c("longitude","latitude"),crs = "EPSG:4326")%>%
  select(gauge_ID,country,calculated_drain_area,longitude,latitude)%>%
  st_as_sf(coords = c("longitude","latitude"),crs = "EPSG:4326")%>%
  
  filter(gauge_ID %in% gauge_IDs)


# function to safely read streamflow data 
readNC <- function(stnID) {
  tryCatch(
    {
      
      ncdf4::nc_open(paste0(pth,"model_data/google/dual_lstm/hydrologically_separated/",stnID,".nc"))
    },
    error = function(cond) {
      message(paste("Couldn't load station file:", stnID))
      message("Here's the original error message:")
      message(conditionMessage(cond))
      # Choose a return value in case of error
      NULL
    }
    
  )
}

# initialize list to store performance metrics
stnPerf<-vector(mode = "list", length = nrow(stns))

# initialize list to store seasonality metrics
stnSeas<-vector(mode = "list", length = nrow(stns))

# initialize dataframe to store variance components
stn_var<-data.frame(gauge_ID = stns$gauge_ID,
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
  
  # The observed data are in the 'lstm' file
  lstm_dat<-readNC(paste0("GRDC_",stns$gauge_ID[it]))
  if(is.null(lstm_dat)){next}
  
  timeUnits<-lstm_dat$dim$time$units
  
  dat_obs<-data.frame(
    
    QObs = ncvar_get(lstm_dat,"observation",start = c(1,1),count = c(1,-1)),
    dt = ymd(timeUnits%>%str_remove("days since "))+ncvar_get(lstm_dat,"time")-1
  )
  
  # read and join the GLOFAS simulated data
  dat<-data.frame(
    QSim = ncvar_get(dat_nc,"dis",start = c(it,1),count = c(1,-1))*86.4/stns$calculated_drain_area[it], # convert to mm/d
    
    dt = ymd("1978-01-02")+ncvar_get(dat_nc,"time")
  )
  dat<-left_join(dat,dat_obs)
  
  if(any(is.na(dat$QSim))){
    print(sprintf("THERE ARE NA VALUES IN STATION %s (number %d), SKIPPING",stns$gauge_ID[it],it))
    {next}
  }
  
  
  # calculate goodness of fit statistics
  perfMetrics<-calcAllGOF(dat,waterYearStart = 10)
  perfMetrics$gauge_ID<-(stns$gauge_ID[it])
  stnPerf[[it]]<-data.frame(perfMetrics)
  
  # calculate seasonality statistics
  seas<-calcSeasonality(dat)
  stnSeas[[it]] = data.frame(gauge_ID = stns$gauge_ID[it],
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
  
  
  
  nc_close(lstm_dat)
  
  toc()
}


x<-bind_rows(stnPerf)%>%
  left_join(bind_rows(stnSeas))%>%
  left_join(stn_var)

stns<-left_join(stns,x)%>%st_drop_geometry()

stns$mdl<-"nearing-glofas"

write.csv(stns%>%st_drop_geometry(),"2.data/highLowBenchmarkGOF/GOF_Nearing_glofas.csv")


