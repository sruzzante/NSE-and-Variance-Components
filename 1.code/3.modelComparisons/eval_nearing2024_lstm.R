# climatological benchmarks and variance components
# Author: Sacha Ruzzante
# sachawruzzante@gmail.com
# Last Update: 2025-06-18

# Evaluate LSTM from Nearing et al (2024)
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

stns<-read.csv(paste0(pth,"metadata/basin_attributes.csv"))

# get available gauge IDs
gauge_IDs<-list.files("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/Nearing2024/model_data/google/dual_lstm/hydrologically_separated/")%>%
  str_remove(".nc")

# 33 gauges seem to be missing from basin_attributes - that's OK
stns_missing<-read.csv(paste0(pth,"metadata/basin_attributes.csv"))%>%
  rename(ID = X)%>%
  mutate(gauge_ID = str_remove(ID,"GRDC_"))%>%
  st_as_sf(coords = c("longitude","latitude"),crs = "EPSG:4326")%>%
  select(ID, gauge_ID,calculated_drain_area)%>%
  filter(!ID %in% gauge_IDs)

# read in GRDC station spreadsheet
stns_grdc<-readxl::read_xlsx("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/GRDC_Stations.xlsx")

# join station IDs with GRDC attributes
stns<-read.csv(paste0(pth,"metadata/basin_attributes.csv"))%>%
  rename(ID = X)%>%
  
  mutate(gauge_ID = str_remove(ID,"GRDC_"))%>%
  left_join(stns_grdc%>%mutate(gauge_ID = as.character(grdc_no))%>%select(gauge_ID,country),
            by = c("gauge_ID"))%>%
  st_as_sf(coords = c("longitude","latitude"),crs = "EPSG:4326")%>%
  select(ID, gauge_ID,country,calculated_drain_area)%>%
  filter(ID %in% gauge_IDs)

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
  
  
  lstm_dat<-readNC(stns$ID[it])
  if(is.null(lstm_dat)){next}
  
  
  timeUnits<-lstm_dat$dim$time$units
  
  
  dat<-data.frame(
    QSim = ncvar_get(lstm_dat,"google_prediction",start = c(1,1),count = c(1,-1)),
    QObs = ncvar_get(lstm_dat,"observation",start = c(1,1),count = c(1,-1)),
    dt = ymd(timeUnits%>%str_remove("days since "))+ncvar_get(lstm_dat,"time")-1
  )
  
  if(any(is.na(dat$QSim))){
    print(sprintf("THERE ARE NA VALUES IN STATION %s (number %d), SKIPPING",stns$ID[it],it))
    {next}
  }
  
  
  # calculate goodness of fit statistics
  perfMetrics<-calcAllGOF(dat,waterYearStart = 10)
  perfMetrics$gauge_ID<-stns$gauge_ID[it]
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


stns<-left_join(stns,x)%>%
  st_drop_geometry()

stns$mdl<-"nearing-lstm"

write.csv(stns,"2.data/highLowBenchmarkGOF/GOF_Nearing_lstm.csv")

stn_result<-compareMetricsHighLow(stns)

write.csv(stn_result,"2.data/highLowBenchmarkGOF/highLowBenchmarkGOF_Nearing_lstm.csv")
