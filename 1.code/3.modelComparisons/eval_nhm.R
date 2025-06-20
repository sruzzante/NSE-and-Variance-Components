# climatological benchmarks and variance components
# Author: Sacha Ruzzante
# sachawruzzante@gmail.com
# Last Update: 2025-06-18

# Evaluate National Hydrologic Model
# Regan, R. S., Juracek, K. E., Hay, L. E., Markstrom, S. L., Viger, R. J., Driscoll, J. M., LaFontaine, J. H., & Norton, P. A. (2019). The U. S. Geological Survey National Hydrologic Model infrastructure: Rationale, description, and application of a watershed-scale model for the conterminous United States. Environmental Modelling & Software, 111, 192–203. https://doi.org/10.1016/j.envsoft.2018.09.023

# We will only evaluate the performance at 'reference' gages as defined by Falcone (2011)
# Falcone, J. A. (2011). GAGES-II: Geospatial Attributes of Gages for Evaluating Streamflow. U.S. Geological Survey. https://doi.org/10.3133/70046617


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

source('1.code/utils.R')

# read station metadata
stns<-read.csv("../DATA/1.Spatial_data/regional/USA/sw_surfacewater_streamflow_runoff_river_network_waterstress/sw1_USGS_NHM/GAGE_ID_crosswalk/poi_gage_id.csv",colClasses = "character")%>%
  cbind(read.csv("../DATA/1.Spatial_data/regional/USA/sw_surfacewater_streamflow_runoff_river_network_waterstress/sw1_USGS_NHM/GAGE_ID_crosswalk/poi_gage_segment.csv")%>%select(poi_gage_segment))

# read GAGES data (Falcone, 2011)
GAGES<-read.csv("../DATA/1.Spatial_data/regional/North_America/sw_surfacewater_streamflow_runoff_river_network_waterstress/sw1_GAGES-II/basinchar_and_report_sept_2011/conterm_bas_classif.txt",
                colClasses = "character"
)


sum(stns$poi_gage_id %in% GAGES$STAID)

# filter to reference gages
stns<-stns%>%
  filter(poi_gage_id %in% GAGES$STAID[GAGES$CLASS == "Ref"])

# load ncdf data
nc_data<-ncdf4::nc_open("../DATA/1.Spatial_data/regional/USA/sw_surfacewater_streamflow_runoff_river_network_waterstress/sw1_USGS_NHM/netcdf/seg_outflow.nc")

t<-ymd("1980-10-01")+ncvar_get(nc_data,"time")
seg_ID<-ncvar_get(nc_data,"segment")



# initialize list to store performance metrics
stnPerf<-vector(mode = "list", length = nrow(stns))

# initialize list to store seasonality metrics
stnSeas<-vector(mode = "list", length = nrow(stns))


# initialize dataframe to store variance components

stn_var<-data.frame(poi_gage_id = stns$poi_gage_id,
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
  if(stns$poi_gage_segment[it]==0){next}
  
  # Download data if necessary, or read from file
  if(!file.exists(sprintf("2.data/NHM_Q/%s.RDS",stns$poi_gage_id[it]))){
    streamflowData<-dataRetrieval::readNWISdv(stns$poi_gage_id_0[it],
                                              parameterCd = "00060",
                                              startDate = "1980-10-01",
                                              endDate = "2016-12-31")
    saveRDS(streamflowData, paste0("2.data/NHM_Q/",stns$poi_gage_id_0[it],".RDS"))
  }else{
    # read observed data 
    OBSDATA<-readRDS(sprintf("2.data/NHM_Q/%s.RDS",stns$poi_gage_id[it]))
  }
  
  
  
  
  # if <10 years, skip
  if(nrow(OBSDATA)<=3652){next}
  
  # some columns are named differently
  if("X_00060_00003"%in% names(OBSDATA)){
    "No need to change column names"
  }else if("X_PUBLISHED_00060_00003" %in% names(OBSDATA)){
    OBSDATA$X_00060_00003<-OBSDATA$X_PUBLISHED_00060_00003
  }else{
    next
  }
  
  OBSDATA<-OBSDATA%>%
    mutate(dt = ymd(Date),
           QObs = X_00060_00003)%>%
    select(dt,QObs)
  
  # join simulated and observed data
  
  dat<-data.frame(
    QSim = ncvar_get(nc_data,"seg_outflow",start = c(stns$poi_gage_segment[it],1),count = c(1,-1)),
    
    dt = t
  )%>%
    left_join(OBSDATA)
  
  
  if(any(is.na(dat$QSim))){
    print(sprintf("THERE ARE NA VALUES IN STATION %s (number %d), SKIPPING",stns$ID[it],it))
    {next}
  }
  
  
  # calculate goodness of fit statistics
  perfMetrics<-calcAllGOF(dat,waterYearStart = 10)
  perfMetrics$poi_gage_id<-stns$poi_gage_id[it]
  stnPerf[[it]]<-data.frame(perfMetrics)
  
  # calculate seasonality statistics
  seas<-calcSeasonality(dat)
  stnSeas[[it]] = data.frame(poi_gage_id = stns$poi_gage_id[it],
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



stns$mdl<-"nhm-ref"


x<-bind_rows(stnPerf)%>%
  left_join(bind_rows(stnSeas))%>%
  left_join(stn_var)


stns<-left_join(stns,x)
write.csv(stns,"2.data/highLowBenchmarkGOF/GOF_NHM_ref.csv")

