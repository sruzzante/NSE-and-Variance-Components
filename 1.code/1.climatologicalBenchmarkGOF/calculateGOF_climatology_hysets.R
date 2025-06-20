# climatological benchmarks and variance components
# Author: Sacha Ruzzante
# sachawruzzante@gmail.com
# Last Update: 2025-06-18


# Arsenault, R., Martel, J.-L., Brunet, F., Brissette, F., & Mai, J. (2023). Continuous streamflow prediction in ungauged basins: Long short-term memory neural networks clearly outperform traditional hydrological models. Hydrology and Earth System Sciences, 27(1), 139–157. https://doi.org/10.5194/hess-27-139-2023


library(ncdf4)
library(hydroGOF)
library(dplyr)
library(ggplot2)

#library(sf)
#library(tmap)
library(tictoc)
library(lubridate)
setwd("/home/ruzzante/projects/def-tgleeson/ruzzante/climatological_benchmarks/")
source("1.code/utils.R")

## HYSETS ########

nc_dat<-ncdf4::nc_open("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/hysets/HYSETS_2023_update_QC_stations.nc")

stns<-read.csv("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/hysets/HYSETS_watershed_properties.txt")%>%
  select(Watershed_ID,Source,Name,Official_ID,Hydrometric_station_latitude,Hydrometric_station_longitude)
# st_as_sf(coords = c("Hydrometric_station_longitude","Hydrometric_station_latitude"),remove = FALSE)#%>%
# filter((!Source == "USGS")|
#          Hydrometric_station_latitude>50)
unique(stns$Source)

# initialize dataframe for GOF statistics
stns$NSE<-NA
stns$KGE<-NA
stns$KGE_r<-NA
stns$KGE_a<-NA
stns$KGE_b<-NA
stns$N<-NA
# Q<-ncvar_get(nc_dat,"discharge",start = c(14200,1))
time_length <- nc_dat$dim$time$len
dt<-as.Date("1950-01-01")+
  days(ncvar_get(nc = nc_dat,
                 varid = "time"))

yr = year(dt)
mnth = month(dt)
dy = day(dt)
ydy = pmin(yday(dt),365)


# initialize dataframe for variance components

stn_var<-data.frame(Official_ID = stns$Official_ID,
                    varSeas_stl = NA,
                    varInterannual_stl = NA,
                    varRem_stl = NA,
                    varSeas_clas = NA,
                    varInterannual_clas = NA,
                    varRem_clas = NA,
                    varSeas_fourier = NA,
                    varInterannual_fourier = NA,
                    varRem_fourier = NA)


#Loop through stations and calculate statistics
for(it in 1:nrow(stns)){
  tic(sprintf("station %d",it))
  # gauge_id<-stns$gauge_id[it]
  
  
  # read in data
  dat<-data.frame(q= ncvar_get(nc = nc_dat,
                               varid = "discharge",
                               start = c(1, stns$Watershed_ID[it]), # Start at the first time step and the desired catchment
                               count = c(time_length, 1)),     # Read all time steps for the catchment)
                  dt = dt,
                  year = yr,
                  month = mnth,
                  day = dy,
                  yday = ydy
  )
  
  allNA<-
    dat%>%group_by(year)%>%
    summarize(allNA = sum(is.na(q))>=365)%>%
    filter(allNA)
  dat<-dat%>%filter(!year %in% allNA$year)
  
  
  
  # nc
  
  
  # initial check for time series length >=10 years
  if(nrow(dat)<(365*10)){next}
  if(sum(!is.na(dat$q))<=(365*10)){next}
  
  
  # ensure each day of the year has at least 10 years of data
  numYearsFull<-
    dat%>%
    group_by(yday)%>%
    summarise(notNA = sum(!is.na(q)))%>%
    ungroup()%>%
    summarize(allnotNA = min(notNA))
  
  if(numYearsFull<10){next}
  
  
  yrs<-unique(dat$year)
  
  # initialize empty data frame for joined test and climatology data
  dat_test2<-data.frame()
  
  # LOOCV routine for GOF
  for(it_yr in yrs){
    
    dat_test<-dat%>%filter(year  %in%  it_yr)
    dat_train = dat%>%filter(!(year  %in%  it_yr))
    
    # calculate climatology on training years
    dat_sum<-dat_train%>%
      group_by(yday)%>%
      summarize(q = mean(q,na.rm = TRUE))
    
    
    # join climatology to test year
    dat_test<-left_join(dat_test,dat_sum,by = "yday",suffix = c(".obs",".avg"))
    
    # append joined test & climatology data to data_test2
    dat_test2<-rbind(dat_test2,dat_test)
    
  }
  stns$NSE[it]<- hydroGOF::NSE(dat_test2$q.avg,dat_test2$q.obs,na.rm = TRUE)
  stns$KGE[it]<- hydroGOF::KGE(dat_test2$q.avg,dat_test2$q.obs,na.rm = TRUE)
  
  stns$KGE_r[it]<- hydroGOF::KGE(dat_test2$q.avg,dat_test2$q.obs,na.rm = TRUE, s  =c(1,0,0))
  stns$KGE_a[it]<- hydroGOF::KGE(dat_test2$q.avg,dat_test2$q.obs,na.rm = TRUE, s  =c(0,1,0))
  stns$KGE_b[it]<- hydroGOF::KGE(dat_test2$q.avg,dat_test2$q.obs,na.rm = TRUE, s  =c(0,0,1))
  
  stns$N[it]<-length(yrs)
  
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

write.csv(stn_var,"2.data/varComponents/hysets.csv",row.names = FALSE)

st_write(stns,"2.data/GOF/hysets.gpkg",append = FALSE)
