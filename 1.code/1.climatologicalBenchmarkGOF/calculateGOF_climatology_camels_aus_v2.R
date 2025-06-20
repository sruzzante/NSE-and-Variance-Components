# climatological benchmarks and variance components
# Author: Sacha Ruzzante
# sachawruzzante@gmail.com
# Last Update: 2025-06-18

# Camels-AUS-v2 data:
# Fowler, K. J. A., Zhang, Z., & Hou, X. (2024). CAMELS-AUS v2: Updated hydrometeorological timeseries and landscape attributes for an enlarged set of catchments in Australia. Earth System Science Data Discussions, 1–21. https://doi.org/10.5194/essd-2024-263


library(ncdf4)
library(hydroGOF)
library(dplyr)
library(ggplot2)
library(hydroGOF)
library(sf)
#library(tmap)
library(tictoc)
library(lubridate)
setwd("/home/ruzzante/projects/def-tgleeson/ruzzante/climatological_benchmarks/")
source("1.code/utils.R")


## CAMELS-AUS-V2 ########

stns<-st_read("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-aus-v2/02_location_boundary_area/shp/CAMELS_AUS_v2_BasinOutlets_adopted.shp")


# initialize dataframe for GOF statistics
stns$NSE<-NA
stns$KGE<-NA
stns$KGE_r<-NA
stns$KGE_a<-NA
stns$KGE_b<-NA

stns$N<-NA


# read in data

streamflow_dat<-read.csv("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-aus-v2/streamflow_MLd.csv")
stns$station_id[stns$station_id=="227225A"] = "227225"



# initialize dataframe for variance components
stn_var<-data.frame(station_id = stns$station_id,
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
  
  
  ID_str<-stns$station_id[it]
  if((substr(ID_str,1,1))%in% as.character(0:9)){ID_str = paste0("X",ID_str)}
  dat<-streamflow_dat%>%
    select(year,month,day,all_of(ID_str))%>%
    mutate(dt = ymd(paste(year,month,day)),
           yday = pmin(yday(dt),365))
  
  names(dat)[4]<-"q"
  
  dat<-filter(dat,!is.na(dat$q))
  
  
  dat$q[dat$q== -99.99]<-NA
  
  
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
  if(numYearsFull$allnotNA<10){next}
  
  
  yrs<-unique(dat$year)
  dat_test2<-data.frame()
  
  # LOOCV routine for GOF
  for(it_yr in yrs){
    
    dat_test<-dat%>%filter(year  %in%  it_yr)
    dat_train = dat%>%filter(!(year  %in%  it_yr))
    
    dat_sum<-dat_train%>%
      group_by(yday)%>%
      summarize(q = mean(q,na.rm = TRUE))
    
    
    
    dat_test<-left_join(dat_test,dat_sum,by = "yday",suffix = c(".obs",".avg"))
    
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

write.csv(stn_var,"2.data/varComponents/camels_AUS.csv",row.names = FALSE)



ggplot(stns,aes(x= NSE))+stat_ecdf()+
  scale_x_continuous(limits = c(-1,1))


ggplot(stns,aes(x= KGE))+stat_ecdf()+
  scale_x_continuous(limits = c(-1,1))
median(stns$NSE,na.rm = TRUE)

tmap_mode("view")
tm_shape(stns)+tm_dots(col = "NSE",breaks = c(-Inf,seq(-1,1,0.2)),
                       palette = "RdBu")

stns<-stns%>%st_transform("EPSG:4326")

st_write(stns,"2.data/GOF/camels_AUS.gpkg",append = FALSE)
