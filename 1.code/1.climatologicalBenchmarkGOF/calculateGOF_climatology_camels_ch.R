# climatological benchmarks and variance components
# Author: Sacha Ruzzante
# sachawruzzante@gmail.com
# Last Update: 2025-06-18


# Höge, M., Kauzlaric, M., Siber, R., Schönenberger, U., Horton, P., Schwanbeck, J., Floriancic, M. G., Viviroli, D., Wilhelm, S., Sikorska-Senoner, A. E., Addor, N., Brunner, M., Pool, S., Zappa, M., & Fenicia, F. (2023). CAMELS-CH: Hydro-meteorological time series and landscape attributes for 331 catchments in hydrologic Switzerland. Earth System Science Data, 15(12), 5755–5784. 




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
source("1.code/utils.R")


## CAMELS-CH ########

stns<-st_read("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-ch/camels_ch/catchment_delineations/CAMELS_CH_gauging_stations.shp")%>%
  filter(type == "stream")


# initialize dataframe for GOF statistics
stns$NSE<-NA
stns$KGE<-NA
stns$KGE_r<-NA
stns$KGE_a<-NA
stns$KGE_b<-NA

stns$N<-NA

# initialize dataframe for variance components
stn_var<-data.frame(gauge_id = stns$gauge_id,
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
  
  dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-ch/camels_ch/timeseries/observation_based/CAMELS_CH_obs_based_",
                         stns$gauge_id[it],".csv"),
                  # skip  =7,
                  sep  =",")%>%
    mutate(dt = ymd(date),
           q = discharge_vol.m3.s.,
           yday = pmin(yday(dt),365),
           year = year(dt))
  

  dat<-filter(dat,!is.na(dat$q))

  
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
write.csv(stn_var,"2.data/varComponents/camels_ch.csv",row.names = FALSE)

ggplot(stns,aes(x= NSE))+stat_ecdf()+
  scale_x_continuous(limits = c(-1,1))


ggplot(stns,aes(x= KGE))+stat_ecdf()+
  scale_x_continuous(limits = c(-1,1))
median(stns$NSE,na.rm = TRUE)


tmap_mode("view")
tm_shape(stns%>%filter(!is.na(NSE)))+tm_dots(col = "NSE",breaks = c(-Inf,seq(-1,1,0.2)),
                                             palette = "RdBu")


st_write(stns,"2.data/GOF/camels_ch.gpkg",append = FALSE)
