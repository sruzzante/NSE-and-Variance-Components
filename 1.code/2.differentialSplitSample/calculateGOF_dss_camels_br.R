# climatological benchmarks and variance components

# Author: Sacha Ruzzante
# sachawruzzante@gmail.com
# Last Update: 2025-06-18


# calculated the benchmark KGE and NSE based on differential split samples, for Camels -BR

# Camels-BR-v2 data:
# Chagas, V. B. P., Chaffe, P. L. B., Addor, N., Fan, F. M., Fleischmann, A. S., Paiva, R. C. D., & Siqueira, V. A. (2020). CAMELS-BR: Hydrometeorological time series and landscape attributes for 897 catchments in Brazil. Earth System Science Data, 12(3), 2075–2096. https://doi.org/10.5194/essd-12-2075-2020

# Updated v1.2 dataset: 
# Chagas, V. B. P., Chaffe, P. L. B., Addor, N., Fan, F. M., Fleischmann, A. S., Paiva, R. C. D., & Siqueira, V. A. (2025). CAMELS-BR: Hydrometeorological time series and landscape attributes for 897 catchments in Brazil - link to files. (Version 1.2) [Dataset]. Zenodo. https://doi.org/10.5281/zenodo.15025488



library(ncdf4)
library(hydroGOF)
library(dplyr)
library(ggplot2)
library(hydroGOF)
library(sf)
library(tmap)
library(tictoc)
library(lubridate)
library(stringr)
setwd("/home/ruzzante/projects/def-tgleeson/ruzzante/climatological_benchmarks/")


## CAMELS-BR ########

stns<-st_read("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-br/13_CAMELS_BR_gauge_location/location_gauges_streamflow.gpkg")

precip_IDs<-list.files("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-br/05_CAMELS_BR_precipitation/")%>%
  str_remove("_precipitation.txt")

# filter to only include 897 gauges with meteorological data
stns<-filter(stns, gauge_id %in% precip_IDs)

# initialize empty columns
stns$NSE_P<-NA
stns$KGE_P<-NA

stns$NSE_T<-NA
stns$KGE_T<-NA

stns$NSE_random<-NA
stns$KGE_random<-NA


for(it in 1:nrow(stns)){
  tic(sprintf("station %d",it))
  
  # read streamflow data
  dat<-read.delim(sprintf("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-br/02_CAMELS_BR_streamflow_all_catchments/%s_streamflow.txt",stns$gauge_id[it]),
                  sep = "")%>%
    mutate(dt = ymd(paste(year,month,day)))%>%
    filter(qual_flag == 1)
  
  names(dat)[4]<-"q"
  # We will need at least 20 years of data; If less, then skip to next iteration
  if(nrow(dat)<(365*20)){next}
  
  #empty dataframe
  dat_dt<-data.frame(dt = seq.Date(from = dat$dt[1],to = dat$dt[nrow(dat)], by = "1 day"))%>%
    mutate(year = year(dt),
           yday = pmin(yday(dt),365))
  
  # water year begins october 1
  dat_dt$wateryear<-dat_dt$year
  dat_dt$wateryear[dat_dt$month %in% (10:12)]<-
    dat_dt$wateryear[dat_dt$month %in% (10:12)]+1
  
  # load P data (ERA5-lAND)
  dat_p<-read.delim(sprintf("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-br/05_CAMELS_BR_precipitation/%s_precipitation.txt",stns$gauge_id[it]),
                    sep = "")%>%
    mutate(dt = ymd(paste(year,month,day)),
           precip = p_era5land)
  
  # load T data (ERA5-lAND)
  dat_t<-read.delim(sprintf("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-br/09_CAMELS_BR_temperature/%s_temperature.txt",stns$gauge_id[it]),
                    sep = "")%>%
    mutate(dt = ymd(paste(year,month,day)),
           temp = tmean_era5land)
  
  # join all data
  dat<-left_join(dat_dt,
                 dat%>%select(dt,q))%>%
    left_join(dat_p%>%select(dt,precip))%>%
    left_join(dat_t%>%select(dt,temp))
  
  # summarize T and P data by year
  yearlySummary<-
    dat%>%
    group_by(wateryear)%>%
    summarize(nNA = sum(is.na(q)| is.na(precip)),
              N = n(),
              meanT = mean(temp,na.rm = TRUE),
              sumP = sum(precip,na.rm = TRUE))%>%
    filter((N-nNA)>350)%>% # require at least 350 days of data
    
    mutate(period_P = as.numeric(sumP<median(sumP))+1,   #label periods (1 or 2)
           period_T = as.numeric(meanT<median(meanT))+1)
  
  
  dat<-filter(dat,
              wateryear %in% yearlySummary$wateryear)
  
  
  dat<-left_join(dat,yearlySummary%>%select(wateryear,period_P,period_T))
  
  # require at least 20 years of (mostly) complete data
  if(nrow(yearlySummary)<(20)){next}
  
  
  NSE_P<-c()
  KGE_P<-c()
  # loop over periods 1 and 2
  for(it_P in 1:2){
    
    dat_test<-dat%>%filter(period_P  == it_P)
    dat_train = dat%>%filter(period_P != it_P)
    
    dat_sum<-dat_train%>%
      
      group_by(yday)%>%
      summarize(q = mean(q,na.rm = TRUE))
    
    
    dat_test<-left_join(dat_test,dat_sum,by = "yday",suffix = c(".obs",".avg"))
    
    NSE_P[it_P] <- hydroGOF::NSE(dat_test$q.avg,dat_test$q.obs,na.rm = TRUE)
    KGE_P[it_P] <- hydroGOF::KGE(dat_test$q.avg,dat_test$q.obs,na.rm = TRUE)
    
  }
  
  # take median NSE/KGE of the two periods
  stns$NSE_P[it]<- median(NSE_P)
  stns$KGE_P[it]<- median(KGE_P)
  
  
  NSE_T<-c()
  KGE_T<-c()
  # loop over periods 1 and 2
  for(it_T in 1:2){
    
    dat_test<-dat%>%filter(period_T  == it_T)
    dat_train = dat%>%filter(period_T != it_T)
    
    dat_sum<-dat_train%>%
      group_by(yday)%>%
      summarize(q = mean(q,na.rm = TRUE))
    
    dat_test<-left_join(dat_test,dat_sum,by = "yday",suffix = c(".obs",".avg"))
    
    NSE_T[it_T] <- hydroGOF::NSE(dat_test$q.avg,dat_test$q.obs,na.rm = TRUE)
    KGE_T[it_T] <- hydroGOF::KGE(dat_test$q.avg,dat_test$q.obs,na.rm = TRUE)
    
  }
  # take median NSE/KGE of the two periods
  stns$NSE_T[it]<- median(NSE_T)
  stns$KGE_T[it]<- median(KGE_T)
  
  
  # 10 random splits of data
  NSE_random<-c()
  KGE_random<-c()
  set.seed(1)
  # loop over 10 random splits
  for(it_random in 1:10){
    # randomize periods
    yearlySummary$period_rand<-sample(yearlySummary$period_P,
                                      size = nrow(yearlySummary),
                                      replace = FALSE)
    dat_x<-left_join(dat,yearlySummary%>%select(wateryear,period_rand))
    
    # loop over 2 periods
    for(it_T in 1:2){
      
      dat_test<-dat_x%>%filter(period_rand  == it_T)
      dat_train = dat_x%>%filter(period_rand != it_T)
      
      dat_sum<-dat_train%>%
        group_by(yday)%>%
        summarize(q = mean(q,na.rm = TRUE))
      
      
      
      dat_test<-left_join(dat_test,dat_sum,by = "yday",suffix = c(".obs",".avg"))
      
      NSE_random <- c(NSE_random,
                      hydroGOF::NSE(dat_test$q.avg,dat_test$q.obs,na.rm = TRUE))
      KGE_random <- c(KGE_random,
                      hydroGOF::KGE(dat_test$q.avg,dat_test$q.obs,na.rm = TRUE))
      
    }
    
    
    
  }
  # take median NSE/KGE over 20 splits & periods
  stns$NSE_random  [it]<- NSE_random%>%group_by(it_split)%>%summarise(NSE = mean(NSE))%>%pull(NSE)%>%median()
  
  stns$KGE_random  [it]<- KGE_random%>%group_by(it_split)%>%summarise(KGE = mean(KGE))%>%pull(KGE)%>%median()
  
  toc()
}

st_write(stns,"2.data/GOF/camels_br_dss.gpkg",append = FALSE)


