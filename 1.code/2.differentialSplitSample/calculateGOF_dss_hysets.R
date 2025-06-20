# climatological benchmarks and variance components

# Author: Sacha Ruzzante
# sachawruzzante@gmail.com
# Last Update: 2025-06-18


# calculated the benchmark KGE and NSE based on differential split samples, for HYSETS
# Arsenault, R., Brissette, F., Martel, J.-L., Troin, M., Lévesque, G., Davidson-Chaput, J., Gonzalez, M. C., Ameli, A., & Poulin, A. (2020). A comprehensive, multisource database for hydrometeorological modeling of 14,425 North American watersheds. Scientific Data, 7(1), 243. https://doi.org/10.1038/s41597-020-00583-2


# Since this is a big dataset, it is set up to accept 2 arguments: the number of splits into which to partition the dataset, and the particular split to evaluate. This can be called from a batch script


args = commandArgs(trailingOnly=TRUE)
it_split = 1
numsplit = 1

# test if there is at least two arguments: if not, continue without splitting
if (length(args)<2) {
  print("No arguments supplied, continuing without splitting data")
}else{
  it_split = as.integer(args[1])
  numsplit = as.integer(args[2])
}




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
library(rnaturalearth)
setwd("/home/ruzzante/projects/def-tgleeson/ruzzante/climatological_benchmarks/")


## hysets ########
########

# load streamflow data
nc_dat<-ncdf4::nc_open("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/hysets/HYSETS_2023_update_QC_stations.nc")

# load meteorological data
nc_dat_forcing<-ncdf4::nc_open("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/hysets/HYSETS_2023_update_monthly_meteorological_data.nc")

# station metadata
stns<-read.csv("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/hysets/HYSETS_watershed_properties.txt")%>%
  select(Watershed_ID,Source,Name,Official_ID,Hydrometric_station_latitude,Hydrometric_station_longitude)%>%
  st_as_sf(coords = c("Hydrometric_station_longitude","Hydrometric_station_latitude"),remove = FALSE,crs  = "EPSG:4326")#%>%


# partition the stations
if(numsplit>1){
  x<-split(1:nrow(stns),cut(seq_along(1:nrow(stns)), breaks = numsplit, labels = FALSE))
  stns<-stns[x[[it_split]],]
  
}


time_length <- nc_dat$dim$time$len
dt<-as.Date("1950-01-01")+
  days(ncvar_get(nc = nc_dat,
                 varid = "time"))

yr = year(dt)
mnth = month(dt)
dy = day(dt)
ydy = pmin(yday(dt),365)

# water year
wyr = yr
wyr[mnth %in% 10:12]<-wyr[mnth %in% 10:12]+1


time_length_for <- nc_dat_forcing$dim$time$len

dt_for<-as.Date("1950-01-01")+
  months(ncvar_get(nc = nc_dat_forcing,
                   varid = "time")-1)
year_for<-year(dt_for)
month_for = month(dt_for)

stns$NSE_P<-NA
stns$KGE_P<-NA

stns$NSE_T<-NA
stns$KGE_T<-NA

stns$NSE_random<-NA
stns$KGE_random<-NA

# loop over stations
for(it in 1:nrow(stns)){
  tic(sprintf("station %d",it))
 
  
  dat<-data.frame(q= ncvar_get(nc = nc_dat,
                               varid = "discharge",
                               start = c(1, stns$Watershed_ID[it]), # Start at the first time step and the desired catchment
                               count = c(time_length, 1)),     # Read all time steps for the catchment)
                  dt = dt,
                  year = yr,
                  wateryear = wyr,
                  month = mnth,
                  day = dy,
                  yday = ydy
  )
  dat_forcing<-data.frame(pr = ncvar_get(nc = nc_dat_forcing,
                                         varid = "ERA5Land_pr_month",
                                         start = c(1, stns$Watershed_ID[it]), # Start at the first time step and the desired catchment
                                         count = c(time_length_for, 1)),
                          tasmax = ncvar_get(nc = nc_dat_forcing,
                                             varid = "ERA5Land_tasmax_month",
                                             start = c(1, stns$Watershed_ID[it]), # Start at the first time step and the desired catchment
                                             count = c(time_length_for, 1)),
                          tasmin = ncvar_get(nc = nc_dat_forcing,
                                             varid = "ERA5Land_tasmin_month",
                                             start = c(1, stns$Watershed_ID[it]), # Start at the first time step and the desired catchment
                                             count = c(time_length_for, 1)),
                          dt = dt_for,
                          year = year_for,
                          month = month_for)%>%
    mutate(tasmean = (tasmax+tasmin)/2, # mean t is the average of max and min temperatures
           days_in_month = days_in_month(dt),
           pr = pr/days_in_month)%>%
    filter(year>=1979)
  
  # join streamflow and meteorological data
  dat<-left_join(dat,dat_forcing%>%
                   select(year,month,pr, tasmean))
  
  # We will need at least 20 years of data; If less, then skip to next iteration
  if(nrow(dat)<(365*20)){next}
  
  # summarize T and P data by year
  yearlySummary<-
    dat%>%
    group_by(wateryear)%>%
    summarize(nNA = sum(is.na(q)| is.na(pr)),
              N = n(),
              meanT = mean(tasmean,na.rm = TRUE),
              sumP = sum(pr,na.rm = TRUE))%>%
    filter((N-nNA)>350)%>%
    mutate(period_P = as.numeric(sumP<median(sumP))+1,
           period_T = as.numeric(meanT<median(meanT))+1)
  
  # keep only complete years
  dat<-filter(dat,
              wateryear %in% yearlySummary$wateryear)
  
  
  dat<-left_join(dat,yearlySummary%>%select(wateryear,period_P,period_T))
  
  
  # require at least 20 years of (mostly) complete data
  if(nrow(yearlySummary)<(20)){next}
  
  
  NSE_P<-c()
  KGE_P<-c()
  
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
  stns$NSE_P[it]<- median(NSE_P)
  stns$KGE_P[it]<- median(KGE_P)
  
  NSE_T<-c()
  KGE_T<-c()
  
  for(it_T in 1:2){
    
    dat_test<-dat%>%filter(period_T  == it_T)
    dat_train = dat%>%filter(period_T != it_T)
    
    dat_sum<-dat_train%>%
      # mutate(q[year == it_year] = NA)%>%
      group_by(yday)%>%
      summarize(q = mean(q,na.rm = TRUE))
    
    
    
    dat_test<-left_join(dat_test,dat_sum,by = "yday",suffix = c(".obs",".avg"))
    
    NSE_T[it_T] <- hydroGOF::NSE(dat_test$q.avg,dat_test$q.obs,na.rm = TRUE)
    KGE_T[it_T] <- hydroGOF::KGE(dat_test$q.avg,dat_test$q.obs,na.rm = TRUE)
    
  }
  stns$NSE_T[it]<- median(NSE_T)
  stns$KGE_T[it]<- median(KGE_T)
  
  
  # 10 random splits
  NSE_random<-data.frame()
  KGE_random<-data.frame()
  set.seed(1)
  
  for(it_random in 1:10){
    yearlySummary$period_rand<-sample(yearlySummary$period_P,
                                      size = nrow(yearlySummary),
                                      replace = FALSE)
    dat_x<-left_join(dat,yearlySummary%>%select(wateryear,period_rand))
    
    
    for(it_T in 1:2){
      
      dat_test<-dat_x%>%filter(period_rand  == it_T)
      dat_train = dat_x%>%filter(period_rand != it_T)
      
      dat_sum<-dat_train%>%
        group_by(yday)%>%
        summarize(q = mean(q,na.rm = TRUE))
      
      
      
      dat_test<-left_join(dat_test,dat_sum,by = "yday",suffix = c(".obs",".avg"))
      
      NSE_random <- rbind(NSE_random,
                          data.frame(
                            it_split = it_random,
                            NSE = hydroGOF::NSE(dat_test$q.avg,dat_test$q.obs,na.rm = TRUE)
                          )
      )    
      KGE_random <- rbind(KGE_random,
                          data.frame(
                            it_split = it_random,
                            KGE = hydroGOF::KGE(dat_test$q.avg,dat_test$q.obs,na.rm = TRUE)
                          )
      )
      
      
    }
    
    
    
  }
  
  stns$NSE_random  [it]<- NSE_random%>%group_by(it_split)%>%summarise(NSE = mean(NSE))%>%pull(NSE)%>%median()
  
  stns$KGE_random  [it]<- KGE_random%>%group_by(it_split)%>%summarise(KGE = mean(KGE))%>%pull(KGE)%>%median()
  
  toc()
  
}

st_write(stns,sprintf("2.data/GOF/hysets_dss_%s.gpkg",it_split),append = FALSE)