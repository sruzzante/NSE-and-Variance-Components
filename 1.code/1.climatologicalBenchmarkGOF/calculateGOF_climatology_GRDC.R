# climatological benchmarks and variance components
# Author: Sacha Ruzzante
# sachawruzzante@gmail.com
# Last Update: 2025-06-18


# Global Runoff Data Centre, The (2025) 56068 Koblenz, Germany. https://grdc.bafg.de/

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
library(readxl)
setwd("/home/ruzzante/projects/def-tgleeson/ruzzante/climatological_benchmarks/")

source("1.code/utils.R")


# GRDC

fls<-list.files("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/asia/",
                pattern  ="*_Q_Day.Cmd.txt")

stnIDs<-str_remove(fls,"_Q_Day.Cmd.txt")

stns<-read_xlsx("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/GRDC_Stations.xlsx",sheet = 1)

stns<-filter(stns,
             grdc_no %in% stnIDs&
               d_yrs>=10)

# initialize dataframe for GOF statistics
stns$NSE<-NA
stns$KGE<-NA
stns$KGE_r<-NA
stns$KGE_a<-NA
stns$KGE_b<-NA
stns$QCI <-NA
stns$N<-NA


# initialize dataframe for variance components
stn_var<-data.frame(grdc_no = stns$grdc_no,
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
  
  dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/asia/",
                         stns$grdc_no,
                         "_Q_Day.Cmd.txt")[it],
                  sep = ";",skip = 36,header =  TRUE)%>%
    mutate(dt = ymd(YYYY.MM.DD ),
           year = year(dt),
           month = month(dt),
           yday = pmin(yday(dt),365),
           q = Value)
  
  
  # initial check for time series length >=10 years
  dat$q[dat$q==-999]<-NA
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
  
  dat_monthly = dat%>%
    group_by(year,month)%>%
    summarise(q = mean(q))%>%
    group_by(month)%>%
    summarise(q= mean(q,na.rm = TRUE))
  
  stns$QCI[it]<-sum(dat_monthly$q^2)/(sum(dat_monthly$q))^2*100
  
  
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



stn_var<-stn_var%>%
  filter(!grdc_no %in% c(2909100,
                         2911502,
                         2910450,
                         2907200,
                         2903200,
                         2997200,
                         2909250,
                         2910500,
                         2909700,
                         2903100,
                         2903960,
                         2998800,
                         2998702,
                         2903920,
                         2903400,
                         2903500)) # remove duplicates from arcticnet
write.csv(stn_var,"2.data/varComponents/GRDC_asia.csv",row.names = FALSE)

ggplot(stns,aes(x = QCI,y = NSE))+geom_point()+
  geom_smooth(method = "lm",se = FALSE)+
  scale_x_log10()

ggplot(stns,aes(x= NSE))+stat_ecdf()+
  scale_x_continuous(limits = c(-1,1))


ggplot(stns,aes(x= KGE))+stat_ecdf()+
  scale_x_continuous(limits = c(-1,1))

stns_sf<-st_as_sf(stns,coords = c("long","lat"),crs = "EPSG:4326")
tmap_mode("view")
tm_shape(stns_sf%>%filter(!is.na(NSE)))+tm_dots(col = "NSE",breaks = c(-Inf,seq(-1,1,0.25)),
                                                palette = "RdBu")

stns_sf<-stns_sf%>%
  filter(!grdc_no %in% c(2909100,
                         2911502,
                         2910450,
                         2907200,
                         2903200,
                         2997200,
                         2909250,
                         2910500,
                         2909700,
                         2903100,
                         2903960,
                         2998800,
                         2998702,
                         2903920,
                         2903400,
                         2903500)) # remove duplicates from arcticnet
st_write(stns_sf,"2.data/GOF/GRDC_asia.gpkg",append = FALSE)

median(stns$NSE,na.rm = TRUE)

## Africa


fls<-list.files("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/africa/",
                pattern  ="*_Q_Day.Cmd.txt")

stnIDs<-str_remove(fls,"_Q_Day.Cmd.txt")

stns<-read_xlsx("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/GRDC_Stations.xlsx",sheet = 1)

stns<-filter(stns,
             grdc_no %in% stnIDs&
               d_yrs>=10)
# initialize dataframe for GOF statistics
stns$NSE<-NA
stns$KGE<-NA
stns$KGE_r<-NA
stns$KGE_a<-NA
stns$KGE_b<-NA
stns$QCI <-NA
stns$N<-NA

# initialize dataframe for variance components
stn_var<-data.frame(grdc_no = stns$grdc_no,
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
  
  dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/africa/",
                         stns$grdc_no,
                         "_Q_Day.Cmd.txt")[it],
                  sep = ";",skip = 36,header =  TRUE)%>%
    mutate(dt = ymd(YYYY.MM.DD ),
           year = year(dt),
           month = month(dt),
           yday = pmin(yday(dt),365),
           q = Value)
  
  dat$q[dat$q==-999]<-NA
  if(nrow(dat)<(365*10)){next}
  if(sum(!is.na(dat$q))<=(365*10)){next}
  
  numYearsFull<-
    dat%>%
    group_by(yday)%>%
    summarise(notNA = sum(!is.na(q)))%>%
    ungroup()%>%
    summarize(allnotNA = min(notNA))
  if(numYearsFull$allnotNA<10){next}
  
  dat_monthly = dat%>%
    group_by(year,month)%>%
    summarise(q = mean(q))%>%
    group_by(month)%>%
    summarise(q= mean(q,na.rm = TRUE))
  
  stns$QCI[it]<-sum(dat_monthly$q^2)/(sum(dat_monthly$q))^2*100
  
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

write.csv(stn_var,"2.data/varComponents/GRDC_africa.csv",row.names = FALSE)

ggplot(stns,aes(x = QCI,y = NSE))+geom_point()+
  geom_smooth(method = "lm",se = FALSE)+
  scale_x_log10()

ggplot(stns,aes(x= NSE))+stat_ecdf()+
  scale_x_continuous(limits = c(-1,1))


ggplot(stns,aes(x= KGE))+stat_ecdf()+
  scale_x_continuous(limits = c(-1,1))

stns_sf<-st_as_sf(stns,coords = c("long","lat"),crs = "EPSG:4326")
tmap_mode("view")
tm_shape(stns_sf%>%filter(!is.na(NSE)))+tm_dots(col = "NSE",breaks = c(-Inf,seq(-1,1,0.25)),
                                                palette = "RdBu")
st_write(stns_sf,"2.data/GOF/GRDC_africa.gpkg",append = FALSE)

median(stns$NSE,na.rm = TRUE)







## south_west_pacific


fls<-list.files("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/south_west_pacific/",
                pattern  ="*_Q_Day.Cmd.txt")

stnIDs<-str_remove(fls,"_Q_Day.Cmd.txt")

stns<-read_xlsx("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/GRDC_Stations.xlsx",sheet = 1)

stns<-filter(stns,
             grdc_no %in% stnIDs&
               d_yrs>=10)
unique(stns$country)
stns<-filter(stns,!country %in% c("AU"))
# initialize dataframe for GOF statistics
stns$NSE<-NA
stns$KGE<-NA
stns$KGE_r<-NA
stns$KGE_a<-NA
stns$KGE_b<-NA
stns$QCI <-NA
stns$N<-NA

# initialize dataframe for variance components
stn_var<-data.frame(grdc_no = stns$grdc_no,
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
  
  dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/south_west_pacific/",
                         stns$grdc_no,
                         "_Q_Day.Cmd.txt")[it],
                  sep = ";",skip = 36,header =  TRUE)%>%
    mutate(dt = ymd(YYYY.MM.DD ),
           year = year(dt),
           month = month(dt),
           yday = pmin(yday(dt),365),
           q = Value)
  
  dat$q[dat$q==-999]<-NA
  if(nrow(dat)<(365*10)){next}
  if(sum(!is.na(dat$q))<=(365*10)){next}
  
  numYearsFull<-
    dat%>%
    group_by(yday)%>%
    summarise(notNA = sum(!is.na(q)))%>%
    ungroup()%>%
    summarize(allnotNA = min(notNA))
  if(numYearsFull$allnotNA<10){next}
  
  dat_monthly = dat%>%
    group_by(year,month)%>%
    summarise(q = mean(q))%>%
    group_by(month)%>%
    summarise(q= mean(q,na.rm = TRUE))
  
  stns$QCI[it]<-sum(dat_monthly$q^2)/(sum(dat_monthly$q))^2*100
  
  
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

write.csv(stn_var,"2.data/varComponents/GRDC_south_west_pacific.csv",row.names = FALSE)


ggplot(stns,aes(x = QCI,y = NSE))+geom_point()+
  geom_smooth(method = "lm",se = FALSE)+
  scale_x_log10()

ggplot(stns,aes(x= NSE))+stat_ecdf()+
  scale_x_continuous(limits = c(-1,1))


ggplot(stns,aes(x= KGE))+stat_ecdf()+
  scale_x_continuous(limits = c(-1,1))

stns_sf<-st_as_sf(stns,coords = c("long","lat"),crs = "EPSG:4326")
tmap_mode("view")
tm_shape(stns_sf%>%filter(!is.na(NSE)))+tm_dots(col = "NSE",breaks = c(-Inf,seq(-1,1,0.25)),
                                                palette = "RdBu")
st_write(stns_sf,"2.data/GOF/GRDC_south_west_pacific.gpkg",append = FALSE)

median(stns$NSE,na.rm = TRUE)





## south_america


fls<-list.files("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/south_america/",
                pattern  ="*_Q_Day.Cmd.txt")

stnIDs<-str_remove(fls,"_Q_Day.Cmd.txt")

stns<-read_xlsx("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/GRDC_Stations.xlsx",sheet = 1)

stns<-filter(stns,
             grdc_no %in% stnIDs&
               d_yrs>=10)
unique(stns$country)
stns<-filter(stns,!country %in% c("CL","BR"))
# initialize dataframe for GOF statistics
stns$NSE<-NA
stns$KGE<-NA
stns$KGE_r<-NA
stns$KGE_a<-NA
stns$KGE_b<-NA
stns$QCI <-NA
stns$N<-NA

# initialize dataframe for variance components
stn_var<-data.frame(grdc_no = stns$grdc_no,
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
  
  dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/south_america/",
                         stns$grdc_no,
                         "_Q_Day.Cmd.txt")[it],
                  sep = ";",skip = 36,header =  TRUE)%>%
    mutate(dt = ymd(YYYY.MM.DD ),
           year = year(dt),
           month = month(dt),
           yday = pmin(yday(dt),365),
           q = Value)
  
  dat$q[dat$q==-999]<-NA
  if(nrow(dat)<(365*10)){next}
  if(sum(!is.na(dat$q))<=(365*10)){next}
  
  numYearsFull<-
    dat%>%
    group_by(yday)%>%
    summarise(notNA = sum(!is.na(q)))%>%
    ungroup()%>%
    summarize(allnotNA = min(notNA))
  if(numYearsFull$allnotNA<10){next}
  
  ggplot(dat,aes(yday,q,group = factor(year)))+geom_line()
  
  dat_monthly = dat%>%
    group_by(year,month)%>%
    summarise(q = mean(q))%>%
    group_by(month)%>%
    summarise(q= mean(q,na.rm = TRUE))
  
  stns$QCI[it]<-sum(dat_monthly$q^2)/(sum(dat_monthly$q))^2*100
  
  
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

write.csv(stn_var,"2.data/varComponents/GRDC_south_america.csv",row.names = FALSE)

ggplot(stns,aes(x = QCI,y = NSE))+geom_point()+
  geom_smooth(method = "lm",se = FALSE)+
  scale_x_log10()

ggplot(stns,aes(x= NSE))+stat_ecdf()+
  scale_x_continuous(limits = c(-1,1))


ggplot(stns,aes(x= KGE))+stat_ecdf()+
  scale_x_continuous(limits = c(-1,1))

stns_sf<-st_as_sf(stns,coords = c("long","lat"),crs = "EPSG:4326")
tmap_mode("view")
tm_shape(stns_sf%>%filter(!is.na(NSE)))+tm_dots(col = "NSE",breaks = c(-Inf,seq(-1,1,0.25)),
                                                palette = "RdBu")
st_write(stns_sf,"2.data/GOF/GRDC_south_america.gpkg",append = FALSE)

median(stns$NSE,na.rm = TRUE)


## north_america


fls<-list.files("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/north_america/",
                pattern  ="*_Q_Day.Cmd.txt")

stnIDs<-str_remove(fls,"_Q_Day.Cmd.txt")

stns<-read_xlsx("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/GRDC_Stations.xlsx",sheet = 1)

stns<-filter(stns,
             grdc_no %in% stnIDs&
               d_yrs>=10)
unique(stns$country)
stns<-filter(stns,!country %in% c("US","CA","MX"))
# initialize dataframe for GOF statistics
stns$NSE<-NA
stns$KGE<-NA
stns$KGE_r<-NA
stns$KGE_a<-NA
stns$KGE_b<-NA
stns$N<-NA

# initialize dataframe for variance components
stn_var<-data.frame(grdc_no = stns$grdc_no,
                    varSeas_stl = NA,
                    varInterannual_stl = NA,
                    varRem_stl = NA,
                    varSeas_clas = NA,
                    varInterannual_clas = NA,
                    varRem_clas = NA,
                    varSeas_fourier = NA,
                    varInterannual_fourier = NA,
                    varRem_fourier = NA)

stns$QCI <-NA

#Loop through stations and calculate statistics
for(it in 1:nrow(stns)){
  tic(sprintf("station %d",it))
  
  dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/north_america/",
                         stns$grdc_no,
                         "_Q_Day.Cmd.txt")[it],
                  sep = ";",skip = 36,header =  TRUE)%>%
    mutate(dt = ymd(YYYY.MM.DD ),
           year = year(dt),
           month = month(dt),
           yday = pmin(yday(dt),365),
           q = Value)
  
  dat$q[dat$q==-999]<-NA
  if(nrow(dat)<(365*10)){next}
  if(sum(!is.na(dat$q))<=(365*10)){next}
  
  numYearsFull<-
    dat%>%
    group_by(yday)%>%
    summarise(notNA = sum(!is.na(q)))%>%
    ungroup()%>%
    summarize(allnotNA = min(notNA))
  if(numYearsFull$allnotNA<10){next}
  

  dat_monthly = dat%>%
    group_by(year,month)%>%
    summarise(q = mean(q))%>%
    group_by(month)%>%
    summarise(q= mean(q,na.rm = TRUE))
  
  stns$QCI[it]<-sum(dat_monthly$q^2)/(sum(dat_monthly$q))^2*100
  
  
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

write.csv(stn_var,"2.data/varComponents/GRDC_north_america.csv",row.names = FALSE)


ggplot(stns,aes(x = QCI,y = NSE))+geom_point()+
  geom_smooth(method = "lm",se = FALSE)+
  scale_x_log10()

ggplot(stns,aes(x= NSE))+stat_ecdf()+
  scale_x_continuous(limits = c(-1,1))


ggplot(stns,aes(x= KGE))+stat_ecdf()+
  scale_x_continuous(limits = c(-1,1))

stns_sf<-st_as_sf(stns,coords = c("long","lat"),crs = "EPSG:4326")
tmap_mode("view")
tm_shape(stns_sf%>%filter(!is.na(NSE)))+tm_dots(col = "NSE",breaks = c(-Inf,seq(-1,1,0.25)),
                                                palette = "RdBu")
st_write(stns_sf,"2.data/GOF/GRDC_north_america.gpkg",append = FALSE)

median(stns$NSE,na.rm = TRUE)




## europe_aq


fls<-list.files("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/europe_a_q/",
                pattern  ="*_Q_Day.Cmd.txt")

stnIDs<-str_remove(fls,"_Q_Day.Cmd.txt")

stns<-read_xlsx("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/GRDC_Stations.xlsx",sheet = 1)

stns<-filter(stns,
             grdc_no %in% stnIDs&
               d_yrs>=10)
unique(stns$country)%>%sort()

stns<-filter(stns,!country %in% c("DE","DK","MX","GB","LU","FR","AT","ES","IL","CH","IS")|
               grdc_no %in% c("6603330","6603510","6603500","6603320","6603350","6603310","6603100",
                              "6603110","6603120","6603300","6603340")) # 11 GRDC stations in Northern Ireland are incorrectly labelled GB 


stns<-filter(stns,!grdc_no %in% (c("6142170","6142250",
                                   "6142110","6142150","6142120","6142300",
                                   "6142260","6142160","6142100",
                                   "6142190","6142101","6142270"))) # remove stations that are in Lamah-ce (Czechia)
stns_sf<-st_as_sf(stns,coords = c("long","lat"),crs = "EPSG:4326")
tmap_mode("view")
tm_shape(stns_sf)+tm_dots(col = "country")


# initialize dataframe for GOF statistics
stns$NSE<-NA
stns$KGE<-NA
stns$KGE_r<-NA
stns$KGE_a<-NA
stns$KGE_b<-NA
stns$QCI <-NA
stns$N<-NA

# initialize dataframe for variance components
stn_var<-data.frame(grdc_no = stns$grdc_no,
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
  
  dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/europe_a_q/",
                         stns$grdc_no,
                         "_Q_Day.Cmd.txt")[it],
                  sep = ";",skip = 36,header =  TRUE)%>%
    mutate(dt = ymd(YYYY.MM.DD ),
           year = year(dt),
           month = month(dt),
           yday = pmin(yday(dt),365),
           q = Value)
  
  dat$q[dat$q==-999]<-NA
  if(nrow(dat)<(365*10)){next}
  if(sum(!is.na(dat$q))<=(365*10)){next}
  
  numYearsFull<-
    dat%>%
    group_by(yday)%>%
    summarise(notNA = sum(!is.na(q)))%>%
    ungroup()%>%
    summarize(allnotNA = min(notNA))
  if(numYearsFull$allnotNA<10){next}
  

  dat_monthly = dat%>%
    group_by(year,month)%>%
    summarise(q = mean(q))%>%
    group_by(month)%>%
    summarise(q= mean(q,na.rm = TRUE))
  
  stns$QCI[it]<-sum(dat_monthly$q^2)/(sum(dat_monthly$q))^2*100
  
  
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

write.csv(stn_var,"2.data/varComponents/GRDC_europe_a_q.csv",row.names = FALSE)


ggplot(stns,aes(x = QCI,y = NSE))+geom_point()+
  geom_smooth(method = "lm",se = FALSE)+
  scale_x_log10()

ggplot(stns,aes(x= NSE))+stat_ecdf()+
  scale_x_continuous(limits = c(-1,1))


ggplot(stns,aes(x= KGE))+stat_ecdf()+
  scale_x_continuous(limits = c(-1,1))

stns_sf<-st_as_sf(stns,coords = c("long","lat"),crs = "EPSG:4326")
tmap_mode("view")
tm_shape(stns_sf%>%filter(!is.na(NSE)))+tm_dots(col = "NSE",breaks = c(-Inf,seq(-1,1,0.25)),
                                                palette = "RdBu")
st_write(stns_sf,"2.data/GOF/GRDC_europe_a_q.gpkg",append = FALSE)

median(stns$NSE,na.rm = TRUE)




#europe_r_z


fls<-list.files("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/europe_r_z/",
                pattern  ="*_Q_Day.Cmd.txt")

stnIDs<-str_remove(fls,"_Q_Day.Cmd.txt")

stns<-read_xlsx("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/GRDC_Stations.xlsx",sheet = 1)

stns<-filter(stns,
             grdc_no %in% stnIDs&
               d_yrs>=10)
unique(stns$country)%>%sort()
stns<-filter(stns,!country %in% c("DE","DK","MX","GB","LU","FR","AT","ES","IL","CH","IS"))
stns_sf<-st_as_sf(stns,coords = c("long","lat"),crs = "EPSG:4326")
tmap_mode("view")
tm_shape(stns_sf)+tm_dots(col = "country")




# initialize dataframe for GOF statistics
stns$NSE<-NA
stns$KGE<-NA
stns$KGE_r<-NA
stns$KGE_a<-NA
stns$KGE_b<-NA
stns$QCI <-NA
stns$N<-NA

# initialize dataframe for variance components
stn_var<-data.frame(grdc_no = stns$grdc_no,
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
  
  dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/europe_r_z/",
                         stns$grdc_no,
                         "_Q_Day.Cmd.txt")[it],
                  sep = ";",skip = 36,header =  TRUE)%>%
    mutate(dt = ymd(YYYY.MM.DD ),
           year = year(dt),
           month = month(dt),
           yday = pmin(yday(dt),365),
           q = Value)
  
  dat$q[dat$q==-999]<-NA
  if(nrow(dat)<(365*10)){next}
  if(sum(!is.na(dat$q))<=(365*10)){next}
  
  numYearsFull<-
    dat%>%
    group_by(yday)%>%
    summarise(notNA = sum(!is.na(q)))%>%
    ungroup()%>%
    summarize(allnotNA = min(notNA))
  if(numYearsFull$allnotNA<10){next}
  

  dat_monthly = dat%>%
    group_by(year,month)%>%
    summarise(q = mean(q))%>%
    group_by(month)%>%
    summarise(q= mean(q,na.rm = TRUE))
  
  stns$QCI[it]<-sum(dat_monthly$q^2)/(sum(dat_monthly$q))^2*100
  
  
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

write.csv(stn_var,"2.data/varComponents/GRDC_europe_r_z.csv",row.names = FALSE)


ggplot(stns,aes(x = QCI,y = NSE))+geom_point()+
  geom_smooth(method = "lm",se = FALSE)+
  scale_x_log10()

ggplot(stns,aes(x= NSE))+stat_ecdf()+
  scale_x_continuous(limits = c(-1,1))


ggplot(stns,aes(x= KGE))+stat_ecdf()+
  scale_x_continuous(limits = c(-1,1))

stns_sf<-st_as_sf(stns,coords = c("long","lat"),crs = "EPSG:4326")
tmap_mode("view")
tm_shape(stns_sf%>%filter(!is.na(NSE)))+tm_dots(col = "NSE",breaks = c(-Inf,seq(-1,1,0.25)),
                                                palette = "RdBu")
st_write(stns_sf,"2.data/GOF/GRDC_europe_r_z.gpkg",append = FALSE)

median(stns$NSE,na.rm = TRUE)






quit()





