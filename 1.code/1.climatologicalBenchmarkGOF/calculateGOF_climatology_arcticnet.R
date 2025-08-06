# climatological benchmarks and variance components
# Author: Sacha Ruzzante
# sachawruzzante@gmail.com
# Last Update: 2025-08-06

# Arcticnet data:
#Lammers, R. B., & Shiklomanov, A. I. (2000). R-ArcticNet, A Regional Hydrographic Data Network for the Pan-Arctic Region. [Dataset]. https://www.r-arcticnet.sr.unh.edu/v4.0/AllData/index.html


library(ncdf4)
library(hydroGOF)
library(dplyr)
library(ggplot2)
library(hydroGOF)
library(sf)
library(tmap)
library(tictoc)
library(lubridate)
library(tidyr)
setwd("/home/ruzzante/projects/def-tgleeson/ruzzante/climatological_benchmarks/")
source("1.code/5.utils/utils.R")


## arcticnet ########

stns<-read.delim("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/arcticnet/Daily_SiteAttributes.txt",
                 sep = '\t')%>%
  filter(!is.na(Lat))%>%
  st_as_sf(coords = c("Long","Lat"),crs = "EPSG:4326")


# initialize dataframe for GOF statistics
stns$NSE<-NA
stns$KGE<-NA
stns$KGE_r<-NA
stns$KGE_a<-NA
stns$KGE_b<-NA

stns$N<-NA

# read in data

dat_discharge<-read.delim("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/arcticnet/discharge_m3_s_UNH-UCLA.txt")%>%
  pivot_longer(cols = Jan:Dec,names_to = "month",values_to = "q")%>%
  mutate(dt = ymd(paste(Year, month, Day)))%>%
  arrange(Code,dt)%>%
  filter(!is.na(dt))%>%
  mutate( yday =pmin(yday(dt),365))


# initialize dataframe for variance components
stn_var<-data.frame(Code = stns$Code,
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
  
  
  dat<-filter(dat_discharge,Code %in% stns$Code[it])
  
  
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
  
  hydroGOF::KGE(dat_test2$q.avg,dat_test2$q.obs,na.rm = TRUE,out.type	= "full") 
  
  
  stns$KGE_r[it]<- hydroGOF::KGE(dat_test2$q.avg,dat_test2$q.obs,na.rm = TRUE, s  =c(1,0,0))
  stns$KGE_a[it]<- hydroGOF::KGE(dat_test2$q.avg,dat_test2$q.obs,na.rm = TRUE, s  =c(0,1,0))
  stns$KGE_b[it]<- hydroGOF::KGE(dat_test2$q.avg,dat_test2$q.obs,na.rm = TRUE, s  =c(0,0,1))
  stns$N[it]<-length(yrs)
  
  dat$year<-dat$Year
  
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
write.csv(stn_var,"2.data/varComponents/arcticnet.csv",row.names = FALSE)


ggplot(stns,aes(x= NSE))+stat_ecdf()+
  scale_x_continuous(limits = c(-1,1))


ggplot(stns,aes(x= KGE))+stat_ecdf()+
  scale_x_continuous(limits = c(-1,1))
median(stns$NSE,na.rm = TRUE)


#visualize data
tmap_mode("view")
tm_shape(stns%>%filter(!is.na(NSE)))+tm_dots(col = "NSE",breaks = c(-Inf,seq(-1,1,0.2)),
                                             palette = "RdBu")


st_write(stns,"2.data/GOF/arcticnet.gpkg",append = FALSE)

# remove duplicates from GRDC
stns<-st_read("2.data/GOF/arcticnet.gpkg")

stns_grdc<-st_read("2.data/GOF/GRDC_asia.gpkg")

tmap_mode("view")


# manually look for stations in GRDC that are duplicates
stn_buffer<-st_buffer(stns,20000)%>%
  st_union()%>%
  st_make_valid()
stns_grdc_dup<-st_filter(stns_grdc,stn_buffer)

tm_shape(stns%>%filter(!is.na(NSE)))+tm_symbols(col = "NSE",breaks = c(-Inf,seq(-1,1,0.2)),
                                                palette = "RdBu", shape = 21,size = 0.4)+
  tm_shape(stns_grdc_dup)+tm_symbols(col = "NSE",breaks = c(-Inf,seq(-1,1,0.2)),
                                     palette = "RdBu",shape = 22,size = 0.4)


stns_grdc<-stns_grdc%>%
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
                         2903500)) # remove duplicates
st_write(stns_grdc,"2.data/GOF/GRDC_asia.gpkg",append = FALSE)
