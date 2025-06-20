
library(dplyr)
library(stringr)
library(lubridate)
library(tictoc)
setwd("/home/ruzzante/projects/def-tgleeson/ruzzante/climatological_benchmarks/")
fls<-list.files("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-br/preprocessed/",full.names = TRUE)



for(it_fl in 1:length(fls)){
  tic(sprintf("catchment %d",it_fl))
  dat<-read.csv(fls[it_fl])
  
  tmin_cpc_fill<-c(NA, dat$tmin_cpc[1:(nrow(dat)-1)])+
    c( dat$tmin_cpc[2:(nrow(dat))],NA)
  dat$tmin_cpc[is.na(dat$tmin_cpc)]<-tmin_cpc_fill[is.na(dat$tmin_cpc)]
  
  tmax_cpc_fill<-c(NA, dat$tmax_cpc[1:(nrow(dat)-1)])+
    c( dat$tmax_cpc[2:(nrow(dat))],NA)
  dat$tmax_cpc[is.na(dat$tmax_cpc)]<-tmax_cpc_fill[is.na(dat$tmax_cpc)]
  
  write.csv(dat,fls[it_fl],
            row.names = FALSE,
            na = "")
  toc()
}


for(it_fl in 386:length(fls)){
  stnID<-str_split_fixed(fls[it_fl],"/|\\.",13)[,12]
  # tic(sprintf("catchment %d",it_fl))
  
  dat<-read.csv(fls[it_fl])%>%
    mutate(yd = yday(ymd(date)))
  
  dat_med<-dat%>%
    group_by(yd)%>%
    summarize(
      p_cpc = median(p_cpc,na.rm = TRUE),
      tmax_cpc = median(tmax_cpc,na.rm = TRUE),
      tmin_cpc = median(tmin_cpc,na.rm = TRUE))
  
  dat<-left_join(dat,dat_med,by = "yd",suffix = c("",".med"))
  
  mask_na<-is.na(dat$p_cpc & year(ymd(dat$date))<2024)
  dat$p_cpc[mask_na]<-dat$p_cpc.med[mask_na]
  num_p<-sum(mask_na)
  
  mask_na<-is.na(dat$tmin_cpc & year(ymd(dat$date))<2024)
  dat$tmin_cpc[mask_na]<-dat$tmin_cpc.med[mask_na]
  num_tmin<-sum(mask_na)
  
  mask_na<-is.na(dat$tmax_cpc & year(ymd(dat$date))<2024)
  dat$tmax_cpc[mask_na]<-dat$tmax_cpc.med[mask_na]
  num_tmax<-sum(mask_na)
  
  write.csv(dat,fls[it_fl],
            row.names = FALSE,
            na = "")
  
  print(sprintf("Overwrote %d p, %d tmin, and %d tmax values for station %s",num_p,num_tmin,num_tmax,stnID))
  # toc()
}
