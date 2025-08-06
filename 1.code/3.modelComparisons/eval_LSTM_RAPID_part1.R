# climatological benchmarks and variance components
# Author: Sacha Ruzzante
# sachawruzzante@gmail.com
# Last Update: 2025-08-06

# Evaluate GRID-LSTM-RAPID for dataset (part 1)
# Yang, Y., Feng, D., Beck, H. E., Hu, W., Abbas, A., Sengupta, A., Delle Monache, L., Hartman, R., Lin, P., Shen, C., & Pan, M. (2025). Global Daily Discharge Estimation Based on Grid Long Short-Term Memory (LSTM) Model and River Routing. Water Resources Research, 61(6), e2024WR039764. https://doi.org/10.1029/2024WR039764

# This script collects and resaves observed streamflow data from various sources for catchments that are modelled by Yang et al (2025)


library(ncdf4)
library(hydroGOF)
library(dplyr)
library(ggplot2)
library(hydroGOF)
library(sf)
# library(tmap)
library(tictoc)
library(lubridate)
library(stringr)



setwd("/home/ruzzante/projects/def-tgleeson/ruzzante/climatological_benchmarks/")
source("1.code/5.utils/utils.R")

# stations<-read.csv("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRID-LSTM-RAPID/input/training_gauge_attributes.csv")
stations<-read.csv("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRID-LSTM-RAPID/input/training_gauge_information.csv")
src_id<-str_split_fixed(stations$name,"_",n =3)

stations$src<-paste0(src_id[,1],src_id[,2])
stations$gauge_id<-str_remove(src_id[,3],"_")

summary(factor(stations$src))


summary(factor(stations$src))

streamflow_dat_aus<-read.csv("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-aus-v2/streamflow_MLd.csv")

fls_done<-list.files("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRID-LSTM-RAPID/observations/")%>%
  str_remove(".rds")

stations<-stations%>%
  filter(!name %in% fls_done)

summary(factor(stations$src))

#preload polish data
fls<-list.files("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/poland_IMGW/data/",
                pattern = ".csv",
                full.names = TRUE)
fls<-fls[!fls == "../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/poland_IMGW/data//codz_2023.csv" ]

dat_poland<-lapply(fls,function(fl){
  dat<-read.csv(fl,header = FALSE)
  names(dat)<-c("gauge_id","station_nm","river.lake_nm",'hydro_year',"hydro_month","day","water_level","discharge","water_temp","month")
  dat$year<-dat$hydro_year
  dat$year[dat$month %in% c(11,12)]<-dat$hydro_year[dat$month %in% c(11,12)]-1
  dat$dt<-ymd(paste(dat$year,dat$month,dat$day))
  dat$discharge[dat$discharge>=99999. &dat$discharge<1e5]<-NA
  dat%>%select(gauge_id,dt,discharge)
})%>%bind_rows()
# add 2023
dat<-read.delim("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/poland_IMGW/data//codz_2023.csv" ,
                    sep = ";",header = FALSE)
names(dat)<-c("gauge_id","station_nm","river.lake_nm",'hydro_year',"hydro_month","day","water_level","discharge","water_temp","month")
dat$year<-dat$hydro_year
dat$year[dat$month %in% c(11,12)]<-dat$hydro_year[dat$month %in% c(11,12)]-1
dat$dt<-ymd(paste(dat$year,dat$month,dat$day))
dat$discharge[dat$discharge>=99999. &dat$discharge<1e5]<-NA
dat<-dat%>%select(gauge_id,dt,discharge)

dat_poland<-rbind(dat_poland,dat)
dat_poland$gauge_id<-as.numeric(dat_poland$gauge_id)%>%as.character()

readOBS<-function(stations, ID, src){
  
  dat<-data.frame()
  
  
  if(src %in%  c("HYSETS","USGS")){
    
    
    nc_dat<-ncdf4::nc_open("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/hysets/HYSETS_2023_update_QC_stations.nc")
    
    time_length <- nc_dat$dim$time$len
    dt<-as.Date("1950-01-01")+
      days(ncvar_get(nc = nc_dat,
                     varid = "time"))
    
    yr = year(dt)
    mnth = month(dt)
    dy = day(dt)
    ydy = pmin(yday(dt),365)
    
    
    stns<-read.csv("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/hysets/HYSETS_watershed_properties.txt")%>%
      select(Watershed_ID,Source,Name,Official_ID,Hydrometric_station_latitude,Hydrometric_station_longitude)
    
    it_id<-which(stns$Official_ID == ID)
    if(length(it_id)==1){
      dat<-data.frame(q= ncvar_get(nc = nc_dat,
                                   varid = "discharge",
                                   start = c(1, stns$Watershed_ID[it_id]), # Start at the first time step and the desired catchment
                                   count = c(time_length, 1)),     # Read all time steps for the catchment)
                      dt = dt
      )
    }
    
    
    ncdf4::nc_close(nc_dat)
    
  }
  
  if(src == "BOMAustralia"){
    streamflow_dat<-streamflow_dat_aus
    stns<-st_read("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-aus-v2/02_location_boundary_area/shp/CAMELS_AUS_v2_BasinOutlets_adopted.shp")
    stns$station_id[stns$station_id=="227225A"] = "227225"
    
    if((substr(ID,1,1))%in% as.character(0:9)){ID = paste0("X",ID)}
    if(ID %in% names(streamflow_dat)){
      
      dat<-streamflow_dat%>%
        select(year,month,day,all_of(ID))%>%
        mutate(dt = ymd(paste(year,month,day)),
               yday = pmin(yday(dt),365))
      
      names(dat)[4]<-"q"
      
      # dat<-filter(dat,!is.na(dat$q))
      
      
      dat$q[dat$q== -99.99]<-NA
      dat$q<-dat$q/86.4 # MLd to m3s
      dat<-select(dat,q,dt)
    }
    
    
  }
  
  if(src == "CAMELSBR"){
    
    dat<-read.delim(sprintf("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-br/02_CAMELS_BR_streamflow_all_catchments/%s_streamflow.txt",ID),
                    sep = " ")%>%
      filter(qual_flag == 1)%>%
      mutate(dt = ymd(paste(year,month,day)),
             q = streamflow_m3s)%>%
      select(q,dt)
    
  }
  
  if(src =="CAMELSCL"){
    streamflow_dat<-read.csv("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-cl/CAMELS_CL_v202201/q_m3s_day.csv")
    
    dat<-streamflow_dat%>%
      select(date,all_of(paste0("X",ID)))%>%
      mutate(dt = ymd(date))
    
    names(dat)[2]<-"q"
    
    dat<-dplyr::select(dat,c(q,dt))
    
    
    
  }
  
  if(src  == "CAMELSGB"){
    
    dat<-read.csv(paste0(
      "../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-gb/data/timeseries/CAMELS_GB_hydromet_timeseries_",
      ID,"_19701001-20150930.csv"
    ))%>%
      select(date,discharge_vol)%>%
      mutate(q = discharge_vol,
             dt = ymd(date),
      )%>%
      select(q,dt)
  }
  # if(src == "Denmark"){
  #   stns<-st_read("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-dk/Shapefile/CAMELS_DK_304_gauging_stations.shp")%>%
  #     st_transform("EPSG:4326")%>%
  #     select(!c(X,Y))
  #   stns<-cbind(stns,st_coordinates(stns))
  #   stns<-mutate(stns,
  #                lon= round(X,3),
  #                lat = round(Y,3))
  #   it_row<-which(stations$src==src&
  #                   stations$gauge_id == ID)
  #   
  #   it_row_stns<-which(stns$lon==stations$lon[it_row]&
  #                        stns$lat == stations$lat[it_row])
  #   
  #   dat<-read.csv(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-dk/Dynamics/Gauged_catchments/CAMELS_DK_obs_based_",
  #                        ID,".csv"))%>%
  #     select(time,Qobs)%>%
  #     mutate(dt = ymd(time),
  #            q = Qobs)%>%
  #     select(q,dt)
  #   
  # }
  
  if(src == "France"){
    if(file.exists(paste0(
      "../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-fr/CAMELS_FR_time_series/daily/CAMELS_FR_tsd_",
      ID,".csv"
    ))){
      dat<-read.delim(paste0(
        "../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-fr/CAMELS_FR_time_series/daily/CAMELS_FR_tsd_",
        ID,".csv"
      ),
      skip  =7,sep = ";")%>%
        select(tsd_date,tsd_q_l)%>%
        mutate(q = tsd_q_l/1000,
               dt = ymd(tsd_date))%>%
        select(q,dt)
      
    }
    
  }
  
  if(src == "Germany"){
    
    
    
    stns<-read.csv("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-de/CAMELS_DE_topographic_attributes.csv")
    
    it_row = which(stns$provider_id==ID)
    if(length(it_row)==1){
      
      dat<-read.csv(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-de/timeseries/CAMELS_DE_hydromet_timeseries_",
                           stns$gauge_id[it_row],
                           ".csv"))%>%
        select(date,discharge_vol_obs)%>%
        mutate(dt = ymd(date),
               q = discharge_vol_obs)%>%
        select(q,dt)
    }
  }
  
  if(src =="GRDC"){
    fl<-list.files("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/",
                   pattern = paste0(ID,"_Q_Day.Cmd.txt"),
                   recursive = TRUE,
                   full.names = TRUE)
    if(length(fl)==1){
      dat<-read.delim(fl,
                      sep = ";",skip = 36,header =  TRUE)%>%
        mutate(dt = ymd(YYYY.MM.DD ),
               
               q = Value)%>%
        select(q,dt)
      
    }
    
    
  }
  if(src =="Spain"){
    
    if(file.exists(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-es/timeseries/netcdf/camelses/camelses_",
                          ID,".nc"))){
      nc_dat<-nc_open(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-es/timeseries/netcdf/camelses/camelses_",
                             ID,".nc"))
      
      dat<-data.frame(
        q = ncvar_get(nc_dat,"streamflow"),
        dt = as.Date("1990-01-01") + days(ncvar_get(nc_dat,"date"))
      )
    }
    
    
  }
  if(src == "LamaHIce"){
    
    dat<-read.delim(paste0(
      "../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/lamah_ice/lamah_ice/D_gauges/2_timeseries/daily/ID_",
      ID,".csv"
    ),
    sep = ";")%>%
      
      mutate(q = qobs,
             dt = ymd(paste(YYYY,MM,DD)))%>%
      select(dt,q)
    
  }
  if(src == "Thailand"){
    fls<-list.files("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/thailand_RID/",
                    pattern = paste0("disc_d_\\d{4}_RID",substr(ID,1,1),"\\.",substr(ID,2,100),"_m3s-1.txt"),
                    recursive = TRUE,
                    full.names = TRUE)
    dat<-lapply(fls,function(fl){read.delim(fl,sep = "\t",header = FALSE)%>%
        dplyr::mutate(dt = ymd(V1),q = V2)%>%
        select(q,dt)
        })%>%
      bind_rows()%>%
      filter(!is.na(dt))
    dat$q[dat$q<0]<-NA
  }
  
  if(src == "Sweden"){
    fl<-list.files("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-se/2023-173-1/data/catchment time series/",
                   pattern = paste0("catchment_id_",ID,"_*"),
                   full.names = TRUE)
    if(length(fl)==1){
      dat<-read.csv(fl)%>%
        mutate(dt = ymd(paste(Year,Month,Day)),
               q = Qobs_m3s)%>%
        select(q,dt)
    }
    
  }
  if(src == "Poland"){
    dat<-dat_poland%>%
      filter(gauge_id == ID)%>%
      mutate(q=discharge)%>%
      select(q,dt)
    
  }
  
  return(dat)
  
}

for(it in 1:nrow(stations)){
  tic(sprintf("row %d, station %s, %s",it,stations$gauge_id[it],stations$src[it]))
  dat<-readOBS(stations,stations$gauge_id[it],src = stations$src[it])
  if(nrow(dat)>0){
    saveRDS(dat,
            paste0(
              "../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRID-LSTM-RAPID/observations/",
              stations$name[it],".rds"
            )
    )
  }
  
  toc()
}
