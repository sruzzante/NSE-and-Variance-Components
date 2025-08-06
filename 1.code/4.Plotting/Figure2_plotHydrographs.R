# climatological benchmarks for hydrology

library(ncdf4)
library(hydroGOF)
library(dplyr)
library(ggplot2)
library(hydroGOF)
library(sf)
library(tmap)
library(tictoc)
library(lubridate)
library(plotly)
setwd("/home/ruzzante/projects/def-tgleeson/ruzzante/climatological_benchmarks/")

#Brazil
stns<-st_read("2.data/GOF/camels_br.gpkg")
it = which.max(stns$NSE)
it = which(stns$gauge_id=="13880000")
dat<-read.delim(sprintf("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-br/02_CAMELS_BR_streamflow_all_catchments/%s_streamflow.txt",stns$gauge_id[it]),
                sep = "")%>%
  mutate(dt = ymd(paste(year,month,day)),
         yday = pmin(yday(dt),365))%>%
  filter(qual_flag == 1)

names(dat)[4]<-"q"
dat<-left_join(data.frame(dt = seq.Date(dat$dt[1],dat$dt[nrow(dat)],by = "1 day"))%>%
                 mutate(year = year(dt)),
               dat)
dat$fkDate<-ymd("2001-12-31")+days(dat$yday)

dat_yearly<-dat%>%
  group_by(year)%>%
  summarize(full = sum(is.na(q))==0)

dat<-dat%>%filter(year %in% dat_yearly$year[dat_yearly$full])

dat_sum<-dat%>%
  group_by(fkDate)%>%
  summarise(q = mean(q,na.rm = TRUE))

ggplot(dat,aes(x = dt,y = q))+
  geom_line(aes(group = year),col = "grey",alpha = 0.5)

# m_labels<-c("J","F","M","A","M","J","J","A","S","O","N","D")
m_labels<-c("Jan ","Apr ","Jul ","Oct ","Jan  ")

p<-ggplot(dat,aes(x = fkDate,y = q))+
  geom_line(aes(group = year),col = "grey",alpha = 0.5,linewidth = 0.25)+
  geom_line(data = dat_sum,col = "black")+
  
  scale_x_date(name = "Date",
               # date_minor_breaks = "10 year",
               limits = ymd(c(20020101,20030101)),
               expand = c(0,0),
               # date_breaks = "3 month",
               # breaks
               breaks  = seq(from=as.Date("2002-01-01"),to=as.Date("2003-01-01"),by="3 month"),
               labels = m_labels)+
  scale_y_continuous(name = expression(Q~(m^3*s^-1)),
                     limits = c(0,14000),
                     expand = c(0,0),
                     breaks = c(0,10000),
                     # labels =c("0","10^5"),
                     labels = parse(text = c("0","10^4"))
                     )+
  theme_bw(base_size = 6)+
  # ggtitle("Rio Purus at Canutama (13880000)")+
  theme(panel.grid = element_blank())
p

ggsave(p,filename= "3.figures/Canutama.png",width = 2,height = 2,dpi = 600)
ggplotly(p)

x<-left_join(dat,dat_sum,by = "fkDate")


#Takhini River Near Whitehorse

stns<-st_read("2.data/GOF/hysets.gpkg")
it = which.max(stns$NSE)

it=which(stns$Official_ID == "09AB009")

nc_dat<-ncdf4::nc_open("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/hysets/HYSETS_2023_update_QC_stations.nc")
dt<-as.Date("1950-01-01")+
  days(ncvar_get(nc = nc_dat,
                 varid = "time"))
time_length <- nc_dat$dim$time$len
yr = year(dt)
mnth = month(dt)
dy = day(dt)
ydy = pmin(yday(dt),365)

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

dat<-left_join(data.frame(dt = seq.Date(dat$dt[1],dat$dt[nrow(dat)],by = "1 day"))%>%
                 mutate(year = year(dt)),
               dat)
dat$fkDate<-ymd("2001-12-31")+days(dat$yday)

dat_yearly<-dat%>%
  group_by(year)%>%
  summarize(full = sum(is.na(q))==0)

dat<-dat%>%filter(year %in% dat_yearly$year[dat_yearly$full])

dat_sum<-dat%>%
  group_by(fkDate)%>%
  summarise(q = mean(q,na.rm = TRUE))
m_labels<-c("Jan ","Apr ","Jul ","Oct ","Jan  ")

p<-ggplot(dat,aes(x = fkDate,y = q))+
  geom_line(aes(group = year),col = "grey",alpha = 0.5,linewidth = 0.25)+
  geom_line(data = dat_sum,col = "black")+
  
  scale_x_date(name = "Date",
               # date_minor_breaks = "10 year",
               limits = ymd(c(20020101,20030101)),
               expand = c(0,0),
               # date_breaks = "3 month",
               # breaks
               breaks  = seq(from=as.Date("2002-01-01"),to=as.Date("2003-01-01"),by="3 month"),
               labels = m_labels)+
  scale_y_continuous(name = expression(Q~(m^3*s^-1)),
                     # limits = c(0,350),
                     expand = c(0,0)
                     # breaks = c(0,100,200,300)
                     # labels =c("0","10^5"),
                     # labels = parse(text = c("0","10^4"))
  )+
  theme_bw(base_size = 6)+
  # ggtitle("Rio Purus at Canutama (13880000)")+
  theme(panel.grid = element_blank())
p

ggsave(p,filename= "3.figures/Takhini.png",width = 2,height = 2,dpi = 600)

stns<-stns%>%st_set_crs(st_crs("EPSG:4326"))
stns<-stns%>%filter(Hydrometric_station_latitude>56&Hydrometric_station_longitude< -130)

stns<-st_transform(stns,"+proj=wintri +lon_0=0 +lat_1=50.467 +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs +type=crs")
tmap_mode("plot")
cst<- rnaturalearth::ne_countries(country = c("Canada","United States of America"),scale = 50)

tm_shape(stns%>%filter(!is.na(NSE)))+tm_dots()+tm_shape(stns%>%filter(Official_ID == "09AB009"))+tm_dots(fill  = "red")+
  tm_shape(cst)+tm_borders()


# Niger
stns<-st_read("2.data/GOF/GRDC_africa.gpkg")

dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/africa/",
                       "1134460",
                       "_Q_Day.Cmd.txt"),
                sep = ";",skip = 36,header =  TRUE)%>%
  mutate(dt = ymd(YYYY.MM.DD ),
         year = year(dt),
         month = month(dt),
         yday = pmin(yday(dt),365),
         q = Value)

dat$q[dat$q==-999]<-NA

dat<-left_join(data.frame(dt = seq.Date(dat$dt[1],dat$dt[nrow(dat)],by = "1 day"))%>%
                 mutate(year = year(dt)),
               dat)
dat$fkDate<-ymd("2001-12-31")+days(dat$yday)

dat_yearly<-dat%>%
  group_by(year)%>%
  summarize(full = sum(is.na(q))==0)

dat<-dat%>%filter(year %in% dat_yearly$year[dat_yearly$full])

dat_sum<-dat%>%
  group_by(fkDate)%>%
  summarise(q = mean(q,na.rm = TRUE))
m_labels<-c("Jan ","Apr ","Jul ","Oct ","Jan  ")

p<-ggplot(dat,aes(x = fkDate,y = q))+
  geom_line(aes(group = year),col = "grey",alpha = 0.5,linewidth = 0.25)+
  geom_line(data = dat_sum,col = "black")+
  
  scale_x_date(name = "Date",
               # date_minor_breaks = "10 year",
               limits = ymd(c(20020101,20030101)),
               expand = c(0,0),
               # date_breaks = "3 month",
               # breaks
               breaks  = seq(from=as.Date("2002-01-01"),to=as.Date("2003-01-01"),by="3 month"),
               labels = m_labels)+
  scale_y_continuous(name = expression(Q~(m^3*s^-1)),
                     # limits = c(0,350),
                     expand = c(0,0)
                     # breaks = c(0,100,200,300)
                     # labels =c("0","10^5"),
                     # labels = parse(text = c("0","10^4"))
  )+
  theme_bw(base_size = 6)+
  # ggtitle("Rio Purus at Canutama (13880000)")+
  theme(panel.grid = element_blank())
p

ggsave(p,filename= "3.figures/Niger.png",width = 2,height = 2,dpi = 600)


stns<-st_transform(stns,"+proj=wintri +lon_0=0 +lat_1=50.467 +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs +type=crs")
tmap_mode("plot")
cst<- rnaturalearth::ne_countries(scale = 50,country = c("Niger","Mali"))
tm_shape(cst)+tm_borders()+
tm_shape(stns%>%filter(!is.na(NSE)))+tm_dots()+tm_shape(stns%>%filter(grdc_no == "1134460"))+tm_dots(fill  = "red")+
  tm_shape(cst)+tm_borders()

#Switzerland
stns<-st_read("2.data/GOF/camels_ch.gpkg")
it= which.max(stns$NSE)
dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-ch/camels_ch/timeseries/observation_based/CAMELS_CH_obs_based_",
                       stns$gauge_id[it],".csv"),
                # skip  =7,
                sep  =",")%>%
  mutate(dt = ymd(date),
         q = discharge_vol.m3.s.,
         yday = pmin(yday(dt),365),
         year = year(dt))

dat<-left_join(data.frame(dt = seq.Date(dat$dt[1],dat$dt[nrow(dat)],by = "1 day"))%>%
                 mutate(year = year(dt)),
               dat)
dat$fkDate<-ymd("2001-12-31")+days(dat$yday)

dat_yearly<-dat%>%
  group_by(year)%>%
  summarize(full = sum(is.na(q))==0)

dat<-dat%>%filter(year %in% dat_yearly$year[dat_yearly$full])

dat_sum<-dat%>%
  group_by(fkDate)%>%
  summarise(q = mean(q,na.rm = TRUE))
m_labels<-c("Jan ","Apr ","Jul ","Oct ","Jan  ")

p<-ggplot(dat,aes(x = fkDate,y = q))+
  geom_line(aes(group = year),col = "grey",alpha = 0.5,linewidth = 0.25)+
  geom_line(data = dat_sum,col = "black")+
  
  scale_x_date(name = "Date",
               # date_minor_breaks = "10 year",
               limits = ymd(c(20020101,20030101)),
               expand = c(0,0),
               # date_breaks = "3 month",
               # breaks
               breaks  = seq(from=as.Date("2002-01-01"),to=as.Date("2003-01-01"),by="3 month"),
               labels = m_labels)+
  scale_y_continuous(name = expression(Q~(m^3*s^-1)),
                     limits = c(0,99),
                     expand = c(0,0)
                     # breaks = c(0,100,200,300)
                     # labels =c("0","10^5"),
                     # labels = parse(text = c("0","10^4"))
  )+
  theme_bw(base_size = 6)+
  # ggtitle("Rio Purus at Canutama (13880000)")+
  theme(panel.grid = element_blank())
p

ggsave(p,filename= "3.figures/Massa.png",width = 2,height = 2,dpi = 600)




# Asia
#mekong
stns<-st_read("2.data/GOF/GRDC_asia.gpkg")

dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/asia/",
                       "2569002",
                       "_Q_Day.Cmd.txt"),
                sep = ";",skip = 36,header =  TRUE)%>%
  mutate(dt = ymd(YYYY.MM.DD ),
         year = year(dt),
         month = month(dt),
         yday = pmin(yday(dt),365),
         q = Value)

dat$q[dat$q==-999]<-NA

dat<-left_join(data.frame(dt = seq.Date(dat$dt[1],dat$dt[nrow(dat)],by = "1 day"))%>%
                 mutate(year = year(dt)),
               dat)
dat$fkDate<-ymd("2001-12-31")+days(dat$yday)

dat_yearly<-dat%>%
  group_by(year)%>%
  summarize(full = sum(is.na(q))==0)

dat<-dat%>%filter(year %in% dat_yearly$year[dat_yearly$full])

dat_sum<-dat%>%
  group_by(fkDate)%>%
  summarise(q = mean(q,na.rm = TRUE))
m_labels<-c("Jan ","Apr ","Jul ","Oct ","Jan  ")

p<-ggplot(dat,aes(x = fkDate,y = q))+
  geom_line(aes(group = year),col = "grey",alpha = 0.5,linewidth = 0.25)+
  geom_line(data = dat_sum,col = "black")+
  
  scale_x_date(name = "Date",
               # date_minor_breaks = "10 year",
               limits = ymd(c(20020101,20030101)),
               expand = c(0,0),
               # date_breaks = "3 month",
               # breaks
               breaks  = seq(from=as.Date("2002-01-01"),to=as.Date("2003-01-01"),by="3 month"),
               labels = m_labels)+
  scale_y_continuous(name = expression(Q~(m^3*s^-1)),
                     # limits = c(0,350),
                     limits = c(0,50000),
                     expand = c(0,0),
                     breaks = c(0,10000,20000,30000,40000,50000),
                     # labels =c("0","10^5"),
                     labels = parse(text = c("0","1%*%10^4","2%*%10^4","3%*%10^4","4%*%10^4","5%*%10^4"))
  )+
  theme_bw(base_size = 6)+
  # ggtitle("Rio Purus at Canutama (13880000)")+
  theme(panel.grid = element_blank())
p

ggsave(p,filename= "3.figures/Mekong.png",width = 2,height = 2,dpi = 600)

#Irrawaddy
dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/asia/",
                       "2260700",
                       "_Q_Day.Cmd.txt"),
                sep = ";",skip = 36,header =  TRUE)%>%
  mutate(dt = ymd(YYYY.MM.DD ),
         year = year(dt),
         month = month(dt),
         yday = pmin(yday(dt),365),
         q = Value)

dat$q[dat$q==-999]<-NA

dat<-left_join(data.frame(dt = seq.Date(dat$dt[1],dat$dt[nrow(dat)],by = "1 day"))%>%
                 mutate(year = year(dt)),
               dat)
dat$fkDate<-ymd("2001-12-31")+days(dat$yday)

dat_yearly<-dat%>%
  group_by(year)%>%
  summarize(full = sum(is.na(q))==0)

dat<-dat%>%filter(year %in% dat_yearly$year[dat_yearly$full])

dat_sum<-dat%>%
  group_by(fkDate)%>%
  summarise(q = mean(q,na.rm = TRUE))
m_labels<-c("Jan ","Apr ","Jul ","Oct ","Jan  ")

p<-ggplot(dat,aes(x = fkDate,y = q))+
  geom_line(aes(group = year),col = "grey",alpha = 0.5,linewidth = 0.25)+
  geom_line(data = dat_sum,col = "black")+
  
  scale_x_date(name = "Date",
               # date_minor_breaks = "10 year",
               limits = ymd(c(20020101,20030101)),
               expand = c(0,0),
               # date_breaks = "3 month",
               # breaks
               breaks  = seq(from=as.Date("2002-01-01"),to=as.Date("2003-01-01"),by="3 month"),
               labels = m_labels)+
  scale_y_continuous(name = expression(Q~(m^3*s^-1)),
                     limits = c(0,49000),
                     expand = c(0,0),
                     breaks = c(0,10000,20000,30000,40000),
                     # labels =c("0","10^5"),
                     labels = parse(text = c("0","1%*%10^4","2%*%10^4","3%*%10^4","4%*%10^4"))
  )+
  theme_bw(base_size = 6)+
  # ggtitle("Rio Purus at Canutama (13880000)")+
  theme(panel.grid = element_blank())
p

ggsave(p,filename= "3.figures/Irrawaddy.png",width = 2,height = 2,dpi = 600)


#Sokh
dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/asia/",
                       "2903430",
                       "_Q_Day.Cmd.txt"),
                sep = ";",skip = 36,header =  TRUE)%>%
  mutate(dt = ymd(YYYY.MM.DD ),
         year = year(dt),
         month = month(dt),
         yday = pmin(yday(dt),365),
         q = Value)

dat$q[dat$q==-999]<-NA

dat<-left_join(data.frame(dt = seq.Date(dat$dt[1],dat$dt[nrow(dat)],by = "1 day"))%>%
                 mutate(year = year(dt)),
               dat)
dat$fkDate<-ymd("2001-12-31")+days(dat$yday)

dat_yearly<-dat%>%
  group_by(year)%>%
  summarize(full = sum(is.na(q))==0)

dat<-dat%>%filter(year %in% dat_yearly$year[dat_yearly$full])

dat_sum<-dat%>%
  group_by(fkDate)%>%
  summarise(q = mean(q,na.rm = TRUE))
m_labels<-c("Jan ","Apr ","Jul ","Oct ","Jan  ")

p<-ggplot(dat,aes(x = fkDate,y = q))+
  geom_line(aes(group = year),col = "grey",alpha = 0.5,linewidth = 0.25)+
  geom_line(data = dat_sum,col = "black")+
  
  scale_x_date(name = "Date",
               # date_minor_breaks = "10 year",
               limits = ymd(c(20020101,20030101)),
               expand = c(0,0),
               # date_breaks = "3 month",
               # breaks
               breaks  = seq(from=as.Date("2002-01-01"),to=as.Date("2003-01-01"),by="3 month"),
               labels = m_labels)+
  scale_y_continuous(name = expression(Q~(m^3*s^-1)),
                     limits = c(0,190000),
                     expand = c(0,0),
                     breaks = c(0,100000,200000),
                     # labels =c("0","10^5"),
                     labels = parse(text = c("0","10^5","2%*%10^5"))
  )+
  theme_bw(base_size = 6)+
  # ggtitle("Rio Purus at Canutama (13880000)")+
  theme(panel.grid = element_blank())
p

ggsave(p,filename= "3.figures/Lena.png",width = 2,height = 2,dpi = 600)


#VAKHSH


dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/asia/",
                       "2517920",
                       "_Q_Day.Cmd.txt"),
                sep = ";",skip = 36,header =  TRUE)%>%
  mutate(dt = ymd(YYYY.MM.DD ),
         year = year(dt),
         month = month(dt),
         yday = pmin(yday(dt),365),
         q = Value)

dat$q[dat$q==-999]<-NA

dat<-left_join(data.frame(dt = seq.Date(dat$dt[1],dat$dt[nrow(dat)],by = "1 day"))%>%
                 mutate(year = year(dt)),
               dat)
dat$fkDate<-ymd("2001-12-31")+days(dat$yday)

dat_yearly<-dat%>%
  group_by(year)%>%
  summarize(full = sum(is.na(q))==0)

dat<-dat%>%filter(year %in% dat_yearly$year[dat_yearly$full])

dat_sum<-dat%>%
  group_by(fkDate)%>%
  summarise(q = mean(q,na.rm = TRUE))
m_labels<-c("Jan ","Apr ","Jul ","Oct ","Jan  ")

p<-ggplot(dat,aes(x = fkDate,y = q))+
  geom_line(aes(group = year),col = "grey",alpha = 0.5,linewidth = 0.25)+
  geom_line(data = dat_sum,col = "black")+
  
  scale_x_date(name = "Date",
               # date_minor_breaks = "10 year",
               limits = ymd(c(20020101,20030101)),
               expand = c(0,0),
               # date_breaks = "3 month",
               # breaks
               breaks  = seq(from=as.Date("2002-01-01"),to=as.Date("2003-01-01"),by="3 month"),
               labels = m_labels)+
  scale_y_continuous(name = expression(Q~(m^3*s^-1)),
                     limits = c(0,2450),
                     expand = c(0,0),
                     breaks = c(0,1000,2000),
                     # labels =c("0","10^5"),
                     # labels = parse(text = c("0","10^5","2%*%10^5"))
  )+
  theme_bw(base_size = 6)+
  # ggtitle("Rio Purus at Canutama (13880000)")+
  theme(panel.grid = element_blank())
p

ggsave(p,filename= "3.figures/Vakhsh.png",width = 2,height = 2,dpi = 600)


# australia
stns<-st_read("2.data/GOF/camels_AUS.gpkg")
stns_meta<-read.csv("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-aus-v2/01_id_name_metadata/id_name_metadata.csv")

streamflow_dat<-read.csv("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-aus-v2/streamflow_MLd.csv")

ID_str<-"221201"
if((substr(ID_str,1,1))%in% as.character(0:9)){ID_str = paste0("X",ID_str)}
dat<-streamflow_dat%>%
  select(year,month,day,all_of(ID_str))%>%
  mutate(dt = ymd(paste(year,month,day)),
         yday = pmin(yday(dt),365))

names(dat)[4]<-"q"

dat$q[dat$q== -99.99]<-NA
dat$q<-dat$q/1000
# dat$q[dat$q==-999]<-NA

dat<-left_join(data.frame(dt = seq.Date(dat$dt[1],dat$dt[nrow(dat)],by = "1 day"))%>%
                 mutate(year = year(dt)),
               dat)
dat$fkDate<-ymd("2001-12-31")+days(dat$yday)

dat_yearly<-dat%>%
  group_by(year)%>%
  summarize(full = sum(is.na(q))==0)

dat<-dat%>%filter(year %in% dat_yearly$year[dat_yearly$full])

dat_sum<-dat%>%
  group_by(fkDate)%>%
  summarise(q = mean(q,na.rm = TRUE))
m_labels<-c("Jan ","Apr ","Jul ","Oct ","Jan  ")

p<-ggplot(dat,aes(x = fkDate,y = q))+
  geom_line(aes(group = year),col = "grey",alpha = 0.5,linewidth = 0.25)+
  geom_line(data = dat_sum,col = "black")+
  
  scale_x_date(name = "Date",
               # date_minor_breaks = "10 year",
               limits = ymd(c(20020101,20030101)),
               expand = c(0,0),
               # date_breaks = "3 month",
               # breaks
               breaks  = seq(from=as.Date("2002-01-01"),to=as.Date("2003-01-01"),by="3 month"),
               labels = m_labels)+
  scale_y_continuous(name = expression(Q~(m^3*s^-1)),
                     limits = c(0,2450)/1000,
                     expand = c(0,0),
                     breaks = c(0,1000,2000)/1000,
                     # labels =c("0","10^5"),
                     # labels = parse(text = c("0","10^5","2%*%10^5"))
  )+
  theme_bw(base_size = 6)+
  # ggtitle("Rio Purus at Canutama (13880000)")+
  theme(panel.grid = element_blank())
p

ggsave(p,filename= "3.figures/West Branch River at Weeragua.png",width = 2,height = 2,dpi = 600)

