# climatological benchmarks and variance components
# Author: Sacha Ruzzante
# sachawruzzante@gmail.com
# Last Update: 2025-08-06

# This script plots decomposed time series for quite a few example catchments.

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
#install.packages("tricolore")
# library(tricolore)


setwd("/home/ruzzante/projects/def-tgleeson/ruzzante/climatological_benchmarks/")
source('1.code/utils.R')



stns<-rbind(
  st_read("2.data/GOF/camels_AUS.gpkg")%>%select(station_id,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%
    left_join(read.csv("2.data/varComponents/camels_AUS.csv"))%>%
    dplyr::rename(gauge_id = station_id)%>%
    mutate(src = "camels_aus_v2"),
  # 
  st_read("2.data/GOF/camels_br.gpkg")%>%select(gauge_id,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%
    left_join(read.csv("2.data/varComponents/camels_br.csv"))%>%
    mutate(src = "camels_br"),
  
  st_read("2.data/GOF/camels_ch.gpkg")%>%select(gauge_id,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%
    left_join(read.csv("2.data/varComponents/camels_ch.csv"))%>%
    mutate(src = "camels_ch"),
  
  st_read("2.data/GOF/camels_cl.gpkg")%>%select(gauge_id,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%
    
    left_join(read.csv("2.data/varComponents/camels_cl.csv"))%>%
    mutate(src = "camels_cl"),
  
  st_read("2.data/GOF/caravan_de.gpkg")%>%select(gauge_id,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%
    left_join(read.csv("2.data/varComponents/caravan_de.csv"))%>%
    mutate(src = "caravan_de"),
  
  st_read("2.data/GOF/camels_dk.gpkg")%>%select(catch_id,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%
    left_join(read.csv("2.data/varComponents/camels_dk.csv"))%>%
    dplyr::rename(gauge_id = catch_id)%>%
    mutate(src = "camels_dk"),
  
  st_read("2.data/GOF/camels_es.gpkg")%>%select(gauge_id,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%st_make_valid()%>%
    left_join(read.csv("2.data/varComponents/camels_es.csv"))%>%
    mutate(src = "camels_es"),
  
  st_read("2.data/GOF/camels_US.gpkg")%>%select(gauge_id,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%
    left_join(read.csv("2.data/varComponents/camels_US.csv")%>%
                mutate(gauge_id = str_pad(gauge_id, width = 8,pad = "0",side = "left")))%>%
    mutate(src = "camels_us"),
  
  st_read("2.data/GOF/camels_fr.gpkg")%>%select(sta_code_h3,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%
    left_join(read.csv("2.data/varComponents/camels_fr.csv"))%>%
    dplyr::rename(gauge_id = sta_code_h3)%>%
    mutate(src = "camels_fr"),
  
  st_read("2.data/GOF/camels_gb.gpkg")%>%select(gauge_id,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%
    left_join(read.csv("2.data/varComponents/camels_gb.csv"))%>%
    mutate(src = "camels_gb"),
  
  st_read("2.data/GOF/camels_il.gpkg")%>%select(gauge_id,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%
    left_join(read.csv("2.data/varComponents/camels_il.csv")%>%
                mutate(src = "caravan_il")),
  
  st_read("2.data/GOF/camels_ind.gpkg")%>%select(gauge_id,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%
    left_join(read.csv("2.data/varComponents/camels_ind.csv")%>%
                mutate(gauge_id = str_pad(gauge_id,width = 5,pad = "0",side = "left")))%>%
    mutate(src = "camels_ind"),
  
  st_read("2.data/GOF/camels_lux.gpkg")%>%select(gauge_id, NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%
    left_join(read.csv("2.data/varComponents/camels_lux.csv"))%>%
    mutate(src = "camels_lux"),
  # 
  st_read("2.data/GOF/GRDC_africa.gpkg")%>%select(grdc_no,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%
    left_join(read.csv("2.data/varComponents/GRDC_africa.csv"))%>%
    dplyr::rename(gauge_id = grdc_no)%>%
    mutate(src = "grdc_af"),
  
  st_read("2.data/GOF/GRDC_asia.gpkg")%>%select(grdc_no,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%st_transform("EPSG:4326")%>%
    left_join(read.csv("2.data/varComponents/GRDC_asia.csv"))%>%
    dplyr::rename(gauge_id = grdc_no)%>%
    mutate(src = "grdc_as"),
  
  st_read("2.data/GOF/GRDC_europe_a_q.gpkg")%>%select(grdc_no,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%st_transform("EPSG:4326")%>%
    left_join(read.csv("2.data/varComponents/GRDC_europe_a_q.csv"))%>%
    dplyr::rename(gauge_id = grdc_no)%>%
    mutate(src = "grdc_eur_a_q"),
  
  st_read("2.data/GOF/GRDC_europe_r_z.gpkg")%>%select(grdc_no,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%st_transform("EPSG:4326")%>%
    left_join(read.csv("2.data/varComponents/GRDC_europe_r_z.csv"))%>%
    dplyr::rename(gauge_id = grdc_no)%>%
    mutate(src = "grdc_eur_r_z"),
  # 
  st_read("2.data/GOF/GRDC_north_america.gpkg")%>%select(grdc_no,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%st_transform("EPSG:4326")%>%
    left_join(read.csv("2.data/varComponents/GRDC_north_america.csv"))%>%
    dplyr::rename(gauge_id = grdc_no)%>%
    mutate(src = "grdc_na"),
  
  st_read("2.data/GOF/GRDC_south_america.gpkg")%>%select(grdc_no,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%st_transform("EPSG:4326")%>%
    left_join(read.csv("2.data/varComponents/GRDC_south_america.csv"))%>%
    dplyr::rename(gauge_id = grdc_no)%>%
    mutate(src = "grdc_sa"),
  
  st_read("2.data/GOF/GRDC_south_west_pacific.gpkg")%>%select(grdc_no,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%st_transform("EPSG:4326")%>%
    left_join(read.csv("2.data/varComponents/GRDC_south_west_pacific.csv"))%>%
    dplyr::rename(gauge_id = grdc_no)%>%
    mutate(src = "grdc_swp"),
  
  st_read("2.data/GOF/lamah_ce.gpkg")%>%select(govnr,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%
    left_join(read.csv("2.data/varComponents/lamah_ce.csv"))%>%
    dplyr::rename(gauge_id = govnr)%>%
    mutate(src = "lamah_ce"),
  
  st_read("2.data/GOF/lamah_ice.gpkg")%>%select(id,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%
    left_join(read.csv("2.data/varComponents/lamah_ice.csv"))%>%
    dplyr::rename(gauge_id = id)%>%
    mutate(src = "lamah_ice"),
  
  
  st_read("2.data/GOF/hysets.gpkg")%>%select(Official_ID,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_set_crs("EPSG:4326")%>%
    left_join(read.csv("2.data/varComponents/hysets.csv"))%>%
    dplyr::rename(gauge_id = Official_ID)%>%
    mutate(src = "hysets")
)
st_write(stns,"2.data/varComponents/stns_all.gpkg",append = FALSE)
stns<-st_read("2.data/varComponents/stns_all.gpkg")
stns<-filter(stns,!is.na(varSeas_fourier))
Tricolore(stns,  'varRem_fourier', 'varSeas_fourier','varInterannual_fourier',input_validation = TRUE)
colors_and_legend <- Tricolore(stns,  'varRem_fourier', 'varSeas_fourier','varInterannual_fourier')

# histograms ####

stns%>%
  st_drop_geometry()%>%
  select('varRem_fourier', 'varSeas_fourier','varInterannual_fourier')%>%
  pivot_longer(cols = c('varRem_fourier', 'varSeas_fourier','varInterannual_fourier'))%>%
  mutate(lab = factor(name,levels = c('varInterannual_fourier', 'varSeas_fourier','varRem_fourier'),
                      labels = c(expression(sigma[interannual] ^2 / sigma[o]^2),expression(sigma[seasonal] ^2 /sigma[o]^2),expression(sigma[irregular]^2/sigma[o]^2)))
  )%>%
  ggplot(aes(x = value))+
  geom_histogram(breaks = seq(0,1,0.02))+
  scale_x_continuous(name = "Variance Fraction")+
  facet_wrap("lab",ncol = 1, labeller = label_parsed)+
  theme_bw()

ggsave("3.figures/varianceFractionHists.png",width = 3,height = 4)


head(colors_and_legend$rgb)

stns$rgb<-colors_and_legend$rgb
tmap_mode("view")
tm_shape(stns%>%filter(!is.na(varSeas_fourier)))+tm_dots(fill = "rgb",size = 1)

tm_shape(stns)+tm_dots(fill = "varInterannual_stl",col = "varInterannual_stl",shape = 16,size = .75)

tm_shape(stns)+tm_dots(fill = "varInterannual_fourier",col = "varInterannual_fourier",shape = 16,size = .75)

tm_shape(stns)+tm_dots(fill = "varRem_fourier",shape = 16,size = .75, fill.scale = tm_scale(breaks = seq(0,1,0.1)))


tm_shape(stns)+tm_dots(fill = "varSeas_fourier",shape = 16,size = .75, fill.scale = tm_scale(breaks = seq(0,1,0.1)))

mean(stns$varSeas_fourier>0.5)

stns%>%
  mutate(highSeas = varSeas_fourier>0.8)%>%
  filter(highSeas)%>%
  tm_shape()+tm_dots(fill = "highSeas",
                     fill.scale = tm_scale_categorical(values = c("red","black")))


stns%>%
  mutate(highInterannual = varInterannual_fourier>0.8)%>%
  filter(highInterannual)%>%
  tm_shape()+tm_dots(fill = "highInterannual",
                     fill.scale = tm_scale_categorical(values = c("red","black")))

ggplot(stns,aes(x = varSeas_stl,NSE))+geom_point()



# Sturgeon Wier
nc_dat<-ncdf4::nc_open("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/hysets/HYSETS_2023_update_QC_stations.nc")
time_length <- nc_dat$dim$time$len
dt<-as.Date("1950-01-01")+
  days(ncvar_get(nc = nc_dat,
                 varid = "time"))

yr = year(dt)
mnth = month(dt)
dy = day(dt)
ydy = pmin(yday(dt),365)
stns_hysets<-read.csv("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/hysets/HYSETS_watershed_properties.txt")
it = which(stns_hysets$Official_ID=='05KG002')
dat<-data.frame(q= ncvar_get(nc = nc_dat,
                             varid = "discharge",
                             start = c(1, stns_hysets$Watershed_ID[it]), # Start at the first time step and the desired catchment
                             count = c(time_length, 1)),     # Read all time steps for the catchment)
                dt = dt,
                year = yr,
                month = mnth,
                day = dy,
                yday = ydy
)
# ggplot(dat,aes(dt,y=q))+geom_line()+scale_y_log10()
p<- stl_plot(dat,7,365)
p+theme(text = element_text(size = 7))

(p+theme(
  strip.background = element_blank(),
  strip.text.y = element_blank(),
  text = element_text(size = 7),
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  
))%>%
  ggsave(filename = "3.figures/stl_Sturgeon_Weir_05KG002.png",width = 2.23,height = 3)


p<- stl_plot(dat,7,365,include.text = FALSE)


(p+theme(
  strip.background = element_blank(),
  strip.text.y = element_blank(),
  text = element_text(size = 7),
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
))%>%
  ggsave(filename = "3.figures/stl_Sturgeon_Weir_05KG002_noText.png",width = 2.23,height = 3,dpi = 1200)


p<- fourier_plot(dat)


(p+theme(
  strip.background = element_blank(),
  strip.text.y = element_blank(),
  text = element_text(size = 7),
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  
))%>%
  ggsave(filename = "3.figures/fourier_Sturgeon_Weir_05KG002.png",width = 2.23,height = 3)


p<- fourier_plot(dat,include.text = FALSE)


(p+theme(
  strip.background = element_blank(),
  strip.text.y = element_blank(),
  text = element_text(size = 7),
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  
))%>%
  ggsave(filename = "3.figures/fourier_Sturgeon_Weir_05KG002_noText.png",width = 2.23,height = 3,dpi = 1200)






dat<-data.frame(q= ncvar_get(nc = nc_dat,
                             varid = "discharge",
                             start = c(1, stns_hysets$Watershed_ID[(stns_hysets$Official_ID == "06408700")]), # Start at the first time step and the desired catchment
                             count = c(time_length, 1)),     # Read all time steps for the catchment)
                dt = dt,
                year = yr,
                month = mnth,
                day = dy,
                yday = ydy
)


p<-fourier_plot(dat,round.dec = 5)

p

(p+ggtitle("Rhoads Fork near Rochford, SD")+
    # scale_y_continuous(name = "Q (cms)")+
    theme(
      # strip.background = element_blank(),
      # strip.text.y = element_blank(),
      text = element_text(size = 7),
      plot.background = element_rect(fill='white', color=NA), #transparent plot bg
      
    ))%>%
  ggsave(filename = "3.figures/fourier_Rhoads-Fork.png",width = 4,height = 3)



dat<-data.frame(q= ncvar_get(nc = nc_dat,
                             varid = "discharge",
                             start = c(1, stns_hysets$Watershed_ID[(stns_hysets$Official_ID == "02ZK003")]), # Start at the first time step and the desired catchment
                             count = c(time_length, 1)),     # Read all time steps for the catchment)
                dt = dt,
                year = yr,
                month = mnth,
                day = dy,
                yday = ydy
)


p<-fourier_plot(dat,round.dec = 2)

p

(p+ggtitle("Little Barachois River Near Placentia")+
    # scale_y_continuous(name = "Q (cms)")+
    theme(
      # strip.background = element_blank(),
      # strip.text.y = element_blank(),
      text = element_text(size = 7),
      plot.background = element_rect(fill='white', color=NA), #transparent plot bg
      
    ))%>%
  ggsave(filename = "3.figures/fourier_Little Barachois_02ZK003.png",width = 4,height = 3)


dat<-data.frame(q= ncvar_get(nc = nc_dat,
                             varid = "discharge",
                             start = c(1, stns_hysets$Watershed_ID[(stns_hysets$Official_ID == "08405150")]), # Start at the first time step and the desired catchment
                             count = c(time_length, 1)),     # Read all time steps for the catchment)
                dt = dt,
                year = yr,
                month = mnth,
                day = dy,
                yday = ydy
)


p<-fourier_plot(dat,round.dec = 2)

p

(p+ggtitle("Dark Canyon at Carlsbad, NM")+
    # scale_y_continuous(name = "Q (cms)")+
    theme(
      # strip.background = element_blank(),
      # strip.text.y = element_blank(),
      text = element_text(size = 7),
      plot.background = element_rect(fill='white', color=NA), #transparent plot bg
      
    ))%>%
  ggsave(filename = "3.figures/fourier_Dark Canyon at Carlsbad, NM 08405150 .png",width = 4,height = 3)





dat<-data.frame(q= ncvar_get(nc = nc_dat,
                             varid = "discharge",
                             start = c(1, stns_hysets$Watershed_ID[(stns_hysets$Official_ID == "09AC001")]), # Start at the first time step and the desired catchment
                             count = c(time_length, 1)),     # Read all time steps for the catchment)
                dt = dt,
                year = yr,
                month = mnth,
                day = dy,
                yday = ydy
)


p<-fourier_plot(dat,round.dec = 2)

p

(p+ggtitle("Takhini River near Whitehorse")+
    # scale_y_continuous(name = "Q (cms)")+
    theme(
      # strip.background = element_blank(),
      # strip.text.y = element_blank(),
      text = element_text(size = 7),
      plot.background = element_rect(fill='white', color=NA), #transparent plot bg
      
    ))%>%
  ggsave(filename = "3.figures/fourier_Takhini River near Whitehorse 09AC001.png",width = 4,height = 3)


# 06409000 - good option for display, 58% interannual
# also 06476000 - interannual change in variance
# or 13171620
# or 05KG007 (Amisk lake)




# arid, highly irregular
library(readxl)
stns_grdc<-read_xlsx("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/GRDC_Stations.xlsx",sheet = 1)

# 1304800 - OUED KERT

dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/africa/",
                       "1304800",
                       "_Q_Day.Cmd.txt"),
                sep = ";",skip = 36,header =  TRUE)%>%
  mutate(dt = ymd(YYYY.MM.DD ),
         year = year(dt),
         month = month(dt),
         yday = pmin(yday(dt),365),
         q = Value)

dat$q[dat$q==-999]<-NA

p<-fourier_plot(dat)

(p+theme(
  strip.background = element_blank(),
  strip.text.y = element_blank(),
  text = element_text(size = 7),
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
))%>%
  ggsave(filename = "3.figures/fourier_Oued_Kert.png",width = 2.23,height = 3)


p<-fourier_plot(dat,include.text = FALSE)

(p+theme(
  strip.background = element_blank(),
  strip.text.y = element_blank(),
  text = element_text(size = 7),
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
))%>%
  ggsave(filename = "3.figures/fourier_Oued_Kert_noText.png",width = 2.23,height = 3,dpi = 1200)



p<-stl_plot(dat,7,365)

(p+theme(
  strip.background = element_blank(),
  strip.text.y = element_blank(),
  text = element_text(size = 7),
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
))%>%
  ggsave(filename = "3.figures/stl_Oued_Kert.png",width = 2.23,height = 3)


p<-stl_plot(dat,7,365,include.text = FALSE)

(p+theme(
  strip.background = element_blank(),
  strip.text.y = element_blank(),
  text = element_text(size = 7),
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
))%>%
  ggsave(filename = "3.figures/stl_Oued_Kert_noText.png",width = 2.23,height = 3,dpi = 1200)



#camels-br

dat<-read.delim(sprintf("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-br/02_CAMELS_BR_streamflow_all_catchments/%s_streamflow.txt",
                        "13880000"),
                sep = "")%>%
  mutate(dt = ymd(paste(year,month,day)),
         yday = pmin(yday(dt),365))%>%
  filter(qual_flag == 1)

names(dat)[4]<-"q"
dat<-filter(dat,year>=1990&year<=2010)
p<-stl_plot(dat,7,365)
p


(p+theme(
  strip.background = element_blank(),
  strip.text.y = element_blank(),
  text = element_text(size = 7),
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
))%>%
  ggsave(filename = "3.figures/stl_Canutama.png",width = 2.23,height = 3)

p<-stl_plot(dat,7,365,include.text = FALSE)

(p+theme(
  strip.background = element_blank(),
  strip.text.y = element_blank(),
  text = element_text(size = 7),
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
))%>%
  ggsave(filename = "3.figures/stl_Canutama_noText.png",width = 2.23,height = 3,dpi = 1200)



p<-fourier_plot(dat)
p
(p+theme(
  strip.background = element_blank(),
  strip.text.y = element_blank(),
  text = element_text(size = 7),
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
))%>%
  ggsave(filename = "3.figures/fourier_Canutama.png",width = 2.23,height = 3)

p<-fourier_plot(dat,include.text = FALSE)

(p+theme(
  strip.background = element_blank(),
  strip.text.y = element_blank(),
  text = element_text(size = 7),
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
))%>%
  ggsave(filename = "3.figures/fourier_Canutama_noText.png",width = 2.23,height = 3,dpi = 1200)


# 15550000 - santa_isabel


dat<-read.delim(sprintf("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-br/02_CAMELS_BR_streamflow_all_catchments/%s_streamflow.txt",
                        "15550000"),
                sep = "")%>%
  mutate(dt = ymd(paste(year,month,day)),
         yday = pmin(yday(dt),365))%>%
  filter(qual_flag == 1)
names(dat)[4]<-"q"


p<-fourier_plot(dat)
p
(p+theme(
  strip.background = element_blank(),
  strip.text.y = element_blank(),
  text = element_text(size = 7),
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
))%>%
  ggsave(filename = "3.figures/fourier_santa_isabela.png",width = 2.23,height = 3)


p<-fourier_plot(dat,include.text = FALSE)

(p+theme(
  strip.background = element_blank(),
  strip.text.y = element_blank(),
  text = element_text(size = 7),
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
))%>%
  ggsave(filename = "3.figures/fourier_santa_isabela_noText.png",width = 2.23,height = 3,dpi = 1200)


stns_br<-read.csv("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-br/01_CAMELS_BR_attributes/camels_br_location.txt" ,
                  sep = "")

# porto_murtinho
dat<-read.delim(sprintf("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-br/02_CAMELS_BR_streamflow_all_catchments/%s_streamflow.txt",
                        "67100000"),
                sep = "")%>%
  mutate(dt = ymd(paste(year,month,day)),
         yday = pmin(yday(dt),365))%>%
  filter(qual_flag == 1)

names(dat)[4]<-"q"
fourier_plot(dat)


## hysets
it = which(stns_hysets$Official_ID=='14054500')
dat<-data.frame(q= ncvar_get(nc = nc_dat,
                             varid = "discharge",
                             start = c(1, stns_hysets$Watershed_ID[it]), # Start at the first time step and the desired catchment
                             count = c(time_length, 1)),     # Read all time steps for the catchment)
                dt = dt,
                year = yr,
                month = mnth,
                day = dy,
                yday = ydy
)
# ggplot(dat,aes(dt,y=q))+geom_line()+scale_y_log10()
fourier_plot(dat)



#05409890 - good for large irregular
# 09486300 - also - AZ

dat<-read.delim(paste0(
  "../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-fr/CAMELS_FR_time_series/daily/CAMELS_FR_tsd_",
  "G600061010",".csv"
),
skip  =7,sep = ";")%>%
  select(tsd_date,tsd_q_l)%>%
  mutate(q = tsd_q_l/1000,
         date = ymd(tsd_date),
         year = year(date),
         yday = pmin(yday(date),365))%>%
  mutate(dt = date)
# ggplot(dat,aes(date,y=q))+geom_line()
fourier_plot(dat)



p<-fourier_plot(dat,round.dec = 3)
p
(p+
    ggtitle("La Durdent à Vittefleur")+
    theme(
      # strip.background = element_blank(),
      # strip.text.y = element_blank(),
      text = element_text(size = 7),
      plot.background = element_rect(fill='white', color=NA), #transparent plot bg
    ))%>%
  ggsave(filename = "3.figures/fourier_G600061010.png",width = 4,height = 3)



# grdc
stns_grdc<-readxl::read_excel("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/GRDC_Stations.xlsx")

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
ggplot(dat,aes(dt,y=q))+geom_line()

p<-fourier_plot(dat,round.dec = -6)
p
(p+
    ggtitle("Irrawaddy River at Pyay")+
    theme(
      # strip.background = element_blank(),
      # strip.text.y = element_blank(),
      text = element_text(size = 7),
      plot.background = element_rect(fill='white', color=NA), #transparent plot bg
    ))%>%
  ggsave(filename = "3.figures/fourier_Irrawaddy - 2260700.png",width = 4,height = 3)



dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/asia/",
                       "2314400",
                       "_Q_Day.Cmd.txt"),
                sep = ";",skip = 36,header =  TRUE)%>%
  mutate(dt = ymd(YYYY.MM.DD ),
         year = year(dt),
         month = month(dt),
         yday = pmin(yday(dt),365),
         q = Value)

dat$q[dat$q==-999]<-NA

p<-fourier_plot(dat,round.dec = 1)
p

(p+
    ggtitle("Talgar River, Kazakhstan")+
    theme(
      # strip.background = element_blank(),
      # strip.text.y = element_blank(),
      text = element_text(size = 7),
      plot.background = element_rect(fill='white', color=NA), #transparent plot bg
    ))%>%
  ggsave(filename = "3.figures/fourier_Talgar- 2314400.png",width = 4,height = 3)


dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/south_west_pacific//",
                       "5863120",
                       "_Q_Day.Cmd.txt"),
                sep = ";",skip = 36,header =  TRUE)%>%
  mutate(dt = ymd(YYYY.MM.DD ),
         year = year(dt),
         month = month(dt),
         yday = pmin(yday(dt),365),
         q = Value)

dat$q[dat$q==-999]<-NA
ggplot(dat,aes(dt,y=q))+geom_line()
p<-fourier_plot(dat)

(p+
    ggtitle("Rangitaiki River")+
    theme(
      # strip.background = element_blank(),
      # strip.text.y = element_blank(),
      text = element_text(size = 7),
      plot.background = element_rect(fill='white', color=NA), #transparent plot bg
    ))%>%
  ggsave(filename = "3.figures/fourier_Rangitaiki-5863120.png",width = 4,height = 3)


dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/south_america/",
                       "3368100",
                       "_Q_Day.Cmd.txt"),
                sep = ";",skip = 36,header =  TRUE)%>%
  mutate(dt = ymd(YYYY.MM.DD ),
         year = year(dt),
         month = month(dt),
         yday = pmin(yday(dt),365),
         q = Value)

dat$q[dat$q==-999]<-NA
p<-fourier_plot(dat,round.dec = -3)

(p+ggtitle("Ysyry Paraguái")+
    theme(
      # strip.background = element_blank(),
      # strip.text.y = element_blank(),
      text = element_text(size = 7),
      plot.background = element_rect(fill='white', color=NA), #transparent plot bg
    ))%>%
  ggsave(filename = "3.figures/fourier_Paraguái-3368100.png",width = 4,height = 3,dpi = 1200)


dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/south_america/",
                       "3206720",
                       "_Q_Day.Cmd.txt"),
                sep = ";",skip = 36,header =  TRUE)%>%
  mutate(dt = ymd(YYYY.MM.DD ),
         year = year(dt),
         month = month(dt),
         yday = pmin(yday(dt),365),
         q = Value)

dat$q[dat$q==-999]<-NA
p<-fourier_plot(dat,round.dec = -7)
p
(p+ggtitle("Rio Orinoco at Puente Angostura")+
    theme(
      # strip.background = element_blank(),
      # strip.text.y = element_blank(),
      text = element_text(size = 7),
      plot.background = element_rect(fill='white', color=NA), #transparent plot bg
    ))%>%
  ggsave(filename = "3.figures/fourier_Orinoco-3206720.png",width = 4,height = 3,dpi = 1200)




dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/south_america/",
                       "3276800",
                       "_Q_Day.Cmd.txt"),
                sep = ";",skip = 36,header =  TRUE)%>%
  mutate(dt = ymd(YYYY.MM.DD ),
         year = year(dt),
         month = month(dt),
         yday = pmin(yday(dt),365),
         q = Value)

dat$q[dat$q==-999]<-NA
p<-fourier_plot(dat,round.dec = -5)
p
(p+ggtitle("Rio Santa Cruz at Charles Fuhr")+
    theme(
      # strip.background = element_blank(),
      # strip.text.y = element_blank(),
      text = element_text(size = 7),
      plot.background = element_rect(fill='white', color=NA), #transparent plot bg
    ))%>%
  ggsave(filename = "3.figures/fourier-Santa Cruz-3276800.png",width = 4,height = 3,dpi = 1200)


dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/africa/",
                       "1258202",
                       "_Q_Day.Cmd.txt"),
                sep = ";",skip = 36,header =  TRUE)%>%
  mutate(dt = ymd(YYYY.MM.DD ),
         year = year(dt),
         month = month(dt),
         yday = pmin(yday(dt),365),
         q = Value)

dat$q[dat$q==-999]<-NA
p<-fourier_plot(dat,round.dec = 3)
p
(p+ggtitle("Ugab River")+
    theme(
      # strip.background = element_blank(),
      # strip.text.y = element_blank(),
      text = element_text(size = 7),
      plot.background = element_rect(fill='white', color=NA), #transparent plot bg
    ))%>%
  ggsave(filename = "3.figures/fourier_Ugab River Namibia-1258202.png",width = 4,height = 3,dpi = 1200)


dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/africa/",
                       "1134700",
                       "_Q_Day.Cmd.txt"),
                sep = ";",skip = 36,header =  TRUE)%>%
  mutate(dt = ymd(YYYY.MM.DD ),
         year = year(dt),
         month = month(dt),
         yday = pmin(yday(dt),365),
         q = Value)

dat$q[dat$q==-999]<-NA
p<-fourier_plot(dat,round.dec = -5)
p
(p+ggtitle("Niger River at Dire")+
    theme(
      # strip.background = element_blank(),
      # strip.text.y = element_blank(),
      text = element_text(size = 7),
      plot.background = element_rect(fill='white', color=NA), #transparent plot bg
    ))%>%
  ggsave(filename = "3.figures/fourier_Niger River 1134700.png",width = 4,height = 3,dpi = 1200)

dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/africa/",
                       "1749100",
                       "_Q_Day.Cmd.txt"),
                sep = ";",skip = 36,header =  TRUE)%>%
  mutate(dt = ymd(YYYY.MM.DD ),
         year = year(dt),
         month = month(dt),
         yday = pmin(yday(dt),365),
         q = Value)

dat$q[dat$q==-999]<-NA
p<-fourier_plot(dat,round.dec = -6)
p
(p+ggtitle("Ubangi River at Bangui")+
    theme(
      # strip.background = element_blank(),
      # strip.text.y = element_blank(),
      text = element_text(size = 7),
      plot.background = element_rect(fill='white', color=NA), #transparent plot bg
    ))%>%
  ggsave(filename = "3.figures/fourier_Ubangi River 1749100.png",width = 4,height = 3,dpi = 1200)



dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/europe_a_q/",
                       "6731750",
                       "_Q_Day.Cmd.txt"),
                sep = ";",skip = 36,header =  TRUE)%>%
  mutate(dt = ymd(YYYY.MM.DD ),
         year = year(dt),
         month = month(dt),
         yday = pmin(yday(dt),365),
         q = Value)

dat$q[dat$q==-999]<-NA
p<-fourier_plot(dat,round.dec = 2)
p
(p+ggtitle("Gieddejohka River, Norway")+
    theme(
      # strip.background = element_blank(),
      # strip.text.y = element_blank(),
      text = element_text(size = 7),
      plot.background = element_rect(fill='white', color=NA), #transparent plot bg
    ))%>%
  ggsave(filename = "3.figures/fourier_Gieddejohka-6731750.png",width = 4,height = 3,dpi = 1200)
#camels-cl

stns_cl<-read.csv("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-cl/CAMELS_CL_v202201/catchment_attributes.csv",
                  quote = '')

streamflow_dat<-read.csv("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-cl/CAMELS_CL_v202201/q_m3s_day.csv")

dat<-streamflow_dat%>%
  select(year,date,all_of(paste0("X","4313001")))%>%
  mutate(dt = ymd(date),
         yday = pmin(yday(dt),365))

names(dat)[3]<-"q"

p<-fourier_plot(dat,round.dec = 3)

(p+ggtitle("Rio Cochiguaz En El Peñon")+
    theme(
      # strip.background = element_blank(),
      # strip.text.y = element_blank(),
      text = element_text(size = 7),
      plot.background = element_rect(fill='white', color=NA), #transparent plot bg
    ))%>%
  ggsave(filename = "3.figures/fourier_Rio Cochiguaz En El Peñon-4313001.png",width = 4,height = 3,dpi = 1200)



#lamah-ce
stns_lamah_ce<-st_read("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/lamah-ce/D_gauges/3_shapefiles/Gauges.shp")%>%
  filter(!country %in% c("DEU","CHE"))

dat<-read.delim(paste0(
  "../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/lamah-ce/D_gauges/2_timeseries/daily/ID_",
  stns_lamah_ce$ID[stns_lamah_ce$govnr==211722],".csv"
),
sep = ";")%>%
  
  mutate(q = qobs,
         dt = ymd(paste(YYYY,MM,DD)),
         year = year(dt),
         yday = pmin(yday(dt),365))
dat$q[dat$q<0]<-NA

fourier_plot(dat)

stl_plot(dat,7,365)
#caravan_de
stns_de<-read.csv("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/caravan-de/attributes/camelsde/attributes_other_camelsde.csv")%>%
  st_as_sf(coords = c("gauge_lon","gauge_lat"),crs = "EPSG:4326")




nc_dat<-nc_open(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/caravan-de/timeseries/netcdf/camelsde/",
                       "camelsde_DE812110",".nc"))
# gauge_id<-stns$gauge_id[it]
dat<-data.frame(
  q = ncvar_get(nc_dat,"streamflow"),
  dt = as.Date("1951-01-01") + days(ncvar_get(nc_dat,"date"))
)%>%
  mutate(yday = pmin(yday(dt),365),
         year = year(dt))

fourier_plot(dat)
stl_plot(dat,7,365)


#camels-gb


dat<-read.csv(paste0(
  "../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-gb/data/timeseries/CAMELS_GB_hydromet_timeseries_",
  38017,"_19701001-20150930.csv"
))%>%
  select(date,discharge_vol)%>%
  mutate(q = discharge_vol,
         dt = ymd(date),
         year = year(date),
         yday = pmin(yday(date),365))

fourier_plot(dat)

# grdc-africa

dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/africa/",
                       "1160523",
                       "_Q_Day.Cmd.txt"),
                sep = ";",skip = 36,header =  TRUE)%>%
  mutate(dt = ymd(YYYY.MM.DD ),
         year = year(dt),
         month = month(dt),
         yday = pmin(yday(dt),365),
         q = Value)
dat$q[dat$q<0]<-NA
fourier_plot(dat)
dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/africa/",
                       "1591401",
                       "_Q_Day.Cmd.txt"),
                sep = ";",skip = 36,header =  TRUE)%>%
  mutate(dt = ymd(YYYY.MM.DD ),
         year = year(dt),
         month = month(dt),
         yday = pmin(yday(dt),365),
         q = Value)
dat$q[dat$q<0]<-NA
fourier_plot(dat)




dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/asia/",
                       "2316200",
                       "_Q_Day.Cmd.txt"),
                sep = ";",skip = 36,header =  TRUE)%>%
  mutate(dt = ymd(YYYY.MM.DD ),
         year = year(dt),
         month = month(dt),
         yday = pmin(yday(dt),365),
         q = Value)
dat$q[dat$q<0]<-NA
fourier_plot(dat)


p<-fourier_plot(dat,round.dec = -3)

(p+ggtitle("Syr Darya at Tyumen-Aryk ")+
    theme(
      # strip.background = element_blank(),
      # strip.text.y = element_blank(),
      text = element_text(size = 7),
      plot.background = element_rect(fill='white', color=NA), #transparent plot bg
    ))%>%
  ggsave(filename = "3.figures/fourier_Syr Darya-2316200.png",width = 4,height = 3,dpi = 1200)




dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/asia/",
                       "2180712",
                       "_Q_Day.Cmd.txt"),
                sep = ";",skip = 36,header =  TRUE)%>%
  mutate(dt = ymd(YYYY.MM.DD ),
         year = year(dt),
         month = month(dt),
         yday = pmin(yday(dt),365),
         q = Value)
dat$q[dat$q<0]<-NA
fourier_plot(dat)


dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/europe_a_q/",
                       "6854203",
                       "_Q_Day.Cmd.txt"),
                sep = ";",skip = 36,header =  TRUE)%>%
  mutate(dt = ymd(YYYY.MM.DD ),
         year = year(dt),
         month = month(dt),
         yday = pmin(yday(dt),365),
         q = Value)
dat$q[dat$q<0]<-NA
fourier_plot(dat)



dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/south_west_pacific//",
                       "5867710",
                       "_Q_Day.Cmd.txt"),
                sep = ";",skip = 36,header =  TRUE)%>%
  mutate(dt = ymd(YYYY.MM.DD ),
         year = year(dt),
         month = month(dt),
         yday = pmin(yday(dt),365),
         q = Value)

dat$q[dat$q==-999]<-NA
p<-fourier_plot(dat,round.dec = 0)
p
(p+
    ggtitle("Mawheranui (Grey River) at New Waipuna")+
    theme(
      # strip.background = element_blank(),
      # strip.text.y = element_blank(),
      text = element_text(size = 7),
      plot.background = element_rect(fill='white', color=NA), #transparent plot bg
    ))%>%
  ggsave(filename = "3.figures/fourier_Mawheranui-5867710.png",width = 4,height = 3)


# Blatten

dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-ch/camels_ch/timeseries/observation_based/CAMELS_CH_obs_based_",
                       2268,".csv"),
                # skip  =7,
                sep  =",")%>%
  mutate(dt = ymd(date),
         q = discharge_vol.m3.s.,
         yday = pmin(yday(dt),365),
         year = year(dt))

fourier_plot(dat)



# aus
streamflow_dat<-read.csv("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-aus-v2/streamflow_MLd.csv")

dat<-streamflow_dat%>%
  select(year,month,day,all_of("X707002"))%>%
  mutate(dt = ymd(paste(year,month,day)),
         yday = pmin(yday(dt),365))

names(dat)[4]<-"q"

fourier_plot(dat)


# hysets

nc_dat<-ncdf4::nc_open("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/hysets/HYSETS_2023_update_QC_stations.nc")

stns_hysets<-read.csv("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/hysets/HYSETS_watershed_properties.txt")%>%
  select(Watershed_ID,Source,Name,Official_ID,Hydrometric_station_latitude,Hydrometric_station_longitude)
time_length <- nc_dat$dim$time$len
dt<-as.Date("1950-01-01")+
  days(ncvar_get(nc = nc_dat,
                 varid = "time"))

yr = year(dt)
mnth = month(dt)
dy = day(dt)
ydy = pmin(yday(dt),365)


dat<-data.frame(q= ncvar_get(nc = nc_dat,
                             varid = "discharge",
                             start = c(1, stns_hysets$Watershed_ID[(stns_hysets$Official_ID == "03HA009")]), # Start at the first time step and the desired catchment
                             count = c(time_length, 1)),     # Read all time steps for the catchment)
                dt = dt,
                year = yr,
                month = mnth,
                day = dy,
                yday = ydy
)


fourier_plot(dat)



#camels-ind # could be best one!
stns_ind<-read.csv("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-ind/attributes_csv/camels_ind_topo.csv")
stns_ind_info<-read.csv("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-ind/attributes_csv/camels_ind_name.csv")
streamflow<-read.csv("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-ind/streamflow_timeseries/streamflow_observed.csv")
streamflow$dt<-ymd(paste(streamflow$year,streamflow$month,streamflow$day))

q<-streamflow[,paste0("X",as.numeric(04012))]

dt<-streamflow$dt
# nc_close(nc)
dat<-data.frame(dt,q)%>%
  mutate(year = year(dt),
         month = month(dt),
         yday = pmin(yday(dt),365))


p<-fourier_plot(dat, round.dec = 0)
p
(p+ggtitle("Haliya River")+
    theme(
      # strip.background = element_blank(),
      # strip.text.y = element_blank(),
      text = element_text(size = 7),
      plot.background = element_rect(fill='white', color=NA), #transparent plot bg
    ))%>%
  ggsave(filename = "3.figures/fourier_Haliya River 04012.png",width = 4,height = 3,dpi = 1200)



# plots for flowchart
streamflow_dat<-read.csv("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-aus-v2/streamflow_MLd.csv")

dat<-streamflow_dat%>%
  select(year,month,day,all_of("X226218"))%>%
  mutate(dt = ymd(paste(year,month,day)),
         yday = pmin(yday(dt),365))

names(dat)[4]<-"q"

fourier_plot(dat)


# 	gauge_name
# 223	Little Ouse at Abbey Heath

dat<-read.csv(paste0(
  "../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-gb/data/timeseries/CAMELS_GB_hydromet_timeseries_",
  33034,"_19701001-20150930.csv"
))%>%
  select(date,discharge_vol)%>%
  mutate(q = discharge_vol,
         dt = ymd(date),
         year = year(date),
         yday = pmin(yday(date),365))

fourier_plot(dat)



dat$t_diff<-c(NA,diff.Date(dat$dt))

RunLengths<-rle(dat$t_diff)
if(max(RunLengths$lengths[RunLengths$values==1],na.rm = TRUE)<3650){
  return(c(
    varSeas = NA,
    varRem =  NA,
    varInterannual = NA
  ))
}
# take longest continuous run of data 
RunLengths$start<-cumsum(RunLengths$lengths)
maxRun<-which.max(RunLengths$lengths)


dat<-dat[(RunLengths$start[maxRun-1]):
           (RunLengths$start[(maxRun)]),]

# calculate climatological mean of observed data
dat_clim<-dat%>%
  group_by(yday)%>%
  summarise(seasonal.avg=mean(q),
            N = n())
dat<-left_join(dat,dat_clim,by = "yday")

# calculate anomalies
dat$q.anomaly<-dat$q-dat$seasonal.avg

ggplot(data = dat,aes(dt,q.anomaly))+geom_line()+ scale_x_date(name = NULL)+
  theme_bw()+
  scale_y_continuous(name =expression(Anomaly~(m^3/s)),limits = c(-7,23) )

ggsave("3.figures/flowchart/Q.anomaly.png",width = 3,height = 2,dpi = 1200)

ggplot(data = dat,aes(dt,seasonal.avg))+geom_line()+ scale_x_date(name = NULL)+
  theme_bw()+
  scale_y_continuous(name =expression(Seasonal~(m^3/s)),limits = c(-10,20) )

ggsave("3.figures/flowchart/Q.seasonal.png",width = 3,height = 2,dpi = 1200)

# 
# 



# run
FFT<-fft(dat$q.anomaly)

# Need to centre the frequencies 
freq <- (0:(nrow(dat)- 1)) / nrow(dat)
freq <- ifelse(freq > 0.5, freq - 1, freq) * 365  # centered frequency for symmetry

plot(abs(FFT))

df_FFT = data.frame(freq,FFT_Mag = Mod(FFT))%>%
  arrange(freq)%>%
  
  filter(freq>=0)%>%
  mutate(PSD = FFT_Mag^2)

conversion_factor = var(dat$q.anomaly)/(sum(df_FFT$PSD)*diff(df_FFT$freq)[1])

df_FFT$PSD = df_FFT$PSD*conversion_factor
# 
# sp <- spec.pgram(dat$q.anomaly, taper = 0, fast = TRUE, plot = FALSE, demean = FALSE, detrend = TRUE,)
# plot(sp)
# sum(sp$spec) * diff(sp$freq)[1] * 365
# df_FFT$PSD = sp$spec[1:nrow(df_FFT)]

ggplot(df_FFT,aes(freq,PSD))+
  geom_line()+
  # geom_line(aes(y = FFT_mag_smooth),col = "blue")+
  # geom_smooth()+
  geom_vline(xintercept = 2,color = "red")+
  scale_x_sqrt(name = expression(Frequency~(year^-1)),
               breaks = c(0,2,50,100,200))+
  scale_y_continuous(expression(Power~Spectral~Density~(m^2*s^-2*year)))+
  
  theme_bw(base_size = 7)+
  theme(panel.grid.minor = element_blank())
ggsave("3.figures/flowchart/Periodogram.png",width = 3,height = 2,dpi = 1200)

sum(df_FFT$PSD[df_FFT$freq>2])*diff(df_FFT$freq)[1]
sum(df_FFT$PSD[df_FFT$freq<=2])*diff(df_FFT$freq)[1]

var(dat$q.anomaly)
sum(df_FFT$PSD)*diff(df_FFT$freq)[1]

var(dat$q.anomaly)/(sum(df_FFT$PSD)*diff(df_FFT$freq)[1])

FFT_1<-FFT
FFT_1[c(1:50,52:length(FFT_1))]<-0

dat_FFT1 = Re(fft(FFT_1,inverse = TRUE))/nrow(dat)
plot(dat_FFT1)
var(dat_FFT1)
sum(Mod(FFT_1)^2)/var(dat_FFT1)*diff(df_FFT$freq)[1]/86400

# Define cutoff
cutoff <- 2  # frequency cutoff ( cycles/year)

# interannual component has frequencies below or at cutoff
FFT_interannual<-FFT
FFT_interannual[abs(freq) > cutoff]<-0
# irregular component has frequencies above cutoff
FFT_irregular<-FFT
FFT_irregular[abs(freq) <= cutoff]<-0


# compute inverse FFT of the two components
dat$remainder<-Re(fft(FFT_irregular,inverse = TRUE))/nrow(dat)
dat$interannual<-Re(fft(FFT_interannual,inverse = TRUE))/nrow(dat)




ggplot(data = dat,aes(dt,remainder))+geom_line()+ scale_x_date(name = NULL)+
  theme_bw()+
  scale_y_continuous(name =expression(Irregular~(m^3/s)),limits = c(-10,20) )

ggsave("3.figures/flowchart/Q.Irregular.png",width = 3,height = 2,dpi = 1200)


ggplot(data = dat,aes(dt,interannual))+geom_line()+ scale_x_date(name = NULL)+
  theme_bw()+
  scale_y_continuous(name =expression(Interannual~(m^3/s)),limits = c(-10,20) )

ggsave("3.figures/flowchart/Q.Interannual.png",width = 3,height = 2,dpi = 1200)
var(dat$remainder)/var(dat$q  )

var(dat$interannual)/var(dat$q  )







ggplot(data = dat,aes(dt,q))+geom_line()+ scale_x_date(name = NULL)+
  theme_bw()+
  scale_y_continuous(name =expression(Q~(m^3/s)),limits = c(0,30) )
ggsave("3.figures/flowchart/Q.png",width = 3,height = 2,dpi = 1200)

# classical and STL plots
# sturgeon weir

nc_dat<-ncdf4::nc_open("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/hysets/HYSETS_2023_update_QC_stations.nc")
time_length <- nc_dat$dim$time$len
dt<-as.Date("1950-01-01")+
  days(ncvar_get(nc = nc_dat,
                 varid = "time"))

yr = year(dt)
mnth = month(dt)
dy = day(dt)
ydy = pmin(yday(dt),365)
stns_hysets<-read.csv("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/hysets/HYSETS_watershed_properties.txt")
it = which(stns_hysets$Official_ID=='05KG002')
dat<-data.frame(q= ncvar_get(nc = nc_dat,
                             varid = "discharge",
                             start = c(1, stns_hysets$Watershed_ID[it]), # Start at the first time step and the desired catchment
                             count = c(time_length, 1)),     # Read all time steps for the catchment)
                dt = dt,
                year = yr,
                month = mnth,
                day = dy,
                yday = ydy
)
# ggplot(dat,aes(dt,y=q))+geom_line()+scale_y_log10()
p<- stl_plot(dat,7,365)
p+theme(text = element_text(size = 7))

(p+theme(
  strip.background = element_blank(),
  strip.text.y = element_blank(),
  text = element_text(size = 7),
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  
))%>%
  ggsave(filename = "3.figures/stl_Sturgeon_Weir_05KG002.png",width = 2.23,height = 3)



p<- stl_plot(dat,7,365,include.text = FALSE)


(p+theme(
  strip.background = element_blank(),
  strip.text.y = element_blank(),
  text = element_text(size = 7),
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
))%>%
  ggsave(filename = "3.figures/stl_Sturgeon_Weir_05KG002_noText.png",width = 2.23,height = 3,dpi = 1200)

# Sturgeon Weir classical

p<- classical_plot  (dat,include.text = FALSE)


(p+theme(
  strip.background = element_blank(),
  strip.text.y = element_blank(),
  text = element_text(size = 7),
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
))%>%
  ggsave(filename = "3.figures/clas_Sturgeon_Weir_05KG002_noText.png",width = 2.23,height = 3,dpi = 1200)


p<- classical_plot  (dat,include.text = T)


(p+theme(
  strip.background = element_blank(),
  strip.text.y = element_blank(),
  text = element_text(size = 7),
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
))%>%
  ggsave(filename = "3.figures/clas_Sturgeon_Weir_05KG002.png",width = 2.23,height = 3,dpi = 1200)


# Candeias


dat<-read.delim(sprintf("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-br/02_CAMELS_BR_streamflow_all_catchments/%s_streamflow.txt",
                        "15550000"),
                sep = "")%>%
  mutate(dt = ymd(paste(year,month,day)),
         yday = pmin(yday(dt),365))%>%
  filter(qual_flag == 1)
names(dat)[4]<-"q"

p<- stl_plot(dat,7,365)
p+theme(text = element_text(size = 7))

(p+theme(
  strip.background = element_blank(),
  strip.text.y = element_blank(),
  text = element_text(size = 7),
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  
))%>%
  ggsave(filename = "3.figures/stl_Santa_Isabel.png",width = 2.23,height = 3)



p<- stl_plot(dat,7,365,include.text = FALSE)


(p+theme(
  strip.background = element_blank(),
  strip.text.y = element_blank(),
  text = element_text(size = 7),
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
))%>%
  ggsave(filename = "3.figures/stl_Santa_Isabel_noText.png",width = 2.23,height = 3,dpi = 1200)

# Candeias classical

p<- classical_plot  (dat,include.text = FALSE)


(p+theme(
  strip.background = element_blank(),
  strip.text.y = element_blank(),
  text = element_text(size = 7),
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
))%>%
  ggsave(filename = "3.figures/clas_Santa_Isabel_noText.png",width = 2.23,height = 3,dpi = 1200)


p<- classical_plot  (dat,include.text = T)


(p+theme(
  strip.background = element_blank(),
  strip.text.y = element_blank(),
  text = element_text(size = 7),
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
))%>%
  ggsave(filename = "3.figures/clas_Santa_Isabel.png",width = 2.23,height = 3,dpi = 1200)



# Oued Kert


dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/africa/",
                       "1304800",
                       "_Q_Day.Cmd.txt"),
                sep = ";",skip = 36,header =  TRUE)%>%
  mutate(dt = ymd(YYYY.MM.DD ),
         year = year(dt),
         month = month(dt),
         yday = pmin(yday(dt),365),
         q = Value)



dat$q[dat$q==-999]<-NA

p<- stl_plot(dat,7,365)
p+theme(text = element_text(size = 7))

(p+theme(
  strip.background = element_blank(),
  strip.text.y = element_blank(),
  text = element_text(size = 7),
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  
))%>%
  ggsave(filename = "3.figures/stl_Oued_Kert.png",width = 2.23,height = 3)



p<- stl_plot(dat,7,365,include.text = FALSE)


(p+theme(
  strip.background = element_blank(),
  strip.text.y = element_blank(),
  text = element_text(size = 7),
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
))%>%
  ggsave(filename = "3.figures/stl_Oued_Kert_noText.png",width = 2.23,height = 3,dpi = 1200)

# Oued Kert classical

p<- classical_plot  (dat,include.text = FALSE)


(p+theme(
  strip.background = element_blank(),
  strip.text.y = element_blank(),
  text = element_text(size = 7),
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
))%>%
  ggsave(filename = "3.figures/clas_Oued_Kert_noText.png",width = 2.23,height = 3,dpi = 1200)


p<- classical_plot  (dat,include.text = T)


(p+theme(
  strip.background = element_blank(),
  strip.text.y = element_blank(),
  text = element_text(size = 7),
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
))%>%
  ggsave(filename = "3.figures/clas_Oued_Kert.png",width = 2.23,height = 3,dpi = 1200)





# grdc - swp

dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/south_west_pacific//",
                       "5863120",
                       "_Q_Day.Cmd.txt"),
                sep = ";",skip = 36,header =  TRUE)%>%
  mutate(dt = ymd(YYYY.MM.DD ),
         year = year(dt),
         month = month(dt),
         yday = pmin(yday(dt),365),
         q = Value)

fourier_plot(dat)



dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/north_america/",
                       "4185053",
                       "_Q_Day.Cmd.txt"),
                sep = ";",skip = 36,header =  TRUE)%>%
  mutate(dt = ymd(YYYY.MM.DD ),
         year = year(dt),
         month = month(dt),
         yday = pmin(yday(dt),365),
         q = Value)

dat$q[dat$q==-999]<-NA
fourier_plot(dat)


dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/asia/",
                       "2124300",
                       "_Q_Day.Cmd.txt"),
                sep = ";",skip = 36,header =  TRUE)%>%
  mutate(dt = ymd(YYYY.MM.DD ),
         year = year(dt),
         month = month(dt),
         yday = pmin(yday(dt),365),
         q = Value)

dat$q[dat$q==-999]<-NA
fourier_plot(dat)
dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/africa/",
                       "1396210",
                       "_Q_Day.Cmd.txt"),
                sep = ";",skip = 36,header =  TRUE)%>%
  mutate(dt = ymd(YYYY.MM.DD ),
         year = year(dt),
         month = month(dt),
         yday = pmin(yday(dt),365),
         q = Value)

dat$q[dat$q==-999]<-NA
fourier_plot(dat)

stns$totVar_STL = stns$varSeas_stl+stns$varInterannual_stl+stns$varRem_stl
stns$totVar_CLAS = stns$varSeas_clas+stns$varInterannual_clas+stns$varRem_clas
# 09083000
# THOMPSON CREEK NEAR CARBONDALE, CO.
nc_dat<-ncdf4::nc_open("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/hysets/HYSETS_2023_update_QC_stations.nc")
time_length <- nc_dat$dim$time$len
dt<-as.Date("1950-01-01")+
  days(ncvar_get(nc = nc_dat,
                 varid = "time"))

yr = year(dt)
mnth = month(dt)
dy = day(dt)
ydy = pmin(yday(dt),365)
stns_hysets<-read.csv("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/hysets/HYSETS_watershed_properties.txt")
it = which(stns_hysets$Official_ID=='10080000')
dat<-data.frame(q= ncvar_get(nc = nc_dat,
                             varid = "discharge",
                             start = c(1, stns_hysets$Watershed_ID[it]), # Start at the first time step and the desired catchment
                             count = c(time_length, 1)),     # Read all time steps for the catchment)
                dt = dt,
                year = yr,
                month = mnth,
                day = dy,
                yday = ydy
)
# ggplot(dat,aes(dt,y=q))+geom_line()+scale_y_log10()
p<- stl_plot(dat,7,365)
p+theme(text = element_text(size = 7))

(p+theme(
  strip.background = element_blank(),
  strip.text.y = element_blank(),
  text = element_text(size = 7),
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  
))%>%
  ggsave(filename = "3.figures/stl_THOMPSON_CREEK_09083000.png",width = 2.23,height = 3)

p<- fourier_plot(dat)
p+theme(text = element_text(size = 7))

(p+theme(
  strip.background = element_blank(),
  strip.text.y = element_blank(),
  text = element_text(size = 7),
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  
))%>%
  ggsave(filename = "3.figures/stl_THOMPSON_CREEK_09083000.png",width = 2.23,height = 3)



#compare variance fractions by three methods

stns%>%
st_drop_geometry()%>%
  filter(!is.na(varSeas_clas))%>%
           summarize(across(varSeas_stl:varRem_fourier,~median(.x)),
                     N=n())%>%
  pivot_longer(cols = varSeas_stl:varRem_fourier)%>%
  arrange(name)
