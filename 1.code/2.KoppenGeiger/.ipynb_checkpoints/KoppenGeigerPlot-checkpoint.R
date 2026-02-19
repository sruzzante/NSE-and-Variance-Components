# Create the violin plot in figure 2 (b)
# the NSE distribution by Koppen-Geiger zone



setwd("/home/ruzzante/projects/def-tgleeson/ruzzante/climatological_benchmarks/")
library(sf)
library(dplyr)
library(tmap)
library(ggplot2)
library(stringr)
library(terra)


stns<-rbind(
  st_read("2.data/GOF/camels_AUS.gpkg")%>%select(station_id,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%
    dplyr::rename(gauge_id = station_id)%>%
    left_join(read.csv("2.data/worldclim/camels_aus_v2_KG_zone.csv"),by = c("gauge_id" = "ID"))%>%
    mutate(src = "camels_aus_v2"),
  # 
  st_read("2.data/GOF/camels_br.gpkg")%>%select(gauge_id,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%
    left_join(read.csv("2.data/worldclim/camels_br_KG_zone.csv"),by = c("gauge_id" = "ID"))%>%
    mutate(src = "camels_br"),
  
  st_read("2.data/GOF/camels_ch.gpkg")%>%select(gauge_id,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%
    left_join(read.csv("2.data/worldclim/camels_ch_KG_zone.csv"),by = c("gauge_id" = "ID"))%>%
    mutate(src = "camels_ch"),
  
  st_read("2.data/GOF/camels_cl.gpkg")%>%select(gauge_id,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%
    # mutate(gauge_id = (gauge_id))%>%
    left_join(read.csv("2.data/worldclim/camels_cl_KG_zone.csv"),by = c("gauge_id" = "ID"))%>%
    mutate(src = "camels_cl"),
  
  st_read("2.data/GOF/caravan_de.gpkg")%>%select(gauge_id,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%
    left_join(read.csv("2.data/worldclim/caravan_de_KG_zone.csv"),by = c("gauge_id" = "ID"))%>%
    mutate(src = "caravan_de"),
  # 
  st_read("2.data/GOF/camels_dk.gpkg")%>%select(catch_id,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%
    dplyr::rename(gauge_id = catch_id)%>%
    left_join(read.csv("2.data/worldclim/camels_dk_KG_zone.csv"),by = c("gauge_id" = "ID"))%>%
    mutate(src = "camels_dk"),
  
  st_read("2.data/GOF/camels_es.gpkg")%>%select(gauge_id,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%st_make_valid()%>%
    left_join(read.csv("2.data/worldclim/camels_es_KG_zone.csv"),by = c("gauge_id" = "ID"))%>%
    mutate(src = "camels_es"),
  
  st_read("2.data/GOF/camels_US.gpkg")%>%select(gauge_id,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%
    left_join(read.csv("2.data/worldclim/camels_us_KG_zone.csv")%>%
                mutate(ID = str_pad(ID, width = 8,pad = "0",side = "left")),by = c("gauge_id" = "ID"))%>%
    mutate(src = "camels_us"),
  
  st_read("2.data/GOF/camels_fr.gpkg")%>%select(sta_code_h3,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%
    dplyr::rename(gauge_id = sta_code_h3)%>%
    left_join(read.csv("2.data/worldclim/camels_fr_KG_zone.csv"),by = c("gauge_id" = "ID"))%>%
    mutate(src = "camels_fr"),
  
  st_read("2.data/GOF/camels_gb.gpkg")%>%select(gauge_id,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%
    left_join(read.csv("2.data/worldclim/camels_gb_KG_zone.csv"),by = c("gauge_id" = "ID"))%>%
    mutate(src = "camels_gb"),
  
  st_read("2.data/GOF/camels_il.gpkg")%>%select(gauge_id,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%
    left_join(read.csv("2.data/worldclim/caravan_il_KG_zone.csv"),
              by = c("gauge_id" = "ID"))%>%
    mutate(src = "caravan_il"),
  
  st_read("2.data/GOF/camels_ind.gpkg")%>%select(gauge_id,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%
    left_join(read.csv("2.data/worldclim/camels_ind_KG_zone.csv")%>%
                mutate(ID = str_pad(ID,width = 5,pad = "0",side = "left")),
              by = c("gauge_id" = "ID"))%>%
    mutate(src = "camels_ind"),
  
  st_read("2.data/GOF/camels_lux.gpkg")%>%select(gauge_id, NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%
    left_join(read.csv("2.data/worldclim/camels_lux_KG_zone.csv")%>%
                mutate(gauge_id = paste0("ID_",str_pad(ID,width = 2,side = "left",pad = "0")))%>%
                select(!ID),
              by = c("gauge_id" = "gauge_id"))%>%
    mutate(src = "camels_lux"),
  # 
  st_read("2.data/GOF/GRDC_africa.gpkg")%>%select(grdc_no,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%
    dplyr::rename(gauge_id = grdc_no)%>%
    left_join(read.csv("2.data/worldclim/grdc_KG_zone.csv"),
              by = c("gauge_id" = "ID"))%>%
    mutate(src = "grdc_af"),
  
  st_read("2.data/GOF/GRDC_asia.gpkg")%>%select(grdc_no,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%st_transform("EPSG:4326")%>%
    dplyr::rename(gauge_id = grdc_no)%>%
    left_join(read.csv("2.data/worldclim/grdc_KG_zone.csv"),by = c("gauge_id" = "ID"))%>%
    mutate(src = "grdc_as"),
  
  st_read("2.data/GOF/GRDC_europe_a_q.gpkg")%>%select(grdc_no,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%st_transform("EPSG:4326")%>%
    dplyr::rename(gauge_id = grdc_no)%>%
    left_join(read.csv("2.data/worldclim/grdc_KG_zone.csv"),by = c("gauge_id" = "ID"))%>%
    mutate(src = "grdc_eur_a_q"),
  
  st_read("2.data/GOF/GRDC_europe_r_z.gpkg")%>%select(grdc_no,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%st_transform("EPSG:4326")%>%
    dplyr::rename(gauge_id = grdc_no)%>%
    left_join(read.csv("2.data/worldclim/grdc_KG_zone.csv"),by = c("gauge_id" = "ID"))%>%
    mutate(src = "grdc_eur_r_z"),
  # 
  st_read("2.data/GOF/GRDC_north_america.gpkg")%>%select(grdc_no,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%st_transform("EPSG:4326")%>%
    dplyr::rename(gauge_id = grdc_no)%>%
    left_join(read.csv("2.data/worldclim/grdc_KG_zone.csv"),by = c("gauge_id" = "ID"))%>%
    mutate(src = "grdc_na"),
  
  st_read("2.data/GOF/GRDC_south_america.gpkg")%>%select(grdc_no,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%st_transform("EPSG:4326")%>%
    dplyr::rename(gauge_id = grdc_no)%>%
    left_join(read.csv("2.data/worldclim/grdc_KG_zone.csv"),by = c("gauge_id" = "ID"))%>%
    mutate(src = "grdc_sa"),
  
  st_read("2.data/GOF/GRDC_south_west_pacific.gpkg")%>%select(grdc_no,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%st_transform("EPSG:4326")%>%
    dplyr::rename(gauge_id = grdc_no)%>%
    left_join(read.csv("2.data/worldclim/grdc_KG_zone.csv"),by = c("gauge_id" = "ID"))%>%
    mutate(src = "grdc_swp"),
  
  st_read("2.data/GOF/lamah_ce.gpkg")%>%select(govnr,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%
    dplyr::rename(gauge_id = govnr)%>%
    left_join(read.csv("2.data/worldclim/lamah_ce_KG_zone.csv"),by = c("gauge_id" = "ID"))%>%
    mutate(src = "lamah_ce"),
  
  st_read("2.data/GOF/lamah_ice.gpkg")%>%select(id,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_transform("EPSG:4326")%>%
    dplyr::rename(gauge_id = id)%>%
    left_join(read.csv("2.data/worldclim/lamah_ice_KG_zone.csv"),by = c("gauge_id" = "ID"))%>%
    mutate(src = "lamah_ice"),
  
  # 
    st_read("2.data/GOF/hysets.gpkg")%>%select(Official_ID,NSE,KGE,KGE_r,KGE_a,KGE_b)%>%st_set_crs("EPSG:4326")%>%
      dplyr::rename(gauge_id = Official_ID)%>%
      left_join(read.csv("2.data/worldclim/hysets_KG_zone.csv"),by = c("gauge_id" = "ID"))%>%
      mutate(src = "hysets")
)

stns<-filter(stns,!is.na(NSE))
summary(factor(stns$KG_zone))
summary(factor(stns$src[is.na(stns$KG_zone)]))
summary(factor(stns$src[!is.na(stns$KG_zone)]))
sum(!is.na(stns$KG_zone))
sum(is.na(stns$KG_zone))
tmap_mode("view")
tm_shape(stns%>%filter(is.na(KG_zone)))+tm_dots()


KG_rast<-terra::rast("../DATA/1.Spatial_data/global/clim_climate_precip_aridity_permafrost/clim1_Koppen_Geiger/Beck_KG_V1_present_0p0083.tif")
colorTable<-terra::coltab(KG_rast)%>%
  data.frame()%>%
  mutate(hex_code = rgb(red/255,green/255,blue/255))

leg<-read.delim("../DATA/1.Spatial_data/global/clim_climate_precip_aridity_permafrost/clim1_Koppen_Geiger/legend.txt",
                sep = "",
                skip = 2,
                header = FALSE,
                row.names = NULL,
                # col.names = 
)

x<-readLines("../DATA/1.Spatial_data/global/clim_climate_precip_aridity_permafrost/clim1_Koppen_Geiger/legend.txt",
)
x<-x[4:33]
leg<-str_split_fixed(x,regex("\\s+"),10)[,c(2,3,4,5,6,7,8,9)]%>%
  data.frame()%>%
  mutate(num = str_remove(X1,":")%>%as.numeric(),
         code = X2
  )
leg$X5[str_detect(leg$X5,regex("[0-9]"))]<-""
leg$X6[str_detect(leg$X6,regex("[0-9]"))]<-""
leg$X7[str_detect(leg$X7,regex("[0-9]"))]<-""
leg$X8[str_detect(leg$X8,regex("[0-9]"))]<-""
leg$desc<-paste(leg$X3,leg$X4,leg$X5,leg$X6,leg$X7,leg$X8)
leg<-select(leg,num,code,desc)

leg<-left_join(leg,colorTable,c("num" = "value"))

stns<-left_join(stns,leg%>%select(code,hex_code),
                by = c("KG_zone"="code"))

stns2<-stns%>%st_drop_geometry()%>%
  filter(!is.na(KG_zone))

stns2$KG_zone_desc<-(stns2$KG_zone)%>%
  plyr::mapvalues(from = leg$code,
                  to = leg$desc)


stn_med<-stns2%>%
  
  group_by(KG_zone)%>%
  summarize(NSE = mean(NSE),
            N = n())%>%
  filter(N>=50)
# arrange(-NSE)

stns2%>%
  filter(KG_zone%in% stn_med$KG_zone)%>%
  mutate(KG_zone = factor(KG_zone,levels = stn_med$KG_zone))%>%
  ggplot(aes(x = KG_zone,y = NSE))+geom_violin(scale = "width",aes(fill = KG_zone),
                                               bw=0.04,
                                               draw_quantiles=0.5
  )+
  scale_fill_manual(values = leg$hex_code,breaks = leg$code)+
  theme_bw()


# Try combining zones
# leg2<-leg%>%

stns2<-
  stns2%>%
  mutate(KG_zone_agg = dplyr::case_match(
    KG_zone,
    "Af" ~ "Af",
    "Aw" ~ "Aw",
    "Am" ~ "Am",
    "BSh" ~ "B_h",
    "BSk" ~ "B_k",
    "BWh"~ "B_h",
    "BWk" ~ "B_k",
    "Cfa" ~ "C_a",
    "Csa" ~ "C_a",
    "Cwa" ~ "C_a",
    "Cfb" ~ "C_b",
    "Csb" ~ "C_b",
    "Cwb" ~ "C_b",
    "Cfc" ~ "C_c",
    "Csc" ~ "C_c",
    "Cwc" ~ "C_c",
    
    
    "Dfa" ~ "D_a",
    "Dsa" ~ "D_a",
    "Dwa" ~ "D_a",
    "Dfb" ~ "D_b",
    "Dsb" ~ "D_b",
    "Dwb" ~ "D_b",
    "Dfc" ~ "D_cd",
    "Dsc" ~ "D_cd",
    "Dwc" ~ "D_cd",
    "Dfd" ~ "D_cd",
    "Dsd" ~ "D_cd",
    "Dwd" ~ "D_cd",
    "ET" ~ "ET"
    
  ))

leg_agg<-leg%>%
  
  mutate(code_agg = dplyr::case_match(
    code,
    "Af" ~ "Af",
    "Aw" ~ "Aw",
    "Am" ~ "Am",
    "BSh" ~ "B_h",
    "BSk" ~ "B_k",
    "BWh"~ "B_h",
    "BWk" ~ "B_k",
    "Cfa" ~ "C_a",
    "Csa" ~ "C_a",
    "Cwa" ~ "C_a",
    "Cfb" ~ "C_b",
    "Csb" ~ "C_b",
    "Cwb" ~ "C_b",
    "Cfc" ~ "C_c",
    "Csc" ~ "C_c",
    "Cwc" ~ "C_c",
    
    
    "Dfa" ~ "D_a",
    "Dsa" ~ "D_a",
    "Dwa" ~ "D_a",
    "Dfb" ~ "D_b",
    "Dsb" ~ "D_b",
    "Dwb" ~ "D_b",
    "Dfc" ~ "D_cd",
    "Dsc" ~ "D_cd",
    "Dwc" ~ "D_cd",
    "Dfd" ~ "D_cd",
    "Dsd" ~ "D_cd",
    "Dwd" ~ "D_cd",
    "ET" ~ "ET",
    "EF" ~ "EF"
    
  ))%>%
  group_by(code_agg)%>%
  summarize(across(red:blue, ~mean(.x)))%>%
  mutate(hex_code = rgb(red/255,green/255,blue/255))


stn_med_agg<-stns2%>%
  
  group_by(KG_zone_agg)%>%
  summarize(NSE = mean(NSE),
            N = n())

gg<-
  stns2%>%
  
  mutate(KG_zone_agg = factor(KG_zone_agg,levels = stn_med_agg$KG_zone_agg))%>%
  ggplot(aes(x = KG_zone_agg,y = NSE))+
  geom_violin(scale = "width",
              #aes(fill = KG_zone_agg),
              bw=0.04,
              draw_quantiles=0.5,
              linewidth = 0.1
  )+
  # scale_x_discrete(labels = NULL,name = NULL)+
  scale_y_continuous(breaks= c(0,0.5,1), name = "Benchmark NSE",
                     limits = c(-0.2,1),
                     expand = c(0,0))+
  # scale_fill_manual(values = leg_agg$hex_code,breaks = leg_agg$code_agg)+
  theme_bw(base_size = 6)+
  theme(legend.position = "none",axis.ticks = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

gg
ggsave(gg+scale_x_discrete(labels = NULL,name = NULL),
       filename = "3.figures/figure1_KoppenGeiger.png",
       width = 2.5, height = 0.75,
       dpi = 1200)





gg<-
  stns2%>%
  
  mutate(KG_zone_agg = factor(KG_zone_agg,levels = stn_med_agg$KG_zone_agg))%>%
  ggplot(aes(x = KG_zone_agg,y = KGE))+
  geom_violin(scale = "width",
              #aes(fill = KG_zone_agg),
              bw=0.04,
              draw_quantiles=0.5,
              linewidth = 0.1
  )+
  # scale_x_discrete(labels = NULL,name = NULL)+
  scale_y_continuous(breaks= c(0,0.5,1), name = "Benchmark KSE",
                     limits = c(-0.5,1),
                     expand = c(0,0),
                     oob = scales::oob_keep)+
  # scale_fill_manual(values = leg_agg$hex_code,breaks = leg_agg$code_agg)+
  theme_bw(base_size = 6)+
  theme(legend.position = "none",axis.ticks = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

gg
ggsave(gg+scale_x_discrete(labels = NULL,name = NULL),
       filename = "3.figures/figure1_KoppenGeiger_KGE.png",
       width = 2.5, height = 0.75,
       dpi = 1200)


ggplot(stns,aes(x = NSE,y = KGE,col = KGE_a-KGE_r))+geom_point()+
  scale_color_continuous()

ggplot(stns,aes(x = NSE,y = KGE,col = KGE_r))+geom_point()+
  scale_color_continuous()
ggplot(stns,aes(x = NSE,y = KGE,col = KGE_a))+geom_point()+
  scale_color_continuous()

ggplot(stns,aes(x = NSE,y = KGE_r,col = KGE_r-KGE_a))+geom_point()+
  scale_color_continuous()

ggplot(stns,aes(x = NSE,y = KGE_a,col = KGE_r-KGE_a))+geom_point()+
  scale_color_continuous()

ggplot(stns,aes(x = KGE_r,y = KGE_a,col = KGE_r-KGE_a))+geom_point()+
  scale_color_continuous()+
  geom_abline()

ggplot(stns2,aes(x = NSE,y = KGE,col = KGE_r-KGE_a))+geom_point()+
  scale_color_continuous()

ggplot(stns2%>%filter(src == "camels_br"),aes(x = NSE,y = KGE,col = MAT))+geom_point()


ggplot(stns2%>%filter(src == "camels_ch"),aes(x = NSE,y = KGE,col = MAT))+geom_point()

stns$KGE_est<-1-sqrt(stns$KGE_r^2-stns$NSE+2*stns$KGE_a*stns$KGE_r-2*stns$KGE_r-2*stns$KGE_a+2)


ggplot(stns,aes(KGE,KGE_est))+geom_point()
plot(stns$KGE,stns$KGE_est)

x<-stns%>%filter(NSE<0& KGE>0.1)
# hysets 11361000 #########


nc_dat<-ncdf4::nc_open("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/hysets/HYSETS_2023_update_QC_stations.nc")

stns<-read.csv("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/hysets/HYSETS_watershed_properties.txt")%>%
  select(Watershed_ID,Source,Name,Official_ID,Hydrometric_station_latitude,Hydrometric_station_longitude)%>%
  st_as_sf(coords = c("Hydrometric_station_longitude","Hydrometric_station_latitude"),remove = FALSE)
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
                             start = c(1, stns$Watershed_ID[stns$Official_ID == "11361000"]), # Start at the first time step and the desired catchment
                             count = c(time_length, 1)),     # Read all time steps for the catchment)
                dt = dt,
                year = yr,
                month = mnth,
                day = dy,
                yday = ydy
)

ggplot(dat,aes(dt,q))+geom_line()
###########

stns_tropical<-stns%>%
  filter(KG_zone %in% c("Af","Aw","Am"))

stns2%>%
  
  mutate(KG_zone_agg = factor(KG_zone_agg,levels = stn_med_agg$KG_zone_agg))%>%
  
  ggplot(aes(col = KG_zone_agg,x = NSE))+stat_ecdf()+
  # scale_x_discrete(labels = NULL,name = NULL)+
  # scale_y_continuous(breaks= c(0,0.5,1), name = "Benchmark NSE")+
  scale_color_manual(values = leg_agg$hex_code,breaks = leg_agg$code_agg)+
  theme_bw(base_size = 6)
# theme(legend.position = "none",axis.ticks = element_blank())

stns_veryHigh<-stns%>%filter(NSE>0.8)
tmap_mode("view")
tm_shape(stns_veryHigh)+tm_dots(fill ="KG_zone",
                                fill.scale = tm_scale(values = leg$hex_code,
                                                      levels  =leg$code),
                                size = 1)

stns_High<-stns%>%filter(NSE>0.75)


KG_rast<-terra::rast("../DATA/1.Spatial_data/global/clim_climate_precip_aridity_permafrost/clim1_Koppen_Geiger/Beck_KG_V1_present_0p083.tif")

tmap_mode("view")

watersheds_grdc<-st_read("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/shapefiles/2025-02-11 13-54-59.shp")
watersheds_grdc<-filter(watersheds_grdc,grdc_no %in% stns_veryHigh$gauge_id)

tm_shape(KG_rast)+tm_raster(col.scale = tm_scale(n.max = 31))+
  tm_shape(watersheds_grdc)+tm_borders()+
  tm_shape(stns_High)+tm_dots(fill ="KG_zone",
                              fill.scale = tm_scale(values = leg$hex_code,
                                                    levels  =leg$code),
                              size = 1)





KG_rast<-terra::rast("../DATA/1.Spatial_data/global/clim_climate_precip_aridity_permafrost/clim1_Koppen_Geiger/Beck_KG_V1_present_0p083.tif")

KG_rast_df<-data.frame(num_code = values(KG_rast))%>%
  filter(!(Beck_KG_V1_present_0p083%in% c(0,30)))

KG_rast_df$code<-plyr::mapvalues(KG_rast_df$Beck_KG_V1_present_0p083,from = leg$num,
                                 to = leg$code)
KG_rast_df<-KG_rast_df%>%
  mutate(code_agg = dplyr::case_match(
    code,
    "Af" ~ "Af",
    "Aw" ~ "Aw",
    "Am" ~ "Am",
    "BSh" ~ "B_h",
    "BSk" ~ "B_k",
    "BWh"~ "B_h",
    "BWk" ~ "B_k",
    "Cfa" ~ "C_a",
    "Csa" ~ "C_a",
    "Cwa" ~ "C_a",
    "Cfb" ~ "C_b",
    "Csb" ~ "C_b",
    "Cwb" ~ "C_b",
    "Cfc" ~ "C_c",
    "Csc" ~ "C_c",
    "Cwc" ~ "C_c",
    
    
    "Dfa" ~ "D_a",
    "Dsa" ~ "D_a",
    "Dwa" ~ "D_a",
    "Dfb" ~ "D_b",
    "Dsb" ~ "D_b",
    "Dwb" ~ "D_b",
    "Dfc" ~ "D_cd",
    "Dsc" ~ "D_cd",
    "Dwc" ~ "D_cd",
    "Dfd" ~ "D_cd",
    "Dsd" ~ "D_cd",
    "Dwd" ~ "D_cd",
    "ET" ~ "ET",
    "EF" ~ "EF"
    
  ))
summary(factor(KG_rast_df$code_agg))/nrow(KG_rast_df)
# plot(KG_rast)

summary(factor(stns2$KG_zone_agg))/nrow(stns2)

(summary(factor(stns2$KG_zone_agg))/nrow(stns2))/(summary(factor(KG_rast_df$code_agg))/nrow(KG_rast_df))
