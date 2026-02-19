# get climatology for camels-us


library(dplyr)
library(ggplot2)
library(sf)
# library(tmap)
library(tictoc)
library(stringr)
library(tidyr)
library(lubridate)
library(nngeo)
setwd("/home/ruzzante/projects/def-tgleeson/ruzzante/climatological_benchmarks/")

sf_use_s2(FALSE)

stns<-read.delim("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-us/camels_topo.txt",sep = ";",colClasses = c("character","numeric","numeric","numeric","numeric","numeric"))%>%
  st_as_sf(coords = c("gauge_lon","gauge_lat"),crs = "EPSG:4326")%>%
  left_join(read.delim("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-us/camels_name.txt",sep = ";",colClasses = c("character","character","character")))

hysets<-read.csv("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/hysets/HYSETS_watershed_properties.txt")
stns<-stns%>%filter(!gauge_id%in%hysets$Official_ID)


watersheds<-rbind(
  st_read("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-us/basin_dataset_public_v1p2/shapefiles/merge/Region_01_nhru_simplify_100.shp")%>%select(GAGEID),
  st_read("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-us/basin_dataset_public_v1p2/shapefiles/merge/Region_02_nhru_simplify_100.shp")%>%select(GAGEID),
  st_read("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-us/basin_dataset_public_v1p2/shapefiles/merge/Region_03_nhru_simplify_100.shp")%>%select(GAGEID),
  st_read("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-us/basin_dataset_public_v1p2/shapefiles/merge/Region_04_nhru_simplify_100.shp")%>%select(GAGEID),
  st_read("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-us/basin_dataset_public_v1p2/shapefiles/merge/Region_08_nhru_simplify_100.shp")%>%select(GAGEID),
  st_read("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-us/basin_dataset_public_v1p2/shapefiles/merge/Region_09_nhru_simplify_100.shp")%>%select(GAGEID),
  st_read("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-us/basin_dataset_public_v1p2/shapefiles/merge/Region_10U_nhru_.shp")%>%select(GAGEID),
  st_read("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-us/basin_dataset_public_v1p2/shapefiles/merge/Region_10L_nhru_.shp")%>%select(GAGEID),
  st_read("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-us/basin_dataset_public_v1p2/shapefiles/merge/Region_11_nhru_simplify_100.shp")%>%select(GAGEID),
  st_read("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-us/basin_dataset_public_v1p2/shapefiles/merge/Region_12_nhru_simplify_100.shp")%>%select(GAGEID),
  st_read("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-us/basin_dataset_public_v1p2/shapefiles/merge/Region_13_nhru_simplify_100.shp")%>%select(GAGEID),
  st_read("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-us/basin_dataset_public_v1p2/shapefiles/merge/Region_14_nhru_simplify_100.shp")%>%select(GAGEID),
  st_read("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-us/basin_dataset_public_v1p2/shapefiles/merge/Region_15_nhru_simplify_100.shp")%>%select(GAGEID),
  st_read("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-us/basin_dataset_public_v1p2/shapefiles/merge/Region_17_nhru_simplify_100.shp")%>%select(GAGEID)
  
  )%>%
  filter(GAGEID %in% stns$gauge_id)
watersheds<-watersheds%>%
  group_by(GAGEID)%>%
  summarize(geometry = sf::st_union(geometry))

plot(watersheds)
watersheds<-st_simplify(watersheds,preserveTopology = FALSE,dTolerance = 200)
watersheds<-nngeo::st_remove_holes(watersheds)

plot(watersheds)
watersheds<- st_transform(watersheds,st_crs("EPSG:4326"))
st_write(watersheds,"2.data/camels_US.gpkg",append = FALSE)


watersheds$gauge_id<-watersheds$GAGEID

watersheds$area<-st_area(watersheds)%>%as.numeric()

summary(watersheds$area)
watersheds_small<-watersheds%>%filter(area<(1000*10^6))

watersheds_med<-watersheds%>%filter(area>=(1000*10^6)&
                                      area<(5000*10^6))

watersheds_large<-watersheds%>%filter(area>=(5000*10^6)&
                                        area<17000*10^6)

watersheds_huge<-watersheds%>%filter(area>=(17000*10^6))



## small
dat_small<-expand.grid(var = c("tas","pcp"),
                 month = str_pad(1:12,width = 2,pad = "0",side = "left"),
                 ID = watersheds_small$gauge_id,
                 val = NA)%>%
  pivot_wider(names_from =c(var,month),
              values_from = val)
all(watersheds_small$gauge_id==dat_small$ID)


for(it_month in 1:12){
  tic(it_month)
  it_month_st<-str_pad(it_month,width = 2,pad = "0",side = "left")
  T_rast<-terra::rast(sprintf("../DATA/1.Spatial_data/global/clim_climate_precip_aridity_permafrost/worldclim_2.1/wc2.1_30s_tavg_%s.tif",it_month_st))
  dat_small[,paste0("tas_",it_month_st)] <-terra::extract(T_rast,watersheds_small,fun = mean,weights = FALSE,na.rm = TRUE)[,2]
  
  P_rast<-terra::rast(sprintf("../DATA/1.Spatial_data/global/clim_climate_precip_aridity_permafrost/worldclim_2.1/wc2.1_30s_prec_%s.tif",it_month_st))
  dat_small[,paste0("pcp_",it_month_st)] <-terra::extract(P_rast,watersheds_small,fun = mean,weights = FALSE,na.rm = TRUE)[,2]
  toc()
}


# medium
dat_med<-expand.grid(var = c("tas","pcp"),
                     month = str_pad(1:12,width = 2,pad = "0",side = "left"),
                     ID = watersheds_med$gauge_id,
                     val = NA)%>%
  pivot_wider(names_from =c(var,month),
              values_from = val)
all(watersheds_med$gauge_id==dat_med$ID)


for(it_month in 1:12){
  tic(it_month)
  it_month_st<-str_pad(it_month,width = 2,pad = "0",side = "left")
  T_rast<-terra::rast(sprintf("../DATA/1.Spatial_data/global/clim_climate_precip_aridity_permafrost/worldclim_2.1/wc2.1_2.5m_tavg_%s.tif",it_month_st))
  dat_med[,paste0("tas_",it_month_st)] <-terra::extract(T_rast,watersheds_med,fun = mean,weights = FALSE,na.rm = TRUE)[,2]
  
  P_rast<-terra::rast(sprintf("../DATA/1.Spatial_data/global/clim_climate_precip_aridity_permafrost/worldclim_2.1/wc2.1_2.5m_prec_%s.tif",it_month_st))
  dat_med[,paste0("pcp_",it_month_st)] <-terra::extract(P_rast,watersheds_med,fun = mean,weights = FALSE,na.rm = TRUE)[,2]
  toc()
}

# large
dat_large<-expand.grid(var = c("tas","pcp"),
                       month = str_pad(1:12,width = 2,pad = "0",side = "left"),
                       ID = watersheds_large$gauge_id,
                       val = NA)%>%
  pivot_wider(names_from =c(var,month),
              values_from = val)
all(watersheds_large$gauge_id==dat_large$ID)


for(it_month in 1:12){
  tic(it_month)
  it_month_st<-str_pad(it_month,width = 2,pad = "0",side = "left")
  T_rast<-terra::rast(sprintf("../DATA/1.Spatial_data/global/clim_climate_precip_aridity_permafrost/worldclim_2.1/wc2.1_5m_tavg_%s.tif",it_month_st))
  dat_large[,paste0("tas_",it_month_st)] <-terra::extract(T_rast,watersheds_large,fun = mean,weights = FALSE,na.rm = TRUE)[,2]
  
  P_rast<-terra::rast(sprintf("../DATA/1.Spatial_data/global/clim_climate_precip_aridity_permafrost/worldclim_2.1/wc2.1_5m_prec_%s.tif",it_month_st))
  dat_large[,paste0("pcp_",it_month_st)] <-terra::extract(P_rast,watersheds_large,fun = mean,weights = FALSE,na.rm = TRUE)[,2]
  toc()
}

# # huge
# dat_huge<-expand.grid(var = c("tas","pcp"),
#                        month = str_pad(1:12,width = 2,pad = "0",side = "left"),
#                        ID = watersheds_huge$gauge_id,
#                        val = NA)%>%
#   pivot_wider(names_from =c(var,month),
#               values_from = val)
# all(watersheds_huge$gauge_id==dat_huge$ID)
# 
# 
# for(it_month in 1:12){
#   tic(it_month)
#   it_month_st<-str_pad(it_month,width = 2,pad = "0",side = "left")
#   T_rast<-terra::rast(sprintf("../DATA/1.Spatial_data/global/clim_climate_precip_aridity_permafrost/worldclim_2.1/wc2.1_10m_tavg_%s.tif",it_month_st))
#   dat_huge[,paste0("tas_",it_month_st)] <-terra::extract(T_rast,watersheds_huge,fun = mean,weights = FALSE,na.rm = TRUE)[,2]
#   
#   P_rast<-terra::rast(sprintf("../DATA/1.Spatial_data/global/clim_climate_precip_aridity_permafrost/worldclim_2.1/wc2.1_10m_prec_%s.tif",it_month_st))
#   dat_huge[,paste0("pcp_",it_month_st)] <-terra::extract(P_rast,watersheds_huge,fun = mean,weights = FALSE,na.rm = TRUE)[,2]
#   toc()
# }

dat<-rbind(dat_small,
           dat_med,
           dat_large
           # dat_huge
           )

write.csv(dat,"2.data/worldclim/camels_us_clim.csv")
