# get climatology for camels-br


library(dplyr)
library(ggplot2)
library(sf)
library(tmap)
library(tictoc)
library(stringr)
library(tidyr)
library(readxl)
library(lubridate)
setwd("/home/ruzzante/projects/def-tgleeson/ruzzante/climatological_benchmarks/")

sf_use_s2(FALSE)

# stns<-st_read("2.data/GOF/camels_br_fromv1.gpkg")
# %>%
#   filter(!is.na(NSE))
if(FALSE){
  
  # gsim_meta1<-read.csv("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/gsim/GSIM_metadata/GSIM_catalog/GSIM_metadata.csv")%>%
  #   filter(country == "BR")%>%
  #   mutate(reference.no = as.numeric(reference.no))%>%
  #   filter(reference.no %in% stns$gauge_id)
  watersheds<-st_read("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-br/12_CAMELS_BR_catchment_boundaries/camels_br_catchments.gpkg")
    
  st_geometry(watersheds)<-"geometry"
  stns<-st_read("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-br/13_CAMELS_BR_gauge_location/location_gauges_streamflow.gpkg")

  # 
  # grdc_watersheds<-st_read("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/shapefiles/2025-02-11 13-54-59.shp")%>%
  #   st_transform(st_crs("EPSG:4326"))
  # grdc_stns<-read_xlsx("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRDC/GRDC_Stations.xlsx")%>%
  #   fi
  grdc_stns<-read.csv("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/caravan_ext_grdc/GRDC_Caravan_extension_csv/attributes/grdc/attributes_additional_grdc.csv")%>%
    filter(country == "BR")%>%
    filter(!nat_id %in% watersheds$gauge_id)
    # 
  watersheds_grdc<-st_read("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/caravan_ext_grdc/GRDC_Caravan_extension_csv/shapefiles/grdc/grdc_basin_shapes.shp")%>%
    
    mutate(gauge_id = plyr::mapvalues(gauge_id ,from = grdc_stns$gauge_id,grdc_stns$nat_id))%>%
    filter(gauge_id %in% grdc_stns$nat_id)
  
  gsim_meta<-read.csv("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/gsim/GSIM_metadata/GSIM_catalog/GSIM_metadata.csv")%>%
    filter(country == "BR")%>%
    mutate(reference.no = as.numeric(reference.no))%>%
    filter(reference.no %in% stns$gauge_id&
             !reference.no %in% watersheds$gauge_id
             # !reference.no %in% watersheds_grdc$gauge_id
           )
  

  
  watersheds_gsim<-lapply(sprintf("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/gsim/GSIM_metadata/GSIM_catchments/%s.shp",tolower(gsim_meta$gsim.no)),
                     # MARGIN = 1,
                     FUN = st_read)%>%
    bind_rows()
  
  
  watersheds_gsim$gauge_id<-plyr::mapvalues(watersheds_gsim$FILENAME,
                                       from = tolower(gsim_meta$gsim.no),to = gsim_meta$reference.no)
  
  watersheds<-rbind(
    watersheds%>%select(gauge_id),
    watersheds_grdc%>%select(gauge_id),
    watersheds_gsim%>%select(gauge_id)
  )
  watersheds<-st_make_valid(watersheds)
  watersheds<-watersheds[!duplicated(watersheds$gauge_id),]
  st_write(watersheds,"../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-br/combined_boundaries.gpkg",append = FALSE)
  
}
# watersheds<-
watersheds<-st_read("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-br/combined_boundaries.gpkg")
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
  dat_small[,paste0("tas_",it_month_st)] <-terra::extract(T_rast,watersheds_small,fun = mean,weights = FALSE)[,2]
  
  P_rast<-terra::rast(sprintf("../DATA/1.Spatial_data/global/clim_climate_precip_aridity_permafrost/worldclim_2.1/wc2.1_30s_prec_%s.tif",it_month_st))
  dat_small[,paste0("pcp_",it_month_st)] <-terra::extract(P_rast,watersheds_small,fun = mean,weights = FALSE)[,2]
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
  dat_med[,paste0("tas_",it_month_st)] <-terra::extract(T_rast,watersheds_med,fun = mean,weights = FALSE)[,2]
  
  P_rast<-terra::rast(sprintf("../DATA/1.Spatial_data/global/clim_climate_precip_aridity_permafrost/worldclim_2.1/wc2.1_2.5m_prec_%s.tif",it_month_st))
  dat_med[,paste0("pcp_",it_month_st)] <-terra::extract(P_rast,watersheds_med,fun = mean,weights = FALSE)[,2]
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
  dat_large[,paste0("tas_",it_month_st)] <-terra::extract(T_rast,watersheds_large,fun = mean,weights = FALSE)[,2]
  
  P_rast<-terra::rast(sprintf("../DATA/1.Spatial_data/global/clim_climate_precip_aridity_permafrost/worldclim_2.1/wc2.1_5m_prec_%s.tif",it_month_st))
  dat_large[,paste0("pcp_",it_month_st)] <-terra::extract(P_rast,watersheds_large,fun = mean,weights = FALSE)[,2]
  toc()
}

# huge
dat_huge<-expand.grid(var = c("tas","pcp"),
                      month = str_pad(1:12,width = 2,pad = "0",side = "left"),
                      ID = watersheds_huge$gauge_id,
                      val = NA)%>%
  pivot_wider(names_from =c(var,month),
              values_from = val)
all(watersheds_huge$gauge_id==dat_huge$ID)


for(it_month in 1:12){
  tic(it_month)
  it_month_st<-str_pad(it_month,width = 2,pad = "0",side = "left")
  T_rast<-terra::rast(sprintf("../DATA/1.Spatial_data/global/clim_climate_precip_aridity_permafrost/worldclim_2.1/wc2.1_10m_tavg_%s.tif",it_month_st))
  dat_huge[,paste0("tas_",it_month_st)] <-terra::extract(T_rast,watersheds_huge,fun = mean,weights = FALSE)[,2]
  
  P_rast<-terra::rast(sprintf("../DATA/1.Spatial_data/global/clim_climate_precip_aridity_permafrost/worldclim_2.1/wc2.1_10m_prec_%s.tif",it_month_st))
  dat_huge[,paste0("pcp_",it_month_st)] <-terra::extract(P_rast,watersheds_huge,fun = mean,weights = FALSE)[,2]
  toc()
}

dat<-rbind(dat_small,
           dat_med,
           dat_large,
           dat_huge)

write.csv(dat,"2.data/worldclim/camels_br_clim.csv")

