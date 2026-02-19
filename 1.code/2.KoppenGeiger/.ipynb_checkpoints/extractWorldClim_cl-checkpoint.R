# get climatology for camels-cl


library(dplyr)
library(ggplot2)
library(sf)
library(tmap)
library(tictoc)
library(stringr)
library(tidyr)
library(lubridate)
setwd("/home/ruzzante/projects/def-tgleeson/ruzzante/climatological_benchmarks/")

sf_use_s2(FALSE)
watersheds<-st_read("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-cl/CAMELS_CL_v202201/camels_cl_boundaries/camels_cl_boundaries.shp")%>%
  st_transform(st_crs("EPSG:4326"))

# watersheds<-  st_set_crs(watersheds,st_crs("EPSG:4326"))


watersheds$gauge_id<-watersheds$gauge_id

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

write.csv(dat,"2.data/worldclim/camels_cl_clim.csv")
