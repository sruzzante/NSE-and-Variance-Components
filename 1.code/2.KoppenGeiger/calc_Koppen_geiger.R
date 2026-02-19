# caluclate the Koppen-Geiger zone for each catchment, given climatology data from WorldClim

library(dplyr)
library(ggplot2)
library(sf)
# library(tmap)
library(tictoc)
library(stringr)
library(tidyr)
library(lubridate)
library(terra)
setwd("/home/ruzzante/projects/def-tgleeson/ruzzante/climatological_benchmarks/")

fls<-list.files("2.data/worldclim/",
                pattern = "*clim.csv")
for(it in 1:length(fls)){
  
  
  dat<-read.csv(paste0("2.data/worldclim/",fls[it]))%>%
    pivot_longer(cols = tas_01:pcp_12)%>%
    mutate(var = substr(name,1,3),
           month = substr(name,5,6))%>%
    select(!name)%>%
    pivot_wider(names_from = var,values_from = value)
  
  summer<-dat%>%
    group_by(ID)%>%
    summarize(tas_apr_sep = mean(tas[month%in% c("04","05","06","07","08","09")]),
              tas_oct_mar = mean(tas[month%in% c("10","11","12","01","02","03")]))%>%
    mutate(summer= (as.numeric(tas_apr_sep>tas_oct_mar))%>%
             plyr::mapvalues(from = c(0,1), to = c("oct_mar","apr_sep")))%>%
    select(ID,summer)
  
  
  
  dat_summer<-dat%>%left_join(summer)%>%
    filter((month %in% c("04","05","06","07","08","09") & summer == "apr_sep")|
             month %in% c("10","11","12","01","02","03") & summer == "oct_mar")%>%
    select(ID,month,tas,pcp)
  
  
  dat_winter<-dat%>%left_join(summer)%>%
    filter((month %in% c("10","11","12","01","02","03") & summer == "apr_sep")|
             month %in%  c("04","05","06","07","08","09")& summer == "oct_mar")%>%
    select(ID,month,tas,pcp)
  
  dat<-left_join(dat,dat_summer,by = c("ID","month"), suffix = c("",".summer"))%>%
    left_join(dat_winter,by = c("ID","month"), suffix = c("",".winter"))
  
  
  dat_sum<-
    dat%>%
    group_by(ID)%>%
    summarize(MAT = mean(tas),
              Tcold = min(tas),
              Thot = max(tas),
              Tmon10 = sum(tas>10),
              MAP = sum(pcp),
              Pdry = min(pcp),
              Psdry = min(pcp.summer,na.rm = TRUE),
              Pwdry = min(pcp.winter,na.rm = TRUE),
              Pswet = max(pcp.summer,na.rm = TRUE),
              Pwwet = max(pcp.winter,na.rm = TRUE),
              totalSummerP = sum(pcp.summer,na.rm = TRUE),
              totalWinterP = sum(pcp.winter,na.rm = TRUE),
              totalP = sum(pcp))
  
  dat_sum$Pthreshold = NA
  
  dat_sum$Pthreshold[(dat_sum$totalSummerP/dat_sum$totalP)>0.7]<-2*dat_sum$MAT[(dat_sum$totalSummerP/dat_sum$totalP)>0.7]
  
  dat_sum$Pthreshold[(dat_sum$totalWinterP/dat_sum$totalP)>0.7]<-2*dat_sum$MAT[(dat_sum$totalWinterP/dat_sum$totalP)>0.7]
  
  dat_sum$Pthreshold[(dat_sum$totalWinterP/dat_sum$totalP)<=0.7&
                       (dat_sum$totalSummerP/dat_sum$totalP)<=0.7]<-2*dat_sum$MAT[(dat_sum$totalWinterP/dat_sum$totalP)<=0.7&
                                                                                    (dat_sum$totalSummerP/dat_sum$totalP)<=0.7]+14
  
  dat_sum$KG_zone1<-""
  dat_sum$KG_zone2<-""
  dat_sum$KG_zone3<-""
  #A
  
  dat_sum$KG_zone1[dat_sum$Tcold>=18]<-"A"
  mask_A = dat_sum$KG_zone1 == "A"
  dat_sum$KG_zone2[mask_A&
                     dat_sum$Pdry>=(100-dat_sum$MAP/25)]<-"m"
  
  dat_sum$KG_zone2[mask_A&
                     dat_sum$Pdry<(100-dat_sum$MAP/25)]<-"w"
  
  dat_sum$KG_zone2[mask_A&
                     dat_sum$Pdry>=60]<-"f"
  
  #C
  dat_sum$KG_zone1[dat_sum$Thot>10&
                     dat_sum$Tcold>0&dat_sum$Tcold<18]<-"C"
  mask_C = dat_sum$KG_zone1 == "C"
  dat_sum$KG_zone2[mask_C&
                     dat_sum$Psdry<40&dat_sum$Psdry<(dat_sum$Pwwet/3)]<-"s"
  dat_sum$KG_zone2[mask_C&
                     dat_sum$Pwdry<(dat_sum$Pswet/10)]<-"w"
  dat_sum$KG_zone2[mask_C&
                     (dat_sum$KG_zone2 == "")]<-"f"
  dat_sum$KG_zone3[mask_C&
                     dat_sum$Tmon10>=1&dat_sum$Tmon10<4]<-"c"
  dat_sum$KG_zone3[mask_C&
                     dat_sum$Tmon10>=4]<-"b"
  dat_sum$KG_zone3[mask_C&
                     dat_sum$Thot>=22]<-"a"
  
  
  #D
  mask_D<-dat_sum$Thot>10&
    dat_sum$Tcold<=0
  dat_sum$KG_zone1[mask_D]<-"D"
  
  dat_sum$KG_zone2[mask_D & dat_sum$Psdry<40 & dat_sum$Psdry<(dat_sum$Pwwet/3)]<-"s"
  
  dat_sum$KG_zone2[mask_D&
                     dat_sum$Pwdry<(dat_sum$Pswet/10)]<-"w"
  dat_sum$KG_zone2[mask_D&
                     (dat_sum$KG_zone2 == "")]<-"f"
  
  
  dat_sum$KG_zone3[mask_D]<-"c"
  dat_sum$KG_zone3[mask_D&
                     dat_sum$Tcold< -38]<-"d"
  dat_sum$KG_zone3[mask_D&
                     dat_sum$Tmon10>=4]<-"b"
  
  dat_sum$KG_zone3[mask_D&
                     dat_sum$Thot>=22]<-"a"
  
  
  #E
  mask_E<-dat_sum$Thot<=10
  dat_sum$KG_zone1[mask_E]<-"E"
  dat_sum$KG_zone2[mask_E&dat_sum$Thot>0]<-"T"
  dat_sum$KG_zone2[mask_E&dat_sum$Thot<=0]<-"F"
  
  
  #B
  mask_B = dat_sum$MAP<(10*dat_sum$Pthreshold)
  dat_sum$KG_zone1[mask_B]<-"B"
  dat_sum$KG_zone2[mask_B & dat_sum$MAP<(5*dat_sum$Pthreshold)]<-"W"
  dat_sum$KG_zone2[mask_B & dat_sum$MAP>=(5*dat_sum$Pthreshold)]<-"S"
  
  dat_sum$KG_zone3[mask_B & dat_sum$MAT>=18]<-"h"
  dat_sum$KG_zone3[mask_B & dat_sum$MAT<18]<-"k"
  
  summary(dat_sum$KG_zone1%>%factor())
  
  dat
  dat_sum<-
    mutate(
      
      
      dat_sum,KG_zone = paste0(KG_zone1,KG_zone2,KG_zone3))
  summary(dat_sum$KG_zone%>%factor())
  
  write.csv(dat_sum,
            paste0("2.data/worldclim/",
                   str_replace(fls[it],"clim","KG_zone")))
  
}


#sanity check
sf_use_s2(FALSE)

watersheds<-st_read("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-aus-v2/02_location_boundary_area/shp/CAMELS_AUS_v2_Boundaries_adopted.shp")%>%
  st_transform(st_crs("EPSG:4326"))

watersheds<-st_read("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-ch/camels_ch/catchment_delineations/CAMELS_CH_catchments.shp")%>%
  st_zm()%>%
  st_transform(st_crs("EPSG:4326"))%>%
  st_make_valid()

watersheds$gauge_id<-watersheds$gauge_id

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


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

KG_zone_Beck<-terra::extract(terra::rast("../DATA/1.Spatial_data/global/clim_climate_precip_aridity_permafrost/clim1_Koppen_Geiger/Beck_KG_V1_present_0p083.tif"),
                             watersheds,
                             fun = Mode)%>%
  mutate(code = plyr::mapvalues(Beck_KG_V1_present_0p083,from = leg$num,to = leg$code),
         desc = plyr::mapvalues(Beck_KG_V1_present_0p083,from = leg$num,to = leg$desc),
         ID = watersheds$gauge_id)

x<-inner_join(dat_sum%>%select(ID,KG_zone)%>%
                mutate(desc_byMe = plyr::mapvalues(KG_zone,from = leg$code,to = leg$desc)),
              KG_zone_Beck)
x_<-x[x$KG_zone!=x$code,]

KG_rast<-terra::rast("../DATA/1.Spatial_data/global/clim_climate_precip_aridity_permafrost/clim1_Koppen_Geiger/Beck_KG_V1_present_0p0083.tif")
watershed_x<-watersheds%>%filter(gauge_id == "2026")
KG_rast_x = terra::crop(KG_rast,watershed_x%>%st_buffer(0.1))
KG_rast_x<-as.factor(KG_rast_x)
lvls<-unique(KG_rast_x)

levels(KG_rast_x) <- data.frame(
  ID = leg$num[leg$num %in% lvls$Beck_KG_V1_present_0p0083],
  class =leg$code[leg$num %in% lvls$Beck_KG_V1_present_0p0083]
)

tm_shape(KG_rast_x) +
  tm_raster(title = "Köppen-Geiger Classes", legend.show = TRUE,style = "cat") +
  tm_shape(watershed_x) +
  tm_borders()

plot(KG_rast_x)
