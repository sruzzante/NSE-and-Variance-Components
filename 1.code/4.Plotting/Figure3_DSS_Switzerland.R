# climatological benchmarks and variance components

# Author: Sacha Ruzzante
# sachawruzzante@gmail.com
# Last Update: 2025-06-18

# Plot the NSE from the differential split samples for Switzerland


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
library(rnaturalearth)
setwd("/home/ruzzante/projects/def-tgleeson/ruzzante/climatological_benchmarks/")

stns<-st_read("2.data/GOF/camels_ch_dss.gpkg")

stns<-filter(stns,(!is.na(NSE_T) & !is.na(NSE_P) & !is.na(NSE_random)))


# ECDF of NSE
ggplot(stns%>%filter(!is.na(NSE_P)))+
  stat_ecdf(aes(x= NSE_P,col = "precipitation"))+
  stat_ecdf(aes(x= NSE_T,col = "temperature"))+
  stat_ecdf(aes(x= NSE_random,col = "random"))+
  scale_color_manual(name = "Sample Split Scheme",
                     values = c("#7570B3","#D95F02","#1B9E77"),
                     labels = c( "Random","Temperature","Precipitation"),
                     breaks =c( "random","temperature","precipitation"))+
  scale_x_continuous(name = "NSE",
                     limits = c(-1,1),
                     oob = scales::squish,
                     expand = c(0,0))+
  scale_y_continuous("Cumulative Density")+
  
  theme_bw(base_size = 7.5)+
  theme(legend.key.spacing.y = unit(-2, "mm"),
        legend.background = element_rect(colour = "black"),
        legend.position= c(0.75,0.22))
ggsave("3.figures/NSE_ecdf_dss_ch.png",width = 3,height = 2.6,dpi = 600)


# maps
coastline<-ne_countries(scale = 10)%>%
  filter(admin %in% c("Switzerland","Austria","France","Italy","Germany","Liechtenstein"))
plot(coastline)
tmap_mode("plot")

# load a hillshade for Switzerland
hs<-terra::rast("2.data/basemaps/hillshade_ch_1.tif")

# Random split NSE
tm1<-tm_shape(coastline%>%filter(admin == "Switzerland"))+tm_borders()+
  tm_shape(hs,raster.downsample=FALSE )+tm_raster(outliers.trunc = TRUE, 
                                                  col.scale = tm_scale_continuous(values = grey.colors(100),
                                                                                  limits = c(0,225)),
                                                  options = opt_tm_raster(interpolate = FALSE),
                                                  col.legend = tm_legend_hide(),
                                                  col_alpha = 0.5
                                                  
  )+
  
  tm_shape(stns%>%filter(!is.na(NSE_random)))+
  tm_dots(col = "NSE_random",breaks = c(-Inf,seq(-1,1,0.2)),
          size = 0.4,
          palette = scico::scico(n = 11,palette = "roma"),
          legend.show = FALSE)+
  tm_shape(coastline)+tm_borders()

tm1


# Temperature split NSE
tm2<-tm_shape(coastline%>%filter(admin == "Switzerland"))+tm_borders()+
  tm_shape(hs,raster.downsample=FALSE )+tm_raster(outliers.trunc = TRUE, 
                                                  col.scale = tm_scale_continuous(values = grey.colors(100),
                                                                                  limits = c(0,225)),
                                                  options = opt_tm_raster(interpolate = FALSE),
                                                  col.legend = tm_legend_hide(),
                                                  col_alpha = 0.5
                                                  
  )+
  
  tm_shape(stns%>%filter(!is.na(NSE_random)))+
  tm_dots(col = "NSE_T",breaks = c(-Inf,seq(-1,1,0.2)),
          size = 0.4,
          palette = scico::scico(n = 11,palette = "roma"),
          legend.show = FALSE)+
  tm_shape(coastline)+tm_borders()

tm2

# Precipitation split NSE
tm3<-tm_shape(coastline%>%filter(admin == "Switzerland"))+tm_borders()+
  tm_shape(hs,raster.downsample=FALSE )+tm_raster(outliers.trunc = TRUE, 
                                                  col.scale = tm_scale_continuous(values = grey.colors(100),
                                                                                  limits = c(0,225)),
                                                  options = opt_tm_raster(interpolate = FALSE),
                                                  col.legend = tm_legend_hide(),
                                                  col_alpha = 0.5
                                                  
  )+
 
  tm_shape(stns%>%filter(!is.na(NSE_random)))+
  tm_dots(col = "NSE_P",breaks = c(-Inf,seq(-1,1,0.2)),
          size = 0.4,
          palette = scico::scico(n = 11,palette = "roma"),
          legend.show = FALSE)+
  tm_shape(coastline)+tm_borders()





tm4<-tmap_arrange(tm1+tm_layout(legend.show=FALSE),
                  tm2+tm_layout(legend.show=FALSE),
                  tm3+tm_layout(legend.show=FALSE),
                  ncol = 3)
tmap_save(tm4,filename = "3.figures/NSE_dss_ch.png",height = 2.6, dpi = 600,width = 8)
