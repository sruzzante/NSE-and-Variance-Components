# climatological benchmarks and variance components

# Author: Sacha Ruzzante
# sachawruzzante@gmail.com
# Last Update: 2025-06-18

# Plot the NSE from the differential split samples for HYSETS



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


fls<-list.files(path = "2.data/GOF/",
                pattern = "hysets_dss_.*.gpkg",
                full.names = TRUE)
stns<-lapply(fls,st_read)%>%
  bind_rows()

stns<-filter(stns,(!is.na(NSE_T) & !is.na(NSE_P) & !is.na(NSE_random)))

st_write(stns,"2.data/GOF/hysets_dss.gpkg",append = FALSE)


stns<-st_read("2.data/GOF/hysets_dss.gpkg")
ggplot(stns%>%filter(!is.na(NSE_P)))+
  stat_ecdf(aes(x= NSE_P,col = "precipitation"))+
  stat_ecdf(aes(x= NSE_T,col = "temperature"))+
  stat_ecdf(aes(x= NSE_random,col = "random"))+
  scale_color_manual(name = "Sample Split Scheme",
                     values = c("#7570B3","#D95F02","#1B9E77"),
                     labels = c( "Random","Temperature","Precipitation"),
                     breaks =c( "random","temperature","precipitation"))+
  scale_x_continuous(name = "Benchmark NSE",
                     limits = c(-1,1),
                     oob = scales::squish,
                     expand = c(0,0))+
  scale_y_continuous("Cumulative Density")+
  
  theme_bw(base_size = 7.5)+
  theme(legend.key.spacing.y = unit(-2, "mm"),
        legend.background = element_rect(colour = "black"),
        legend.position= c(0.75,0.22))
ggsave("3.figures/NSE_ecdf_dss_hysets.png",width = 3,height = 2.6,dpi = 600)



coastline<-ne_countries(scale = 50)%>%
  filter(admin %in% c("Canada","United States of America","Mexico"))

coastline<-sf::st_crop(coastline,
                       xmin = -180,ymin = 0,ymax =90,xmax = 0
)


plot(coastline)
tmap_mode("plot")

hs<-terra::rast("2.data/basemaps/hillshade_NA.tif")

coastline<-st_transform(coastline,st_crs(hs))
coastline<-sf::st_crop(coastline,
                       xmin = -4202148,ymin= -2734939, xmax= 2976925, ymax= 5252158
)
hs<-terra::mask(hs,coastline)
stns<-st_transform(stns,st_crs(hs))

tm1<-tm_shape(coastline)+tm_borders(lwd = 0.5)+
  tm_shape(hs)+tm_raster(
    col.scale = tm_scale_continuous(values = grey.colors(100),
                                    limits = c(0,255)),
    options = opt_tm_raster(interpolate = FALSE),
    col.legend = tm_legend_hide(),
    col_alpha = 0.5
    
  )+
  
  tm_shape(stns)+
  tm_dots(
    fill = "NSE_random",
    
    size = 0.05,
    
    fill.legend = tm_legend_hide(),
    fill.scale = tm_scale(breaks = c(-Inf,seq(-1,1,0.2)),
                          values = scico::scico(n = 11,palette = "roma")))+
  tm_shape(coastline)+tm_borders(lwd = 0.5)+
  tm_add_legend(type = "symbols", 
                fill =  scico::scico(n=11,palette = "roma"),
                border.lwd = 0,
                size = 0.5,
                shape = 16,
                labels = c("Less than -1.0","-1.0 to -0.8","-0.8 to -0.6","-0.6 to -0.4","-0.4 to -0.2","-0.2 to 0.0",
                           "0.0 to 0.2","0.2 to 0.4","0.4 to 0.6","0.6 to 0.8","0.8 to 1.0"),
                title = "NSE")


tmap_save(tm1,filename = "3.figures/NSE_dss_NA_a.png",height = 2.6, dpi = 600,width = 5)


tm2<-tm_shape(coastline)+tm_borders(lwd = 0.5)+
  tm_shape(hs)+tm_raster(
    col.scale = tm_scale_continuous(values = grey.colors(100),
                                    limits = c(0,255)),
    options = opt_tm_raster(interpolate = FALSE),
    col.legend = tm_legend_hide(),
    col_alpha = 0.5
    
  )+
  
  tm_shape(stns)+
  tm_dots(
    fill = "NSE_T",
    
    size = 0.05,
    
    fill.legend = tm_legend_hide(),
    fill.scale = tm_scale(breaks = c(-Inf,seq(-1,1,0.2)),
                          values = scico::scico(n = 11,palette = "roma")))+
  tm_shape(coastline)+tm_borders(lwd = 0.5)


tm3<-tm_shape(coastline)+tm_borders(lwd = 0.5)+
  tm_shape(hs)+tm_raster(
    col.scale = tm_scale_continuous(values = grey.colors(100),
                                    limits = c(0,255)),
    options = opt_tm_raster(interpolate = FALSE),
    col.legend = tm_legend_hide(),
    col_alpha = 0.5
    
  )+
  
  tm_shape(stns)+
  tm_dots(
    fill = "NSE_P",
    
    size = 0.05,
    
    fill.legend = tm_legend_hide(),
    fill.scale = tm_scale(breaks = c(-Inf,seq(-1,1,0.2)),
                          values = scico::scico(n = 11,palette = "roma")))+
  tm_shape(coastline)+tm_borders(lwd = 0.5)



tm4<-tmap_arrange(tm1+tm_layout(legend.show=FALSE),
                  tm2+tm_layout(legend.show=FALSE),
                  tm3+tm_layout(legend.show=FALSE),
                  ncol = 3)
tmap_save(tm4,filename = "3.figures/NSE_dss_NA.png",height = 2.6, dpi = 600,width = 8)

