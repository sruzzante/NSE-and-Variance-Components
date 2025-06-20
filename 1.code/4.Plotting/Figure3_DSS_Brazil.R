# climatological benchmarks and variance components

# Author: Sacha Ruzzante
# sachawruzzante@gmail.com
# Last Update: 2025-06-18

# Plot the NSE from the differential split samples for Brazil


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


stns<-st_read("2.data/GOF/camels_br_dss.gpkg")
#plot ECDF of NSE
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
ggsave("3.figures/NSE_ecdf_dss_br.png",width = 3,height = 2.6,dpi = 600)

# maps
coastline<-ne_countries(scale = 10)%>%
  filter(admin %in% c("Brazil","Argentina","Uruguay","Bolivia","Chile","Falkland Islands","Paraguay","Suriname","Guyana",
                      "Venezuela","Colombia","Peru","Ecuador","France"))
plot(coastline)
tmap_mode("plot")

# load a hillshade for Brazil
hs<-terra::rast("2.data/basemaps/hillshade_br_2.tif")
hs[hs == 0]<-NA

# Random split NSE
tm1<-tm_shape(coastline%>%filter(admin == "Brazil"))+tm_borders()+

  tm_shape(hs,raster.downsample=FALSE )+tm_raster(outliers.trunc = TRUE, 
                                                  col.scale = tm_scale_continuous(values = grey.colors(100),
                                                                                  limits = c(0,474)),
                                                  options = opt_tm_raster(interpolate = FALSE),
                                                  col.legend = tm_legend_hide(),
                                                  col_alpha = 0.7
                                                  # colorNA = "white"
                                                  
  )+
  tm_shape(stns%>%filter(!is.na(NSE_random)))+
  tm_dots(col = "NSE_random",breaks = c(-Inf,seq(-1,1,0.2)),
          size = 0.4,
          palette = scico::scico(n = 11,palette = "roma"),
          legend.show = FALSE)+
  tm_shape(coastline)+tm_borders()+
  tm_add_legend(type = "symbol", 
                col =  scico::scico(n=11,palette = "roma"),
                border.lwd = 0,
                size = 0.5,
                shape = 16,
                labels = c("Less than -1.0","-1.0 to -0.8","-0.8 to -0.6","-0.6 to -0.4","-0.4 to -0.2","-0.2 to 0.0",
                           "0.0 to 0.2","0.2 to 0.4","0.4 to 0.6","0.6 to 0.8","0.8 to 1.0"),
                title = "NSE")
tm1

# Temperature split NSE
tm2<-tm_shape(coastline%>%filter(admin == "Brazil"))+tm_borders()+
  tm_shape(hs,raster.downsample=FALSE )+tm_raster(outliers.trunc = TRUE, 
                                                  col.scale = tm_scale_continuous(values = grey.colors(100),
                                                                                  limits = c(0,474)),
                                                  options = opt_tm_raster(interpolate = FALSE),
                                                  col.legend = tm_legend_hide(),
                                                  col_alpha = 0.7
                                                  # colorNA = "white"
                                                  
  )+
  tm_shape(stns%>%filter(!is.na(NSE_random)))+
  tm_dots(col = "NSE_T",breaks = c(-Inf,seq(-1,1,0.2)),
          size = 0.4,
          
          palette = scico::scico(n = 11,palette = "roma"),
          legend.show = FALSE)+
  tm_shape(coastline)+tm_borders()

# Precipitation split NSE
tm3<-tm_shape(coastline%>%filter(admin == "Brazil"))+tm_borders()+
  tm_shape(hs,raster.downsample=FALSE )+tm_raster(outliers.trunc = TRUE, 
                                                  col.scale = tm_scale_continuous(values = grey.colors(100),
                                                                                  limits = c(0,474)),
                                                  options = opt_tm_raster(interpolate = FALSE),
                                                  col.legend = tm_legend_hide(),
                                                  col_alpha = 0.7
                                                  # colorNA = "white"
                                                  
  )+
  tm_shape(stns%>%filter(!is.na(NSE_P)))+
  tm_dots(col = "NSE_P",breaks = c(-Inf,seq(-1,1,0.2)),
          size = 0.4,
          palette = scico::scico(n = 11,palette = "roma"))+
  tm_shape(coastline)+tm_borders()





tm4<-tmap_arrange(tm1+tm_layout(legend.show=FALSE),
                  tm2+tm_layout(legend.show=FALSE),
                  tm3+tm_layout(legend.show=FALSE),
                  ncol = 3)
tmap_save(tm4,filename = "3.figures/NSE_dss_br.png",height = 2.6, dpi = 600,width = 8)
tm1
tmap_save(tm1,filename = "3.figures/NSE_dss_br_a.png",height = 2.6, dpi = 600,width = 5)