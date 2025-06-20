# climatological benchmarks and variance components
# Author: Sacha Ruzzante
# sachawruzzante@gmail.com
# Last Update: 2025-06-18

# This script plots the benchmark NSE values (Figure 2).


setwd(paste0(Sys.getenv("USERPROFILE"), "/OneDrive - University of Victoria/climatological_benchmarks/")) #Set the working directory
library(sf)
library(dplyr)
library(tmap)
library(ggplot2)
library(terra)
library(rnaturalearth)
library(stars)

# read in data produced by codes in 1.code/1.climatologicalBenchmarkGOF/
stns<-rbind(
  st_read("2.data/GOF/camels_AUS.gpkg")%>%select(NSE,KGE)%>%st_transform("EPSG:4326"),
  st_read("2.data/GOF/camels_br.gpkg")%>%select(NSE,KGE)%>%st_transform("EPSG:4326"),
  st_read("2.data/GOF/camels_ch.gpkg")%>%select(NSE,KGE)%>%st_transform("EPSG:4326"),
  st_read("2.data/GOF/camels_cl.gpkg")%>%select(NSE,KGE)%>%st_transform("EPSG:4326"),
  # st_read("2.data/GOF/camels_de.gpkg")%>%select(NSE,KGE)%>%st_transform("EPSG:4326"), # Not using because caravan-DE has more stations
  st_read("2.data/GOF/caravan_de.gpkg")%>%select(NSE,KGE)%>%st_transform("EPSG:4326"),
  st_read("2.data/GOF/camels_dk.gpkg")%>%select(NSE,KGE)%>%st_transform("EPSG:4326"),
  st_read("2.data/GOF/camels_es.gpkg")%>%select(NSE,KGE)%>%st_transform("EPSG:4326")%>%st_make_valid(),
  st_read("2.data/GOF/camels_US.gpkg")%>%select(NSE,KGE)%>%st_transform("EPSG:4326"),
  st_read("2.data/GOF/camels_fr.gpkg")%>%select(NSE,KGE)%>%st_transform("EPSG:4326"),
  st_read("2.data/GOF/camels_gb.gpkg")%>%select(NSE,KGE)%>%st_transform("EPSG:4326"),
  st_read("2.data/GOF/camels_il.gpkg")%>%select(NSE,KGE)%>%st_transform("EPSG:4326"),
  st_read("2.data/GOF/camels_ind.gpkg")%>%select(NSE,KGE)%>%st_transform("EPSG:4326"),
  st_read("2.data/GOF/camels_lux.gpkg")%>%select(NSE,KGE)%>%st_transform("EPSG:4326"),
  st_read("2.data/GOF/GRDC_africa.gpkg")%>%select(NSE,KGE)%>%st_transform("EPSG:4326"),
  st_read("2.data/GOF/GRDC_asia.gpkg")%>%select(NSE,KGE)%>%st_transform("EPSG:4326"),
  st_read("2.data/GOF/GRDC_europe_a_q.gpkg")%>%select(NSE,KGE)%>%st_transform("EPSG:4326"),
  st_read("2.data/GOF/GRDC_europe_r_z.gpkg")%>%select(NSE,KGE)%>%st_transform("EPSG:4326"),
  st_read("2.data/GOF/GRDC_north_america.gpkg")%>%select(NSE,KGE)%>%st_transform("EPSG:4326"),
  st_read("2.data/GOF/GRDC_south_america.gpkg")%>%select(NSE,KGE)%>%st_transform("EPSG:4326"),
  st_read("2.data/GOF/GRDC_south_west_pacific.gpkg")%>%select(NSE,KGE)%>%st_transform("EPSG:4326"),
  st_read("2.data/GOF/lamah_ce.gpkg")%>%select(NSE,KGE)%>%st_transform("EPSG:4326"),
  st_read("2.data/GOF/lamah_ice.gpkg")%>%select(NSE,KGE)%>%st_transform("EPSG:4326"),
  st_read("2.data/GOF/arcticnet.gpkg")%>%select(NSE,KGE)%>%st_transform("EPSG:4326")

)

hysets<-st_read("2.data/GOF/hysets.gpkg")%>%select(NSE,KGE)
hysets<-st_set_crs(hysets,"EPSG:4326")

stns<-rbind(stns,hysets)
stns<-filter(stns,!is.na(NSE))


coastline<-rnaturalearth::ne_countries(scale = 50)
coastline<-st_transform(coastline,"+proj=wintri +lon_0=0 +lat_1=50.467 +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs +type=crs") # winket tripel projection
coastline<-st_union(coastline%>%st_make_valid())

plot(coastline)

stns<-st_transform(stns,"+proj=wintri +lon_0=0 +lat_1=50.467 +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs +type=crs")

median(stns$NSE)
sum(stns$NSE>=0.8)

ggplot(stns,aes(x= NSE))+stat_ecdf()+
  scale_x_continuous(limits = c(-.2,1),
                     expand = c(0,0),
                     breaks = seq(-1,1,0.5),
                     labels = c(-1,-0.5,0,0.5,1)
  )+
  scale_y_continuous(expand = c(0.01,0.01),
                     name = "Cumulative Density")+
  scale_x_continuous(name = "Benchmark NSE")+
  theme_bw()+
  theme(plot.background = element_rect())

ggsave("3.figures/ecdf_all.png",width = 3,height=3, dpi = 600 )

summary(stns$NSE)
quantile(stns$NSE,seq(0,1,0.1))
mean(stns$NSE>0.8)
stns_highNSE<-stns%>%filter(NSE>0.8)

tm_shape(stns_highNSE)+tm_dots() # visualize the high-NSE stations


# downsample the stations to 300X600 cell grid
rast_grid<-terra::rast(stns,nrows = 300,ncols = 600,vals = 1)

stns$cell<-(terra::extract(rast_grid,stns,cells = TRUE))[,"cell"]
length(unique(stns$cell))

# downsample by median NSE and KGE, mean lat/long
stns_ds<-stns%>%
  cbind(st_coordinates(stns))%>%
  st_drop_geometry()%>%
  filter(!is.na(NSE))%>%
  group_by(cell)%>%
  summarise(NSE = median(NSE),
            KGE = median(KGE),
            X = mean(X),
            Y = mean(Y)
  )%>%
  st_as_sf(coords = c("X","Y"),crs = st_crs(stns))

quantile(stns_ds$NSE,seq(0,1,0.1))


stns<-stns%>%filter(!is.na(NSE))

earthOutline<-
  st_bbox(st_as_stars()) |> 
  st_as_sfc() |> 
  st_set_crs(NA) |> 
  st_segmentize(1) |> 
  st_set_crs("OGC:CRS84") |> 
  st_transform(st_crs(stns)) 

# Hydrorivers V10 - big dataset: https://www.hydrosheds.org/products/hydrorivers
hydroRIVERS<-st_read("D:/DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/sw_HYDRORIVERS/HydroRIVERS_v10_shp/HydroRIVERS_v10_shp/HydroRIVERS_v10.shp")

hydroRIVERS_large<- hydroRIVERS%>%filter(DIS_AV_CMS>6) # select only larger rivers (6 cms was a visually good cutoff)
head(hydroRIVERS_large)
hydroRIVERS_large<-arrange(hydroRIVERS_large,DIS_AV_CMS) # plot larger rivers on top
hydroRIVERS_large$lwd.scale<-log10(hydroRIVERS_large$DIS_AV_CMS)-log10(6)+0.1 # plot linewidth scale 
summary(hydroRIVERS_large$lwd.scale)
hydroRIVERS_large$col.scale<-hydroRIVERS_large$lwd.scale/max(hydroRIVERS_large$lwd.scale)
# hydroRIVERS_large$col.scale<-hydroRIVERS_large$lwd.scale/4.7
summary(hydroRIVERS_large$col.scale)



hydroRIVERS_large<-st_transform(hydroRIVERS_large,st_crs(coastline))

#Two of my highlighted stations got downsampled. Add them back in for the plot
stns_sel<-rbind(
  st_read("2.data/GOF/hysets.gpkg")%>%filter(Official_ID == "09AB009")%>%select(NSE,KGE)%>%st_set_crs("EPSG:4326"),
  st_read("2.data/GOF/camels_ch.gpkg")%>%filter(gauge_id == 2161)%>%select(NSE,KGE)%>%st_transform("EPSG:4326")
  # st_read("2.data/GOF/camels_br.gpkg")%>%filter(gauge_id == 13880000)%>%select(NSE,KGE)%>%st_transform("EPSG:4326"),
  # st_read("2.data/GOF/GRDC_africa.gpkg")%>%filter(grdc_no == 1134460)%>%select(NSE,KGE)%>%st_transform("EPSG:4326"),
  # st_read("2.data/GOF/GRDC_asia.gpkg")%>%filter(grdc_no %in% c())%>%select(NSE,KGE)%>%st_transform("EPSG:4326")
)%>%
  st_transform("+proj=wintri +lon_0=0 +lat_1=50.467 +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs +type=crs")




tm1<-
  
  tm_shape(coastline)+
  tm_polygons(fill = "grey92",
              lwd = 0.1)+
  tm_shape(hydroRIVERS_large)+
  tm_lines(col = "col.scale",
           col.scale = tm_scale_continuous(values = c("grey92","grey80")),
           col.legend = tm_legend_hide(),
           lwd.legend = tm_legend_hide(),
           # col_alpha = "alpha.scale",
           # col_alpha.scale = tm_scale_asis(values.scale =1),
           lwd = "lwd.scale",
           lwd.scale = tm_scale_asis(values.scale = 0.3
           )
  )+
  tm_shape(earthOutline)+tm_borders()+
  tm_shape(stns_ds)+tm_symbols(
    fill.scale = tm_scale(breaks = c(-Inf,seq(-1,1,0.2)),
                          values = scico::scico(n=11,palette = "roma"),
                          midpoint=0),
    fill.legend = tm_legend_hide(),
    col = "grey50",
    lwd  = 0.0001,
    # col_alpha = 1,
    shape = 21,
    fill = "NSE",
    size = 0.1)+
  tm_shape(stns_sel)+tm_symbols(
    fill.scale = tm_scale(breaks = c(-Inf,seq(-1,1,0.2)),
                          values = scico::scico(n=11,palette = "roma"),
                          midpoint=0),
    fill.legend = tm_legend_hide(),
    col = "grey50",
    lwd  = 0.0001,
    # col_alpha = 1,
    shape = 21,
    fill = "NSE",
    size = 0.1)+
  tm_shape(coastline)+tm_borders(lwd = 0.1)+
  tm_add_legend(type = "symbols", 
                fill =  scico::scico(n=11,palette = "roma")[6:11],
                lwd = 0.01,
                size = 0.4,
                shape = 21,
                labels = c("Less than -1.0","-1.0 to -0.8","-0.8 to -0.6","-0.6 to -0.4","-0.4 to -0.2","-0.2 to 0.0",
                           "0.0 to 0.2","0.2 to 0.4","0.4 to 0.6","0.6 to 0.8","0.8 to 1.0")[6:11],
                title = "Benchmark NSE")+
  tm_layout(legend.frame = TRUE,
            legend.position = tm_pos_in(pos.h = "left",pos.v = "bottom"),
            legend.bg.color = "white",
            legend.text.size = 0.7,
            legend.title.size = 0.8,
            earth.boundary = TRUE)
# tm1
# tmap_save(tm1,"3.figures/NSE_alt2.pdf",width = 7,dpi = 1200)
tmap_save(tm1,"3.figures/Figure2_Panel_A.png",width = 7,dpi = 1200)


# same figure for KGE

tm1<-
  
  tm_shape(coastline)+
  tm_polygons(fill = "grey92",
              lwd = 0.1)+
  tm_shape(hydroRIVERS_large)+
  tm_lines(col = "col.scale",
           col.scale = tm_scale_continuous(values = c("grey92","grey80")),
           col.legend = tm_legend_hide(),
           lwd.legend = tm_legend_hide(),
           # col_alpha = "alpha.scale",
           # col_alpha.scale = tm_scale_asis(values.scale =1),
           lwd = "lwd.scale",
           lwd.scale = tm_scale_asis(values.scale = 0.3
           )
  )+
  tm_shape(earthOutline)+tm_borders()+
  tm_shape(stns_ds)+tm_symbols(
    fill.scale = tm_scale(breaks = c(-Inf,seq(-1,1,0.2)),
                          values = scico::scico(n=11,palette = "roma"),
                          midpoint=0),
    fill.legend = tm_legend_hide(),
    col = "grey50",
    lwd  = 0.0001,
    # col_alpha = 1,
    shape = 21,
    fill = "KGE",
    size = 0.1)+
  tm_shape(stns_sel)+tm_symbols(
    fill.scale = tm_scale(breaks = c(-Inf,seq(-1,1,0.2)),
                          values = scico::scico(n=11,palette = "roma"),
                          midpoint=0),
    fill.legend = tm_legend_hide(),
    col = "grey50",
    lwd  = 0.0001,
    # col_alpha = 1,
    shape = 21,
    fill = "KGE",
    size = 0.1)+
  tm_shape(coastline)+tm_borders(lwd = 0.1)+
  tm_add_legend(type = "symbols", 
                fill =  scico::scico(n=11,palette = "roma")[4:11],
                lwd = 0.01,
                size = 0.4,
                shape = 21,
                labels = c("Less than -1.0","-1.0 to -0.8","-0.8 to -0.6","-0.6 to -0.4","-0.4 to -0.2","-0.2 to 0.0",
                           "0.0 to 0.2","0.2 to 0.4","0.4 to 0.6","0.6 to 0.8","0.8 to 1.0")[4:11],
                title = "Benchmark KGE")+
  tm_layout(legend.frame = TRUE,
            legend.position = tm_pos_in(pos.h = "left",pos.v = "bottom"),
            legend.bg.color = "white",
            legend.text.size = 0.7,
            legend.title.size = 0.8,
            earth.boundary = TRUE)

tmap_save(tm1,"3.figures/FigureS1_Panel_A.png",width = 7,dpi = 1200)

