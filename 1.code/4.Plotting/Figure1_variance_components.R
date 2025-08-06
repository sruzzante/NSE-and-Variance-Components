# climatological benchmarks and variance components
# Author: Sacha Ruzzante
# sachawruzzante@gmail.com
# Last Update: 2025-08-06

# This script plots the variance components (Figure 1 ).

setwd(paste0(Sys.getenv("USERPROFILE"), "/OneDrive - University of Victoria/climatological_benchmarks/")) #Set the working directory
library(sf)

library(dplyr)
library(tmap)
library(ggplot2)
library(terra)
library(rnaturalearth)


# reach data from all sources
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

stns<-filter(stns, !is.na(varRem_fourier))

coastline<-rnaturalearth::ne_countries(scale = 50)
coastline<-filter(coastline,sovereignt!="Antarctica")
coastline<-st_transform(coastline,"+proj=wintri +lon_0=0 +lat_1=50.467 +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs +type=crs")

coastline<-st_union(coastline%>%st_make_valid())
plot(coastline)

stns<-st_transform(stns,"+proj=wintri +lon_0=0 +lat_1=50.467 +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs +type=crs")


rast_grid<-terra::rast(stns,nrows = 300,ncols = 600,vals = 1)

# rast_grid<-terra::rasterize(stns,rast_grid,field = "NSE",fun = "median")
plot(rast_grid)
stns$cell<-(terra::extract(rast_grid,stns,cells = TRUE))[,"cell"]
length(unique(stns$cell))
stns$varSum<-stns$varInterannual_stl+stns$varRem_stl+stns$varSeas_stl


ggplot(stns,aes(x= varSeas_stl+varInterannual_stl+varRem_stl))+geom_histogram()

ggplot(stns,aes(x= varSeas_fourier+varInterannual_fourier+varRem_fourier))+geom_histogram()

quantile(stns$varSum,probs = c(0.05,0.95))
quantile(stns$varSum)

stns$varInterannual_stl<-stns$varInterannual_stl/stns$varSum
stns$varRem_stl<-stns$varRem_stl/stns$varSum
stns$varSeas_stl<-stns$varSeas_stl/stns$varSum

stns_ds<-stns%>%
  cbind(st_coordinates(stns))%>%
  st_drop_geometry()%>%
  filter(!is.na(NSE))%>%
  group_by(cell)%>%
  summarise(varSeas_stl=mean(varSeas_stl),
            varInterannual_stl = mean(varInterannual_stl),
            varRem_stl = mean(varRem_stl),
            varSeas_fourier=mean(varSeas_fourier),
            varInterannual_fourier = mean(varInterannual_fourier),
            varRem_fourier = mean(varRem_fourier),
            X = mean(X),
            Y = mean(Y)
  )%>%
  st_as_sf(coords = c("X","Y"),crs = st_crs(stns))


library(stars)
earthOutline<-
  st_bbox(st_as_stars()) |> 
  st_as_sfc() |> 
  st_set_crs(NA) |> 
  st_segmentize(1) |> 
  st_set_crs("OGC:CRS84") |> 
  st_transform(st_crs(stns)) 

hydroRIVERS<-st_read("D:/DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/sw_HYDRORIVERS/HydroRIVERS_v10_shp/HydroRIVERS_v10_shp/HydroRIVERS_v10.shp")

hydroRIVERS_large<- hydroRIVERS%>%filter(DIS_AV_CMS>6)
head(hydroRIVERS_large)
# hydroRIVERS_large_ex<-hydroRIVERS_large%>%filter(MAIN_RIV == "60443230")
hydroRIVERS_large<-arrange(hydroRIVERS_large,DIS_AV_CMS)
hydroRIVERS_large$lwd.scale<-log10(hydroRIVERS_large$DIS_AV_CMS)-log10(6)+0.1
summary(hydroRIVERS_large$lwd.scale)
hydroRIVERS_large$col.scale<-hydroRIVERS_large$lwd.scale/max(hydroRIVERS_large$lwd.scale)
# hydroRIVERS_large$col.scale<-hydroRIVERS_large$lwd.scale/4.7
summary(hydroRIVERS_large$col.scale)



hydroRIVERS_large<-st_transform(hydroRIVERS_large,st_crs(coastline))

library(tricolore)
colors_and_legend <- tricolore::Tricolore(stns_ds, 'varRem_stl', 'varSeas_stl', 'varInterannual_stl')
# the color-coded compositions
colors_and_legend$key

head(colors_and_legend$rgb)

stns_ds$rgb<-colors_and_legend$rgb

p_key<-colors_and_legend$key
p_key<-
  p_key+theme(text = element_blank(),
              panel.background = element_rect(fill='transparent'), #transparent panel bg
              plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
              panel.grid.major = element_blank(), #remove major gridlines
              panel.grid.minor = element_blank(), #remove minor gridlines
              legend.background = element_rect(fill='transparent'), #transparent legend bg
              legend.box.background = element_rect(fill='transparent') #transparent legend pane)
  )
ggsave(plot = p_key,"3.figures/var_comp_key.png",dpi = 600)

ggsave(plot = p_key,"3.figures/var_comp_key.svg")

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
    fill = "rgb",
    # fill.scale = tm_scale(breaks = c(-Inf,seq(-1,1,0.2)),
    #                       values = scico::scico(n=11,palette = "roma"),
    #                       midpoint=0),
    # fill.legend = tm_legend_hide(),
    col = "grey50",
    lwd  = 0.0001,
    # col_alpha = 1,
    shape = 21,
    # fill = "NSE",
    size = 0.1)+
  
  tm_shape(coastline)+tm_borders(lwd = 0.1)+
  
  tm_layout(earth.boundary = TRUE)
# tm1
# tmap_save(tm1,"3.figures/NSE_alt2.pdf",width = 7,dpi = 1200)
tmap_save(tm1,"3.figures/VarComponents.png",width = 7,dpi = 1200)



# fourier


library(tricolore)

tricolore::Tricolore(stns_ds, 'varRem_fourier', 'varSeas_fourier', 'varInterannual_fourier')
tricolore::Tricolore(stns, 'varRem_fourier', 'varSeas_fourier', 'varInterannual_fourier')


colors_and_legend <- tricolore::Tricolore(stns_ds, 'varRem_fourier', 'varSeas_fourier', 'varInterannual_fourier')
# the color-coded compositions
colors_and_legend$key

head(colors_and_legend$rgb)

stns_ds$rgb<-colors_and_legend$rgb

p_key<-colors_and_legend$key
p_key<-
  p_key+theme(
    # text = element_blank(),
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend pane)
  )
ggsave(plot = p_key,"3.figures/var_comp_key_fourier.png",
       dpi = 600,width = 7.085,
       height = 6.436667,
       bg = NULL)

ggsave(plot = p_key,"3.figures/var_comp_key_fourier.svg")

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
    fill = "rgb",
    # fill.scale = tm_scale(breaks = c(-Inf,seq(-1,1,0.2)),
    #                       values = scico::scico(n=11,palette = "roma"),
    #                       midpoint=0),
    # fill.legend = tm_legend_hide(),
    col = "grey50",
    lwd  = 0.0001,
    # col_alpha = 1,
    shape = 21,
    # fill = "NSE",
    size = 0.1)+
  
  tm_shape(coastline)+tm_borders(lwd = 0.1)+
  
  tm_layout(earth.boundary = TRUE)
# tm1
# tmap_save(tm1,"3.figures/NSE_alt2.pdf",width = 7,dpi = 1200)
tmap_save(tm1,"3.figures/VarComponents_fourier.png",width = 7,dpi = 1200)
