# climatological benchmarks and variance components

# Author: Sacha Ruzzante
# sachawruzzante@gmail.com
# Last Update: 2025-06-18


# calculated the benchmark KGE and NSE based on differential split samples, for Camels-CH

# Höge, M., Kauzlaric, M., Siber, R., Schönenberger, U., Horton, P., Schwanbeck, J., Floriancic, M. G., Viviroli, D., Wilhelm, S., Sikorska-Senoner, A. E., Addor, N., Brunner, M., Pool, S., Zappa, M., & Fenicia, F. (2023). CAMELS-CH: Hydro-meteorological time series and landscape attributes for 331 catchments in hydrologic Switzerland. Earth System Science Data, 15(12), 5755–5784. 

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
setwd("/home/ruzzante/projects/def-tgleeson/ruzzante/climatological_benchmarks/")


## CAMELS-CH ########

# load station data
stns<-st_read("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-ch/camels_ch/catchment_delineations/CAMELS_CH_gauging_stations.shp")%>%
  filter(type == "stream")

# initialize empty columns
stns$NSE_P<-NA
stns$KGE_P<-NA

stns$NSE_T<-NA
stns$KGE_T<-NA

stns$NSE_random<-NA
stns$KGE_random<-NA


for(it in 1:nrow(stns)){
  tic(sprintf("station %d",it))
  
  # read in streamflow data
  dat<-read.delim(paste0("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-ch/camels_ch/timeseries/observation_based/CAMELS_CH_obs_based_",
                         stns$gauge_id[it],".csv"),
                  # skip  =7,
                  sep  =",")%>%
    mutate(dt = ymd(date),
           q = discharge_vol.m3.s.,
           yday = pmin(yday(dt),365),
           year = year(dt))
  
  
  # We will need at least 20 years of data; If less, then skip to next iteration
  if(nrow(dat)<(365*20)){next}
  
  dat_dt<-data.frame(dt = seq.Date(from = dat$dt[1],to = dat$dt[nrow(dat)], by = "1 day"))%>%
    mutate(year = year(dt),
           month = month(dt),
           yday = pmin(yday(dt),365))
  
  # water year begins october 1
  dat_dt$wateryear<-dat_dt$year
  dat_dt$wateryear[dat_dt$month %in% (10:12)]<-
    dat_dt$wateryear[dat_dt$month %in% (10:12)]+1
  
  
  dat<-left_join(dat_dt,
                 dat%>%select(dt,q,temperature_mean.degC.,precipitation.mm.d.))
  
  # summarize T and P data by year
  yearlySummary<-
    dat%>%
    group_by(wateryear)%>%
    summarize(nNA = sum(is.na(q)| is.na(precipitation.mm.d.)),
              N = n(),
              meanT = mean(temperature_mean.degC.,na.rm = TRUE),
              sumP = sum(precipitation.mm.d.,na.rm = TRUE))%>%
    filter((N-nNA)>350)%>%
    mutate(period_P = as.numeric(sumP<median(sumP))+1,
           period_T = as.numeric(meanT<median(meanT))+1)
  
  dat<-filter(dat,
              wateryear %in% yearlySummary$wateryear)
  
  
  dat<-left_join(dat,yearlySummary%>%select(wateryear,period_P,period_T))
  
  
  # require at least 20 years of (mostly) complete data
  if(nrow(yearlySummary)<(20)){next}
  
  
  # dat_test2<-data.frame()
  NSE_P<-c()
  KGE_P<-c()
  
  # loop over periods 1 and 2
  for(it_P in 1:2){
    
    dat_test<-dat%>%filter(period_P  == it_P)
    dat_train = dat%>%filter(period_P != it_P)
    
    dat_sum<-dat_train%>%
      # mutate(q[year == it_year] = NA)%>%
      group_by(yday)%>%
      summarize(q = mean(q,na.rm = TRUE))
    
    
    
    # ggplot(dat_sum,aes(yday,q))+geom_line()
    # ggplot(dat_sum,aes(yday,q30))+geom_line()
    dat_test<-left_join(dat_test,dat_sum,by = "yday",suffix = c(".obs",".avg"))
    
    NSE_P[it_P] <- hydroGOF::NSE(dat_test$q.avg,dat_test$q.obs,na.rm = TRUE)
    KGE_P[it_P] <- hydroGOF::KGE(dat_test$q.avg,dat_test$q.obs,na.rm = TRUE)
    
  }
  stns$NSE_P[it]<- median(NSE_P)
  stns$KGE_P[it]<- median(KGE_P)
  
  # dat_test2<-data.frame()
  NSE_T<-c()
  KGE_T<-c()
  
  # loop over periods 1 and 2
  for(it_T in 1:2){
    
    dat_test<-dat%>%filter(period_T  == it_T)
    dat_train = dat%>%filter(period_T != it_T)
    
    dat_sum<-dat_train%>%
      # mutate(q[year == it_year] = NA)%>%
      group_by(yday)%>%
      summarize(q = mean(q,na.rm = TRUE))
    
    
    
    # ggplot(dat_sum,aes(yday,q))+geom_line()
    # ggplot(dat_sum,aes(yday,q30))+geom_line()
    dat_test<-left_join(dat_test,dat_sum,by = "yday",suffix = c(".obs",".avg"))
    
    NSE_T[it_T] <- hydroGOF::NSE(dat_test$q.avg,dat_test$q.obs,na.rm = TRUE)
    KGE_T[it_T] <- hydroGOF::KGE(dat_test$q.avg,dat_test$q.obs,na.rm = TRUE)
    
  }
  stns$NSE_T[it]<- median(NSE_T)
  stns$KGE_T[it]<- median(KGE_T)
  
  
  # 10 random splits
  NSE_random<-data.frame()
  KGE_random<-data.frame()
  set.seed(1)
  
  # loop over periods 1 and 2
  for(it_random in 1:10){
    yearlySummary$period_rand<-sample(yearlySummary$period_P,
                                      size = nrow(yearlySummary),
                                      replace = FALSE)
    dat_x<-left_join(dat,yearlySummary%>%select(wateryear,period_rand))
    
    
    for(it_T in 1:2){
      
      dat_test<-dat_x%>%filter(period_rand  == it_T)
      dat_train = dat_x%>%filter(period_rand != it_T)
      
      dat_sum<-dat_train%>%
        # mutate(q[year == it_year] = NA)%>%
        group_by(yday)%>%
        summarize(q = mean(q,na.rm = TRUE))
      
      
      
      # ggplot(dat_sum,aes(yday,q))+geom_line()
      # ggplot(dat_sum,aes(yday,q30))+geom_line()
      dat_test<-left_join(dat_test,dat_sum,by = "yday",suffix = c(".obs",".avg"))
      
      NSE_random <- rbind(NSE_random,
                          data.frame(
                            it_split = it_random,
                            NSE = hydroGOF::NSE(dat_test$q.avg,dat_test$q.obs,na.rm = TRUE)
                          )
      )    
      KGE_random <- rbind(KGE_random,
                          data.frame(
                            it_split = it_random,
                            KGE = hydroGOF::KGE(dat_test$q.avg,dat_test$q.obs,na.rm = TRUE)
                          )
      )
      
      
    }
    
    
    
  }
  
  stns$NSE_random  [it]<- NSE_random%>%group_by(it_split)%>%summarise(NSE = mean(NSE))%>%pull(NSE)%>%median()
  
  stns$KGE_random  [it]<- KGE_random%>%group_by(it_split)%>%summarise(KGE = mean(KGE))%>%pull(KGE)%>%median()
  
  
  toc()
  
}

st_write(stns,"2.data/GOF/camels_ch_dss.gpkg",append = FALSE)
stns<-st_read("2.data/GOF/camels_ch_dss.gpkg")

stns<-filter(stns,(!is.na(NSE_T) & !is.na(NSE_P) & !is.na(NSE_random)))

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

coastline<-ne_countries(scale = 10)%>%
  filter(admin %in% c("Switzerland","Austria","France","Italy","Germany","Liechtenstein"))
plot(coastline)
tmap_mode("plot")

# dem<-terra::rast("2.data/basemaps/dem_ch.tif")
hs<-terra::rast("2.data/basemaps/hillshade_ch_1.tif")
# dem.hs = dem*hs
# plot(dem.hs)


tm1<-tm_shape(coastline%>%filter(admin == "Switzerland"))+tm_borders()+
  tm_shape(hs,raster.downsample=FALSE )+tm_raster(outliers.trunc = TRUE, 
                                                  col.scale = tm_scale_continuous(values = grey.colors(100),
                                                                                  limits = c(0,225)),
                                                  options = opt_tm_raster(interpolate = FALSE),
                                                  col.legend = tm_legend_hide(),
                                                  col_alpha = 0.5
                                                  
  )+
  # tm_shape(coastline)+tm_polygons(fill = 'grey97')+
  tm_shape(stns%>%filter(!is.na(NSE_random)))+
  tm_dots(col = "NSE_random",breaks = c(-Inf,seq(-1,1,0.2)),
          size = 0.4,
          palette = scico::scico(n = 11,palette = "roma"),
          legend.show = FALSE)+
  tm_shape(coastline)+tm_borders()

tm1

tm2<-tm_shape(coastline%>%filter(admin == "Switzerland"))+tm_borders()+
  tm_shape(hs,raster.downsample=FALSE )+tm_raster(outliers.trunc = TRUE, 
                                                  col.scale = tm_scale_continuous(values = grey.colors(100),
                                                                                  limits = c(0,225)),
                                                  options = opt_tm_raster(interpolate = FALSE),
                                                  col.legend = tm_legend_hide(),
                                                  col_alpha = 0.5
                                                  
  )+
  # tm_shape(coastline)+tm_polygons(fill = 'grey97')+
  tm_shape(stns%>%filter(!is.na(NSE_random)))+
  tm_dots(col = "NSE_T",breaks = c(-Inf,seq(-1,1,0.2)),
          size = 0.4,
          palette = scico::scico(n = 11,palette = "roma"),
          legend.show = FALSE)+
  tm_shape(coastline)+tm_borders()

tm2
tm3<-tm_shape(coastline%>%filter(admin == "Switzerland"))+tm_borders()+
  tm_shape(hs,raster.downsample=FALSE )+tm_raster(outliers.trunc = TRUE, 
                                                  col.scale = tm_scale_continuous(values = grey.colors(100),
                                                                                  limits = c(0,225)),
                                                  options = opt_tm_raster(interpolate = FALSE),
                                                  col.legend = tm_legend_hide(),
                                                  col_alpha = 0.5
                                                  
  )+
  # tm_shape(coastline)+tm_polygons(fill = 'grey97')+
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


