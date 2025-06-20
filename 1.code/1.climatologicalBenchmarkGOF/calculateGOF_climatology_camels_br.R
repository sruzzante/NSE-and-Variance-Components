# climatological benchmarks and variance components

# Author: Sacha Ruzzante
# sachawruzzante@gmail.com
# Last Update: 2025-06-18

# Camels-BR-v2 data:
# Chagas, V. B. P., Chaffe, P. L. B., Addor, N., Fan, F. M., Fleischmann, A. S., Paiva, R. C. D., & Siqueira, V. A. (2020). CAMELS-BR: Hydrometeorological time series and landscape attributes for 897 catchments in Brazil. Earth System Science Data, 12(3), 2075–2096. https://doi.org/10.5194/essd-12-2075-2020

# Updated v1.2 dataset: 
# Chagas, V. B. P., Chaffe, P. L. B., Addor, N., Fan, F. M., Fleischmann, A. S., Paiva, R. C. D., & Siqueira, V. A. (2025). CAMELS-BR: Hydrometeorological time series and landscape attributes for 897 catchments in Brazil - link to files. (Version 1.2) [Dataset]. Zenodo. https://doi.org/10.5281/zenodo.15025488


library(ncdf4)
library(hydroGOF)
library(dplyr)
library(ggplot2)
library(hydroGOF)
library(sf)
library(tmap)
library(tictoc)
library(lubridate)
setwd("/home/ruzzante/projects/def-tgleeson/ruzzante/climatological_benchmarks/")
source("1.code/utils.R")


## CAMELS-BR ########

stns<-st_read("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-br/13_CAMELS_BR_gauge_location/location_gauges_streamflow.gpkg")

# initialize dataframe for GOF statistics
stns$NSE<-NA
stns$KGE<-NA
stns$KGE_r<-NA
stns$KGE_a<-NA
stns$KGE_b<-NA
stns$N<-NA


stns$varFracSeas = NA



# initialize dataframe for variance components
stn_var<-data.frame(gauge_id = stns$gauge_id,
                    varSeas_stl = NA,
                    varInterannual_stl = NA,
                    varRem_stl = NA,
                    varSeas_clas = NA,
                    varInterannual_clas = NA,
                    varRem_clas = NA,
                    varSeas_fourier = NA,
                    varInterannual_fourier = NA,
                    varRem_fourier = NA)

#Loop through stations and calculate statistics
for(it in 1:nrow(stns)){
  tic(sprintf("station %d",it))
  # gauge_id<-stns$gauge_id[it]
  
  # read in data
  dat<-read.delim(sprintf("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-br/02_CAMELS_BR_streamflow_all_catchments/%s_streamflow.txt",stns$gauge_id[it]),
                  sep = "")%>%
    mutate(dt = ymd(paste(year,month,day)),
           yday = pmin(yday(dt),365))%>%
    filter(qual_flag == 1)
  
  names(dat)[4]<-"q"
  
  dat<-filter(dat,!is.na(dat$q))
  # nc
  
  
  # initial check for time series length >=10 years
  if(nrow(dat)<(365*10)){next}
  if(sum(!is.na(dat$q))<=(365*10)){next}
  
  # ensure each day of the year has at least 10 years of data
  numYearsFull<-
    dat%>%
    group_by(yday)%>%
    summarise(notNA = sum(!is.na(q)))%>%
    ungroup()%>%
    summarize(allnotNA = min(notNA))
  if(numYearsFull$allnotNA<10){next}
  
  # no cross validation
  dat_sum<-dat%>%
    
    group_by(yday)%>%
    summarize(q = mean(q,na.rm = TRUE))
  dat_test<-left_join(dat,dat_sum,by = "yday",suffix = c(".obs",".avg"))
  
  stns$varFracSeas[it]<- hydroGOF::NSE(dat_test$q.avg,dat_test$q.obs,na.rm = TRUE)
  
  
  if(round(stns$varFracSeas[it],5)!= round(var(dat_test$q.avg)/var(dat_test$q.obs),5)){
    print("variance fraction is not equal to NSE")
    quit()}
  
  yrs<-unique(dat$year)
  dat_test2<-data.frame()
  
  
  # LOOCV routine for GOF
  for(it_yr in yrs){
    
    dat_test<-dat%>%filter(year  %in%  it_yr)
    dat_train = dat%>%filter(!(year  %in%  it_yr))
    
    dat_sum<-dat_train%>%
      group_by(yday)%>%
      summarize(q = mean(q,na.rm = TRUE))
    
    
    
    dat_test<-left_join(dat_test,dat_sum,by = "yday",suffix = c(".obs",".avg"))
    
    dat_test2<-rbind(dat_test2,dat_test)
    
  }
  stns$NSE[it]<- hydroGOF::NSE(dat_test2$q.avg,dat_test2$q.obs,na.rm = TRUE)
  stns$KGE[it]<- hydroGOF::KGE(dat_test2$q.avg,dat_test2$q.obs,na.rm = TRUE)
  stns$KGE_r[it]<- hydroGOF::KGE(dat_test2$q.avg,dat_test2$q.obs,na.rm = TRUE, s  =c(1,0,0))
  stns$KGE_a[it]<- hydroGOF::KGE(dat_test2$q.avg,dat_test2$q.obs,na.rm = TRUE, s  =c(0,1,0))
  stns$KGE_b[it]<- hydroGOF::KGE(dat_test2$q.avg,dat_test2$q.obs,na.rm = TRUE, s  =c(0,0,1))
  
  stns$N[it]<-length(yrs)
  
  
  
  
  # calculate variance components
  
  # STL decomposition
  stl_var_x<-stl_var(dat,s.window = 7,t.window = 365) 
  stn_var$varSeas_stl[it] <-stl_var_x["varSeas"] 
  stn_var$varInterannual_stl[it] <-stl_var_x["varInterannual"] 
  stn_var$varRem_stl[it] <-stl_var_x["varRem"]    
  
  #classical decomposition 
  clas_var_x<-clas_var(dat)    
  stn_var$varSeas_clas[it] <-clas_var_x["varSeas"]  
  stn_var$varInterannual_clas[it] <-clas_var_x["varInterannual"]  
  stn_var$varRem_clas[it] <-clas_var_x["varRem"]
  
  # Fourier decomposition
  fourier_var_x<-fourier_var(dat)    
  stn_var$varSeas_fourier[it] <-fourier_var_x["varSeas"]  
  stn_var$varInterannual_fourier[it] <-fourier_var_x["varInterannual"]  
  stn_var$varRem_fourier[it] <-fourier_var_x["varRem"]
  
  
  toc()
  
}
write.csv(stn_var,"2.data/varComponents/camels_br.csv",row.names = FALSE)



ggplot(stns,aes(x= NSE))+stat_ecdf()+
  scale_x_continuous(limits = c(-1,1))


ggplot(stns,aes(x= KGE))+stat_ecdf()+
  scale_x_continuous(limits = c(-1,1))
median(stns$NSE,na.rm = TRUE)
tmap_mode("view")
tm_shape(stns%>%filter(!is.na(NSE)))+tm_dots(col = "NSE",breaks = c(-Inf,seq(-1,1,0.2)),
                                             palette = "RdBu")


st_write(stns,"2.data/GOF/camels_br.gpkg",append = FALSE)


# relationship between seasonal variance fraction and NSE
plotly::ggplotly(
  stns%>%
    # filter(N==11)%>%
    mutate(N = factor(N))%>%
    ggplot(aes(x = varFracSeas,y = NSE,col = N))+geom_point()+geom_smooth(method = "lm",se = FALSE)
)

stns%>%
  ggplot(aes(x = varFracSeas,y = NSE,col = N))+geom_point()+
  geom_abline()+
  
  scico::scale_color_scico(
    name = "# years",
    palette = "hawaii",
    breaks = c(10, 30, 100),
    limits = c(10, 100),
    oob = scales::oob_squish,
    trans = "log"
  ) +
  scale_y_continuous(name = expression(NSE*"\'"[cb]~" (with cross-validation)"))+
  scale_x_continuous(name = expression(NSE[cb]~"(without cross-validation)"))+
  theme_bw() 
