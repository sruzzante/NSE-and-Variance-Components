# climatological benchmarks and variance components
# Author: Sacha Ruzzante
# sachawruzzante@gmail.com
# Last Update: 2025-08-06

# Evaluate GRID-LSTM-RAPID for dataset (1st, save all data)
#Yang, Y., Feng, D., Beck, H. E., Hu, W., Abbas, A., Sengupta, A., Delle Monache, L., Hartman, R., Lin, P., Shen, C., & Pan, M. (2025). Global Daily Discharge Estimation Based on Grid Long Short-Term Memory (LSTM) Model and River Routing. Water Resources Research, 61(6), e2024WR039764. https://doi.org/10.1029/2024WR039764

# Note that you need to find and save all observation data first, and then run part 1 to save the multt-source data to one folder

library(ncdf4)
library(hydroGOF)
library(dplyr)
library(ggplot2)
library(hydroGOF)
library(sf)
# library(tmap)
library(tictoc)
library(lubridate)
library(stringr)



setwd("/home/ruzzante/projects/def-tgleeson/ruzzante/climatological_benchmarks/")
source("1.code/5.utils/utils.R")

stns<-read.csv("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRID-LSTM-RAPID/input/training_gauge_information.csv")

fls<-list.files('../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRID-LSTM-RAPID/observations/')%>%
  str_remove(".rds")
stns<-filter(stns,name %in% fls)  

test_group_all<-data.frame()
for(it in 0:9){
  test_group<-read.csv(sprintf("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRID-LSTM-RAPID/input/10fold_gauge/test_basin_random_group_0%s.csv",
                               it))
  test_group$test_group<-it
  test_group$test_row<-1:nrow(test_group)
  test_group_all<-rbind(test_group_all,test_group)
  
}
stns<-left_join(stns,test_group_all%>%select(name,test_group,test_row))
# initialize list to store performance metrics
stnPerf<-vector(mode = "list", length = nrow(stns))

# initialize list to store seasonality metrics
stnSeas<-vector(mode = "list", length = nrow(stns))

# initialize list to store goodness of fit for variance components
stnCompNSE<- vector(mode = "list", length = nrow(stns))

stn_var<-data.frame(ID = stns$name,
                    varSeas_stl = NA,
                    varInterannual_stl = NA,
                    varRem_stl = NA,
                    varSeas_clas = NA,
                    varInterannual_clas = NA,
                    varRem_clas = NA,
                    varSeas_fourier = NA,
                    varInterannual_fourier = NA,
                    varRem_fourier = NA
)
# dat_sim<-read.csv("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRID-LSTM-RAPID/simulation/Basin-scale_LSTM/kfold_00_test_basin_19800101_20231001_4seedmean_Streamflow.csv",
#                   header = FALSE)
# dat_sim<-t(dat_sim)%>%
#   data.frame()
# 
# dat_sim$dt<-ymd("1979-12-31")+days(1:nrow(dat_sim))
# dat_sim<-pivot_longer(cols = 1:(ncol(dat_sim)-1),
#                       names_to = )

# loop through stations
for(it in 1:nrow(stns)){
  tic(sprintf("station %d",it))
  dat_obs<-readRDS(sprintf("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRID-LSTM-RAPID/observations/%s.rds",stns$name[it]))
  
  
  dat_sim<-read.csv(sprintf(
    "../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/GRID-LSTM-RAPID/simulation/Basin-scale_LSTM/kfold_0%d_test_basin_19800101_20231001_4seedmean_Streamflow.csv",
    stns$test_group[it]
  ),
  header = FALSE,
  skip = stns$test_row[it]-1,
  nrows = 1)%>%
    t()%>%
    data.frame()
  names(dat_sim)<-"q"
  dat_sim$dt<-ymd("1979-12-31")+days(1:15979)
  
  dat<-left_join(dat_obs,dat_sim,by = "dt")%>%
    mutate(
           yday = pmin(yday(dt),365),
           QObs =q.x,
           QSim = q.y)
  
  if(nrow(dat)<(365)){next}
  dat<-filter(dat,dt %in% seq.Date(from = ymd("1980-01-01"),
                                   to = ymd("2023-09-30"),
              by = "1 day"))
  
  
  if(any(is.na(dat$QSim))){
    print(sprintf("THERE ARE NA VALUES IN STATION %s (number %d), SKIPPING",stns$name[it],it))
    {next}
  }
  
  # calculate goodness of fit statistics
  perfMetrics<-calcAllGOF(dat,waterYearStart = 10,minYears = 7)
  perfMetrics$ID<-stns$name[it]
  stnPerf[[it]]<-data.frame(perfMetrics)
  
  # calculate seasonality statistics
  seas<-calcSeasonality(dat)
  stnSeas[[it]] = data.frame(ID = stns$name[it],
                             CoV = seas["CoV"],
                             QCI = seas["QCI"])
  
  # Goodness of fit of variance components
  NSE_comps<- gof_components(dat)
  NSE_comps$ID<-stns$name[it]
  stnCompNSE[[it]] = data.frame(NSE_comps)
  
  # calculate variance components
  dat<-mutate(dat,
              q = QObs,
              year  =year(dt),
              yday = pmin(yday(dt),365))
  
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

stns$mdl<-"LSTM-RAPID"

x<-bind_rows(stnPerf)%>%
  left_join(bind_rows(stnSeas))%>%
  left_join(bind_rows(stnCompNSE),by = c("ID"),suffix = c("",".maxRun"))%>%
  left_join(stn_var)


stns<-left_join(stns,x,by =c("name" = "ID"))
stns<-st_drop_geometry(stns)
write.csv(stns,"2.data/highLowBenchmarkGOF/GOF_LSTM-RAPID.csv")


