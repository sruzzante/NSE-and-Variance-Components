# benchmark Brasil

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
library(reticulate)
Sys.setenv(UV_OFFLINE=1)
use_python("/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/python/3.11.5/bin/python3", required = TRUE)

reticulate::py_require("pandas")
reticulate::py_require("xarray")
pd <- import("pandas")

setwd("/home/ruzzante/projects/def-tgleeson/ruzzante/climatological_benchmarks/")
source('1.code/utils.R')


yaml_read<- function(pth){
  config<-read_yaml(paste0(pth,"/config.yml"))
  data.frame(
    # epoch = it_epoch,
    batch_size = config$batch_size,
    # forcings = config$forcings,
    head = config$head,
    hidden_size = config$hidden_size,
    initial_forget_bias = config$initial_forget_bias,
    learning_rate_0 = config$learning_rate$`0`,
    # learning_rate_20 = config$learning_rate$`30`,
    # learning_rate_25 = config$learning_rate$`40`,
    # static_attributes= paste(config$static_attributes,collapse = ", "),
    output_dropout = config$output_dropout,
    loss = config$loss
    # dynamic_inputs = paste(config$dynamic_inputs,collapse = ", ")
    )
  
} 
# dat_obs<-read.csv(sprintf("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-br/preprocessed/12500000.csv"))
# dat_obs%>%
#   summarize(across(streamflow_m3s:sm_layer4_era5land,~sum(is.na(.x))))
if(FALSE){
  
  dat_obs%>%filter(is.na(aet_gleam))%>%mutate(yr = year(date))%>%select(yr)%>%unique()
  dat_obs%>%filter(is.na(sm_surface_gleam))%>%mutate(yr = year(date))%>%select(yr)%>%unique()
  dat_obs%>%filter(is.na(p_chirps))%>%mutate(yr = year(date))%>%select(yr)%>%unique()
  dat_obs%>%filter(is.na(p_cpc))%>%mutate(yr = year(date))%>%select(yr)%>%unique()
  # dat_obs%>%filter(is.na(p_cpc))%>%mutate(yr = year(date))%>%select(yr)%>%unique()
  dat_obs%>%filter(is.na(p_mswep))%>%mutate(yr = factor(year(date)))%>%select(yr)%>%summary()
  dat_obs%>%filter(is.na(tmax_cpc))%>%mutate(yr = year(date))%>%select(yr)%>%unique()
  
  lapply(list.files(path = "lstm-camelsbr/runs/test_run_0804_171629/validation/",pattern = "validation_metrics.csv",recursive = TRUE,full.names = TRUE),
         function(fl){read.csv(fl)%>%mutate(epoch = stringr::str_split_fixed(fl,"/",7)[,6])})%>%
    bind_rows()%>%
    group_by(epoch)%>%
    summarize(across(NSE:KGE,~median(.x,na.rm = TRUE)))%>%
    print(n=50)
  
  lapply(list.files(path = "lstm-camelsbr/runs/test_run_0804_173313/validation/",pattern = "validation_metrics.csv",recursive = TRUE,full.names = TRUE),
         function(fl){read.csv(fl)%>%mutate(epoch = stringr::str_split_fixed(fl,"/",7)[,6])})%>%
    bind_rows()%>%
    group_by(epoch)%>%
    summarize(across(NSE:KGE,~median(.x,na.rm = TRUE)))
  
  lapply(list.files(path = "lstm-camelsbr/runs/test_run_0804_192503/validation/",pattern = "validation_metrics.csv",recursive = TRUE,full.names = TRUE),
         function(fl){read.csv(fl)%>%mutate(epoch = stringr::str_split_fixed(fl,"/",7)[,6])})%>%
    bind_rows()%>%
    group_by(epoch)%>%
    summarize(across(NSE:KGE,~median(.x,na.rm = TRUE)))%>%print(n=50)
  
  
  lapply(list.files(path = "lstm-camelsbr/runs/test_run_0804_200655/validation/",pattern = "validation_metrics.csv",recursive = TRUE,full.names = TRUE),
         function(fl){read.csv(fl)%>%mutate(epoch = stringr::str_split_fixed(fl,"/",7)[,6])})%>%
    bind_rows()%>%
    group_by(epoch)%>%
    summarize(across(NSE:KGE,~median(.x,na.rm = TRUE)))%>%print(n=100)
  
  
  lapply(list.files(path = "lstm-camelsbr/runs/test_run_0904_091949/validation/",pattern = "validation_metrics.csv",recursive = TRUE,full.names = TRUE),
         function(fl){read.csv(fl)%>%mutate(epoch = stringr::str_split_fixed(fl,"/",7)[,6])})%>%
    bind_rows()%>%
    group_by(epoch)%>%
    summarize(across(NSE:KGE,~median(.x,na.rm = TRUE)))%>%print(n=100)
  
  lapply(list.files(path = "lstm-camelsbr/runs/test_run_0904_122342/validation/",pattern = "validation_metrics.csv",recursive = TRUE,full.names = TRUE),
         function(fl){read.csv(fl)%>%mutate(epoch = stringr::str_split_fixed(fl,"/",7)[,6])})%>%
    bind_rows()%>%
    group_by(epoch)%>%
    summarize(across(NSE:KGE,~median(.x,na.rm = TRUE)))%>%print(n=100)
  
  lapply(list.files(path = "lstm-camelsbr/runs/test_run_0904_140923/validation/",pattern = "validation_metrics.csv",recursive = TRUE,full.names = TRUE),
         function(fl){read.csv(fl)%>%mutate(epoch = stringr::str_split_fixed(fl,"/",7)[,6])})%>%
    bind_rows()%>%
    group_by(epoch)%>%
    summarize(across(NSE:KGE,~median(.x,na.rm = TRUE)))%>%print(n=100)
  
  lapply(list.files(path = "lstm-camelsbr/runs/test_run_1004_174500/validation/",pattern = "validation_metrics.csv",recursive = TRUE,full.names = TRUE),
         function(fl){read.csv(fl)%>%mutate(epoch = stringr::str_split_fixed(fl,"/",7)[,6])})%>%
    bind_rows()%>%
    group_by(epoch)%>%
    summarize(across(NSE:KGE,~median(.x,na.rm = TRUE)))%>%print(n=100)
  lapply(list.files(path = "lstm-camelsbr/runs/test_run_1104_125627/validation/",pattern = "validation_metrics.csv",recursive = TRUE,full.names = TRUE),
         function(fl){read.csv(fl)%>%mutate(epoch = stringr::str_split_fixed(fl,"/",7)[,6])})%>%
    bind_rows()%>%
    group_by(epoch)%>%
    summarize(across(NSE:KGE,~median(.x,na.rm = TRUE)))%>%print(n=100)
  
  lapply(list.files(path = "lstm-camelsbr/runs/test_run_2605_142113/validation/",pattern = "validation_metrics.csv",recursive = TRUE,full.names = TRUE),
         function(fl){read.csv(fl)%>%mutate(epoch = stringr::str_split_fixed(fl,"/",7)[,6])})%>%
    bind_rows()%>%
    group_by(epoch)%>%
    summarize(across(NSE:KGE,~median(.x,na.rm = TRUE)))%>%print(n=100)
  
  lapply(list.files(path = "lstm-camelsbr/runs/test_run_2605_142818/validation/",pattern = "validation_metrics.csv",recursive = TRUE,full.names = TRUE),
         function(fl){read.csv(fl)%>%mutate(epoch = stringr::str_split_fixed(fl,"/",7)[,6])})%>%
    bind_rows()%>%
    group_by(epoch)%>%
    summarize(across(NSE:KGE,~median(.x,na.rm = TRUE)))%>%print(n=100)
  
  
  lapply(list.files(path = "lstm-camelsbr/runs/test_run_2705_130421/validation/",pattern = "validation_metrics.csv",recursive = TRUE,full.names = TRUE),
         function(fl){read.csv(fl)%>%mutate(epoch = stringr::str_split_fixed(fl,"/",7)[,6])})%>%
    dplyr::bind_rows()%>%
    group_by(epoch)%>%
    dplyr::summarize(across(NSE:KGE,~median(.x,na.rm = TRUE)))%>%print(n=100)
  yaml_read("lstm-camelsbr/runs/test_run_2705_130421")
  
  
  lapply(list.files(path = "lstm-camelsbr/runs/test_run_2705_154358/validation/",pattern = "validation_metrics.csv",recursive = TRUE,full.names = TRUE),
         function(fl){read.csv(fl)%>%mutate(epoch = stringr::str_split_fixed(fl,"/",7)[,6])})%>%
    dplyr::bind_rows()%>%
    group_by(epoch)%>%
    dplyr::summarize(across(NSE:KGE,~median(.x,na.rm = TRUE)))%>%print(n=100)
  yaml_read("lstm-camelsbr/runs/test_run_2705_154358")
  
  
  pth = "lstm-camelsbr/runs/test_run_2705_172137/"
  lapply(list.files(path =paste0(pth,"/validation"), pattern = "validation_metrics.csv",recursive = TRUE,full.names = TRUE),
         function(fl){read.csv(fl)%>%mutate(epoch = stringr::str_split_fixed(fl,"/",7)[,6])})%>%
    dplyr::bind_rows()%>%
    group_by(epoch)%>%
    dplyr::summarize(across(NSE:KGE,~median(.x,na.rm = TRUE)))%>%print(n=100)
  yaml_read(pth)
  
  pth = "lstm-camelsbr/runs/test_run_2705_174150/"
  lapply(list.files(path =paste0(pth,"/validation"), pattern = "validation_metrics.csv",recursive = TRUE,full.names = TRUE),
         function(fl){read.csv(fl)%>%mutate(epoch = stringr::str_split_fixed(fl,"/",7)[,6])})%>%
    dplyr::bind_rows()%>%
    group_by(epoch)%>%
    dplyr::summarize(across(NSE:KGE,~median(.x,na.rm = TRUE)))%>%print(n=100)
  yaml_read(pth)
  
  
  pth = "lstm-camelsbr/runs/test_run_2705_154358/"
  lapply(list.files(path =paste0(pth,"/validation"), pattern = "validation_metrics.csv",recursive = TRUE,full.names = TRUE),
         function(fl){read.csv(fl)%>%mutate(epoch = stringr::str_split_fixed(fl,"/",7)[,6])})%>%
    dplyr::bind_rows()%>%
    group_by(epoch)%>%
    dplyr::summarize(across(NSE:KGE,~median(.x,na.rm = TRUE)))%>%print(n=100)
  yaml_read(pth)
  
  pth = "lstm-camelsbr/runs/test_run_2705_130421/"
  lapply(list.files(path =paste0(pth,"/validation"), pattern = "validation_metrics.csv",recursive = TRUE,full.names = TRUE),
         function(fl){read.csv(fl)%>%mutate(epoch = stringr::str_split_fixed(fl,"/",7)[,6])})%>%
    dplyr::bind_rows()%>%
    group_by(epoch)%>%
    dplyr::summarize(across(NSE:KGE,~median(.x,na.rm = TRUE)))%>%print(n=100)
  yaml_read(pth)
  
  pth = "lstm-camelsbr/runs/test_run_2605_142818/"
  lapply(list.files(path =paste0(pth,"/validation"), pattern = "validation_metrics.csv",recursive = TRUE,full.names = TRUE),
         function(fl){read.csv(fl)%>%mutate(epoch = stringr::str_split_fixed(fl,"/",7)[,6])})%>%
    dplyr::bind_rows()%>%
    group_by(epoch)%>%
    dplyr::summarize(across(NSE:KGE,~median(.x,na.rm = TRUE)))%>%print(n=100)
  yaml_read(pth)
  
  pth = "lstm-camelsbr/runs/test_run_2605_142818/"
  lapply(list.files(path =paste0(pth,"/validation"), pattern = "validation_metrics.csv",recursive = TRUE,full.names = TRUE),
         function(fl){read.csv(fl)%>%mutate(epoch = stringr::str_split_fixed(fl,"/",7)[,6])})%>%
    dplyr::bind_rows()%>%
    group_by(epoch)%>%
    dplyr::summarize(across(NSE:KGE,~median(.x,na.rm = TRUE)))%>%print(n=100)
  yaml_read(pth)
  
  pth = "lstm-camelsbr/runs/test_run_2705_172137/"
  lapply(list.files(path =paste0(pth,"/validation"), pattern = "validation_metrics.csv",recursive = TRUE,full.names = TRUE),
         function(fl){read.csv(fl)%>%mutate(epoch = stringr::str_split_fixed(fl,"/",7)[,6])})%>%
    dplyr::bind_rows()%>%
    group_by(epoch)%>%
    dplyr::summarize(across(NSE:KGE,~median(.x,na.rm = TRUE)))%>%print(n=100)
  yaml_read(pth)
  
  
  pth = "lstm-camelsbr/runs/test_run_2705_174150//"
  lapply(list.files(path =paste0(pth,"/validation"), pattern = "validation_metrics.csv",recursive = TRUE,full.names = TRUE),
         function(fl){read.csv(fl)%>%mutate(epoch = stringr::str_split_fixed(fl,"/",7)[,6])})%>%
    dplyr::bind_rows()%>%
    group_by(epoch)%>%
    dplyr::summarize(across(NSE:KGE,~median(.x,na.rm = TRUE)))%>%print(n=100)
  yaml_read(pth)
  
  
  pth = "lstm-camelsbr/runs/test_run_2705_180944"
  lapply(list.files(path =paste0(pth,"validation"), pattern = "validation_metrics.csv",recursive = TRUE,full.names = TRUE),
         function(fl){read.csv(fl)%>%mutate(epoch = stringr::str_split_fixed(fl,"/",7)[,6])})%>%
    dplyr::bind_rows()%>%
    group_by(epoch)%>%
    dplyr::summarize(across(NSE:KGE,~median(.x,na.rm = TRUE)))%>%print(n=100)
  yaml_read(pth)
  
  
  
  #5-member ensemble
  pth = "lstm-camelsbr/runs/test_run_2605_194845/"
  x<-lapply(list.files(path =paste0(pth,"/validation"),pattern = "validation_metrics.csv",recursive = TRUE,full.names = TRUE),
         function(fl){read.csv(fl)%>%mutate(epoch = stringr::str_split_fixed(fl,"/",7)[,6])})%>%
    dplyr::bind_rows()%>%
    group_by(epoch)%>%
    dplyr::summarize(across(NSE:KGE,~median(.x,na.rm = TRUE)))%>%print(n=100)
  yaml_read(pth)
  which.max(x$NSE)
  
  pth=  "lstm-camelsbr/runs/ens1_2705_130421/"
  x<-lapply(list.files(path =paste0(pth,"/validation"),pattern = "validation_metrics.csv",recursive = TRUE,full.names = TRUE),
         function(fl){read.csv(fl)%>%mutate(epoch = stringr::str_split_fixed(fl,"/",7)[,6])})%>%
    dplyr::bind_rows()%>%
    group_by(epoch)%>%
    dplyr::summarize(across(NSE:KGE,~median(.x,na.rm = TRUE)))%>%print(n=100)
  yaml_read(pth)
  which.max(x$NSE)
  
  # lapply(list.files(path = "lstm-camelsbr/runs/ens2_2705_154103/validation/",pattern = "validation_metrics.csv",recursive = TRUE,full.names = TRUE),
  #        function(fl){read.csv(fl)%>%mutate(epoch = stringr::str_split_fixed(fl,"/",7)[,6])})%>%
  #   dplyr::bind_rows()%>%
  #   group_by(epoch)%>%
  #   dplyr::summarize(across(NSE:KGE,~median(.x,na.rm = TRUE)))%>%print(n=100)
  pth=  "lstm-camelsbr/runs/ens2_2705_164646/"
  x<-lapply(list.files(path =paste0(pth,"/validation"),pattern = "validation_metrics.csv",recursive = TRUE,full.names = TRUE),
            function(fl){read.csv(fl)%>%mutate(epoch = stringr::str_split_fixed(fl,"/",7)[,6])})%>%
    dplyr::bind_rows()%>%
    group_by(epoch)%>%
    dplyr::summarize(across(NSE:KGE,~median(.x,na.rm = TRUE)))%>%print(n=100)
  yaml_read(pth)
  which.max(x$NSE)
  
  pth=  "lstm-camelsbr/runs/ens3_2705_154103/"
  x<-lapply(list.files(path =paste0(pth,"/validation"),pattern = "validation_metrics.csv",recursive = TRUE,full.names = TRUE),
            function(fl){read.csv(fl)%>%mutate(epoch = stringr::str_split_fixed(fl,"/",7)[,6])})%>%
    dplyr::bind_rows()%>%
    group_by(epoch)%>%
    dplyr::summarize(across(NSE:KGE,~median(.x,na.rm = TRUE)))%>%print(n=100)
  yaml_read(pth)
  which.max(x$NSE)
  
  pth=  "lstm-camelsbr/runs/ens4_2705_154103/"
  x<-lapply(list.files(path =paste0(pth,"/validation"),pattern = "validation_metrics.csv",recursive = TRUE,full.names = TRUE),
            function(fl){read.csv(fl)%>%mutate(epoch = stringr::str_split_fixed(fl,"/",7)[,6])})%>%
    dplyr::bind_rows()%>%
    group_by(epoch)%>%
    dplyr::summarize(across(NSE:KGE,~median(.x,na.rm = TRUE)))%>%print(n=100)
  yaml_read(pth)
  which.max(x$NSE)
  
  
}


# average validation performance


pickle_data<-list()
pickle_data[[1]] <- pd$read_pickle("lstm-camelsbr/runs/test_run_2605_194845/validation/model_epoch004/validation_results.p")
pickle_data[[2]] <- pd$read_pickle("lstm-camelsbr/runs/ens1_2705_130421/validation/model_epoch002/validation_results.p")
pickle_data[[3]] <- pd$read_pickle("lstm-camelsbr/runs/ens2_2705_164646/validation/model_epoch006/validation_results.p")
pickle_data[[4]] <- pd$read_pickle("lstm-camelsbr/runs/ens3_2705_154103/validation/model_epoch003/validation_results.p")
pickle_data[[5]] <- pd$read_pickle("lstm-camelsbr/runs/ens4_2705_154103/validation/model_epoch003/validation_results.p")
IDs<-names(pickle_data[[1]])

stns<-data.frame(
  ID =  IDs
)

stnPerf<-vector(mode = "list", length = nrow(stns))

stnSeas<-vector(mode = "list", length = nrow(stns))
stn_var<-data.frame(ID = stns$ID,
                    varSeas_stl = NA,
                    varInterannual_stl = NA,
                    varRem_stl = NA,
                    varSeas_clas = NA,
                    varInterannual_clas = NA,
                    varRem_clas = NA)

for(it in 1:nrow(stns)){
  tic(sprintf("station %d",it))
  # gauge_id<-stns$gauge_id[it]
  dat_obs<-read.csv(sprintf("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-br/preprocessed/%s.csv",stns$ID[it]))
  # select(date,streamflow_mm)
  # dat_sim<-pickle_data
  # gauge_id<-stns$gauge_id[it]
  dat_ls<-list()
  for(it_ens in 1:5){
    xarray_dat<-pickle_data[[it_ens]][[it]]$`1D`$xr
    dat_ls[[it_ens]]<-
      data.frame(
        QObs = xarray_dat$streamflow_mm_obs  $data,
        QSim = xarray_dat$streamflow_mm_sim  $data,
        dt = as.Date(xarray_dat$date$data)
        
      )
    
  }
  dat<-dat_ls%>%
    bind_rows()%>%
    group_by(dt)%>%
    summarize(across(QObs:QSim,~mean(.x)))
  
  # if(any(is.na(dat$QSim))){
  #   print(sprintf("THERE ARE NA VALUES IN STATION %s (number %d), SKIPPING",stns$ID[it],it))
  #   {next}
  # }
  dat<-filter(dat,!is.na(QSim))
  if(nrow(dat)<2){next}
  
  
  perfMetrics<-calcAllGOF(dat,waterYearStart = 10,minYears = 7)
  perfMetrics$ID<-stns$ID[it]
  
  stnPerf[[it]]<-data.frame(perfMetrics)
  
  seas<-calcSeasonality(dat)
  stnSeas[[it]] = data.frame(ID = stns$ID[it],
                             CoV = seas["CoV"],
                             QCI = seas["QCI"])
  
  
  dat<-mutate(dat,
              q = QObs,
              year  =year(dt),
              yday = pmin(yday(dt),365))
  stl_var_x<-stl_var(dat,s.window = 7,t.window = 365) 
  stn_var$varSeas_stl[it] <-stl_var_x["varSeas"] 
  stn_var$varInterannual_stl[it] <-stl_var_x["varInterannual"] 
  stn_var$varRem_stl[it] <-stl_var_x["varRem"]    
  clas_var_x<-clas_var(dat)    
  stn_var$varSeas_clas[it] <-clas_var_x["varSeas"]  
  stn_var$varInterannual_clas[it] <-clas_var_x["varInterannual"]  
  stn_var$varRem_clas[it] <-clas_var_x["varRem"]
  
  
  
  
  toc()
}
# is the benchmark efficiency worse in high-NSE benchmark stations
stns$mdl<-"lstm-camels-br"

x<-bind_rows(stnPerf)%>%
  left_join(bind_rows(stnSeas))%>%
  left_join(stn_var)

stns<-left_join(stns,x)

ggplot(stns,aes(x = NSE))+stat_ecdf(geom = "step")+
  scale_x_continuous(limits = c(-1,1),expand = c(0,0),oob = scales::oob_keep)

stns_ens<-rbind(
  read.csv("lstm-camelsbr/runs/test_run_2605_194845/validation/model_epoch004/validation_metrics.csv")%>%mutate(ens = "1"),
  read.csv("lstm-camelsbr/runs/ens1_2705_130421/validation/model_epoch002/validation_metrics.csv")%>%mutate(ens = "2"),
  read.csv("lstm-camelsbr/runs/ens2_2705_164646/validation/model_epoch006/validation_metrics.csv")%>%mutate(ens = "3"),
  read.csv("lstm-camelsbr/runs/ens3_2705_154103/validation/model_epoch003/validation_metrics.csv")%>%mutate(ens = "4"),
  read.csv("lstm-camelsbr/runs/ens4_2705_154103/validation/model_epoch003/validation_metrics.csv")%>%mutate(ens = "5")
)

ggplot(stns,aes(x = NSE))+stat_ecdf(geom = "step")+
  stat_ecdf(data = stns_ens,aes(x = NSE,col = ens) ,geom = "step")+
  scale_x_continuous(limits = c(-1,1),expand = c(0,0),oob = scales::oob_keep)


stns_valid<-stns


read.csv("lstm-camelsbr/runs/test_run_2605_194845/test/model_epoch004/test_metrics.csv")%>%
  dplyr::summarize(across(NSE:KGE,~median(.x,na.rm = TRUE)))


pickle_data<-list()
pickle_data[[1]] <- pd$read_pickle("lstm-camelsbr/runs/test_run_2605_194845/test/model_epoch004/test_results.p")
pickle_data[[2]] <- pd$read_pickle("lstm-camelsbr/runs/ens1_2705_130421/test/model_epoch002/test_results.p")
pickle_data[[3]] <- pd$read_pickle("lstm-camelsbr/runs/ens2_2705_164646/test/model_epoch006/test_results.p")
pickle_data[[4]] <- pd$read_pickle("lstm-camelsbr/runs/ens3_2705_154103/test/model_epoch003/test_results.p")
pickle_data[[5]] <- pd$read_pickle("lstm-camelsbr/runs/ens4_2705_154103/test/model_epoch003/test_results.p")



stnPerf<-vector(mode = "list", length = nrow(stns))

stnSeas<-vector(mode = "list", length = nrow(stns))
for(it in 1:nrow(stns)){
  tic(sprintf("station %d",it))
  # gauge_id<-stns$gauge_id[it]
  dat_obs<-read.csv(sprintf("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-br/preprocessed/%s.csv",stns$ID[it]))
  # select(date,streamflow_mm)
  # dat_sim<-pickle_data
  # gauge_id<-stns$gauge_id[it]
  dat_ls<-list()
  for(it_ens in 1:5){
    xarray_dat<-pickle_data[[it_ens]][[it]]$`1D`$xr
    dat_ls[[it_ens]]<-
      data.frame(
        QObs = xarray_dat$streamflow_mm_obs  $data,
        QSim = xarray_dat$streamflow_mm_sim  $data,
        dt = as.Date(xarray_dat$date$data)
        
      )
    
  }
  dat<-dat_ls%>%
    bind_rows()%>%
    group_by(dt)%>%
    summarize(across(QObs:QSim,~mean(.x)))
  
  if(any(is.na(dat$QSim))){
    print(sprintf("THERE ARE NA VALUES IN STATION %s (number %d), SKIPPING",stns$ID[it],it))
    {next}
  }
  
  
  perfMetrics<-calcAllGOF(dat,waterYearStart = 10,minYears = 7)
  perfMetrics$ID<-stns$ID[it]
  
  stnPerf[[it]]<-data.frame(perfMetrics)
  
  seas<-calcSeasonality(dat)
  stnSeas[[it]] = data.frame(ID = stns$ID[it],
                             CoV = seas["CoV"],
                             QCI = seas["QCI"])
  
  toc()
}
stns<-data.frame(
  ID =  IDs
)
# is the benchmark efficiency worse in high-NSE benchmark stations
stns$mdl<-"lstm-camels-br"

x<-bind_rows(stnPerf)%>%
  left_join(bind_rows(stnSeas))

stns<-left_join(stns,x)

write.csv(stns,"2.data/highLowBenchmarkGOF/GOF_camels-br.csv")

stns<-read.csv("2.data/highLowBenchmarkGOF/GOF_camels-br.csv")

stns_ens<-rbind(
  read.csv("lstm-camelsbr/runs/test_run_2605_194845/test/model_epoch004/test_metrics.csv")%>%mutate(ens = "1"),
  read.csv("lstm-camelsbr/runs/ens1_2705_130421/test/model_epoch002/test_metrics.csv")%>%mutate(ens = "2"),
  read.csv("lstm-camelsbr/runs/ens2_2705_164646/test/model_epoch006/test_metrics.csv")%>%mutate(ens = "3"),
  read.csv("lstm-camelsbr/runs/ens3_2705_154103/test/model_epoch003/test_metrics.csv")%>%mutate(ens = "4"),
  read.csv("lstm-camelsbr/runs/ens4_2705_154103/test/model_epoch003/test_metrics.csv")%>%mutate(ens = "5")
)

ggplot(stns,aes(x = NSE))+stat_ecdf(geom = "step")+
  stat_ecdf(data = stns_ens,aes(x = NSE,col = ens) ,geom = "step")+
  scale_x_continuous(limits = c(-1,1),expand = c(0,0),oob = scales::oob_keep)

plotly::ggplotly(
  ggplot(stns,aes(x = NSE))+stat_ecdf(geom = "step",aes(col = "test"))+
  stat_ecdf(data = stns_valid,geom = "step",aes(col = "validation"))+
  scale_x_continuous(limits = c(0,1),expand = c(0,0),oob = scales::oob_keep)
)

ggplot(stns,aes(x = NSE))+stat_ecdf(geom = "step",aes(col = "test (1990-2009)"))+
  stat_ecdf(data = stns_valid,geom = "step",aes(col = "validation (1980-1989)"))+
  scale_x_continuous(limits = c(0,1),expand = c(0,0),oob = scales::oob_keep)+
  scale_colour_discrete(name = NULL)+
  scale_y_continuous(name = "CDF")+
  theme_bw()+
  theme(legend.position = c(0.2,0.9),
        legend.box.background = element_rect())
  
ggsave("3.figures/lstm-brazil-NSE.png",width = 4,height = 4,dpi = 600)

stns<-left_join(stns,read.csv("2.data/worldclim/indices/camels_br_indices.csv")%>%select(!X)%>%mutate(ID = as.character(ID)))


stn_result<-compareMetricsHighLow(stns,splitVar = "NSEB",thresh = 0.5)
write.csv(stn_result,"2.data/highLowBenchmarkGOF/NSEB_0.5/highLowBenchmarkGOF_camels-br.csv")
stn_result<-compareMetricsHighLow(stns,splitVar = "NSEB",thresh = 0.4)
write.csv(stn_result,"2.data/highLowBenchmarkGOF/NSEB_0.4/highLowBenchmarkGOF_camels-br.csv")
stn_result<-compareMetricsHighLow(stns,splitVar = "NSEB",thresh = 0.6)
write.csv(stn_result,"2.data/highLowBenchmarkGOF/NSEB_0.6/highLowBenchmarkGOF_camels-br.csv")


stn_result%>%
  pivot_longer(cols = highNSE.FALSE:highNSE.TRUE)%>%
  filter(!metric %in% c("KGEB","NSEB"))%>%
  # mutate(name = recode(name, highNSE.FALSE = ))
  ggplot(aes(y = mdl,x= value,alpha = w.p<0.05))+
  geom_path(arrow = arrow(angle = 30,length = unit(0.1, "inches"),
                          ends = "last", type = "open"))+
  facet_wrap("metric",scales = "free_x")+
  theme_bw()








