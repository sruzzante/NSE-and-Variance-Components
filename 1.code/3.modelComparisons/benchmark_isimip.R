# benchmark the ISIMIP models
#extract ISIMIP
args = commandArgs(trailingOnly=TRUE)
it_mdl = 1
# numsplit = 1

# test if there is at least two arguments: if not, continue without splitting
if (length(args)<1) {
  print("No arguments supplied, continuing without splitting data")
}else{
  it_mdl = as.integer(args[1])
  
}
# numsplit = as.integer(args[2])





# library(ncdf4)
library(hydroGOF)
library(dplyr)
library(tidyr)
library(ggplot2)
library(hydroGOF)
library(sf)
# library(tmap)
library(tictoc)
library(lubridate)
library(stringr)
setwd("/home/ruzzante/projects/def-tgleeson/ruzzante/climatological_benchmarks/")

source('1.code/utils.R')

fldrs<-dir("2.data/isimip/discharge_simulated/",no.. = TRUE)

watersheds<-st_read("2.data/isimip/catchments/watersheds_combined_filtered.gpkg")

stns<-watersheds%>%st_drop_geometry()%>%
  select(gauge_id,gauge_id2,src,area_km2)



dat_obs<-readRDS("2.data/isimip/discharge/combined.RDS")

fls_model<-list.files(sprintf("2.data/isimip/discharge_simulated/%s/",fldrs[it_mdl]),full.names = TRUE)


dat_model<-lapply(fls_model,readRDS )%>%
  bind_rows()


stnPerf<-vector(mode = "list", length = nrow(stns))

for(it in 1:nrow(stns)){
  tic(sprintf("station %d",it))
  # gauge_id<-stns$gauge_id[it]
  
  dat<-
    data.frame(
      
      QSim = dat_model$q.mm.d[dat_model$gauge_id2 == stns$gauge_id2[it]]*86400,
      dt = dat_model$date[dat_model$gauge_id2 == stns$gauge_id2[it]]
      # ind = ncvar_get(dat_nc,"index")
    )%>%
    left_join(dat_obs%>%
                filter(gauge_id2 == stns$gauge_id2[it])%>%
                select(date,q.mm.d)%>%
                dplyr::rename(QObs = q.mm.d),
              by = c("dt" = "date"))
  if(fldrs[it_mdl] == "jules-w2_gswp3-w5e5_obsclim_histsoc_default_qtot_global_daily"){
      dat<-dat%>%filter(dt<ymd("2019-12-31"))
  }

  
  if(any(is.na(dat$QSim))){
    print(sprintf("THERE ARE NA VALUES IN STATION %s (number %d), SKIPPING",stns$gauge_id2[it],it))
    {next}
  }
  
  
  
  perfMetrics<-calcAllGOF(dat,waterYearStart = 10)
  perfMetrics$gauge_id2<-(stns$gauge_id2[it])
  
  stnPerf[[it]]<-data.frame(perfMetrics)
  
  
  toc()
}




x<-bind_rows(stnPerf)

stns<-left_join(stns,x)

stns$mdl<-fldrs[it_mdl]

write.csv(stns,sprintf("2.data/highLowBenchmarkGOF/GOF_isimip_%s.csv",fldrs[it_mdl]))


stn_result<-compareMetricsHighLow(stns)

write.csv(stn_result,sprintf("2.data/highLowBenchmarkGOF/highLowBenchmarkGOF_isimip_%s.csv",fldrs[it_mdl]))

# }
stns<-stns%>%
  # bind_rows()%>%
  # filter(!is.na(NSE))%>%
  filter(src == "camels-br")



stn_result<-compareMetricsHighLow(stns)
write.csv(stn_result,sprintf("2.data/highLowBenchmarkGOF/highLowBenchmarkGOF_isimip_%s_brazil.csv",fldrs[it_mdl]))
