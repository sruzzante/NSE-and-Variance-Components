# climatological benchmarks and variance components

# Author: Sacha Ruzzante
# sachawruzzante@gmail.com
# Last Update: 2025-08-06

# These are utility functions and some labelling list for other scripts. 

require(dplyr)
require(hydroGOF)
require(adc)
require(tidyr)
library(ggh4x)
library(scales)

# Calculate benchmark efficiency (Schaefli and Gupta, 2007)
fun_BE<-function(QSim,QObs.obs,QObs.avg){
  mask_NA<-is.na(QSim)|is.na(QObs.obs)|is.na(QObs.avg)
  QSim = QSim[!mask_NA]
  QObs.obs = QObs.obs[!mask_NA]
  QObs.avg = QObs.avg[!mask_NA]
  
  1-(sum((QObs.obs-QSim)^2)/
       sum((QObs.obs - QObs.avg)^2)
  ) # based on Schaefli & Gupta (2007)
}

# Modified Kendall correlation coefficient - if either x or y is constant, set tau = 0
# This is mostly useful when the observed time series is all 0
modKendallCor = function(x,y){
  tau = cor(x,y,method = "kendall",use = "complete.obs")
  if(length(x)<2){return(NA)}
  if(sd(x)==0 | sd(y) ==0){tau = 0}
  return(tau)
}

# Modified Spearman correlation coefficient - if either x or y is constant, set tau = 0
# also accpepts an argument specifying the minimum number of observations 
rSpearmanMin10<-function(sim,obs,minYears = 10){
  if(is.data.frame(sim)){
    sim<-sim[[1]]
  }
  if(is.data.frame(obs)){
    obs<-obs[[1]]
  }
  rs<-hydroGOF::rSpearman(sim,obs)
  if (sum(!is.na(sim) & !is.na(obs))<minYears){
    rs<-NA
  }else if(sd(sim,na.rm = TRUE)==0&!(sd(obs,na.rm = TRUE)==0)){
    rs<-0
  }
  return(rs)
}


NSEMin10<-function(sim,obs,minYears = 10){
  if(is.data.frame(sim)){
    sim<-sim[[1]]
  }
  if(is.data.frame(obs)){
    obs<-obs[[1]]
  }
  nse<-hydroGOF::NSE(sim,obs)
  if (sum(!is.na(sim) & !is.na(obs))<minYears){
    nse<-NA
  }
  return(nse)
}


KGEMin10<-function(sim,obs,minYears = 10){
  if(is.data.frame(sim)){
    sim<-sim[[1]]
  }
  if(is.data.frame(obs)){
    obs<-obs[[1]]
  }
  kge<-hydroGOF::KGE(sim,obs,out.type = "full",
                     method = "2009")
  kge_comps<-c(kge$KGE.value,kge$KGE.elements)
  names(kge_comps)[1]<-"kge"
  
  if (sum(!is.na(sim) & !is.na(obs))<minYears){
    kge_comps<-c(kge = NA,r = NA, Beta = NA,Alpha =NA)
  }else if(sd(sim,na.rm = TRUE)==0&!(sd(obs,na.rm = TRUE)==0)){
    kge_comps[["r"]]<-0
    kge_comps[["kge"]]<-
      1-sqrt(
        (kge_comps[["r"]]-1)^2 +
          (kge_comps[["Beta"]]-1)^2 +
          (kge_comps[["Alpha"]]-1)^2
      )
    
  }
  return(kge_comps)
}


# Calculate all goodness-of-fit statistics 
calcAllGOF<-function(dat,waterYearStart,minYears = NULL){
  
  # require minimum number of years
  if(is.null(minYears)){
    minYears = 10
  }
  
  # initialize list of GOF metrics
  perfMetrics<-list(
    NSE = NA,
    KGE = NA,
    NSEB = NA,
    KGEB = NA,
    BE = NA,
    
    r_halfFlowDay = NA,
    
    r_QCI = NA,
    r_slope_fdc = NA,
    r_BFI = NA,
    
    rPearson = NA,
    Beta = NA,
    Beta_p = NA,
    Gamma = NA,
    Gamma_p = NA,
    KGElf = NA,
    
    #IHA
    r_mag_m1 = NA,
    
    
    r_mag_m2 = NA,
    
    
    r_mag_m3 = NA,
    
    
    r_mag_m4 = NA,
    
    r_mag_m5 = NA,
    
    
    r_mag_m6 = NA,
    r_mag_m7 = NA,
    
    
    r_mag_m8 = NA,
    
    r_mag_m9 = NA,
    
    r_mag_m10 = NA,
    
    r_mag_m11 = NA,
    
    r_mag_m12 = NA,
    
    r_mag_ann = NA,
    
    r_min_1d = NA,
    r_min_3d = NA,
    r_min_7d = NA,
    r_min_30d = NA,
    r_min_90d = NA,
    
    r_max_1d = NA,
    r_max_3d = NA,
    r_max_7d = NA,
    r_max_30d = NA,
    r_max_90d = NA,
    
    r_min_timing = NA,
    
    r_max_timing = NA,
    
    
    r_no_high_pulses = NA,
    r_no_low_pulses = NA,
    
    r_dur_high_pulses = NA,
    r_dur_low_pulses = NA,
    
    r_numDays_low = NA,
    r_numDays_high = NA,
    
    
    r_mean_pos_diff = NA,
    r_mean_neg_diff = NA,
    
    r_no_rises = NA,
    r_no_falls = NA,
    
    r_rld = NA,
    r_fld = NA,
    
    
    # NSE of indicators
    
    NSE_halfFlowDay = NA,
    NSE_QCI = NA,
    NSE_slope_fdc = NA,
    NSE_BFI = NA,
    
    #IHA
    NSE_mag_m1 = NA,
    NSE_mag_m2 = NA,
    NSE_mag_m3 = NA,
    NSE_mag_m4 = NA,
    NSE_mag_m5 = NA,
    NSE_mag_m6 = NA,
    NSE_mag_m7 = NA,
    NSE_mag_m8 = NA,
    NSE_mag_m9 = NA,
    NSE_mag_m10 = NA,
    NSE_mag_m11 = NA,
    NSE_mag_m12 = NA,
    NSE_mag_ann = NA,
    
    NSE_min_1d = NA,
    NSE_min_3d = NA,
    NSE_min_7d = NA,
    NSE_min_30d = NA,
    NSE_min_90d = NA,
    
    NSE_max_1d = NA,
    NSE_max_3d = NA,
    NSE_max_7d = NA,
    NSE_max_30d = NA,
    NSE_max_90d = NA,
    
    NSE_min_timing = NA,
    NSE_max_timing = NA,
    NSE_no_high_pulses = NA,
    NSE_no_low_pulses = NA,
    NSE_dur_high_pulses = NA,
    NSE_dur_low_pulses = NA,
    NSE_numDays_low = NA,
    NSE_numDays_high = NA,
    NSE_mean_pos_diff = NA,
    NSE_mean_neg_diff = NA,
    NSE_no_rises = NA,
    NSE_no_falls = NA,
    NSE_rld = NA,
    NSE_fld = NA,
    
    
    # rPearson of indicators
    
    rp_halfFlowDay = NA,
    rp_QCI = NA,
    rp_slope_fdc = NA,
    rp_BFI = NA,
    
    #IHA
    rp_mag_m1 = NA,
    rp_mag_m2 = NA,
    rp_mag_m3 = NA,
    rp_mag_m4 = NA,
    rp_mag_m5 = NA,
    rp_mag_m6 = NA,
    rp_mag_m7 = NA,
    rp_mag_m8 = NA,
    rp_mag_m9 = NA,
    rp_mag_m10 = NA,
    rp_mag_m11 = NA,
    rp_mag_m12 = NA,
    rp_mag_ann = NA,
    
    rp_min_1d = NA,
    rp_min_3d = NA,
    rp_min_7d = NA,
    rp_min_30d = NA,
    rp_min_90d = NA,
    
    rp_max_1d = NA,
    rp_max_3d = NA,
    rp_max_7d = NA,
    rp_max_30d = NA,
    rp_max_90d = NA,
    
    rp_min_timing = NA,
    rp_max_timing = NA,
    rp_no_high_pulses = NA,
    rp_no_low_pulses = NA,
    rp_dur_high_pulses = NA,
    rp_dur_low_pulses = NA,
    rp_numDays_low = NA,
    rp_numDays_high = NA,
    rp_mean_pos_diff = NA,
    rp_mean_neg_diff = NA,
    rp_no_rises = NA,
    rp_no_falls = NA,
    rp_rld = NA,
    rp_fld = NA,
    
    
    
    # bias of indicators
    
    bias_halfFlowDay = NA,
    bias_QCI = NA,
    bias_slope_fdc = NA,
    bias_BFI = NA,
    
    #IHA
    bias_mag_m1 = NA,
    bias_mag_m2 = NA,
    bias_mag_m3 = NA,
    bias_mag_m4 = NA,
    bias_mag_m5 = NA,
    bias_mag_m6 = NA,
    bias_mag_m7 = NA,
    bias_mag_m8 = NA,
    bias_mag_m9 = NA,
    bias_mag_m10 = NA,
    bias_mag_m11 = NA,
    bias_mag_m12 = NA,
    bias_mag_ann = NA,
    
    bias_min_1d = NA,
    bias_min_3d = NA,
    bias_min_7d = NA,
    bias_min_30d = NA,
    bias_min_90d = NA,
    
    bias_max_1d = NA,
    bias_max_3d = NA,
    bias_max_7d = NA,
    bias_max_30d = NA,
    bias_max_90d = NA,
    
    bias_min_timing = NA,
    bias_max_timing = NA,
    bias_no_high_pulses = NA,
    bias_no_low_pulses = NA,
    bias_dur_high_pulses = NA,
    bias_dur_low_pulses = NA,
    bias_numDays_low = NA,
    bias_numDays_high = NA,
    bias_mean_pos_diff = NA,
    bias_mean_neg_diff = NA,
    bias_no_rises = NA,
    bias_no_falls = NA,
    bias_rld = NA,
    bias_fld = NA,
    
    
    
    # var ratio of indicators
    
    alpha_halfFlowDay = NA,
    alpha_QCI = NA,
    alpha_slope_fdc = NA,
    alpha_BFI = NA,
    
    #IHA
    alpha_mag_m1 = NA,
    alpha_mag_m2 = NA,
    alpha_mag_m3 = NA,
    alpha_mag_m4 = NA,
    alpha_mag_m5 = NA,
    alpha_mag_m6 = NA,
    alpha_mag_m7 = NA,
    alpha_mag_m8 = NA,
    alpha_mag_m9 = NA,
    alpha_mag_m10 = NA,
    alpha_mag_m11 = NA,
    alpha_mag_m12 = NA,
    alpha_mag_ann = NA,
    
    alpha_min_1d = NA,
    alpha_min_3d = NA,
    alpha_min_7d = NA,
    alpha_min_30d = NA,
    alpha_min_90d = NA,
    
    alpha_max_1d = NA,
    alpha_max_3d = NA,
    alpha_max_7d = NA,
    alpha_max_30d = NA,
    alpha_max_90d = NA,
    
    alpha_min_timing = NA,
    alpha_max_timing = NA,
    alpha_no_high_pulses = NA,
    alpha_no_low_pulses = NA,
    alpha_dur_high_pulses = NA,
    alpha_dur_low_pulses = NA,
    alpha_numDays_low = NA,
    alpha_numDays_high = NA,
    alpha_mean_pos_diff = NA,
    alpha_mean_neg_diff = NA,
    alpha_no_rises = NA,
    alpha_no_falls = NA,
    alpha_rld = NA,
    alpha_fld = NA
    
    
  )
  
  
  
  
  # Initial check for less than minYears, return empty list
  if(sum(!is.na(dat$QObs))<=(365*minYears)){return(perfMetrics)}
  
  # Ensure missing days are included in the dataframe, and that they are marked NA
  dat<-left_join(data.frame(dt = seq.Date(from = min(dat$dt),
                                          to = max(dat$dt),
                                          by = "1 day")),
                 dat)%>%
    mutate(yday = pmin(yday(dt),365),
           year = year(dt),
           month = month(dt),
           wateryear = year +floor(month(dt)/waterYearStart) )
  
  
  numYearsFull<-
    dat%>%
    group_by(yday)%>%
    summarise(notNA = sum(!is.na(QObs)))%>%
    ungroup()%>%
    summarize(allnotNA = min(notNA))
  
  # again check for minimum number of years
  if(numYearsFull$allnotNA<minYears){return(perfMetrics)}
  
  
  
  #NSE
  perfMetrics$NSE<- hydroGOF::NSE(dat$QSim,dat$QObs,na.rm = TRUE)
  # KGE and its elements
  kge<-hydroGOF::KGE(dat$QSim,dat$QObs,na.rm = TRUE,out.type = "full",method="2012")
  perfMetrics$KGE<- kge$KGE.value
  perfMetrics$rPearson<-kge$KGE.elements[["r"]]
  perfMetrics$Beta<-kge$KGE.elements[["Beta"]]
  
  # modified beta that goes from -inf to 1
  perfMetrics$Beta_p<-1-sqrt((kge$KGE.elements[["Beta"]] -1)^2 )
  perfMetrics$Gamma<-kge$KGE.elements[["Gamma"]]
  # modified gamma that goes from -inf to 1
  perfMetrics$Gamma_p<-1-sqrt((kge$KGE.elements[["Gamma"]] -1)^2 )
  # perfMetrics$Beta<-hydroGOF::rPearson(dat$QSim,dat$QObs,na.rm = TRUE)
  
  
  # add 0.01* mean (Qobs ) to all data so 1/Q doesn't result in infinite values (Pushpalatha, 2012)
  eps_val<-0.01*mean(dat$QObs,na.rm = TRUE)
  
  perfMetrics$KGElf<- hydroGOF::KGE(1/(dat$QSim+eps_val),
                                    1/(dat$QObs+eps_val),
                                    epsilon.value = 0,
                                    na.rm = TRUE,
                                    method="2012")
  
  
  # benchmark NSEcb and KGEcb
  yrs<-unique(dat$wateryear)
  dat_test2<-data.frame()
  
  for(it_yr in yrs){
    
    dat_test<-dat%>%filter(wateryear  %in%  it_yr)
    dat_train = dat%>%filter(!(wateryear  %in%  it_yr))
    
    dat_sum<-dat_train%>%
      group_by(yday)%>%
      summarize(QObs = mean(QObs,na.rm = TRUE))
    
    dat_test<-left_join(dat_test,dat_sum,by = "yday",suffix = c(".obs",".avg"))
    
    dat_test2<-rbind(dat_test2,dat_test)
    
  }
  perfMetrics$NSEB<- hydroGOF::NSE(dat_test2$QObs.avg,dat_test2$QObs.obs,na.rm = TRUE)
  
  perfMetrics$KGEB<- hydroGOF::KGE(dat_test2$QObs.avg,dat_test2$QObs.obs,na.rm = TRUE)
  
  # benchmark efficiency
  perfMetrics$BE<-fun_BE(
    dat_test2$QSim,
    dat_test2$QObs.obs,
    dat_test2$QObs.avg
  ) # based on Schaefli & Gupta (2007)
  
  # 3, 7, 30, and 90-day rolling means
  dat<-dat%>%
    mutate(
      QObs.mm.d._obs.3= RcppRoll::roll_meanr(QObs,n = 3,fill = NA),
      QObs.mm.d._sim.3= RcppRoll::roll_meanr(QSim,n = 3,fill = NA),
      QObs.mm.d._obs.7= RcppRoll::roll_meanr(QObs,n = 7,fill = NA),
      QObs.mm.d._sim.7= RcppRoll::roll_meanr(QSim,n = 7,fill = NA),
      QObs.mm.d._obs.30= RcppRoll::roll_meanr(QObs,n = 30,fill = NA),
      QObs.mm.d._sim.30= RcppRoll::roll_meanr(QSim,n = 30,fill = NA),
      QObs.mm.d._obs.90= RcppRoll::roll_meanr(QObs,n = 90,fill = NA),
      QObs.mm.d._sim.90= RcppRoll::roll_meanr(QSim,n = 90,fill = NA)
    )
  
  
  dat_yearly<-dat%>%
    
    group_by(wateryear)%>%
    mutate(Q_sum.obs = cumsum(QObs),
           Q_sum.sim = cumsum(QSim))%>%
    
    summarize(
      # half flow day
      halfFlowDay.obs = (yday(dt[which(Q_sum.obs>(max(Q_sum.obs)/2))[1]])+92)%%365,
      halfFlowDay.sim = (yday(dt[which(Q_sum.sim>(max(Q_sum.sim)/2))[1]])+92)%%365,
      
      # max and min 1-day flows
      maxFlow1.Obs = max(QObs),
      maxFlow1.Sim = max(QSim),
      minFlow1.Obs = min(QObs),
      minFlow1.Sim = min(QSim),
      # max and min 3-day flows
      minFlow3.Obs = min(QObs.mm.d._obs.3),
      minFlow3.Sim = min(QObs.mm.d._sim.3),
      maxFlow3.Obs = max(QObs.mm.d._obs.3),
      maxFlow3.Sim = max(QObs.mm.d._sim.3),
      # max and min 7-day flows
      minFlow7.Obs = min(QObs.mm.d._obs.7),
      minFlow7.Sim = min(QObs.mm.d._sim.7),
      maxFlow7.Obs = max(QObs.mm.d._obs.7),
      maxFlow7.Sim = max(QObs.mm.d._sim.7),
      # max and min 30-day flows
      minFlow30.Obs = min(QObs.mm.d._obs.30),
      minFlow30.Sim = min(QObs.mm.d._sim.30),
      maxFlow30.Obs = max(QObs.mm.d._obs.30),
      maxFlow30.Sim = max(QObs.mm.d._sim.30),
      # max and min 90-day flows
      minFlow90.Obs = min(QObs.mm.d._obs.90),
      minFlow90.Sim = min(QObs.mm.d._sim.90),
      maxFlow90.Obs = max(QObs.mm.d._obs.90),
      maxFlow90.Sim = max(QObs.mm.d._sim.90),
      
      # mean annual flow
      QObs.mean = mean(QObs),
      QSim.mean = mean(QSim),
      
      # slope of the flow duration curve
      slope_fdc.Obs = (quantile(log(QObs+eps_val),0.33,na.rm = TRUE)-quantile(log(QObs+eps_val),0.66,na.rm = TRUE))/0.33,
      slope_fdc.Sim = (quantile(log(QSim+eps_val),0.33,na.rm = TRUE)-quantile(log(QSim+eps_val),0.66,na.rm = TRUE))/0.33,
      
      # number of observations
      N = sum(!is.na(QObs)& !is.na(QSim))
    )%>%
    filter(N>=(365+leap_year(wateryear)))
  
  # NSE of min and max
  # stns$NSE_min7[it] = hydroGOF::NSE(dat_yearly$minFlow.Sim,dat_yearly$minFlow.Obs,na.rm = TRUE)
  # stns$NSE_max[it] = hydroGOF::NSE(dat_yearly$maxFlow.Sim,dat_yearly$maxFlow.Obs,na.rm = TRUE)
  
  if(nrow(dat_yearly)>=minYears){
    
    # half flow day
    perfMetrics$r_halfFlowDay = rSpearmanMin10(dat_yearly$halfFlowDay.sim,dat_yearly$halfFlowDay.obs,minYears = minYears)
    perfMetrics$NSE_halfFlowDay = NSEMin10(dat_yearly$halfFlowDay.sim,dat_yearly$halfFlowDay.obs,minYears = minYears)
    kge_comps<-KGEMin10(dat_yearly$halfFlowDay.sim,dat_yearly$halfFlowDay.obs,minYears = minYears)
    perfMetrics$KGE_halfFlowDay= kge_comps[["kge"]]
    perfMetrics$rp_halfFlowDay = kge_comps[["r"]]
    perfMetrics$bias_halfFlowDay = kge_comps[["Beta"]]
    perfMetrics$alpha_halfFlowDay = kge_comps[["Alpha"]]
    
    
    perfMetrics$r_slope_fdc = rSpearmanMin10(dat_yearly$slope_fdc.Sim,dat_yearly$slope_fdc.Obs,minYears = minYears)
    perfMetrics$NSE_slope_fdc = NSEMin10(dat_yearly$slope_fdc.Sim,dat_yearly$slope_fdc.Obs,minYears = minYears)
    kge_comps<-KGEMin10(dat_yearly$slope_fdc.Sim,dat_yearly$slope_fdc.Obs,minYears = minYears)
    perfMetrics$KGE_slope_fdc= kge_comps[["kge"]]
    perfMetrics$rp_slope_fdc = kge_comps[["r"]]
    perfMetrics$bias_slope_fdc = kge_comps[["Beta"]]
    perfMetrics$alpha_slope_fdc = kge_comps[["Alpha"]]
    
    
    
    # streamflow concentration
    QCI = dat%>%
      group_by(wateryear,month)%>%
      summarize(across(c(QObs,QSim),~mean(.x)),
                comp = n() ==days_in_month(dt[1]))%>%
      filter(comp)%>%
      group_by(wateryear)%>%
      summarise(QCI.obs = 100*sum(QObs^2)/sum(QObs)^2,
                QCI.sim = 100*sum(QSim^2)/sum(QSim)^2,
                N = n())%>%
      filter(N ==12& !is.na(QCI.obs))
    
    # spearman correlation and NSE of QCI
    perfMetrics$r_QCI = rSpearmanMin10(QCI$QCI.sim,QCI$QCI.obs,minYears = minYears)
    perfMetrics$NSE_QCI = NSEMin10(QCI$QCI.sim,QCI$QCI.obs,minYears = minYears)
    kge_comps<-KGEMin10(QCI$QCI.sim,QCI$QCI.obs,minYears = minYears)
    perfMetrics$KGE_QCI = kge_comps[["kge"]]
    perfMetrics$rp_QCI = kge_comps[["r"]]
    perfMetrics$bias_QCI = kge_comps[["Beta"]]
    perfMetrics$alpha_QCI = kge_comps[["Alpha"]]
    
    # spearman correlation  and NSEof annual mean flows
    perfMetrics$r_mag_ann = rSpearmanMin10(dat_yearly$QSim.mean,dat_yearly$QObs.mean,minYears = minYears)
    perfMetrics$NSE_mag_ann = NSEMin10(dat_yearly$QSim.mean,dat_yearly$QObs.mean,minYears = minYears)
    kge_comps<-KGEMin10(dat_yearly$QSim.mean,dat_yearly$QObs.mean,minYears = minYears)
    perfMetrics$KGE_mag_ann = kge_comps[["kge"]]
    perfMetrics$rp_mag_ann = kge_comps[["r"]]
    perfMetrics$bias_mag_ann = kge_comps[["Beta"]]
    perfMetrics$alpha_mag_ann = kge_comps[["Alpha"]]
    
    
    # spearman correlation of max and min flows
    for(it_roll in c(1,3,7,30,90)){
      perfMetrics[paste0("r_min_",it_roll,"d")]<-rSpearmanMin10(dat_yearly[,paste0("minFlow",it_roll,".Sim")],
                                                                dat_yearly[,paste0("minFlow",it_roll,".Obs")],
                                                                minYears = minYears)
      
      perfMetrics[paste0("NSE_min_",it_roll,"d")]<-NSEMin10(dat_yearly[,paste0("minFlow",it_roll,".Sim")],
                                                            dat_yearly[,paste0("minFlow",it_roll,".Obs")],
                                                            minYears = minYears)
      kge_comps<-KGEMin10(dat_yearly[,paste0("minFlow",it_roll,".Sim")],
                          dat_yearly[,paste0("minFlow",it_roll,".Obs")],
                          minYears = minYears)
      perfMetrics[paste0("KGE_min_",it_roll,"d")] <- kge_comps[["kge"]]
      perfMetrics[paste0("rp_min_",it_roll,"d")] = kge_comps[["r"]]
      perfMetrics[paste0("bias_min_",it_roll,"d")]= kge_comps[["Beta"]]
      perfMetrics[paste0("alpha_min_",it_roll,"d")] = kge_comps[["Alpha"]]
      
      
      perfMetrics[paste0("r_max_",it_roll,"d")]<-rSpearmanMin10(dat_yearly[,paste0("maxFlow",it_roll,".Sim")],
                                                                dat_yearly[,paste0("maxFlow",it_roll,".Obs")],
                                                                minYears = minYears)
      
      perfMetrics[paste0("NSE_max_",it_roll,"d")]<-NSEMin10(dat_yearly[,paste0("maxFlow",it_roll,".Sim")],
                                                            dat_yearly[,paste0("maxFlow",it_roll,".Obs")],
                                                            minYears = minYears)
      
      kge_comps<-KGEMin10(dat_yearly[,paste0("maxFlow",it_roll,".Sim")],
                          dat_yearly[,paste0("maxFlow",it_roll,".Obs")],
                          minYears = minYears)
      perfMetrics[paste0("KGE_max_",it_roll,"d")] <- kge_comps[["kge"]]
      perfMetrics[paste0("rp_max_",it_roll,"d")] = kge_comps[["r"]]
      perfMetrics[paste0("bias_max_",it_roll,"d")]= kge_comps[["Beta"]]
      perfMetrics[paste0("alpha_max_",it_roll,"d")] = kge_comps[["Alpha"]]
      
    }
    
  }
  
  # monthly flows
  perf_monthly<-dat%>%
    group_by(wateryear,month)%>%
    summarize(QObs = mean(QObs),
              QSim = mean(QSim))%>%
    group_by(month)%>%
    summarize(
      r_mag = rSpearmanMin10(QSim,QObs,
                             minYears = minYears),
      NSE_mag = NSEMin10(QSim,QObs,
                         minYears = minYears),
      KGE_mag = KGEMin10(QSim,QObs,
                         minYears = minYears)[["kge"]],
      rp_mag = KGEMin10(QSim,QObs,
                        minYears = minYears)[["r"]],
      bias_mag = KGEMin10(QSim,QObs,
                          minYears = minYears)[["Beta"]],
      alpha_mag = KGEMin10(QSim,QObs,
                           minYears = minYears)[["Alpha"]],
      N = sum(!is.na(QSim)&!is.na(QObs)))
  
  perf_monthly$r_mag[perf_monthly$N<minYears] = NA
  
  
  
  for(it_month in 1:12){
    perfMetrics[paste0("r_mag_m",it_month)]<-perf_monthly$r_mag[perf_monthly$month==it_month]
    perfMetrics[paste0("NSE_mag_m",it_month)]<-perf_monthly$NSE_mag[perf_monthly$month==it_month]
    perfMetrics[paste0("KGE_mag_m",it_month)]<-perf_monthly$KGE_mag[perf_monthly$month==it_month]
    perfMetrics[paste0("rp_mag_m",it_month)]<-perf_monthly$rp_mag[perf_monthly$month==it_month]
    perfMetrics[paste0("bias_mag_m",it_month)]<-perf_monthly$bias_mag[perf_monthly$month==it_month]
    perfMetrics[paste0("alpha_mag_m",it_month)]<-perf_monthly$alpha_mag[perf_monthly$month==it_month]
  }
  
  # for high and low flow timing, rotate the water year to be diagonally opposite to the median low/high flow timing
  medFlowDat = dat%>%
    group_by(wateryear)%>%
    dplyr::summarize(lfd = yday[which.min(QObs)],
                     hfd = yday[which.max(QObs)],
    )%>%
    ungroup()%>%
    dplyr::summarize(across(lfd:hfd,~median(.x,na.rm = TRUE)))%>%
    dplyr::mutate(
      wateryearStart4lowflow = (lfd+183)%%365+1,
      wateryearStart4highflow = (hfd+183)%%365+1,
    )
  
  # set water years for low and high flows
  dat$wateryear4lowFlow<-dat$year+as.numeric(dat$yday>=medFlowDat$wateryearStart4lowflow)
  dat$wateryear4highFlow<-dat$year+as.numeric(dat$yday>=medFlowDat$wateryearStart4highflow)
  
  dat$yday4lowflow<-(dat$yday-medFlowDat$wateryearStart4lowflow)%%365+1
  dat$yday4highflow<-(dat$yday-medFlowDat$wateryearStart4highflow)%%365+1
  
  # low flow day
  lowFlowTiming = 
    dat%>%
    group_by(wateryear4lowFlow)%>%
    summarize(lfd.Obs = yday4lowflow[which.min(QObs)],
              lfd.Sim = yday4lowflow[which.min(QSim)],
              N = sum(!is.na(QObs)&!is.na(QSim)))%>%
    filter(N>=365)
  
  if(nrow(lowFlowTiming)>=minYears){
    perfMetrics$r_min_timing <- rSpearmanMin10(lowFlowTiming$lfd.Sim,lowFlowTiming$lfd.Obs,
                                               minYears = minYears)
    perfMetrics$NSE_min_timing <- NSEMin10(lowFlowTiming$lfd.Sim,lowFlowTiming$lfd.Obs,
                                           minYears = minYears)
    kge_comps<- KGEMin10(lowFlowTiming$lfd.Sim,
                         lowFlowTiming$lfd.Obs,
                         minYears = minYears)
    perfMetrics$KGE_min_timing = kge_comps[["kge"]]
    perfMetrics$rp_min_timing = kge_comps[["r"]]
    perfMetrics$bias_min_timing = kge_comps[["Beta"]]
    perfMetrics$alpha_min_timing = kge_comps[["Alpha"]]
    
  }
  # high flow day
  highFlowTiming = 
    dat%>%
    group_by(wateryear4highFlow)%>%
    summarize(hfd.Obs = yday4highflow[which.max(QObs)],
              hfd.Sim = yday4highflow[which.max(QSim)],
              N = sum(!is.na(QObs)&!is.na(QSim)))%>%
    filter(N>=365)
  if(nrow(lowFlowTiming)>=minYears){
    perfMetrics$r_max_timing <- rSpearmanMin10(highFlowTiming$hfd.Sim,highFlowTiming$hfd.Obs,
                                               minYears = minYears)
    perfMetrics$NSE_max_timing <- NSEMin10(highFlowTiming$hfd.Sim,highFlowTiming$hfd.Obs,
                                           minYears = minYears)
    
    kge_comps<- KGEMin10(highFlowTiming$hfd.Sim,
                         highFlowTiming$hfd.Obs,
                         minYears = minYears)
    
    perfMetrics$KGE_max_timing = kge_comps[["kge"]]
    perfMetrics$rp_max_timing = kge_comps[["r"]]
    perfMetrics$bias_max_timing = kge_comps[["Beta"]]
    perfMetrics$alpha_max_timing = kge_comps[["Alpha"]]
    
  }
  
  
  
  # get observed 25th and 75th percentiles
  qntls<-quantile(
    dat_test2%>%
      filter(wateryear%in%dat_yearly$wateryear)%>% # make sure it's a full year of data
      pull(QObs.obs),
    c(0.25,0.75)
  )
  
  # days under 25th or above 75th percentile (called pulses in IHA framework)
  pulseDat<- dat%>%
    select(dt,QObs,QSim, wateryear)%>%
    filter(wateryear%in%dat_yearly$wateryear )%>%
    group_by(wateryear)%>%
    mutate(highPulse.Obs = QObs>=qntls["75%"],
           highPulse.Sim = QSim>=qntls["75%"],
           LowPulse.Obs = QObs<=qntls["25%"],
           LowPulse.Sim = QSim<=qntls["25%"]
    )%>%
    summarize(
      # number of pulses
      numHighPulse.Obs = sum(rle(highPulse.Obs)$values==TRUE),
      numHighPulse.Sim = sum(rle(highPulse.Sim)$values==TRUE),
      numLowPulse.Obs = sum(rle(LowPulse.Obs)$values==TRUE),
      numLowPulse.Sim = sum(rle(LowPulse.Sim)$values==TRUE),
      
      # average length of pulses
      lengthHighPulse.Obs = sum(highPulse.Obs)/numHighPulse.Obs,
      lengthHighPulse.Sim = sum(highPulse.Sim)/numHighPulse.Sim,
      lengthLowPulse.Obs = sum(LowPulse.Obs)/numLowPulse.Obs,
      lengthLowPulse.Sim = sum(LowPulse.Sim)/numLowPulse.Sim,
      
      #total number of days at high/low pulse
      totalHighPulse.Obs = sum(highPulse.Obs),
      totalHighPulse.Sim = sum(highPulse.Sim),
      totalLowPulse.Obs = sum(LowPulse.Obs),
      totalLowPulse.Sim = sum(LowPulse.Sim)
    )
  
  # fix length for years that do not reach threshold (above is divided by 0)
  pulseDat$lengthHighPulse.Obs[pulseDat$numHighPulse.Obs == 0]<-0
  pulseDat$lengthHighPulse.Sim[pulseDat$numHighPulse.Sim == 0]<-0
  pulseDat$lengthLowPulse.Obs[pulseDat$numLowPulse.Obs == 0]<-0
  pulseDat$lengthLowPulse.Sim[pulseDat$numLowPulse.Sim == 0]<-0
  
  pulseDat$totalHighPulse.Obs[pulseDat$numHighPulse.Obs == 0]<-0
  pulseDat$totalHighPulse.Sim[pulseDat$numHighPulse.Sim == 0]<-0
  pulseDat$totalLowPulse.Obs[pulseDat$numLowPulse.Obs == 0]<-0
  pulseDat$totalLowPulse.Sim[pulseDat$numLowPulse.Sim == 0]<-0
  
  
  
  # number, average length, and total days low pulse
  # Number of pulses is usually a small number, with lots of ties, so use Kendall
  if(nrow(pulseDat)>=minYears){
    
    # no_low_pulses
    perfMetrics$r_no_low_pulses = modKendallCor(pulseDat$numLowPulse.Obs,pulseDat$numLowPulse.Sim)
    perfMetrics$NSE_no_low_pulses = NSEMin10(pulseDat$numLowPulse.Obs,pulseDat$numLowPulse.Sim,minYears = minYears)
    kge_comps<- KGEMin10(pulseDat$numLowPulse.Obs,pulseDat$numLowPulse.Sim,
                         minYears = minYears)
    perfMetrics$KGE_no_low_pulses = kge_comps[["kge"]]
    perfMetrics$rp_no_low_pulses = kge_comps[["r"]]
    perfMetrics$bias_no_low_pulses = kge_comps[["Beta"]]
    perfMetrics$alpha_no_low_pulses = kge_comps[["Alpha"]]
    
    
    # dur_low_pulses
    perfMetrics$r_dur_low_pulses = rSpearmanMin10(pulseDat$lengthLowPulse.Sim,pulseDat$lengthLowPulse.Obs,
                                                  minYears = minYears)
    perfMetrics$NSE_dur_low_pulses = NSEMin10(pulseDat$lengthLowPulse.Sim,pulseDat$lengthLowPulse.Obs,
                                              minYears = minYears)
    kge_comps<- KGEMin10(pulseDat$lengthLowPulse.Sim,pulseDat$lengthLowPulse.Obs,
                         minYears = minYears)
    perfMetrics$KGE_dur_low_pulses = kge_comps[["kge"]]
    perfMetrics$rp_dur_low_pulses = kge_comps[["r"]]
    perfMetrics$bias_dur_low_pulses = kge_comps[["Beta"]]
    perfMetrics$alpha_dur_low_pulses = kge_comps[["Alpha"]]
    
    
    
    # numDays_low
    perfMetrics$r_numDays_low = rSpearmanMin10(pulseDat$totalLowPulse.Sim,pulseDat$totalLowPulse.Obs,
                                               minYears = minYears)
    perfMetrics$NSE_numDays_low = NSEMin10(pulseDat$totalLowPulse.Sim,pulseDat$totalLowPulse.Obs,
                                           minYears = minYears)
    
    kge_comps<- KGEMin10(pulseDat$totalLowPulse.Sim,pulseDat$totalLowPulse.Obs,
                         minYears = minYears)
    perfMetrics$KGE_numDays_low = kge_comps[["kge"]]
    perfMetrics$rp_numDays_low = kge_comps[["r"]]
    perfMetrics$bias_numDays_low = kge_comps[["Beta"]]
    perfMetrics$alpha_numDays_low = kge_comps[["Alpha"]]
    
    
    # no_high_pulses
    perfMetrics$r_no_high_pulses = modKendallCor(pulseDat$numHighPulse.Obs,pulseDat$numHighPulse.Sim)
    perfMetrics$NSE_no_high_pulses = NSEMin10(pulseDat$numHighPulse.Obs,pulseDat$numHighPulse.Sim,
                                              minYears = minYears)
    
    kge_comps<- KGEMin10(pulseDat$numHighPulse.Obs,pulseDat$numHighPulse.Sim,
                         minYears = minYears)
    perfMetrics$KGE_no_high_pulses = kge_comps[["kge"]]
    perfMetrics$rp_no_high_pulses = kge_comps[["r"]]
    perfMetrics$bias_no_high_pulses = kge_comps[["Beta"]]
    perfMetrics$alpha_no_high_pulses = kge_comps[["Alpha"]]
    
    # dur_high_pulses
    perfMetrics$r_dur_high_pulses = rSpearmanMin10(pulseDat$lengthHighPulse.Sim,pulseDat$lengthHighPulse.Obs,
                                                   minYears = minYears)
    perfMetrics$NSE_dur_high_pulses = NSEMin10(pulseDat$lengthHighPulse.Sim,pulseDat$lengthHighPulse.Obs,
                                               minYears = minYears)
    kge_comps<- KGEMin10(pulseDat$lengthHighPulse.Sim,pulseDat$lengthHighPulse.Obs,
                         minYears = minYears)
    perfMetrics$KGE_dur_high_pulses = kge_comps[["kge"]]
    perfMetrics$rp_dur_high_pulses = kge_comps[["r"]]
    perfMetrics$bias_dur_high_pulses = kge_comps[["Beta"]]
    perfMetrics$alpha_dur_high_pulses = kge_comps[["Alpha"]]
    
    # r_numDays_high
    perfMetrics$r_numDays_high = rSpearmanMin10(pulseDat$totalHighPulse.Sim,pulseDat$totalHighPulse.Obs,
                                                minYears = minYears)
    perfMetrics$NSE_numDays_high = NSEMin10(pulseDat$totalHighPulse.Sim,pulseDat$totalHighPulse.Obs,
                                            minYears = minYears)
    kge_comps<- KGEMin10(pulseDat$totalHighPulse.Sim,pulseDat$totalHighPulse.Obs,
                         minYears = minYears)
    perfMetrics$KGE_numDays_high = kge_comps[["kge"]]
    perfMetrics$rp_numDays_high  = kge_comps[["r"]]
    perfMetrics$bias_numDays_high  = kge_comps[["Beta"]]
    perfMetrics$alpha_numDays_high  = kge_comps[["Alpha"]]
    
  }
  
  # number of rises and falls
  
  # get daily difference
  dat$diff.Obs<-c(NA, diff(dat$QObs,lag = 1))
  dat$diff.Sim<-c(NA, diff(dat$QSim,lag = 1))
  
  DiffData<-
    dat%>%
    group_by(wateryear)%>%
    summarize(
      # mean rise and fall
      meanRise.Obs = mean(diff.Obs[diff.Obs>0]),
      meanRise.Sim = mean(diff.Obs[diff.Sim>0]),
      meanFall.Obs = mean(diff.Obs[diff.Obs<0]),
      meanFall.Sim = mean(diff.Obs[diff.Sim<0]),
      
      # number of continuous rises and falls
      numRise.Obs = sum(rle(diff.Obs>=0)$values==TRUE),
      numRise.Sim = sum(rle(diff.Sim>=0)$values==TRUE),
      numFall.Obs = sum(rle(diff.Obs<=0)$values==TRUE),
      numFall.Sim = sum(rle(diff.Sim<=0)$values==TRUE),
      
      # rising/falling limb density
      rld.Obs = numRise.Obs/sum(diff.Obs>=0,na.rm = TRUE),
      rld.Sim = numRise.Sim/sum(diff.Sim>=0,na.rm = TRUE),
      fld.Obs = meanFall.Obs/sum(diff.Obs<=0,na.rm = TRUE),
      fld.Sim = meanFall.Sim/sum(diff.Sim<=0,na.rm = TRUE),
    )%>%
    filter(wateryear%in%dat_yearly$wateryear)
  
  # fix entries that were divided by 0 above
  DiffData$meanRise.Sim[DiffData$numRise.Sim==0]<-0
  DiffData$meanFall.Sim[DiffData$numFall.Sim==0]<-0
  DiffData$fld.Sim[DiffData$numFall.Sim==0]<-0
  DiffData$rld.Sim[DiffData$numRise.Sim==0]<-0
  
  DiffData$meanRise.Obs[DiffData$numRise.Obs==0]<-0
  DiffData$meanFall.Sim[DiffData$numFall.Obs==0]<-0
  DiffData$fld.Obs[DiffData$numFall.Obs==0]<-0
  DiffData$rld.Obs[DiffData$numRise.Obs==0]<-0
  
  
  # spearman correlations of rising/falling statistics
  if(nrow(DiffData)>=minYears){
    
    # mean_pos_diff
    perfMetrics$r_mean_pos_diff <-rSpearmanMin10(DiffData$meanRise.Sim,DiffData$meanRise.Obs,
                                                 minYears = minYears)
    perfMetrics$NSE_mean_pos_diff <-NSEMin10(DiffData$meanRise.Sim,DiffData$meanRise.Obs,
                                             minYears = minYears)
    kge_comps<- KGEMin10(DiffData$meanRise.Sim,DiffData$meanRise.Obs,
                         minYears = minYears)
    
    perfMetrics$KGE_mean_pos_diff = kge_comps[["kge"]]
    perfMetrics$rp_mean_pos_diff  = kge_comps[["r"]]
    perfMetrics$bias_mean_pos_diff  = kge_comps[["Beta"]]
    perfMetrics$alpha_mean_pos_diff  = kge_comps[["Alpha"]]
    
    # mean_neg_diff
    perfMetrics$r_mean_neg_diff <-rSpearmanMin10(DiffData$meanFall.Sim,DiffData$meanFall.Obs,
                                                 minYears = minYears)
    perfMetrics$NSE_mean_neg_diff <-NSEMin10(DiffData$meanFall.Sim,DiffData$meanFall.Obs,
                                             minYears = minYears)
    kge_comps<- KGEMin10(DiffData$meanFall.Sim,DiffData$meanFall.Obs,
                         minYears = minYears)
    
    perfMetrics$KGE_mean_neg_diff = kge_comps[["kge"]]
    perfMetrics$rp_mean_neg_diff  = kge_comps[["r"]]
    perfMetrics$bias_mean_neg_diff  = kge_comps[["Beta"]]
    perfMetrics$alpha_mean_neg_diff  = kge_comps[["Alpha"]]
    
    
    #no_rises
    perfMetrics$r_no_rises <-rSpearmanMin10(DiffData$numRise.Sim,DiffData$numRise.Obs,
                                            minYears = minYears)
    perfMetrics$NSE_no_rises <-NSEMin10(DiffData$numRise.Sim,DiffData$numRise.Obs,
                                        minYears = minYears)
    
    kge_comps<- KGEMin10(DiffData$numRise.Sim,DiffData$numRise.Obs,
                         minYears = minYears)
    
    perfMetrics$KGE_no_rises = kge_comps[["kge"]]
    perfMetrics$rp_no_rises  = kge_comps[["r"]]
    perfMetrics$bias_no_rises  = kge_comps[["Beta"]]
    perfMetrics$alpha_no_rises  = kge_comps[["Alpha"]]
    
    #no_falls
    perfMetrics$r_no_falls <-rSpearmanMin10(DiffData$numFall.Sim,DiffData$numFall.Obs,
                                            minYears = minYears)
    perfMetrics$NSE_no_falls <-NSEMin10(DiffData$numFall.Sim,DiffData$numFall.Obs,
                                        minYears = minYears)
    
    kge_comps<- KGEMin10(DiffData$numFall.Sim,DiffData$numFall.Obs,
                         minYears = minYears)
    
    perfMetrics$KGE_no_falls = kge_comps[["kge"]]
    perfMetrics$rp_no_falls  = kge_comps[["r"]]
    perfMetrics$bias_no_falls  = kge_comps[["Beta"]]
    perfMetrics$alpha_no_falls  = kge_comps[["Alpha"]]
    
    
    # rld
    perfMetrics$r_rld <-rSpearmanMin10(DiffData$rld.Sim,DiffData$rld.Obs,
                                       minYears = minYears)
    perfMetrics$NSE_rld <-NSEMin10(DiffData$rld.Sim,DiffData$rld.Obs,
                                   minYears = minYears)
    
    
    kge_comps<- KGEMin10(DiffData$rld.Sim,DiffData$rld.Obs,
                         minYears = minYears)
    perfMetrics$KGE_rld = kge_comps[["kge"]]
    perfMetrics$rp_rld  = kge_comps[["r"]]
    perfMetrics$bias_rld  = kge_comps[["Beta"]]
    perfMetrics$alpha_rld  = kge_comps[["Alpha"]]
    
    
    # fld
    perfMetrics$r_fld <-rSpearmanMin10(DiffData$fld.Sim,DiffData$fld.Obs,
                                       minYears = minYears)
    perfMetrics$NSE_fld <-NSEMin10(DiffData$fld.Sim,DiffData$fld.Obs,
                                   minYears = minYears)
    
    
    kge_comps<- KGEMin10(DiffData$fld.Sim,DiffData$fld.Obs,
                         minYears = minYears)
    perfMetrics$KGE_fld = kge_comps[["kge"]]
    perfMetrics$rp_fld  = kge_comps[["r"]]
    perfMetrics$bias_fld  = kge_comps[["Beta"]]
    perfMetrics$alpha_fld  = kge_comps[["Alpha"]]
    
    
    
  }
  
  # baseflow index
  # require a run (no missing values) of at least 365 days
  dat<- dat%>%
    mutate(
      run_id = cumsum(is.na(lag(QObs)) != is.na(QObs))
    ) %>%
    mutate(
      run_id = ifelse(is.na(QObs), NA, run_id)
    )
  
  
  run_id_length<-dat%>%
    group_by(run_id)%>%
    summarize(N = n())%>%
    filter(N>=365&!is.na(run_id))
  if(nrow(run_id_length)>0){
    
    # loop over continuous runs and calculate baseflow
    for(it_run in run_id_length$run_id){
      dat$BF.Obs[dat$run_id %in% it_run]<-bf_sep_lh( dat$QObs[dat$run_id %in% it_run], a = 0.98, n = 3, reflect = 30)
      dat$BF.Sim[dat$run_id %in% it_run]<-bf_sep_lh( dat$QSim[dat$run_id %in% it_run], a = 0.98, n = 3, reflect = 30)
    }
    
    BFIdat<-
      dat%>%
      group_by(wateryear)%>%
      summarize(BFI.Obs = sum(BF.Obs+eps_val)/sum(QObs+eps_val),
                BFI.Sim = sum(BF.Obs+eps_val)/sum(QSim+eps_val))
    
    
    perfMetrics$r_BFI = rSpearmanMin10(BFIdat$BFI.Sim,BFIdat$BFI.Obs,
                                       minYears = minYears)
    
    
    perfMetrics$NSE_BFI <-NSEMin10(BFIdat$BFI.Sim,BFIdat$BFI.Obs,
                                   minYears = minYears)
    
    
    kge_comps<- KGEMin10(BFIdat$BFI.Sim,BFIdat$BFI.Obs,
                         minYears = minYears)
    perfMetrics$KGE_BFI = kge_comps[["kge"]]
    perfMetrics$rp_BFI  = kge_comps[["r"]]
    perfMetrics$bias_BFI  = kge_comps[["Beta"]]
    perfMetrics$alpha_BFI  = kge_comps[["Alpha"]]
    
  }
  
  
  return(perfMetrics)
}


# This function compares a metric across two subsets of a dataframe.
# The function accepts the (possibly spatial) dataframe stns, and optionally
# accepts the name of a column (splitVar) on which to split the data, and 
# thresh, the threshold at which to perform the split.






compareMetricsHighLow<-function(stns,splitVar = "NSEB",thresh = 0.5){
  stns<-st_drop_geometry(stns) # if the dataframe is spatial, drop the geometry
  
  # Define the two subsets based on splitVar and thresh
  stns$split<-stns[,splitVar]>thresh
  
  
  stns<-stns%>%
    mutate(across(starts_with("alpha_"), ~(1-sqrt((1-.x)^2))))%>%
    mutate(across(starts_with("bias_"), ~(1-sqrt((1-.x)^2))))
    
  # Calculate the median metric value for both subsets, for each model
  stn_median<-
    stns%>%
    
    group_by(mdl,split)%>%
    # calculate median of all columns of the dataframe (from NSE to KGE_BFI)
    reframe(across(NSE:KGE_BFI,~median(.x,na.rm = TRUE))    )%>%
    
    pivot_longer(cols = NSE:KGE_BFI,names_to = "metric",values_to = "value")%>%
    
    pivot_wider(names_from = split,names_prefix = "above.",values_from = value)%>%
    arrange(metric,mdl)
  
  #calculate the number of observations for each subset and model
  stn_N<-
    stns%>%
    
    pivot_longer(cols = NSE:KGE_BFI,names_to = "metric",values_to = "value")%>%
    group_by(mdl,split,metric)%>%
    
    summarize(N = sum(!is.na(value)))%>%
    
    pivot_wider(names_from = split,names_prefix = "N.above.",values_from = N)
  
  if(any(stns$split==TRUE,na.rm = TRUE)){ # as long as there are some observation above the threshold
    
    # calculate Wilcoxon rank-sum test for difference between two subsets
    stn_wilcox<-stns%>%
      pivot_longer(cols = NSE:KGE_BFI,names_to = "metric",values_to = "value")%>%
      group_by(mdl,metric)%>%
      do(w.p = ((wilcox.test(value~split, data=., paired=FALSE))$p.value))%>%
      mutate(w.p = as.numeric(w.p))
    
  }else{
    stn_wilcox = data.frame(mdl = unique(stns$mdl),
                            w.p = NA)
  }
  
  
  # Return the medians, the number of observations, and the p-value of the wilcoxon rank-sum test
  stn_result<-left_join(stn_median,stn_N)%>%
    left_join(stn_wilcox)
  return(stn_result)
}


compareMetricsHighLow_v2<-function(stns,splitVar = "NSEB",thresh = 0.5){
  stns<-st_drop_geometry(stns) # if the dataframe is spatial, drop the geometry
  
  # Define the two subsets based on splitVar and thresh
  stns$split<-stns[,splitVar]>thresh
  
  # Calculate the median metric value for both subsets, for each model
  stn_median<-
    stns%>%
    
    group_by(mdl,split)%>%
    # calculate median of all columns of the dataframe (from NSE to r_fld)
    reframe(across(NSE.maxRun:NSE.Rem,~median(.x,na.rm = TRUE))    )%>%
    
    pivot_longer(cols = NSE.maxRun:NSE.Rem,names_to = "metric",values_to = "value")%>%
    
    pivot_wider(names_from = split,names_prefix = "above.",values_from = value)%>%
    arrange(metric,mdl)
  
  #calculate the number of observations for each subset and model
  stn_N<-
    stns%>%
    
    pivot_longer(cols = NSE.maxRun:NSE.Rem, names_to = "metric",values_to = "value")%>%
    group_by(mdl,split,metric)%>%
    
    summarize(N = sum(!is.na(value)))%>%
    
    pivot_wider(names_from = split,names_prefix = "N.above.",values_from = N)
  
  if(any(stns$split==TRUE,na.rm = TRUE)){ # as long as there are some observation above the threshold
    
    # calculate Wilcoxon rank-sum test for difference between two subsets
    stn_wilcox<-stns%>%
      pivot_longer(cols = NSE.maxRun:NSE.Rem,names_to = "metric",values_to = "value")%>%
      group_by(mdl,metric)%>%
      do(w.p = ((wilcox.test(value~split, data=., paired=FALSE))$p.value))%>%
      mutate(w.p = as.numeric(w.p))
    
  }else{
    stn_wilcox = data.frame(mdl = unique(stns$mdl),
                            w.p = NA)
  }
  
  
  # Return the medians, the number of observations, and the p-value of the wilcoxon rank-sum test
  stn_result<-left_join(stn_median,stn_N)%>%
    left_join(stn_wilcox)
  return(stn_result)
}

compareMetricsHighLow_v3<-function(stns,splitVar = "NSEB",thresh = 0.5){
  stns<-st_drop_geometry(stns)%>% # if the dataframe is spatial, drop the geometry
    
    mutate(NSE.Interannual_v2 = 1- (1-NSE.Interannual)*varInterannual_fourier,
           NSE.Seas_v2 = 1- (1-NSE.Seas)*varSeas_fourier,
           NSE.Rem_v2 = 1- (1-NSE.Rem)*varRem_fourier)
  # Define the two subsets based on splitVar and thresh
  stns$split<-stns[,splitVar]>thresh
  
  
  # Calculate the median metric value for both subsets, for each model
  stn_median<-
    stns%>%
    
    
    group_by(mdl,split)%>%
    # calculate median of all columns of the dataframe (from NSE to r_fld)
    reframe(across(c(NSE.maxRun, NSE.Seas,NSE.Interannual,NSE.Rem, NSE.Seas_v2,NSE.Interannual_v2,NSE.Rem_v2),
                   ~median(.x,na.rm = TRUE)))%>%
    
    pivot_longer(cols = c(NSE.maxRun, NSE.Seas,NSE.Interannual,NSE.Rem, NSE.Seas_v2,NSE.Interannual_v2,NSE.Rem_v2),
                 names_to = "metric",values_to = "value")%>%
    
    pivot_wider(names_from = split,names_prefix = "above.",values_from = value)%>%
    arrange(metric,mdl)
  
  #calculate the number of observations for each subset and model
  stn_N<-
    stns%>%
    
    pivot_longer(cols =c(NSE.maxRun, NSE.Seas,NSE.Interannual,NSE.Rem, NSE.Seas_v2,NSE.Interannual_v2,NSE.Rem_v2),
                 names_to = "metric",values_to = "value")%>%
    group_by(mdl,split,metric)%>%
    
    summarize(N = sum(!is.na(value)))%>%
    
    pivot_wider(names_from = split,names_prefix = "N.above.",values_from = N)
  
  if(any(stns$split==TRUE,na.rm = TRUE)){ # as long as there are some observation above the threshold
    
    # calculate Wilcoxon rank-sum test for difference between two subsets
    stn_wilcox<-stns%>%
      pivot_longer(cols =c(NSE.maxRun, NSE.Seas,NSE.Interannual,NSE.Rem, NSE.Seas_v2,NSE.Interannual_v2,NSE.Rem_v2),
                   names_to = "metric",values_to = "value")%>%
      group_by(mdl,metric)%>%
      do(w.p = ((wilcox.test(value~split, data=., paired=FALSE))$p.value))%>%
      mutate(w.p = as.numeric(w.p))
    
  }else{
    stn_wilcox = data.frame(mdl = unique(stns$mdl),
                            w.p = NA)
  }
  
  
  # Return the medians, the number of observations, and the p-value of the wilcoxon rank-sum test
  stn_result<-left_join(stn_median,stn_N)%>%
    left_join(stn_wilcox)
  return(stn_result)
}

# calculate two seasonality metrics (coefficient of variation abd the streamflow concentration index)
calcSeasonality<-function(dat){
  dat<-left_join(data.frame(dt = seq.Date(from = min(dat$dt),
                                          to = max(dat$dt),
                                          by = "1 day")),
                 dat)%>%
    mutate(yday = pmin(yday(dt),365),
           year = year(dt),
           month = month(dt))
  # coefficient of variation
  CoV<-
    dat%>%
    group_by(yday)%>%
    summarize(QObs = mean(QObs,na.rm = TRUE))%>%
    ungroup()%>%
    summarize(CoV = sd(QObs)/mean(QObs))%>%
    pull(CoV)
  
  
  # streamflow concentration index
  QCI = dat%>%
    group_by(year,month)%>%
    summarize(across(c(QObs),~mean(.x)),
              comp = n() ==days_in_month(dt[1]))%>%
    filter(comp& !is.na(QObs))%>%
    group_by(month)%>%
    summarize(QObs = mean(QObs))%>%
    
    summarise(QCI.obs = 100*sum(QObs^2)/sum(QObs)^2,
              
              N = n())%>%
    filter(N ==12& !is.na(QCI.obs))  %>%
    pull(QCI.obs)
  
  if(length(QCI)==0){
    QCI = NA
  }
  if(length(CoV)==0){
    CoV = NA
  }
  
  return(c(
    CoV = CoV,
    QCI = QCI)
  )
}

# these are plotting labels used for Figure 4, to relabel the models

labData<-data.frame(mdl =  c("nearing-lstm",
                             "LSTM-RAPID",
                             "lstm-camels-br",
                             "lstm-Kraft2025",
                             "lstm-arsenault2022",
                             "lstm-kratzert2024",
                             "dHBV2.0_Song",
                             
                             "nearing-glofas",
                             "MGB-SA-camels-br",
                             "cosero-lamah-ce",
                             "prevah-Kraft2025",
                             # "prevah-camels_ch",
                             "nhm-ref",
                             "q_sim_fuse_900",
                             #"q_sim_fuse_902",
                             #"q_sim_fuse_904",
                             "HBV_ub",
                             "mHm_basin",
                             #"mHm_conus",
                             "SAC_SMA" ,
                             "VIC_basin" ,
                             "vic-gl"
                             # "VIC_conus"
                             
                             # "classic_gswp3-w5e5_obsclim_histsoc_default_qtot_global_daily",
                             # 'cwatm_20crv3_obsclim_histsoc_default_qtot_global_daily',
                             # 'cwatm_20crv3-era5_obsclim_histsoc_default_qtot_global_daily',
                             # "cwatm_20crv3-w5e5_obsclim_histsoc_default_qtot_global_daily",
                             # "cwatm_gswp3-w5e5_obsclim_histsoc_default_qtot_global_daily",
                             # "h08_20crv3_obsclim_histsoc_default_qtot_global_daily",
                             # "h08_20crv3-era5_obsclim_histsoc_default_qtot_global_daily",
                             # "h08_20crv3-w5e5_obsclim_histsoc_default_qtot_global_daily",
                             # "h08_gswp3-w5e5_obsclim_histsoc_default_qtot_global_daily",
                             # "hydropy_gswp3-w5e5_obsclim_histsoc_default_qtot_global_daily",
                             # "lpjml5-7-10-fire_20crv3_obsclim_histsoc_default_qtot_global_daily",
                             # "lpjml5-7-10-fire_20crv3-era5_obsclim_histsoc_default_qtot_global_daily",
                             # "lpjml5-7-10-fire_20crv3-w5e5_obsclim_histsoc_default_qtot_global_daily",
                             # "lpjml5-7-10-fire_gswp3-w5e5_obsclim_histsoc_default_qtot_global_daily",
                             # "miroc-integ-land_20crv3_obsclim_histsoc_default_qtot_global_daily",
                             # "miroc-integ-land_20crv3-era5_obsclim_histsoc_default_qtot_global_daily",
                             # "miroc-integ-land_20crv3-w5e5_obsclim_histsoc_default_qtot_global_daily",
                             # "miroc-integ-land_gswp3-w5e5_obsclim_histsoc_default_qtot_global_daily",
                             # "orchidee-mict_gswp3-w5e5_obsclim_histsoc_default_qtot_global_daily",
                             # "watergap2-2e_20crv3-era5_obsclim_histsoc_default_qtot_global_daily",
                             # "watergap2-2e_20crv3-w5e5_obsclim_histsoc_default_qtot_global_daily",
                             # "watergap2-2e_gswp3-w5e5_obsclim_histsoc_default_qtot_global_daily",
                             # "web-dhm-sg_gswp3-w5e5_obsclim_histsoc_default_qtot_global_daily"
                             
                             
                             
),

labels = c("GLOB-LSTM1",
           "GLOB-LSTM2",
           "BR-LSTM",
           "CH-LSTM",
           "ENA-LSTM",
           "US-LSTM",
           "US-deltaHBV2.0UH",
           
           "GLOB-GLOFAS",
           "BR-MGB-SA",
           "CE-COSERO",
           "CH-PREVAH",
           "US-NHM",
           "US-FUSE",
           #"US-FUSE-902",
           #"US-FUSE-904",
           
           "US-HBV",
           "US-mHM",
           #"US-mHM-regional",
           "US-SAC-SMA",
           
           "US-VIC",
           "WNA-VIC-Gl"
           #"US-VIC-regional"
           # "US-HBV-lb",
           
           # "classic_gswp3-w5e5",
           # "cwatm_20crv3",
           # "cwatm_20crv3-era5",
           # "cwatm_20crv3-w5e5",
           # "cwatm_gswp3-w5e5",
           # "h08_20crv3",
           # "h08_20crv3-era5",
           # "h08_20crv3-w5e5",
           # "h08_gswp3-w5e5",
           # "hydropy_gswp3-w5e5",
           # "lpjml5-7-10-fire_20crv3",
           # "lpjml5-7-10-fire_20crv3-era5",
           # "lpjml5-7-10-fire_20crv3-w5e5",
           # "lpjml5-7-10-fire_gswp3-w5e5",
           # "miroc-integ-land_20crv3",
           # "miroc-integ-land_20crv3-era5",
           # "miroc-integ-land_20crv3-w5e5",
           # "miroc-integ-land_gswp3-w5e5",
           # "orchidee-mict_gswp3-w5e5",
           # "watergap2-2e_20crv3-era5",
           # "watergap2-2e_20crv3-w5e5",
           # "watergap2-2e_gswp3-w5e5",
           # "isimip_web-dhm-sg_gswp3-w5e5"
           
))

labData2<-c(
  "nearing-lstm"="GLOB-LSTM1",
  "LSTM-RAPID"="GLOB-LSTM2",
  "lstm-camels-br"= "BR-LSTM",
  "lstm-Kraft2025"="CH-LSTM",
  "lstm-arsenault2022"="ENA-LSTM",
  "lstm-kratzert2024"="US-LSTM",
  "dHBV2.0_Song"=expression("US-" * delta * "HBV2.0UH"),
  
  "nearing-glofas"="GLOB-GLOFAS",
  "MGB-SA-camels-br"="BR-MGB-SA",
  "cosero-lamah-ce"="CE-COSERO",
  "prevah-Kraft2025"="CH-PREVAH",
  "nhm-ref"="US-NHM",
  "q_sim_fuse_900"="US-FUSE",
  "HBV_ub"="US-HBV",
  "mHm_basin"="US-mHM",
  "SAC_SMA"="US-SAC-SMA" ,
  "VIC_basin"="US-VIC" ,
  "vic-gl"="WNA-VIC-Gl"
  
)

# these are plotting labels used for Figure 4, to relabel the metrics
metricLabs<-
  list(NSE = "NSE",
       KGE = "KGE",
       Beta_p = "1-sqrt((1-beta)^2)",
       Gamma_p = "1-sqrt((1-gamma)^2)",
       rPearson = "Pearson~r",
       KGElf = "KGE(1/Q)",
       b1 = "~x", # these were to create spaces in the plot
       # BE="BE",
       r_mag_m1 = "rho(bar(Q)[Jan])",
       r_mag_m2 = "rho(bar(Q)[Feb])",
       r_mag_m3 = "rho(bar(Q)[Mar])",
       r_mag_m4 = "rho(bar(Q)[Apr])",
       r_mag_m5= "rho(bar(Q)[May])",
       r_mag_m6 = "rho(bar(Q)[Jun])",
       r_mag_m7 = "rho(bar(Q)[Jul])",
       r_mag_m8 = "rho(bar(Q)[Aug])",
       r_mag_m9 = "rho(bar(Q)[Sep])",
       r_mag_m10 = "rho(bar(Q)[Oct])",
       r_mag_m11 = "rho(bar(Q)[Nov])",
       r_mag_m12 = "rho(bar(Q)[Dec])",
       
       b2 = " ~~xx",
       
       r_max_1d = "rho(max~1-day~Q)",
       r_max_3d = "rho(max~3-day~Q)",
       r_max_7d = "rho(max~7-day~Q)",
       r_max_30d = "rho(max~30-day~Q)",
       r_max_90d = "rho(max~90-day~Q)",
       r_max_timing = "rho(day~of~annual~maximum)",
       
       b3 = " ~~~xxx ",
       
       r_min_1d = "rho(min~1-day~Q)",
       r_min_3d = "rho(min~3-day~Q)",
       r_min_7d = "rho(min~7-day~Q)",
       r_min_30d = "rho(min~30-day~Q)",
       r_min_90d = "rho(min~90-day~Q)",
       r_min_timing = "rho(day~of~annual~minimum)",
       
       
       b4 = "~~~~ xxxx  ",
       
       r_no_high_pulses = "tau(No.~high~pulses)",
       r_dur_high_pulses = "rho(high~pulse~duration)",
       r_numDays_high = "rho(No.~high~flow~days)",
       
       r_no_low_pulses = "tau(No.~low~pulses)",
       r_dur_low_pulses = "rho(low~pulse~duration)",
       r_numDays_low = "rho(No.~low~flow~days)",
       
       
       b5 = " ~~~~~xxxxx   ",
       
       r_mean_pos_diff = "rho(mean~daily~rise)",
       r_no_rises = "rho(No.~rises)",
       r_rld = "rho(rising~limb~density)",
       r_mean_neg_diff = "rho(mean~daily~fall)",
       r_no_falls = "rho(No.~falls)",
       r_fld = "rho(falling~limb~density)",
       
       b6 = " ~~~~~~  xxxxxx  ",
       r_QCI = "rho(QCI)",
       r_halfFlowDay = "rho(half~flow~day)",
       r_mag_ann = "rho(bar(Q)[annual])",
       r_slope_fdc = "rho(slope~FDC)",
       r_BFI = "rho(BFI)"
       # BE="BE"
       
  )%>%
  data.frame()%>%
  t()%>%
  data.frame()
names(metricLabs)<-"label"
metricLabs$code<-rownames(metricLabs)


# these are plotting labels used for Figure 4, to relabel the metrics

metricLabs_rp<-
  list(NSE = "NSE",
       KGE = "KGE",
       Beta_p = "1-sqrt((1-beta)^2)",
       Gamma_p = "1-sqrt((1-gamma)^2)",
       rPearson = "Pearson~r",
       KGElf = "KGE(1/Q)",
       b1 = "~x", # these were to create spaces in the plot
       # BE="BE",
       rp_mag_m1 = "r(bar(Q)[Jan])",
       rp_mag_m2 = "r(bar(Q)[Feb])",
       rp_mag_m3 = "r(bar(Q)[Mar])",
       rp_mag_m4 = "r(bar(Q)[Apr])",
       rp_mag_m5= "r(bar(Q)[May])",
       rp_mag_m6 = "r(bar(Q)[Jun])",
       rp_mag_m7 = "r(bar(Q)[Jul])",
       rp_mag_m8 = "r(bar(Q)[Aug])",
       rp_mag_m9 = "r(bar(Q)[Sep])",
       rp_mag_m10 = "r(bar(Q)[Oct])",
       rp_mag_m11 = "r(bar(Q)[Nov])",
       rp_mag_m12 = "r(bar(Q)[Dec])",
       
       b2 = " ~~xx",
       
       rp_max_1d = "r(max~1-day~Q)",
       rp_max_3d = "r(max~3-day~Q)",
       rp_max_7d = "r(max~7-day~Q)",
       rp_max_30d = "r(max~30-day~Q)",
       rp_max_90d = "r(max~90-day~Q)",
       rp_max_timing = "r(day~of~annual~maximum)",
       
       b3 = " ~~~xxx ",
       
       rp_min_1d = "r(min~1-day~Q)",
       rp_min_3d = "r(min~3-day~Q)",
       rp_min_7d = "r(min~7-day~Q)",
       rp_min_30d = "r(min~30-day~Q)",
       rp_min_90d = "r(min~90-day~Q)",
       rp_min_timing = "r(day~of~annual~minimum)",
       
       
       b4 = "~~~~ xxxx  ",
       
       rp_no_high_pulses = "r(No.~high~pulses)",
       rp_dur_high_pulses = "r(high~pulse~duration)",
       rp_numDays_high = "r(No.~high~flow~days)",
       
       rp_no_low_pulses = "r(No.~low~pulses)",
       rp_dur_low_pulses = "r(low~pulse~duration)",
       rp_numDays_low = "r(No.~low~flow~days)",
       
       
       b5 = " ~~~~~xxxxx   ",
       
       rp_mean_pos_diff = "r(mean~daily~rise)",
       rp_no_rises = "r(No.~rises)",
       rp_rld = "r(rising~limb~density)",
       rp_mean_neg_diff = "r(mean~daily~fall)",
       rp_no_falls = "r(No.~falls)",
       rp_fld = "r(falling~limb~density)",
       
       b6 = " ~~~~~~  xxxxxx  ",
       rp_QCI = "r(QCI)",
       rp_halfFlowDay = "r(half~flow~day)",
       rp_mag_ann = "r(bar(Q)[annual])",
       rp_slope_fdc = "r(slope~FDC)",
       rp_BFI = "r(BFI)"
       # BE="BE"
       
  )%>%
  data.frame()%>%
  t()%>%
  data.frame()
names(metricLabs_rp)<-"label"
metricLabs_rp$code<-rownames(metricLabs_rp)



# these are plotting labels used for Figure 4, with NSE, to relabel the metrics

metricLabs_NSE<-
  list(NSE = "NSE",
       KGE = "KGE",
       Beta_p = "1-sqrt((1-beta)^2)",
       Gamma_p = "1-sqrt((1-gamma)^2)",
       rPearson = "Pearson~r",
       KGElf = "KGE(1/Q)",
       b1 = "~x", # these were to create spaces in the plot
       # BE="BE",
       NSE_mag_m1 = "NSE(bar(Q)[Jan])",
       NSE_mag_m2 = "NSE(bar(Q)[Feb])",
       NSE_mag_m3 = "NSE(bar(Q)[Mar])",
       NSE_mag_m4 = "NSE(bar(Q)[Apr])",
       NSE_mag_m5= "NSE(bar(Q)[May])",
       NSE_mag_m6 = "NSE(bar(Q)[Jun])",
       NSE_mag_m7 = "NSE(bar(Q)[Jul])",
       NSE_mag_m8 = "NSE(bar(Q)[Aug])",
       NSE_mag_m9 = "NSE(bar(Q)[Sep])",
       NSE_mag_m10 = "NSE(bar(Q)[Oct])",
       NSE_mag_m11 = "NSE(bar(Q)[Nov])",
       NSE_mag_m12 = "NSE(bar(Q)[Dec])",
       
       b2 = " ~~xx",
       
       NSE_max_1d = "NSE(max~1-day~Q)",
       NSE_max_3d = "NSE(max~3-day~Q)",
       NSE_max_7d = "NSE(max~7-day~Q)",
       NSE_max_30d = "NSE(max~30-day~Q)",
       NSE_max_90d = "NSE(max~90-day~Q)",
       NSE_max_timing = "NSE(day~of~annual~maximum)",
       
       b3 = " ~~~xxx ",
       
       NSE_min_1d = "NSE(min~1-day~Q)",
       NSE_min_3d = "NSE(min~3-day~Q)",
       NSE_min_7d = "NSE(min~7-day~Q)",
       NSE_min_30d = "NSE(min~30-day~Q)",
       NSE_min_90d = "NSE(min~90-day~Q)",
       NSE_min_timing = "NSE(day~of~annual~minimum)",
       
       
       b4 = "~~~~ xxxx  ",
       
       NSE_no_high_pulses = "NSE(No.~high~pulses)",
       NSE_dur_high_pulses = "NSE(high~pulse~duration)",
       NSE_numDays_high = "NSE(No.~high~flow~days)",
       
       NSE_no_low_pulses = "NSE(No.~low~pulses)",
       NSE_dur_low_pulses = "NSE(low~pulse~duration)",
       NSE_numDays_low = "NSE(No.~low~flow~days)",
       
       
       b5 = " ~~~~~xxxxx   ",
       
       NSE_mean_pos_diff = "NSE(mean~daily~rise)",
       NSE_no_rises = "NSE(No.~rises)",
       NSE_rld = "NSE(rising~limb~density)",
       NSE_mean_neg_diff = "NSE(mean~daily~fall)",
       NSE_no_falls = "NSE(No.~falls)",
       NSE_fld = "NSE(falling~limb~density)",
       
       b6 = " ~~~~~~  xxxxxx  ",
       NSE_QCI = "NSE(QCI)",
       NSE_halfFlowDay = "NSE(half~flow~day)",
       NSE_mag_ann = "NSE(bar(Q)[annual])",
       NSE_slope_fdc = "NSE(slope~FDC)",
       NSE_BFI = "NSE(BFI)"
       # BE="BE"
       
  )%>%
  data.frame()%>%
  t()%>%
  data.frame()
names(metricLabs_NSE)<-"label"
metricLabs_NSE$code<-rownames(metricLabs_NSE)


# these are plotting labels used for Figure 4, to relabel the metrics
metricLabs_KGE<-
  list(NSE = "NSE",
       KGE = "KGE",
       Beta_p = "1-sqrt((1-beta)^2)",
       Gamma_p = "1-sqrt((1-gamma)^2)",
       rPearson = "Pearson~r",
       KGElf = "KGE(1/Q)",
       b1 = "~x", # these were to create spaces in the plot
       # BE="BE",
       KGE_mag_m1 = "KGE(bar(Q)[Jan])",
       KGE_mag_m2 = "KGE(bar(Q)[Feb])",
       KGE_mag_m3 = "KGE(bar(Q)[Mar])",
       KGE_mag_m4 = "KGE(bar(Q)[Apr])",
       KGE_mag_m5= "KGE(bar(Q)[May])",
       KGE_mag_m6 = "KGE(bar(Q)[Jun])",
       KGE_mag_m7 = "KGE(bar(Q)[Jul])",
       KGE_mag_m8 = "KGE(bar(Q)[Aug])",
       KGE_mag_m9 = "KGE(bar(Q)[Sep])",
       KGE_mag_m10 = "KGE(bar(Q)[Oct])",
       KGE_mag_m11 = "KGE(bar(Q)[Nov])",
       KGE_mag_m12 = "KGE(bar(Q)[Dec])",
       
       b2 = " ~~xx",
       
       KGE_max_1d = "KGE(max~1-day~Q)",
       KGE_max_3d = "KGE(max~3-day~Q)",
       KGE_max_7d = "KGE(max~7-day~Q)",
       KGE_max_30d = "KGE(max~30-day~Q)",
       KGE_max_90d = "KGE(max~90-day~Q)",
       KGE_max_timing = "KGE(day~of~annual~maximum)",
       
       b3 = " ~~~xxx ",
       
       KGE_min_1d = "KGE(min~1-day~Q)",
       KGE_min_3d = "KGE(min~3-day~Q)",
       KGE_min_7d = "KGE(min~7-day~Q)",
       KGE_min_30d = "KGE(min~30-day~Q)",
       KGE_min_90d = "KGE(min~90-day~Q)",
       KGE_min_timing = "KGE(day~of~annual~minimum)",
       
       
       b4 = "~~~~ xxxx  ",
       
       KGE_no_high_pulses = "KGE(No.~high~pulses)",
       KGE_dur_high_pulses = "KGE(high~pulse~duration)",
       KGE_numDays_high = "KGE(No.~high~flow~days)",
       
       KGE_no_low_pulses = "KGE(No.~low~pulses)",
       KGE_dur_low_pulses = "KGE(low~pulse~duration)",
       KGE_numDays_low = "KGE(No.~low~flow~days)",
       
       
       b5 = " ~~~~~xxxxx   ",
       
       KGE_mean_pos_diff = "KGE(mean~daily~rise)",
       KGE_no_rises = "KGE(No.~rises)",
       KGE_rld = "KGE(rising~limb~density)",
       KGE_mean_neg_diff = "KGE(mean~daily~fall)",
       KGE_no_falls = "KGE(No.~falls)",
       KGE_fld = "KGE(falling~limb~density)",
       
       b6 = " ~~~~~~  xxxxxx  ",
       KGE_QCI = "KGE(QCI)",
       KGE_halfFlowDay = "KGE(half~flow~day)",
       KGE_mag_ann = "KGE(bar(Q)[annual])",
       KGE_slope_fdc = "KGE(slope~FDC)",
       KGE_BFI = "KGE(BFI)"
       # BE="BE"
       
  )%>%
  data.frame()%>%
  t()%>%
  data.frame()
names(metricLabs_KGE)<-"label"
metricLabs_KGE$code<-rownames(metricLabs_KGE)

# these are plotting labels used for Figure 4, to relabel the metrics
metricLabs_bias<-
  list(NSE = "NSE",
       KGE = "KGE",
       Beta_p = "1-sqrt((1-beta)^2)",
       Gamma_p = "1-sqrt((1-gamma)^2)",
       rPearson = "Pearson~r",
       KGElf = "KGE(1/Q)",
       b1 = "~x", # these were to create spaces in the plot
       # BE="BE",
       bias_mag_m1 = "beta*\"'\"*(bar(Q)[Jan])",
       bias_mag_m2 = "beta*\"'\"*(bar(Q)[Feb])",
       bias_mag_m3 = "beta*\"'\"*(bar(Q)[Mar])",
       bias_mag_m4 = "beta*\"'\"*(bar(Q)[Apr])",
       bias_mag_m5= "beta*\"'\"*(bar(Q)[May])",
       bias_mag_m6 = "beta*\"'\"*(bar(Q)[Jun])",
       bias_mag_m7 = "beta*\"'\"*(bar(Q)[Jul])",
       bias_mag_m8 = "beta*\"'\"*(bar(Q)[Aug])",
       bias_mag_m9 = "beta*\"'\"*(bar(Q)[Sep])",
       bias_mag_m10 = "beta*\"'\"*(bar(Q)[Oct])",
       bias_mag_m11 = "beta*\"'\"*(bar(Q)[Nov])",
       bias_mag_m12 = "beta*\"'\"*(bar(Q)[Dec])",
       
       b2 = " ~~xx",
       
       bias_max_1d = "beta*\"'\"*(max~1-day~Q)",
       bias_max_3d = "beta*\"'\"*(max~3-day~Q)",
       bias_max_7d = "beta*\"'\"*(max~7-day~Q)",
       bias_max_30d = "beta*\"'\"*(max~30-day~Q)",
       bias_max_90d = "beta*\"'\"*(max~90-day~Q)",
       bias_max_timing = "beta*\"'\"*(day~of~annual~maximum)",
       
       b3 = " ~~~xxx ",
       
       bias_min_1d = "beta*\"'\"*(min~1-day~Q)",
       bias_min_3d = "beta*\"'\"*(min~3-day~Q)",
       bias_min_7d = "beta*\"'\"*(min~7-day~Q)",
       bias_min_30d = "beta*\"'\"*(min~30-day~Q)",
       bias_min_90d = "beta*\"'\"*(min~90-day~Q)",
       bias_min_timing = "beta*\"'\"*(day~of~annual~minimum)",
       
       
       b4 = "~~~~ xxxx  ",
       
       bias_no_high_pulses = "beta*\"'\"*(No.~high~pulses)",
       bias_dur_high_pulses = "beta*\"'\"*(high~pulse~duration)",
       bias_numDays_high = "beta*\"'\"*(No.~high~flow~days)",
       
       bias_no_low_pulses = "beta*\"'\"*(No.~low~pulses)",
       bias_dur_low_pulses = "beta*\"'\"*(low~pulse~duration)",
       bias_numDays_low = "beta*\"'\"*(No.~low~flow~days)",
       
       
       b5 = " ~~~~~xxxxx   ",
       
       bias_mean_pos_diff = "beta*\"'\"*(mean~daily~rise)",
       bias_no_rises = "beta*\"'\"*(No.~rises)",
       bias_rld = "beta*\"'\"*(rising~limb~density)",
       bias_mean_neg_diff = "beta*\"'\"*(mean~daily~fall)",
       bias_no_falls = "beta*\"'\"*(No.~falls)",
       bias_fld = "beta*\"'\"*(falling~limb~density)",
       
       b6 = " ~~~~~~  xxxxxx  ",
       bias_QCI = "beta*\"'\"*(QCI)",
       bias_halfFlowDay = "beta*\"'\"*(half~flow~day)",
       bias_mag_ann = "beta*\"'\"*(bar(Q)[annual])",
       bias_slope_fdc = "beta*\"'\"*(slope~FDC)",
       bias_BFI = "beta*\"'\"*(BFI)"
       # BE="BE"
       
  )%>%
  data.frame()%>%
  t()%>%
  data.frame()
names(metricLabs_bias)<-"label"
metricLabs_bias$code<-rownames(metricLabs_bias)



# these are plotting labels used for Figure 4, to relabel the metrics
metricLabs_alpha<-
  list(NSE = "NSE",
       KGE = "KGE",
       Beta_p = "1-sqrt((1-beta)^2)",
       Gamma_p = "1-sqrt((1-gamma)^2)",
       rPearson = "Pearson~r",
       KGElf = "KGE(1/Q)",
       b1 = "~x", # these were to create spaces in the plot
       # BE="BE",
       alpha_mag_m1 = "alpha*\"'\"*(bar(Q)[Jan])",
       alpha_mag_m2 = "alpha*\"'\"*(bar(Q)[Feb])",
       alpha_mag_m3 = "alpha*\"'\"*(bar(Q)[Mar])",
       alpha_mag_m4 = "alpha*\"'\"*(bar(Q)[Apr])",
       alpha_mag_m5= "alpha*\"'\"*(bar(Q)[May])",
       alpha_mag_m6 = "alpha*\"'\"*(bar(Q)[Jun])",
       alpha_mag_m7 = "alpha*\"'\"*(bar(Q)[Jul])",
       alpha_mag_m8 = "alpha*\"'\"*(bar(Q)[Aug])",
       alpha_mag_m9 = "alpha*\"'\"*(bar(Q)[Sep])",
       alpha_mag_m10 = "alpha*\"'\"*(bar(Q)[Oct])",
       alpha_mag_m11 = "alpha*\"'\"*(bar(Q)[Nov])",
       alpha_mag_m12 = "alpha*\"'\"*(bar(Q)[Dec])",
       
       b2 = " ~~xx",
       
       alpha_max_1d = "alpha*\"'\"*(max~1-day~Q)",
       alpha_max_3d = "alpha*\"'\"*(max~3-day~Q)",
       alpha_max_7d = "alpha*\"'\"*(max~7-day~Q)",
       alpha_max_30d = "alpha*\"'\"*(max~30-day~Q)",
       alpha_max_90d = "alpha*\"'\"*(max~90-day~Q)",
       alpha_max_timing = "alpha*\"'\"*(day~of~annual~maximum)",
       
       b3 = " ~~~xxx ",
       
       alpha_min_1d = "alpha*\"'\"*(min~1-day~Q)",
       alpha_min_3d = "alpha*\"'\"*(min~3-day~Q)",
       alpha_min_7d = "alpha*\"'\"*(min~7-day~Q)",
       alpha_min_30d = "alpha*\"'\"*(min~30-day~Q)",
       alpha_min_90d = "alpha*\"'\"*(min~90-day~Q)",
       alpha_min_timing = "alpha*\"'\"*(day~of~annual~minimum)",
       
       
       b4 = "~~~~ xxxx  ",
       
       alpha_no_high_pulses = "alpha*\"'\"*(No.~high~pulses)",
       alpha_dur_high_pulses = "alpha*\"'\"*(high~pulse~duration)",
       alpha_numDays_high = "alpha*\"'\"*(No.~high~flow~days)",
       
       alpha_no_low_pulses = "alpha*\"'\"*(No.~low~pulses)",
       alpha_dur_low_pulses = "alpha*\"'\"*(low~pulse~duration)",
       alpha_numDays_low = "alpha*\"'\"*(No.~low~flow~days)",
       
       
       b5 = " ~~~~~xxxxx   ",
       
       alpha_mean_pos_diff = "alpha*\"'\"*(mean~daily~rise)",
       alpha_no_rises = "alpha*\"'\"*(No.~rises)",
       alpha_rld = "alpha*\"'\"*(rising~limb~density)",
       alpha_mean_neg_diff = "alpha*\"'\"*(mean~daily~fall)",
       alpha_no_falls = "alpha*\"'\"*(No.~falls)",
       alpha_fld = "alpha*\"'\"*(falling~limb~density)",
       
       b6 = " ~~~~~~  xxxxxx  ",
       alpha_QCI = "alpha*\"'\"*(QCI)",
       alpha_halfFlowDay = "alpha*\"'\"*(half~flow~day)",
       alpha_mag_ann = "alpha*\"'\"*(bar(Q)[annual])",
       alpha_slope_fdc = "alpha*\"'\"*(slope~FDC)",
       alpha_BFI = "alpha*\"'\"*(BFI)"
       # BE="BE"
       
  )%>%
  data.frame()%>%
  t()%>%
  data.frame()
names(metricLabs_alpha)<-"label"
metricLabs_alpha$code<-rownames(metricLabs_alpha)





# help function to color the arrows in figure 4
func_color<-function(lowNSEVal,highNSEVal,w.p){
  x<-rep(NA,length(w.p))
  x<-factor(x,levels = c("Significantly Higher","Higher","Equal","Lower","Significantly Lower"))
  x[lowNSEVal<highNSEVal& w.p<0.05]<-("Significantly Higher")
  x[lowNSEVal<highNSEVal& w.p>=0.05]<-("Higher")
  x[lowNSEVal>highNSEVal& w.p>=0.05]<-("Lower")
  x[lowNSEVal>highNSEVal& w.p<0.05]<-("Significantly Lower")
  x[lowNSEVal==highNSEVal]<-("Equal")
  return(x)
}

# Main plotting function for Figure 4
# Accepts the following arguments:
# stn_all: a dataframe (produced by compareMetricsHighLow),
# flName: the name of a file to save to,
# groupLabels:legend labels for the two subsets

plotComps<-
  function(stn_all,flName,groupLabels, stat = "r",minXVal = 0,base_size=5,strip.text.y= 7,
           transform = "identity",
           breaks = c(-1,0,1),
           minor_breaks = c(-0.5,0.5)){
    
    if(stat == "r"){
      # get rid of blanks in the  metricLabs dataframe
      metricLabs1<-metricLabs%>%dplyr::filter(!code %in% c("b1","b2","b3","b4","b5","b6"))
      
    }else if(stat == "NSE"){
      # get rid of blanks in the  metricLabs dataframe
      metricLabs1<-metricLabs_NSE%>%dplyr::filter(!code %in% c("b1","b2","b3","b4","b5","b6"))
      
    }else if(stat == "KGE"){
      # get rid of blanks in the  metricLabs dataframe
      metricLabs1<-metricLabs_KGE%>%dplyr::filter(!code %in% c("b1","b2","b3","b4","b5","b6"))
      
    }else if(stat == "bias"){
      # get rid of blanks in the  metricLabs dataframe
      metricLabs1<-metricLabs_bias%>%dplyr::filter(!code %in% c("b1","b2","b3","b4","b5","b6"))
      
    }else if(stat == "alpha"){
      # get rid of blanks in the  metricLabs dataframe
      metricLabs1<-metricLabs_alpha%>%dplyr::filter(!code %in% c("b1","b2","b3","b4","b5","b6"))
      
    }else if(stat == "rp"){
      # get rid of blanks in the  metricLabs dataframe
      metricLabs1<-metricLabs_rp%>%dplyr::filter(!code %in% c("b1","b2","b3","b4","b5","b6"))
      
    }
    
    
    metricSummary<-stn_all%>%
      mutate(leg_col = func_color(above.FALSE,above.TRUE,w.p))%>%
      mutate(highNSE_worse = above.TRUE<above.FALSE)%>%
      filter(mdl %in% labData$mdl)%>%
      
      
      filter( (substr(metric,1,2)=="r_")
      )%>%
      group_by(mdl)%>%
      summarize(SigHigher = sum(leg_col == "Significantly Higher",na.rm = TRUE),
                Higher = sum(leg_col == "Higher"),
                Equal = sum(leg_col == "Equal"),
                Lower = sum(leg_col == "Lower"),
                SigLower = sum(leg_col == "Significantly Lower"),
                
                percLower_or_sigLower = mean(leg_col%in% c("Lower","Significantly Lower")),
                
                N = sum(!is.na(leg_col))
      )
    
    
    stn_all$metric_prefix<-stn_all$metric%>%
      str_split_fixed("_",2)%>%
      .[,1]
    
    # There are two dataframes in order to plot the points and segments separately
    # ggplot had trouble coloring the segments correctly when also using facets
    
    gg_dat1<-
      stn_all%>%
      mutate(leg_col = func_color(above.FALSE,above.TRUE,w.p))%>%
      mutate(above_worse = above.TRUE<above.FALSE)%>%
      pivot_longer(cols = c(above.FALSE,above.TRUE))%>%
      filter(!metric %in% c("KGEB","NSEB"))%>%
      filter(mdl %in% labData$mdl)%>%
      
      mutate(
        metricLab = plyr::mapvalues(metric,from = metricLabs1$code,to = metricLabs1$label),
        metricLab = factor(metricLab,
                           levels = metricLabs1$label,ordered = TRUE),#,levels = c("NSE","KGE","NSE_max","NSE_min7","R2_halfFlowDay","r_halfFlowDay","R2_QCI","r_QCI","R2_lowFlowDays","r_lowFlowDays","BE")),
        mdlLab = factor(mdl,labData$mdl,
                        labels = labData$labels,ordered = TRUE))%>%
      filter(metric %in% c("NSE","KGE","Beta_p","Gamma_p","rPearson","KGElf","b1","b2","b3","b4","b5","b6")|
               metric_prefix==stat)
    
    
    gg_dat2<-
      stn_all%>%
      
      mutate(leg_col2 = func_color(above.FALSE,above.TRUE,w.p))%>%
      # mutate(above = above.TRUE<above.FALSE)%>%
      
      filter(!metric %in% c("KGEB","NSEB"))%>%
      filter(mdl %in% labData$mdl)%>%
      
      mutate(
        metricLab = plyr::mapvalues(metric,from = metricLabs1$code,to = metricLabs1$label),
        metricLab = factor(metricLab,
                           levels = metricLabs1$label,ordered = TRUE),#,levels = c("NSE","KGE","NSE_max","NSE_min7","R2_halfFlowDay","r_halfFlowDay","R2_QCI","r_QCI","R2_lowFlowDays","r_lowFlowDays","BE")),
        mdlLab = factor(mdl,labData$mdl,
                        labels = labData$labels,ordered = TRUE))%>%
      
      filter(metric %in% c("NSE","KGE","Beta_p","Gamma_p","rPearson","KGElf","b1","b2","b3","b4","b5","b6")|
               metric_prefix==stat) %>%
      
      
      mutate(arrowL = plyr::mapvalues(leg_col2,from = c("Significantly Higher","Higher","Equal","Lower","Significantly Lower"),
                                      to = c(0.04,0.02,0.001,0.02,0.04))%>%
               as.character()%>%
               as.numeric()
      ) 
    
    
    
    
    unique(gg_dat1$metric)
    unique(gg_dat2$metric)
    
    
    ordered_levels <- names(labData2)
    
    gg_dat1$mdl <- factor(gg_dat1$mdl, levels = ordered_levels)
    gg_dat2$mdl <- factor(gg_dat2$mdl, levels = ordered_levels)
    # Parse character labels to expressions where needed
    label_list <- unname(labData2)
    names(label_list) <- names(labData2)
    
    p<-
      ggplot(gg_dat1,
             aes(y = mdl,x= value),size = 0.4)+
      lapply(
        split(gg_dat2,gg_dat2$metricLab),
        function(df)
          geom_segment(data = df,
                       aes(y = mdl,yend = mdl,x= above.FALSE,xend = above.TRUE,col = leg_col2),
                       arrow = arrow(length = unit(df$arrowL,  "inches"),
                                     ends = "last", type = "open"))
        
        
      )+
      geom_point(aes(shape = name,col = leg_col),show.legend = c(color = FALSE),size = 0.5)+
      
      facet_wrap("metricLab",scales = "free_x",ncol = 6,strip.position = "right",
                 labeller = label_parsed)+
      
      scale_x_continuous(name = "Median value of performance metric",limits = c(minXVal,1),
                         breaks = breaks,
                         minor_breaks = minor_breaks,
                         oob = scales::oob_keep,
                         transform = transform)+
      scale_y_discrete(name = NULL,
                       limits = rev(ordered_levels),  # Reverse the desired order
                       labels = label_list) +
      scale_color_manual(name = "Significance",
                         breaks = c("Significantly Higher","Higher","Equal","Lower","Significantly Lower"),
                         values = rev(c("#7B3294", "#C2A5CF", "grey", "#A6DBA0", "#008837")))+
      
      scale_shape_manual(name = "Subset",
                         breaks = c("above.FALSE","above.TRUE"),
                         values = c(1,2),
                         labels = groupLabels)+
      theme_bw(base_size = base_size)+
      theme(strip.text.y = element_text(size = strip.text.y))
    
    ggsave(plot = p,filename = flName,width = 7,height = 8,dpi = 600)
    
  }

# Define the custom transform function
reverse_log10_shifted_trans <- function() {
  trans_new(
    name = "reverse_log10_shifted",
    
    transform = function(x) {
      y <- -(x - 1.1)
      ifelse(y > 0, -log10(y), NA_real_)
    },
    
    inverse = function(x) {
      1.1 - 10^(-x)
    },
    
    domain = c(-Inf, 1.1)
  )
}


# v2

metricLabs_v2<-
  list(NSE.maxRun = "NSE",
       NSE.Interannual = "NSE~(Interannual)",
       NSE.Seas = "NSE~(Seasonal)",
       NSE.Rem = "NSE~(Irregular)"
       
       
  )%>%
  data.frame()%>%
  t()%>%
  data.frame()
names(metricLabs_v2)<-"label"
metricLabs_v2$code<-rownames(metricLabs_v2)

plotComps_v2<-
  function(stn_all,flName,groupLabels){
    
    
    
    metricSummary<-stn_all%>%
      mutate(leg_col = func_color(above.FALSE,above.TRUE,w.p))%>%
      mutate(highNSE_worse = above.TRUE<above.FALSE)%>%
      filter(mdl %in% labData$mdl)%>%
      
      
      
      group_by(mdl)%>%
      summarize(SigHigher = sum(leg_col == "Significantly Higher",na.rm = TRUE),
                Higher = sum(leg_col == "Higher"),
                Equal = sum(leg_col == "Equal"),
                Lower = sum(leg_col == "Lower"),
                SigLower = sum(leg_col == "Significantly Lower"),
                
                percLower_or_sigLower = mean(leg_col%in% c("Lower","Significantly Lower")),
                
                N = sum(!is.na(leg_col))
      )
    
    
    # get rid of blanks in the  metricLabs dataframe
    metricLabs1<-metricLabs_v2
    
    
    # There are two dataframes in order to plot the points and segments separately
    # ggplot had trouble coloring the segments correctly when also using facets
    
    gg_dat1<-
      stn_all%>%
      mutate(leg_col = func_color(above.FALSE,above.TRUE,w.p))%>%
      mutate(above_worse = above.TRUE<above.FALSE)%>%
      pivot_longer(cols = c(above.FALSE,above.TRUE))%>%
      filter(!metric %in% c("KGEB","NSEB"))%>%
      filter(mdl %in% labData$mdl)%>%
      
      mutate(
        metricLab = plyr::mapvalues(metric,from = metricLabs1$code,to = metricLabs1$label),
        metricLab = factor(metricLab,
                           levels = metricLabs1$label,ordered = TRUE)#,levels = c("NSE","KGE","NSE_max","NSE_min7","R2_halfFlowDay","r_halfFlowDay","R2_QCI","r_QCI","R2_lowFlowDays","r_lowFlowDays","BE")),
        # mdlLab = factor(mdl,labData$mdl,
        #                 labels = labData$labels,ordered = TRUE)
      )
    
    
    
    
    gg_dat2<-
      stn_all%>%
      mutate(leg_col2 = func_color(above.FALSE,above.TRUE,w.p))%>%
      # mutate(above = above.TRUE<above.FALSE)%>%
      
      filter(!metric %in% c("KGEB","NSEB"))%>%
      filter(mdl %in% labData$mdl)%>%
      
      mutate(
        metricLab = plyr::mapvalues(metric,from = metricLabs1$code,to = metricLabs1$label),
        metricLab = factor(metricLab,
                           levels = metricLabs1$label,ordered = TRUE)#,levels = c("NSE","KGE","NSE_max","NSE_min7","R2_halfFlowDay","r_halfFlowDay","R2_QCI","r_QCI","R2_lowFlowDays","r_lowFlowDays","BE")),
        # mdlLab = factor(mdl,labData$mdl,
        #                 labels = labData$labels,ordered = TRUE)
      )%>%
      
      
      mutate(arrowL = plyr::mapvalues(leg_col2,from = c("Significantly Higher","Higher","Equal","Lower","Significantly Lower"),
                                      to = c(0.04,0.02,0.001,0.02,0.04))%>%
               as.character()%>%
               as.numeric()
      )
    
    unique(gg_dat1$metric)
    
    ordered_levels <- names(labData2)
    
    gg_dat1$mdl <- factor(gg_dat1$mdl, levels = ordered_levels)
    gg_dat2$mdl <- factor(gg_dat2$mdl, levels = ordered_levels)
    # Parse character labels to expressions where needed
    label_list <- unname(labData2)
    names(label_list) <- names(labData2)
    
    
    p<-
      ggplot(gg_dat1,
             aes(y = mdl,x= value),size = 0.4)+
      lapply(
        split(gg_dat2,gg_dat2$metricLab),
        function(df)
          geom_segment(data = df,
                       aes(y = mdl,yend = mdl,x= above.FALSE,xend = above.TRUE,col = leg_col2),
                       arrow = arrow(length = unit(df$arrowL,  "inches"),
                                     ends = "last", type = "open"))
        
        
      )+
      geom_point(aes(shape = name,col = leg_col),show.legend = c(color = FALSE),size = 0.5)+
      
      facet_wrap("metricLab",scales = "free_x",ncol = 1,strip.position = "right",
                 labeller = label_parsed)+
      
      scale_x_continuous(name = "Median NSE",limits = c(0,1),
                         breaks = c(0,1),
                         oob = scales::oob_keep)+
      scale_y_discrete(name = NULL,
                       limits = rev(ordered_levels),  # Reverse the desired order
                       labels = label_list) +
      scale_color_manual(name = "Significance",
                         breaks = c("Significantly Higher","Higher","Equal","Lower","Significantly Lower"),
                         values = rev(c("#7B3294", "#C2A5CF", "grey", "#A6DBA0", "#008837")))+
      
      scale_shape_manual(name = "Subset",
                         breaks = c("above.FALSE","above.TRUE"),
                         values = c(1,2),
                         labels = groupLabels)+
      theme_bw(base_size = 5)
    
    ggsave(plot = p,filename = flName,width = 3,height = 4,dpi = 600)
    
  }



metricLabs_v3<-
  list(NSE.maxRun = "NSE",
       NSE.Interannual = "NSE~(Interannual)",
       NSE.Seas = "NSE~(Seasonal)",
       NSE.Rem = "NSE~(Irregular)",
       NSE.Interannual_v2 = "varNSE~(Interannual)",
       NSE.Seas_v2 = "varNSE~(Seasonal)",
       NSE.Rem_v2 = "varNSE~(Irregular)"
       
       
  )%>%
  data.frame()%>%
  t()%>%
  data.frame()
names(metricLabs_v3)<-"label"
metricLabs_v3$code<-rownames(metricLabs_v3)
plotComps_v3<-
  function(stn_all,flName,groupLabels){
    
    
    
    metricSummary<-stn_all%>%
      mutate(leg_col = func_color(above.FALSE,above.TRUE,w.p))%>%
      mutate(highNSE_worse = above.TRUE<above.FALSE)%>%
      filter(mdl %in% labData$mdl)%>%
      
      
      
      group_by(mdl)%>%
      summarize(SigHigher = sum(leg_col == "Significantly Higher",na.rm = TRUE),
                Higher = sum(leg_col == "Higher"),
                Equal = sum(leg_col == "Equal"),
                Lower = sum(leg_col == "Lower"),
                SigLower = sum(leg_col == "Significantly Lower"),
                
                percLower_or_sigLower = mean(leg_col%in% c("Lower","Significantly Lower")),
                
                N = sum(!is.na(leg_col))
      )
    
    
    # get rid of blanks in the  metricLabs dataframe
    metricLabs1<-metricLabs_v3
    
    
    
    # There are two dataframes in order to plot the points and segments separately
    # ggplot had trouble coloring the segments correctly when also using facets
    
    gg_dat1<-
      stn_all%>%
      mutate(leg_col = func_color(above.FALSE,above.TRUE,w.p))%>%
      mutate(above_worse = above.TRUE<above.FALSE)%>%
      pivot_longer(cols = c(above.FALSE,above.TRUE))%>%
      filter(!metric %in% c("KGEB","NSEB"))%>%
      filter(mdl %in% labData$mdl)%>%
      
      mutate(
        metricLab = plyr::mapvalues(metric,from = metricLabs1$code,to = metricLabs1$label),
        metricLab = factor(metricLab,
                           levels = metricLabs1$label,ordered = TRUE),#,levels = c("NSE","KGE","NSE_max","NSE_min7","R2_halfFlowDay","r_halfFlowDay","R2_QCI","r_QCI","R2_lowFlowDays","r_lowFlowDays","BE")),
        mdlLab = factor(mdl,labData$mdl,
                        labels = labData$labels,ordered = TRUE))
    
    
    
    
    gg_dat2<-
      stn_all%>%
      mutate(leg_col2 = func_color(above.FALSE,above.TRUE,w.p))%>%
      # mutate(above = above.TRUE<above.FALSE)%>%
      
      filter(!metric %in% c("KGEB","NSEB"))%>%
      filter(mdl %in% labData$mdl)%>%
      
      mutate(
        metricLab = plyr::mapvalues(metric,from = metricLabs1$code,to = metricLabs1$label),
        metricLab = factor(metricLab,
                           levels = metricLabs1$label,ordered = TRUE),#,levels = c("NSE","KGE","NSE_max","NSE_min7","R2_halfFlowDay","r_halfFlowDay","R2_QCI","r_QCI","R2_lowFlowDays","r_lowFlowDays","BE")),
        mdlLab = factor(mdl,labData$mdl,
                        labels = labData$labels,ordered = TRUE))%>%
      
      
      mutate(arrowL = plyr::mapvalues(leg_col2,from = c("Significantly Higher","Higher","Equal","Lower","Significantly Lower"),
                                      to = c(0.04,0.02,0.001,0.02,0.04))%>%
               as.character()%>%
               as.numeric()
      )
    
    unique(gg_dat1$metric)
    
    
    ordered_levels <- names(labData2)
    
    gg_dat1$mdl <- factor(gg_dat1$mdl, levels = ordered_levels)
    gg_dat2$mdl <- factor(gg_dat2$mdl, levels = ordered_levels)
    # Parse character labels to expressions where needed
    label_list <- unname(labData2)
    names(label_list) <- names(labData2)
    
    
    p<-
      ggplot(gg_dat1,
             aes(y = mdl,x= value),size = 0.4)+
      lapply(
        split(gg_dat2,gg_dat2$metricLab),
        function(df)
          geom_segment(data = df,
                       aes(y = mdlLab,yend = mdl,x= above.FALSE,xend = above.TRUE,col = leg_col2),
                       arrow = arrow(length = unit(df$arrowL,  "inches"),
                                     ends = "last", type = "open"))
        
        
      )+
      geom_point(aes(shape = name,col = leg_col),show.legend = c(color = FALSE),size = 0.5)+
      
      facet_wrap("metricLab",scales = "free_x",ncol = 1,strip.position = "right",
                 labeller = label_parsed)+
      
      scale_x_continuous(name = "Median NSE",limits = c(0,1),
                         breaks = c(0,1),
                         oob = scales::oob_keep)+
      scale_y_discrete(name = NULL,
                       limits = rev(ordered_levels),  # Reverse the desired order
                       labels = label_list) +
      scale_color_manual(name = "Significance",
                         breaks = c("Significantly Higher","Higher","Equal","Lower","Significantly Lower"),
                         values = rev(c("#7B3294", "#C2A5CF", "grey", "#A6DBA0", "#008837")))+
      
      scale_shape_manual(name = "Subset",
                         breaks = c("above.FALSE","above.TRUE"),
                         values = c(1,2),
                         labels = groupLabels)+
      theme_bw(base_size = 5)
    
    ggsave(plot = p,filename = flName,width = 3,height = 4,dpi = 600)
    
  }


# calculate variance fractions using STL decomposition
stl_var<-function(dat,s.window,t.window){
  
  
  dat<-filter(dat,!is.na(dat$q))
  # ensure minimum 10 years of data
  if(nrow(dat)<3650){ 
    return(c(
      varSeas = NA,
      varRem =  NA,
      varInterannual = NA
    ))}
  
  dat$t_diff<-c(NA,diff.Date(dat$dt))
  
  # ensure continuous run of at least 10 years of data (no missing days)
  RunLengths<-rle(dat$t_diff)
  if(max(RunLengths$lengths[RunLengths$values==1],na.rm = TRUE)<3650){
    return(c(
      varSeas = NA,
      varRem =  NA,
      varInterannual = NA
    ))
  }
  
  # take longest continuous run of data 
  RunLengths$start<-cumsum(RunLengths$lengths)
  maxRun<-which.max(RunLengths$lengths)
  
  
  dat<-dat[(RunLengths$start[maxRun-1]):
             (RunLengths$start[(maxRun)]),]
  
  
  # for leap years this takes the average of December 30 and 31
  dat_stl<-dat%>%
    group_by(year,yday)%>%
    summarize(q=mean(q))
  
  # convert to timeseries object
  ts_stl<-dat_stl$q%>%
    ts(start = c(dat_stl$year[1],dat_stl$yday[1]),frequency = 365)
  
  # perform stl decomposition
  res <- stl(ts_stl,
             s.window = s.window,
             t.window = t.window)
  
  
  res.df<-data.frame(year = dat_stl$year,
                     yday = dat_stl$yday)%>%
    cbind(res$time.series)
  
  # calculate average seasonality
  seas<-res.df%>%
    group_by(yday)%>%
    summarize(seasonal.avg = mean(seasonal))
  res.df<-left_join(res.df,seas,by = "yday")
  
  # add aanomalies in seasonality to interannual series
  res.df$interannual<-res.df$trend+res.df$seasonal-res.df$seasonal.avg
  
  
  varTotal = var(dat_stl$q)
  
  # return variance fractions
  return(c(
    varSeas = var(res.df$seasonal.avg)/varTotal,
    varRem =  var(res.df$remainder)/varTotal,
    varInterannual = var(res.df$interannual)/varTotal
  ))
  
}

# calculate variance fractions using classical decomposition
clas_var<-function(dat){
  
  dat<-filter(dat,!is.na(dat$q))
  # ensure minimum 10 years of data
  if(nrow(dat)<3650){ 
    return(c(
      varSeas = NA,
      varRem =  NA,
      varInterannual = NA
    ))}
  
  dat<-left_join(data.frame(dt = seq.Date(from  = min(dat$dt),to = max(dat$dt),by = "1 day")),dat)
  # calculate 365-day rolling mean 
  dat$q.trend<-RcppRoll::roll_mean(dat$q,n = 365,align = "center",fill = NA)
  dat<-filter(dat,!is.na(q.trend))
  
  # calculate detrended data
  dat$q.detrend = dat$q-dat$q.trend
  
  # Calcualte mean climatology of detrended data
  dat_clim<-dat%>%
    group_by(yday)%>%
    summarise(q=mean(q),
              q.detrend  = mean(q.detrend),
              N = n())
  
  # if <10 years of data, return NA
  if(nrow(dat_clim)<365|any(dat_clim$N<10)){
    return(c(
      varSeas = NA,
      varRem =  NA,
      varInterannual = NA
    ))
  }
  
  dat<-left_join(dat,dat_clim,by ="yday",suffix = c(".obs",".avg"))
  
  # calculate remainder
  dat$q.rem = dat$q.detrend.obs-dat$q.detrend.avg
  
  
  varTotal = var(dat$q.obs)
  varSeas = var(dat$q.detrend.avg)
  varRem = var(dat$q.rem)
  varTrend = var(dat$q.trend)
  
  
  # return variance fractions
  return(c(
    varSeas = varSeas/varTotal,
    varRem = varRem/varTotal,
    varInterannual = varTrend/varTotal
  ))
  
}

# calculate variance fractions using fourier decomposition
fourier_var<-function(dat){
  
  dat<-filter(dat,!is.na(dat$q))
  # ensure minimum 10 years of data
  if(nrow(dat)<3650){ 
    return(c(
      varSeas = NA,
      varRem =  NA,
      varInterannual = NA
    ))}
  
  
  dat$t_diff<-c(NA,diff.Date(dat$dt))
  # ensure continuous run of at least 10 years of data (no missing days)
  
  RunLengths<-rle(dat$t_diff)
  if(max(RunLengths$lengths[RunLengths$values==1],na.rm = TRUE)<3650){
    return(c(
      varSeas = NA,
      varRem =  NA,
      varInterannual = NA
    ))
  }
  
  
  # take longest continuous run of data 
  RunLengths$start<-cumsum(RunLengths$lengths)
  maxRun<-which.max(RunLengths$lengths)
  dat<-dat[(RunLengths$start[maxRun-1]):
             (RunLengths$start[(maxRun)]),]
  
  # calculate climatological mean of observed data
  
  dat_clim<-dat%>%
    group_by(yday)%>%
    summarise(q=mean(q),
              N = n())
  dat<-left_join(dat,dat_clim,by = "yday",suffix = c("",".avg"))
  
  # calculate anomalies
  dat$q.anomaly<-dat$q-dat$q.avg
  
  # take fast fourier transform of anomalies
  FFT<-fft(dat$q.anomaly)
  
  # Need to centre the frequencies 
  freq <- (0:(nrow(dat)- 1)) / nrow(dat)
  freq <- ifelse(freq > 0.5, freq - 1, freq) * 365  # centered frequency for symmetry
  
  
  # Define cutoff
  cutoff <- 2  # frequency cutoff (2 cycles/year)
  
  # interannual component has frequencies below or at cutoff
  FFT_interannual<-FFT
  FFT_interannual[abs(freq) > cutoff]<-0
  # irregular component has frequencies above cutoff
  FFT_irregular<-FFT
  FFT_irregular[abs(freq) <= cutoff]<-0
  
  # compute inverse FFT of the two components
  dat$q.irregular<-Re(fft(FFT_irregular,inverse = TRUE))/nrow(dat)
  dat$q.interannual<-Re(fft(FFT_interannual,inverse = TRUE))/nrow(dat)
  
  
  # 
  var(dat$q)
  var(dat$q.avg)+  var(dat$q.irregular)+  var(dat$q.interannual)
  
  
  varTotal = var(dat$q)
  varSeas = var(dat$q.avg)
  varRem = var(dat$q.irregular)
  varTrend = var(dat$q.interannual)
  
  # return variance fractions
  return(c(
    varSeas = varSeas/varTotal,
    varRem = varRem/varTotal,
    varInterannual = varTrend/varTotal
  ))
}


# plot the decomposed time series using the fourier method
fourier_plot<-function(dat,include.text = TRUE, round.dec = 2){
  
  dat<-dplyr::filter(dat,!is.na(dat$q))
  dat$t_diff<-c(NA,diff.Date(dat$dt))
  
  RunLengths<-rle(dat$t_diff)
  if(max(RunLengths$lengths[RunLengths$values==1],na.rm = TRUE)<3650){
    return(c(
      varSeas = NA,
      varRem =  NA,
      varInterannual = NA
    ))
  }
  # take longest continuous run of data 
  RunLengths$start<-cumsum(RunLengths$lengths)
  maxRun<-which.max(RunLengths$lengths)
  
  
  dat<-dat[(RunLengths$start[maxRun-1]):
             (RunLengths$start[(maxRun)]),]
  
  # calculate climatological mean of observed data
  dat_clim<-dat%>%
    group_by(yday)%>%
    summarise(seasonal.avg=mean(q),
              N = n())
  dat<-left_join(dat,dat_clim,by = "yday")
  
  # calculate anomalies
  dat$q.anomaly<-dat$q-dat$seasonal.avg
  # Define cutoff
  FFT<-fft(dat$q.anomaly)
  
  # Need to centre the frequencies 
  freq <- (0:(nrow(dat)- 1)) / nrow(dat)
  freq <- ifelse(freq > 0.5, freq - 1, freq) * 365  # centered frequency for symmetry
  
  # Define cutoff
  cutoff <- 2  # frequency cutoff (6 cycles/year)
  # interannual component has frequencies below or at cutoff
  FFT_interannual<-FFT
  FFT_interannual[abs(freq) > cutoff]<-0
  # irregular component has frequencies above cutoff
  FFT_irregular<-FFT
  FFT_irregular[abs(freq) <= cutoff]<-0
  
  
  # compute inverse FFT of the two components
  dat$remainder<-Re(fft(FFT_irregular,inverse = TRUE))/nrow(dat)
  dat$interannual<-Re(fft(FFT_interannual,inverse = TRUE))/nrow(dat)
  
  
  
  # create dataframe with observed data and three components
  res.df.2<-dat%>%
    select(dt,q,interannual,seasonal.avg,remainder)%>%
    pivot_longer(cols = q:remainder)%>%
    mutate(name = factor(name,levels = c("q","interannual","seasonal.avg","remainder"),
                         labels = c("Q~(m^3/s)","Interannual","Seasonal","Irregular")))
  varTotal = var(dat$q)
  varSeas = var(dat$seasonal.avg)
  varRem = var(dat$remainder)
  varInterannual = var(dat$interannual)
  
  # calculate all ranges and ensure each facet has the same scale to enable visual comparison
  y_range1<-range(dat$q)
  y_range2<-range(dat$interannual)
  y_range3<-range(dat$seasonal.avg)
  y_range4<-range(dat$remainder)
  
  y_range_max<-max(c(y_range2[2]-y_range2[1],y_range3[2]-y_range3[1],y_range4[2]-y_range4[1]))
  
  y_range_2_ex<-(y_range_max-diff(y_range2))/2
  y_range_3_ex<-(y_range_max-diff(y_range3))/2
  y_range_4_ex<-(y_range_max-diff(y_range4))/2
  
  
  
  #  labels for plot
  lab_df<-data.frame(
    name = factor(c("q","interannual","seasonal.avg","remainder"),
                  levels = c("q","interannual","seasonal.avg","remainder"),
                  labels = c("Q~(m^3/s)","Interannual","Seasonal","Irregular")),
    lab = c((paste0("sigma^{2}*(Q)==",round(varSeas,round.dec))),
            (paste0("sigma^{2} *(interannual)/sigma^2 *(Q)==",round(varInterannual/varTotal,2))),
            (paste0("sigma^{2} *(seasonal)/sigma^2 *(Q)==",round(varSeas/varTotal,2))),
            (paste0("sigma^{2}*(irregular)/sigma^2 *(Q)==",round(varRem/varTotal,2)))
    ),
    yloc = c(
      y_range1[1]*0.2+y_range1[2]*0.8,
      y_range2[1]*0.2+y_range2[2]*0.8+y_range_2_ex*0.6,
      y_range3[1]*0.2+y_range3[2]*0.8+y_range_3_ex*0.6,
      y_range4[1]*0.2+y_range4[2]*0.8+y_range_4_ex*0.6
    )
  )
  
  
  # list of possible breaks
  pos_brks<-c(0.001,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,300,1000,2000,4000,5000,10000,20000,50000)
  # choose most suitable breaks
  brks<-max(pos_brks[pos_brks<   (brks<-(diff(y_range1)/2))])
  #
  
  p<-ggplot(res.df.2,aes(x = dt,y = value))+
    
    geom_line()+
    
    facet_wrap("name",ncol = 1,strip.position = "left",scales = "free_y",
               labeller = label_parsed)+
    scale_y_continuous(name = NULL)+
    scale_x_date(name = NULL)+
    theme_bw()+
    theme(strip.placement = "outside",legend.position = "none")+
    facetted_pos_scales(
      y = list(
        name == "Q~(m^3/s)" ~ scale_y_continuous(limits = y_range1,name = NULL),
        name == "Interannual" ~ scale_y_continuous(limits = c(y_range2[1]-y_range_2_ex,y_range2[2]+y_range_2_ex),
                                                   breaks = seq(floor(-max(res.df.2$value)/brks)*brks,max(res.df.2$value),brks),name = NULL),
        name == "Seasonal" ~ scale_y_continuous(limits = c(y_range3[1]-y_range_3_ex,y_range3[2]+y_range_3_ex),
                                                breaks = seq(floor(-max(res.df.2$value)/brks)*brks,max(res.df.2$value),brks),name = NULL),
        name == "Irregular" ~ scale_y_continuous(limits = c(y_range4[1]-y_range_4_ex,y_range4[2]+y_range_4_ex),
                                                 breaks = seq(floor(-max(res.df.2$value)/brks)*brks,max(res.df.2$value),brks),name = NULL)
      )
    )
  if(include.text){
    p<-p+ geom_label(data = lab_df,aes(x = max(res.df.2$dt)-(max(res.df.2$dt)-min(res.df.2$dt))/2,y=yloc,label = lab),
                    parse = TRUE,
                    size = 2,
                    col = "grey25")
  }
  
  
  return(p)
  
}


# plot the decomposed time series using the STL method
stl_plot<-function(dat,s.window,t.window,include.text = TRUE){
  
  
  dat<-filter(dat,!is.na(dat$q))
  
  
  dat$t_diff<-c(NA,diff.Date(dat$dt))
  
  RunLengths<-rle(dat$t_diff)
  if(max(RunLengths$lengths[RunLengths$values==1],na.rm = TRUE)<3650){
    return(NULL)
  }
  
  RunLengths$start<-cumsum(RunLengths$lengths)
  maxRun<-which.max(RunLengths$lengths)
  
  
  dat<-dat[(RunLengths$start[maxRun-1]):
             (RunLengths$start[(maxRun)]),]
  
  dat_stl<-dat%>%
    group_by(year,yday)%>%
    summarize(q=mean(q))
  
  
  ts_stl<-dat_stl$q%>%
    ts(start = c(dat_stl$year[1],dat_stl$yday[1]),frequency = 365)
  
  res <- stl(ts_stl,
             s.window = s.window,
             t.window = t.window)
  # plot(res)
  
  res.df<-data.frame(year = dat_stl$year,
                     yday = dat_stl$yday)%>%
    cbind(res$time.series)
  
  seas<-res.df%>%
    group_by(yday)%>%
    summarize(seasonal.avg = mean(seasonal))
  res.df<-left_join(res.df,seas,by = "yday")
  
  res.df$interannual<-res.df$trend+res.df$seasonal-res.df$seasonal.avg
  res.df$q<-dat_stl$q
  
  mean(res.df$seasonal.avg)
  mean(res.df$interannual)
  res.df<-
    res.df%>%
    mutate(dt = ymd(paste(year,"-01-01"))+yday-1)
  res.df.2<-res.df%>%
    select(dt,q,interannual,seasonal.avg,remainder)%>%
    pivot_longer(cols = q:remainder)%>%
    mutate(name = factor(name,levels = c("q","interannual","seasonal.avg","remainder"),
                         labels = c("Q~(m^3/s)","Interannual","Seasonal","Irregular")))
  
  y_range1<-range(res.df$q)
  y_range2<-range(res.df$interannual)
  y_range3<-range(res.df$seasonal.avg)
  y_range4<-range(res.df$remainder)
  
  y_range_max<-max(c(y_range2[2]-y_range2[1],y_range3[2]-y_range3[1],y_range4[2]-y_range4[1]))
  
  y_range_2_ex<-(y_range_max-diff(y_range2))/2
  y_range_3_ex<-(y_range_max-diff(y_range3))/2
  y_range_4_ex<-(y_range_max-diff(y_range4))/2
  
  
  
  #
  
  varTotal = var(dat_stl$q)
  varSeas = var(res.df$seasonal.avg)
  varRem =  var(res.df$remainder)
  varInterannual = var(res.df$interannual)
  
  #  
  lab_df<-data.frame(
    name = factor(c("q","interannual","seasonal.avg","remainder"),
                  levels = c("q","interannual","seasonal.avg","remainder"),
                  labels = c("Q~(m^3/s)","Interannual","Seasonal","Irregular")),
    lab = c((paste0("sigma^{2}*(Q)==",round(varSeas,2))),
            (paste0("sigma^{2} *(interannual)/sigma^2 *(Q)==",round(varInterannual/varTotal,2))),
            (paste0("sigma^{2} *(seasonal)/sigma^2 *(Q)==",round(varSeas/varTotal,2))),
            (paste0("sigma^{2}*(irregular)/sigma^2 *(Q)==",round(varRem/varTotal,2)))
    ),
    yloc = c(
      y_range1[1]*0.2+y_range1[2]*0.8,
      y_range2[1]*0.2+y_range2[2]*0.8+y_range_2_ex*0.6,
      y_range3[1]*0.2+y_range3[2]*0.8+y_range_3_ex*0.6,
      y_range4[1]*0.2+y_range4[2]*0.8+y_range_4_ex*0.6
    )
  )
  
  
  
  pos_brks<-c(0.001,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000,2000,4000,5000,10000,20000,50000)
  brks<-max(pos_brks[pos_brks<   (brks<-(diff(y_range1)/2))])
  #
  
  p<-ggplot(res.df.2,aes(x = dt,y = value))+
    # scale_color_manual(
    #   values = c("black","#FF80F7","#00D1D0","#CFB000"))+
    geom_line()+
    
    facet_wrap("name",ncol = 1,strip.position = "left",scales = "free_y",
               labeller = label_parsed)+
    scale_y_continuous(name = NULL)+
    scale_x_date(name = NULL)+
    theme_bw()+
    theme(strip.placement = "outside",legend.position = "none")+
    facetted_pos_scales(
      y = list(
        name == "Q~(m^3/s)" ~ scale_y_continuous(limits = y_range1,name = NULL),
        name == "Interannual" ~ scale_y_continuous(limits = c(y_range2[1]-y_range_2_ex,y_range2[2]+y_range_2_ex),
                                                   breaks = seq(floor(-max(res.df.2$value)/brks)*brks,max(res.df.2$value),brks),name = NULL),
        name == "Seasonal" ~ scale_y_continuous(limits = c(y_range3[1]-y_range_3_ex,y_range3[2]+y_range_3_ex),
                                                breaks = seq(floor(-max(res.df.2$value)/brks)*brks,max(res.df.2$value),brks),name = NULL),
        name == "Irregular" ~ scale_y_continuous(limits = c(y_range4[1]-y_range_4_ex,y_range4[2]+y_range_4_ex),
                                                 breaks = seq(floor(-max(res.df.2$value)/brks)*brks,max(res.df.2$value),brks),name = NULL)
      )
    )
  if(include.text){
    p<-p+ geom_text(data = lab_df,aes(x = max(res.df.2$dt)-(max(res.df.2$dt)-min(res.df.2$dt))/5,y=yloc,label = lab),
                    parse = TRUE,
                    size = 1)
  }
  
  
  return(p)
  
}

gof_components = 
  function(dat){
    dat<-mutate(dat,
                
                year  =year(dt),
                yday = pmin(yday(dt),365))
    
    
    dat<-filter(dat,!is.na(dat$QObs)& !is.na(dat$QSim))
    if(nrow(dat)<3650){
      return(c(
        NSE = NA,
        NSE.Seas = NA,
        NSE.Interannual = NA,
        NSE.Rem = NA,
        r = NA,
        r.Seas = NA,
        r.Interannual = NA,
        r.Rem = NA,
        
        alpha = NA,
        alpha.Seas = NA,
        alpha.Interannual = NA,
        alpha.Rem = NA
        
      ))
    }
    dat$t_diff<-c(NA,diff.Date(dat$dt))
    
    RunLengths<-rle(dat$t_diff)
    if(max(RunLengths$lengths[RunLengths$values==1],na.rm = TRUE)<3650){
      return(c(
        NSE = NA,
        NSE.Seas = NA,
        NSE.Interannual = NA,
        NSE.Rem = NA,
        r = NA,
        r.Seas = NA,
        r.Interannual = NA,
        r.Rem = NA,
        
        alpha = NA,
        alpha.Seas = NA,
        alpha.Interannual = NA,
        alpha.Rem = NA
        
      ))
    }
    # take longest continuous run of data 
    RunLengths$start<-cumsum(RunLengths$lengths)
    maxRun<-which.max(RunLengths$lengths)
    
    
    dat<-dat[(RunLengths$start[maxRun-1]):
               (RunLengths$start[(maxRun)]),]
    
    # calculate climatological mean of observed data
    dat_clim<-dat%>%
      group_by(yday)%>%
      summarise(QObs.seasonal.avg=mean(QObs),
                QSim.seasonal.avg=mean(QSim),
                N = n())
    dat<-left_join(dat,dat_clim,by = "yday")
    
    # calculate anomalies
    dat$QObs.anomaly<-dat$QObs-dat$QObs.seasonal.avg
    dat$QSim.anomaly<-dat$QSim-dat$QSim.seasonal.avg
    # Define cutoff
    FFT.Obs<-fft(dat$QObs.anomaly)
    FFT.Sim<-fft(dat$QSim.anomaly)
    
    # Need to centre the frequencies 
    freq <- (0:(nrow(dat)- 1)) / nrow(dat)
    freq <- ifelse(freq > 0.5, freq - 1, freq) * 365  # centered frequency for symmetry
    
    # Define cutoff
    cutoff <- 2  # frequency cutoff (6 cycles/year)
    # interannual component has frequencies below or at cutoff
    FFT_interannual.Obs<-FFT.Obs
    FFT_interannual.Obs[abs(freq) > cutoff]<-0
    FFT_interannual.Sim<-FFT.Sim
    FFT_interannual.Sim[abs(freq) > cutoff]<-0
    # irregular component has frequencies above cutoff
    FFT_irregular.Obs<-FFT.Obs
    FFT_irregular.Obs[abs(freq) <= cutoff]<-0
    FFT_irregular.Sim<-FFT.Sim
    FFT_irregular.Sim[abs(freq) <= cutoff]<-0
    
    # compute inverse FFT of the two components
    dat$QObs.remainder<-Re(fft(FFT_irregular.Obs,inverse = TRUE))/nrow(dat)
    dat$QObs.interannual<-Re(fft(FFT_interannual.Obs,inverse = TRUE))/nrow(dat)
    
    dat$QSim.remainder<-Re(fft(FFT_irregular.Sim,inverse = TRUE))/nrow(dat)
    dat$QSim.interannual<-Re(fft(FFT_interannual.Sim,inverse = TRUE))/nrow(dat)
    
    return(
      c(
        NSE = NSE(dat$QSim,dat$QObs),
        NSE.Seas = NSE(dat$QSim.seasonal.avg,dat$QObs.seasonal.avg),
        NSE.Interannual = NSE(dat$QSim.interannual,dat$QObs.interannual),
        NSE.Rem = NSE(dat$QSim.remainder,dat$QObs.remainder),
        
        r = rPearson(dat$QSim,dat$QObs),
        r.Seas = rPearson(dat$QSim.seasonal.avg,dat$QObs.seasonal.avg),
        r.Interannual = rPearson(dat$QSim.interannual,dat$QObs.interannual),
        r.Rem = rPearson(dat$QSim.remainder,dat$QObs.remainder),
        
        alpha = rSD(dat$QSim,dat$QObs),
        alpha.Seas = rSD(dat$QSim.seasonal.avg,dat$QObs.seasonal.avg),
        alpha.Interannual = rSD(dat$QSim.interannual,dat$QObs.interannual),
        alpha.Rem = rSD(dat$QSim.remainder,dat$QObs.remainder)
        
      )
    )
    
  }
