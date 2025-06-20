# climatological benchmarks and variance components

# Author: Sacha Ruzzante
# sachawruzzante@gmail.com
# Last Update: 2025-06-19

# These are utility functions and some labelling list for other scripts. 

require(dplyr)
require(hydroGOF)
require(adc)
require(tidyr)
library(ggh4x)

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
    r_fld = NA
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
    perfMetrics$r_slope_fdc = rSpearmanMin10(dat_yearly$slope_fdc.Sim,dat_yearly$slope_fdc.Obs,minYears = minYears)
    
    
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
    
    # spearman correlation of QCI
    perfMetrics$r_QCI = rSpearmanMin10(QCI$QCI.sim,QCI$QCI.obs,minYears = minYears)
    
    # spearman correlation of annual mean flows
    perfMetrics$r_mag_ann = rSpearmanMin10(dat_yearly$QSim.mean,dat_yearly$QObs.mean,minYears = minYears)
    
    
    
    # spearman correlation of max and min flows
    for(it_roll in c(1,3,7,30,90)){
      perfMetrics[paste0("r_min_",it_roll,"d")]<-rSpearmanMin10(dat_yearly[,paste0("minFlow",it_roll,".Sim")],
                                                                dat_yearly[,paste0("minFlow",it_roll,".Obs")],
                                                                minYears = minYears)
      
      perfMetrics[paste0("r_max_",it_roll,"d")]<-rSpearmanMin10(dat_yearly[,paste0("maxFlow",it_roll,".Sim")],
                                                                dat_yearly[,paste0("maxFlow",it_roll,".Obs")],
                                                                minYears = minYears)
      
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
      N = sum(!is.na(QSim)&!is.na(QObs)))
  
  perf_monthly$r_mag[perf_monthly$N<minYears] = NA
  
  
  
  for(it_month in 1:12){
    perfMetrics[paste0("r_mag_m",it_month)]<-perf_monthly$r_mag[perf_monthly$month==it_month]
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
    
    perfMetrics$r_no_low_pulses = modKendallCor(pulseDat$numLowPulse.Obs,pulseDat$numLowPulse.Sim)
    perfMetrics$r_dur_low_pulses = rSpearmanMin10(pulseDat$lengthLowPulse.Sim,pulseDat$lengthLowPulse.Obs,
                                                  minYears = minYears)
    perfMetrics$r_numDays_low = rSpearmanMin10(pulseDat$totalLowPulse.Sim,pulseDat$totalLowPulse.Obs,
                                               minYears = minYears)
    
    perfMetrics$r_no_high_pulses = modKendallCor(pulseDat$numHighPulse.Obs,pulseDat$numHighPulse.Sim)
    perfMetrics$r_dur_high_pulses = rSpearmanMin10(pulseDat$lengthHighPulse.Sim,pulseDat$lengthHighPulse.Obs,
                                                   minYears = minYears)
    perfMetrics$r_numDays_high = rSpearmanMin10(pulseDat$totalHighPulse.Sim,pulseDat$totalHighPulse.Obs,
                                                minYears = minYears)
    
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
      fld.Sim = meanFall.Sim/sum(diff.Sim<00,na.rm = TRUE),
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
    
    
    perfMetrics$r_mean_pos_diff <-rSpearmanMin10(DiffData$meanRise.Sim,DiffData$meanRise.Obs,
                                                 minYears = minYears)
    perfMetrics$r_mean_neg_diff <-rSpearmanMin10(DiffData$meanFall.Sim,DiffData$meanFall.Obs,
                                                 minYears = minYears)
    
    perfMetrics$r_no_rises <-rSpearmanMin10(DiffData$numRise.Sim,DiffData$numRise.Obs,
                                            minYears = minYears)
    perfMetrics$r_no_falls <-rSpearmanMin10(DiffData$numFall.Sim,DiffData$numFall.Obs,
                                            minYears = minYears)
    
    perfMetrics$r_rld <-rSpearmanMin10(DiffData$rld.Sim,DiffData$rld.Obs,
                                       minYears = minYears)
    perfMetrics$r_fld <-rSpearmanMin10(DiffData$fld.Sim,DiffData$fld.Obs,
                                       minYears = minYears)
    
    
    
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
  
  # Calculate the median metric value for both subsets, for each model
  stn_median<-
    stns%>%

    group_by(mdl,split)%>%
    # calculate median of all columns of the dataframe (from NSE to r_fld)
    reframe(across(NSE:r_fld,~median(.x,na.rm = TRUE))    )%>%
    
    pivot_longer(cols = NSE:r_fld,names_to = "metric",values_to = "value")%>%
    
    pivot_wider(names_from = split,names_prefix = "above.",values_from = value)%>%
    arrange(metric,mdl)
  
  #calculate the number of observations for each subset and model
  stn_N<-
    stns%>%

    pivot_longer(cols = NSE:r_fld,names_to = "metric",values_to = "value")%>%
    group_by(mdl,split,metric)%>%
    
    summarize(N = sum(!is.na(value)))%>%
    
    pivot_wider(names_from = split,names_prefix = "N.above.",values_from = N)
  
  if(any(stns$split==TRUE,na.rm = TRUE)){ # as long as there are some observation above the threshold
    
    # calculate Wilcoxon rank-sum test for difference between two subsets
    stn_wilcox<-stns%>%
      pivot_longer(cols = NSE:r_fld,names_to = "metric",values_to = "value")%>%
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
                             "lstm-camels-br",
                             "lstm-Kraft2025",
                             "lstm-arsenault2022",
                             "lstm-kratzert2024",
                             
                             "nearing-glofas",
                             "prevah-Kraft2025", 
                             "nhm-ref",
                             "q_sim_fuse_900",
                             #"q_sim_fuse_902",
                             #"q_sim_fuse_904",
                             "HBV_ub",
                             "mHm_basin",
                             #"mHm_conus",
                             "SAC_SMA" ,
                             "VIC_basin" 
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

labels = c("GLOB-LSTM",
           "BR-LSTM",
           "CH-LSTM",
           "CA.US-LSTM",
           "US-LSTM",
           
           "GLOB-GLOFAS",
           "CH-PREVAH",
           "US-NHM",
           "US-FUSE",
           #"US-FUSE-902",
           #"US-FUSE-904",
           
           "US-HBV",
           "US-mHM",
           #"US-mHM-regional",
           "US-SAC-SMA",
           
           "US-VIC"
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
compareMetricsHighLow
plotComps<-
  function(stn_all,flName,groupLabels){
    
    
    
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
    
    
    # get rid of blanks in the  metricLabs dataframe
    metricLabs1<-metricLabs%>%filter(!code %in% c("b1","b2","b3","b4","b5","b6"))
    
    
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
      filter()%>%
      filter(metric %in% c("NSE","KGE","Beta_p","Gamma_p","rPearson","KGElf","b1","b2","b3","b4","b5","b6")|
               substr(metric,1,2)=="r_")
    
    
    
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
               substr(metric,1,2)=="r_")%>%
      mutate(arrowL = plyr::mapvalues(leg_col2,from = c("Significantly Higher","Higher","Equal","Lower","Significantly Lower"),
                                      to = c(0.04,0.02,0.001,0.02,0.04))%>%
               as.character()%>%
               as.numeric()
      )
    
    unique(gg_dat1$metric)
    
  
    
    p<-
      ggplot(gg_dat1,
             aes(y = mdlLab,x= value),size = 0.4)+
      lapply(
        split(gg_dat2,gg_dat2$metricLab),
        function(df)
          geom_segment(data = df,
                       aes(y = mdlLab,yend = mdlLab,x= above.FALSE,xend = above.TRUE,col = leg_col2),
                       arrow = arrow(length = unit(df$arrowL,  "inches"),
                                     ends = "last", type = "open"))
        
        
      )+
      geom_point(aes(shape = name,col = leg_col),show.legend = c(color = FALSE),size = 0.5)+
      
      facet_wrap("metricLab",scales = "free_x",ncol = 6,strip.position = "right",
                 labeller = label_parsed)+
      
      scale_x_continuous(name = "Median value of performance metric",limits = c(0,1),
                         breaks = c(0,1),
                         oob = scales::oob_keep)+
      scale_y_discrete(name = NULL,limits=rev,
      )+
      scale_color_manual(name = "Significance",
                         breaks = c("Significantly Higher","Higher","Equal","Lower","Significantly Lower"),
                         values = rev(c("#7B3294", "#C2A5CF", "grey", "#A6DBA0", "#008837")))+
      
      scale_shape_manual(name = "Subset",
                         breaks = c("above.FALSE","above.TRUE"),
                         values = c(1,2),
                         labels = groupLabels)+
      theme_bw(base_size = 5)
    
    ggsave(plot = p,filename = flName,width = 7,height = 8,dpi = 600)
    
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
fourier_plot<-function(dat,include.text = TRUE){
  
  dat<-filter(dat,!is.na(dat$q))
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
  
  
  # list of possible breaks
  pos_brks<-c(0.001,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000,2000,4000,5000,10000,20000,50000)
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
    p<-p+ geom_text(data = lab_df,aes(x = max(res.df.2$dt)-(max(res.df.2$dt)-min(res.df.2$dt))/2,y=yloc,label = lab),
                    parse = TRUE,
                    size = 2)
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



