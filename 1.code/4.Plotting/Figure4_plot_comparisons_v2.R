# climatological benchmarks and variance components
# Author: Sacha Ruzzante
# sachawruzzante@gmail.com
# Last Update: 2025-06-19

# This script plots the variance components (Figure 1 ).


library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(sf)

library(scico)
setwd("/home/ruzzante/projects/def-tgleeson/ruzzante/climatological_benchmarks/")

source('1.code/5.utils/utils.R')


# Read in data from 13 models


GOF_ls<-list()
#arsenault
GOF_ls[[1]]<-left_join(read.csv("2.data/highLowBenchmarkGOF/GOF_Arsenault2022.csv"),
                       read.csv("2.data/worldclim/indices/Arsenault_indices.csv"))
#Brazil
GOF_ls[[2]]<-left_join(read.csv("2.data/highLowBenchmarkGOF/GOF_camels-br.csv"),
                       read.csv("2.data/worldclim/indices/camels_br_indices.csv"))%>%
  left_join(read.csv("2.data/varComponents/camels_br.csv"),
            by =  c("ID" = "gauge_id"))

#Kraft
#lstm
GOF<-read.csv("2.data/highLowBenchmarkGOF/GOF_Kraft_lstm.csv")%>%select(!X)
WC_indices<-read.csv("2.data/worldclim/indices/camels_ch_indices.csv")
CH_crossTab<-read.csv("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/Kraft2025/crossTab.csv")
WC_indices$ID<-plyr::mapvalues(WC_indices$ID,CH_crossTab$camels_catchment,CH_crossTab$ID)
GOF_ls[[3]]<-left_join(GOF,WC_indices)

# prevah
GOF<-read.csv("2.data/highLowBenchmarkGOF/GOF_Kraft_Prevah.csv")%>%select(!X)
GOF_ls[[4]]<-left_join(GOF,WC_indices)

# Nearing - lstm
GOF_ls[[5]]<-left_join(read.csv("2.data/highLowBenchmarkGOF/GOF_Nearing_lstm.csv"),
                       read.csv("2.data/worldclim/indices/grdc_indices.csv"),by = c("gauge_ID" = "ID"))

GOF_ls[[13]]<-left_join(read.csv("2.data/highLowBenchmarkGOF/GOF_Nearing_glofas.csv"),
                        read.csv("2.data/worldclim/indices/grdc_indices.csv"),by = c("gauge_ID" = "ID"))


#NHM
GOF_ls[[6]]<-left_join(read.csv("2.data/highLowBenchmarkGOF/GOF_NHM_ref.csv"),
                       rbind(read.csv("2.data/worldclim/indices/camels_us_indices.csv"),
                             read.csv("2.data/worldclim/indices/hysets_indices.csv"))
)

# Kratzert process-based models
# US Indices
inds<-rbind(read.csv("2.data/worldclim/indices/camels_us_indices.csv"),
            read.csv("2.data/worldclim/indices/hysets_indices.csv"))%>%mutate(ID = as.numeric(ID))%>%
  filter(!is.na(ID))

varComps_US<-rbind(read.csv("2.data/varComponents/camels_US.csv"),
                   read.csv("2.data/varComponents/hysets.csv")%>%rename(gauge_id = Official_ID))%>%
  mutate(gauge_id = as.numeric(gauge_id))%>%
  filter(!is.na(gauge_id))

GOF_ls[[7]]<-left_join(read.csv("2.data/highLowBenchmarkGOF/GOF_kratzert_HBV_ub.csv"),
                       inds,
                       by = c("gauge_id" = "ID")
)%>%
  left_join(varComps_US)

GOF_ls[[8]]<-left_join(read.csv("2.data/highLowBenchmarkGOF/GOF_kratzert_mHm_basin.csv"),
                       inds,
                       by = c("gauge_id" = "ID")
)%>%
  left_join(varComps_US)
GOF_ls[[9]]<-left_join(read.csv("2.data/highLowBenchmarkGOF/GOF_kratzert_q_sim_fuse_900.csv"),
                       inds,
                       by = c("gauge_id" = "ID")
)%>%
  left_join(varComps_US)
GOF_ls[[10]]<-left_join(read.csv("2.data/highLowBenchmarkGOF/GOF_kratzert_SAC_SMA.csv"),
                        inds,
                        by = c("gauge_id" = "ID")
)%>%
  left_join(varComps_US)

GOF_ls[[11]]<-left_join(read.csv("2.data/highLowBenchmarkGOF/GOF_kratzert_VIC_basin.csv"),
                        inds,
                        by = c("gauge_id" = "ID")
)%>%
  left_join(varComps_US)

GOF_ls[[12]]<-left_join(read.csv("2.data/highLowBenchmarkGOF/GOF_kratzert_lstm.csv"),
                        inds,
                        by = c("gauge_id" = "ID")
)%>%
  left_join(varComps_US)



#seasonal
stn_comp_ls<-lapply(GOF_ls,function(dat){compareMetricsHighLow(dat,splitVar = "varSeas_fourier",thresh = 0.5)})
stn_all<-bind_rows(stn_comp_ls)
plotComps(stn_all,"3.figures/figure3_SeasonalVar_0.5.pdf",groupLabels = c("SeasVar≤              0.5",
                                                                          "SeasVar>      0.5"))

# NSEB 0.5
stn_comp_ls<-lapply(GOF_ls,function(dat){compareMetricsHighLow(dat,splitVar = "NSEB",thresh = 0.5)})
stn_all<-bind_rows(stn_comp_ls)
plotComps(stn_all,"3.figures/figure3_NSEB_0.5.svg",groupLabels = 
            c("Benchmark NSE≤0.5",
              "Benchmark NSE>0.5"))

stn_all%>%
  filter(!metric %in% c("NSE","KGE","BE","Beta","Beta_p","Gamma","Gamma_p","NSEB","KGEB","rPearson","KGElf"))%>%
  group_by(mdl)%>%
  summarize(mean(above.FALSE>above.TRUE),
            sum(!is.na(above.FALSE)&!is.na(above.TRUE)))

stn_all%>%
  filter(!metric %in% c("NSE","KGE","BE","Beta","Beta_p","Gamma","Gamma_p","NSEB","KGEB","rPearson","KGElf"))%>%
  group_by(metric)%>%
  summarize(mean(above.FALSE>above.TRUE),
            sum(!is.na(above.FALSE)&!is.na(above.TRUE)))%>%
  print(n = 100)


stn_all%>%
  filter(!metric %in% c("NSE","KGE","BE","Beta","Beta_p","Gamma","Gamma_p","NSEB","KGEB","rPearson","KGElf"))%>%
  # group_by(metric)%>%
  summarize(mean(above.FALSE>above.TRUE),
            mean(above.FALSE>above.TRUE&w.p<0.05),
            sum(!is.na(above.FALSE)&!is.na(above.TRUE)))



# NSEB 0.4
stn_comp_ls<-lapply(GOF_ls,function(dat){compareMetricsHighLow(dat,splitVar = "NSEB",thresh = 0.4)})
stn_all<-bind_rows(stn_comp_ls)
plotComps(stn_all,"3.figures/figure3_NSEB_0.4.svg",groupLabels = c("Benchmark NSE≤0.4",
                                                                   "Benchmark NSE>0.5"))
# NSEB 0.6
stn_comp_ls<-lapply(GOF_ls,function(dat){compareMetricsHighLow(dat,splitVar = "NSEB",thresh = 0.6)})
stn_all<-bind_rows(stn_comp_ls)
plotComps(stn_all,"3.figures/figure3_NSEB_0.6.svg",groupLabels = c("Benchmark NSE≤0.6",
                                                                   "Benchmark NSE>0.6"))

# CoV
GOF_ls[[3]]%>%
  
  ggplot(aes(x = CoV))+geom_histogram()

stn_comp_ls<-lapply(GOF_ls,function(dat){compareMetricsHighLow(dat,splitVar = "CoV",thresh = 1)})
stn_all<-bind_rows(stn_comp_ls)
plotComps(stn_all,"3.figures/figure3_CoV_1.svg",groupLabels = c("Benchmark NSE≤9.9",
                                                                "CoV > 1"))

stn_comp_ls<-lapply(GOF_ls,function(dat){compareMetricsHighLow(dat,splitVar = "CoV",thresh = 0.75)})
stn_all<-bind_rows(stn_comp_ls)
plotComps(stn_all,"3.figures/figure3_CoV_0.75.svg",groupLabels = c("Benchmark NSE≤9.9",
                                                                   "CoV > 0.75"))
# I_mr
GOF_ls[[3]]%>%
  
  ggplot(aes(x = I_mr))+geom_histogram()

stn_comp_ls<-lapply(GOF_ls,function(dat){compareMetricsHighLow(dat,splitVar = "I_mr",thresh = 1)})
stn_all<-bind_rows(stn_comp_ls)
plotComps(stn_all,"3.figures/figure3_I_mr_1.svg",groupLabels = c("Benchmark NSE≤9.9",
                                                                 "I_mr > 1"))

stn_comp_ls<-lapply(GOF_ls,function(dat){compareMetricsHighLow(dat,splitVar = "I_mr",thresh = 1.5)})
stn_all<-bind_rows(stn_comp_ls)%>%
  filter(!is.na(w.p))%>%
  select(!c(above.NA,N.above.NA))
plotComps(stn_all,"3.figures/figure3_I_mr_1.5.svg",groupLabels = c("Benchmark NSE≤9.9",
                                                                   "I_mr > 1.5"))
# fs
GOF_ls[[3]]%>%
  
  ggplot(aes(x = fs))+geom_histogram()

stn_comp_ls<-lapply(GOF_ls,function(dat){compareMetricsHighLow(dat,splitVar = "fs",thresh = 0.5)})
stn_all<-bind_rows(stn_comp_ls)%>%
  filter(!is.na(w.p))%>%
  select(!c(above.NA,N.above.NA))
plotComps(stn_all,"3.figures/figure3_fs_0.5.svg",groupLabels = c("Benchmark NSE≤9.9",
                                                                 "sf > 0.5"))

# QCI
GOF_ls[[6]]%>%
  
  ggplot(aes(x = QCI))+geom_histogram()

stn_comp_ls<-lapply(GOF_ls,function(dat){compareMetricsHighLow(dat,splitVar = "QCI",thresh = 20)})
stn_all<-bind_rows(stn_comp_ls)%>%
  filter(!is.na(w.p))%>%
  select(!c(above.NA,N.above.NA))
plotComps(stn_all,"3.figures/figure3_QCI_20.svg",groupLabels = c("Benchmark NSE≤9.9",
                                                                 "QCI > 20"))
GOF_ls[[6]]%>%
  
  ggplot(aes(x = QCI))+geom_histogram()

stn_comp_ls<-lapply(GOF_ls,function(dat){compareMetricsHighLow(dat,splitVar = "QCI",thresh = 15)})
stn_all<-bind_rows(stn_comp_ls)%>%
  filter(!is.na(w.p))%>%
  select(!c(above.NA,N.above.NA))
plotComps(stn_all,"3.figures/figure3_QCI_15.svg",groupLabels = c("Benchmark NSE≤9.9",
                                                                 "QCI > 15"))

stn_comp_ls<-lapply(GOF_ls,function(dat){compareMetricsHighLow(dat,splitVar = "QCI",thresh = 30)})
stn_all<-bind_rows(stn_comp_ls)%>%
  filter(!is.na(w.p))%>%
  select(!c(above.NA,N.above.NA))
plotComps(stn_all,"3.figures/figure3_QCI_30.svg",groupLabels = c("Benchmark NSE≤9.9",
                                                                 "QCI > 30"))



# interannual variance

# Interannual>0.2
stn_comp_ls<-lapply(GOF_ls,function(dat){compareMetricsHighLow(dat,splitVar = "varInterannual_fourier",thresh = 0.2)})
stn_all<-bind_rows(stn_comp_ls)
plotComps(stn_all,"3.figures/figure3_InterannualVar_0.2.pdf",groupLabels = c("InterAnnVar≤0.2",
                                                                             "InterAnnVar>0.2"))
# Interannual>0.2
stn_comp_ls<-lapply(GOF_ls,function(dat){compareMetricsHighLow(dat,splitVar = "varInterannual_fourier",thresh = 0.15)})
stn_all<-bind_rows(stn_comp_ls)
plotComps(stn_all,"3.figures/figure3_InterannualVar_0.15.pdf",groupLabels = c("InterAnnVar≤0.15",
                                                                              "InterAnnVar>0.15"))

# Interannual>0.1
stn_comp_ls<-lapply(GOF_ls,function(dat){compareMetricsHighLow(dat,splitVar = "varInterannual_fourier",thresh = 0.1)})
stn_all<-bind_rows(stn_comp_ls)
plotComps(stn_all,"3.figures/figure3_InterannualVar_0.1.pdf",groupLabels = c("InterAnnVar≤0.1",
                                                                             "InterAnnVar>0.1"))
# Interannual>0.1
stn_comp_ls<-lapply(GOF_ls,function(dat){compareMetricsHighLow(dat,splitVar = "varInterannual_clas",thresh = 0.1)})
stn_all<-bind_rows(stn_comp_ls)
plotComps(stn_all,"3.figures/figure3_InterannualVarClassical_0.1.pdf",groupLabels = c("InterAnnVar≤0.1",
                                                                                      "InterAnnVar>0.1"))



stn_all%>%
  filter(!metric %in% c("NSE","KGE","BE","Beta","Beta_p","Gamma","Gamma_p","NSEB","KGEB","rPearson","KGElf"))%>%
  group_by(mdl)%>%
  summarize(mean(above.FALSE>above.TRUE),
            sum(!is.na(above.FALSE)&!is.na(above.TRUE)))

stn_all%>%
  filter(!metric %in% c("NSE","KGE","BE","Beta","Beta_p","Gamma","Gamma_p","NSEB","KGEB","rPearson","KGElf"))%>%
  group_by(metric)%>%
  summarize(mean(above.FALSE>above.TRUE),
            sum(!is.na(above.FALSE)&!is.na(above.TRUE)))%>%
  print(n = 100)



stn_all%>%
  filter(!metric %in% c("NSE","KGE","BE","Beta","Beta_p","Gamma","Gamma_p","NSEB","KGEB","rPearson","KGElf"))%>%
  # group_by(metric)%>%
  summarize(mean(above.FALSE>above.TRUE),
            mean(above.FALSE>above.TRUE&w.p<0.05),
            sum(!is.na(above.FALSE)&!is.na(above.TRUE)))





# Does the NSE Correlate with each of the metrics?
func_cors<-function(GOF){
  cors<-(GOF%>%select(any_of(metricLabs$code))%>%
           cor(use = "pairwise.complete",method = "spearman"))[,"NSE"]
  cors_df<-data.frame(cors)
  cors_df$var = rownames(cors_df)
  rownames(cors_df)<-NULL
  return(cors_df)
}


GOF<-GOF_ls[[2]]

dat<-list(GOF%>%
            func_cors()%>%
            mutate(sbst = "all"),
          
          GOF%>%
            subset(NSEB>0.5)%>%
            func_cors()%>%
            mutate(sbst = "NSEB>0.5"),
          GOF%>%
            subset(NSEB<=0.5)%>%
            func_cors()%>%
            mutate(sbst = "NSEB<=0.5")
)%>%
  bind_rows()

metricLabs$label_parsed <- lapply(metricLabs$label, function(lbl) parse(text = lbl))

dat$var_lbl<-plyr::mapvalues(dat$var,from = metricLabs$code,to = metricLabs$label)

ggplot(dat,aes(x = sbst,y =var_lbl,fill = cors ))+
  geom_tile()+
  geom_text(aes(label = round(cors,1)))+
  scale_y_discrete(limits = rev(metricLabs$label),
                   labels = scales::label_parse())+
  scale_fill_scico(palette = "roma",midpoint =0)+
  theme_bw()

