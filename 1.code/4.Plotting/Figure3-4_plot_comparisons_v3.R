# climatological benchmarks and variance components
# Author: Sacha Ruzzante
# sachawruzzante@gmail.com
# Last Update: 2025-08-06

# This script plots the comparison figures based on 18 hydrologic models. 


library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(sf)

library(scico)
setwd("/home/ruzzante/projects/def-tgleeson/ruzzante/climatological_benchmarks/")

source('1.code/5.utils/utils.R')


# Read in data from 18 models ######
GOF_ls<-list()

#arsenault
GOF_ls[[1]]<-cbind(read.csv("2.data/highLowBenchmarkGOF/GOF_Arsenault2022.csv"),
                   read.csv("2.data/worldclim/indices/Arsenault_indices.csv")[,c("I_m","I_mr","fs")])
#Brazil
GOF_ls[[2]]<-left_join(read.csv("2.data/highLowBenchmarkGOF/GOF_camels-br.csv"),
                       read.csv("2.data/worldclim/indices/camels_br_indices.csv"))%>%
  select(!c(varSeas_clas,varSeas_stl,varSeas_fourier, # I want to substitute the values from the full time series
            varInterannual_clas,varInterannual_stl,varInterannual_fourier,
            varRem_clas,varRem_stl,varRem_fourier))%>%
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
# WC_indices<-read.csv("2.data/worldclim/indices/camels_ch_indices.csv")%>%select(!X)
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
)%>%  select(!c(varSeas_clas,varSeas_stl,varSeas_fourier,
                varInterannual_clas,varInterannual_stl,varInterannual_fourier,
                varRem_clas,varRem_stl,varRem_fourier))%>%
  left_join(varComps_US)

GOF_ls[[8]]<-left_join(read.csv("2.data/highLowBenchmarkGOF/GOF_kratzert_mHm_basin.csv"),
                       inds,
                       by = c("gauge_id" = "ID")
)%>%  select(!c(varSeas_clas,varSeas_stl,varSeas_fourier,
                varInterannual_clas,varInterannual_stl,varInterannual_fourier,
                varRem_clas,varRem_stl,varRem_fourier))%>%
  left_join(varComps_US)
GOF_ls[[9]]<-left_join(read.csv("2.data/highLowBenchmarkGOF/GOF_kratzert_q_sim_fuse_900.csv"),
                       inds,
                       by = c("gauge_id" = "ID")
)%>% select(!c(varSeas_clas,varSeas_stl,varSeas_fourier,
               varInterannual_clas,varInterannual_stl,varInterannual_fourier,
               varRem_clas,varRem_stl,varRem_fourier))%>%
  left_join(varComps_US)

GOF_ls[[10]]<-left_join(read.csv("2.data/highLowBenchmarkGOF/GOF_kratzert_SAC_SMA.csv"),
                        inds,
                        by = c("gauge_id" = "ID")
)%>%  select(!c(varSeas_clas,varSeas_stl,varSeas_fourier,
                varInterannual_clas,varInterannual_stl,varInterannual_fourier,
                varRem_clas,varRem_stl,varRem_fourier))%>%
  left_join(varComps_US)

GOF_ls[[11]]<-left_join(read.csv("2.data/highLowBenchmarkGOF/GOF_kratzert_VIC_basin.csv"),
                        inds,
                        by = c("gauge_id" = "ID")
)%>%  select(!c(varSeas_clas,varSeas_stl,varSeas_fourier,
                varInterannual_clas,varInterannual_stl,varInterannual_fourier,
                varRem_clas,varRem_stl,varRem_fourier))%>%
  left_join(varComps_US)

GOF_ls[[12]]<-left_join(read.csv("2.data/highLowBenchmarkGOF/GOF_kratzert_lstm.csv"),
                        inds,
                        by = c("gauge_id" = "ID")
)%>%  select(!c(varSeas_clas,varSeas_stl,varSeas_fourier,
                varInterannual_clas,varInterannual_stl,varInterannual_fourier,
                varRem_clas,varRem_stl,varRem_fourier))%>%
  left_join(varComps_US)

GOF_ls[[14]]<-left_join(read.csv("2.data/highLowBenchmarkGOF/GOF_lamah_ce_cosero.csv"),
                        read.csv("2.data/worldclim/indices/lamah_ce_indices.csv"),
                        by = c( "ID")
)


GOF_ls[[15]]<-left_join(read.csv("2.data/highLowBenchmarkGOF/GOF_camels-br-MGB-SA.csv"),
                        read.csv("2.data/worldclim/indices/camels_br_indices.csv"),
                        by = c("gauge_id"= "ID")
)




GOF_ls[[16]]<-left_join(read.csv("2.data/highLowBenchmarkGOF/GOF_vic_gl.csv"),
                        read.csv("2.data/worldclim/indices/hysets_indices.csv"),
                        by = c("gauge_id"= "ID")
)


GOF_ls[[17]]<-cbind(read.csv("2.data/highLowBenchmarkGOF/GOF_LSTM-RAPID.csv"),
                    data.frame(I_m=NA,I_mr = NA, fs = NA)
)


GOF_ls[[18]]<-left_join(read.csv("2.data/highLowBenchmarkGOF/GOF_Song_dHBV.csv"),
                        rbind(
                          read.csv("2.data/worldclim/indices/hysets_indices.csv"),
                          read.csv("2.data/worldclim/indices/camels_us_indices.csv")
                        )%>%mutate(gauge_id = as.numeric(ID)),
                        by = c("gauge_id")
)


# Filter to near-natural where posssible #######
# camels-BR filter to more-natural


human_int_br<-left_join(read.delim("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-br/01_CAMELS_BR_attributes/camels_br_human_intervention.txt",
                                   sep = " "),
                        read.delim("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/camels-br/01_CAMELS_BR_attributes/camels_br_land_cover.txt",
                                   sep = " "))

human_int_br$nearNatural<-human_int_br$consumptive_use_perc<5&
  human_int_br$regulation_degree==0&
  human_int_br$imperv_perc<5

GOF_ls[[2]]<-GOF_ls[[2]]%>%
  filter(ID %in% human_int_br$gauge_id[human_int_br$nearNatural])
GOF_ls[[15]]<-GOF_ls[[15]]%>%
  filter(gauge_id %in% human_int_br$gauge_id[human_int_br$nearNatural])
# near-natural IDs for LamaH-CE
natural_ids_lamah_ce<-read.csv("../DATA/1.Spatial_data/global/sw_surfacewater_streamflow_runoff_river_network_waterstress/lamah-ce/C_basins_intermediate_lowimp/1_attributes/Catchment_attributes.csv",
                               sep = ";")

GOF_ls[[14]]<-GOF_ls[[14]]%>%
  filter(ID %in% natural_ids_lamah_ce$ID)

# For Arsenault paper:
gages_II<-read.csv("../DATA/1.Spatial_data/regional/North_America/sw_surfacewater_streamflow_runoff_river_network_waterstress/sw1_GAGES-II/basinchar_and_report_sept_2011/conterm_bas_classif.txt")
RHBN<-readxl::read_excel("../DATA/1.Spatial_data/regional/North_America/sw_surfacewater_streamflow_runoff_river_network_waterstress/RHBN/RHBN_Metadata.xlsx",skip = 2)
RHBN<-RHBN%>%filter(STATION_NUMBER %in% GOF_ls[[1]]$ID)%>%
  filter(!duplicated(STATION_NUMBER))
gages_II<-gages_II%>%
  mutate(STAID = as.character(STAID%>%str_pad(width = 8,pad = "0")))%>%
  filter(STAID %in%  GOF_ls[[1]]$ID,
         CLASS=="Ref")

GOF_ls[[1]]<-GOF_ls[[1]]%>%
  filter(ID %in% c(gages_II$STAID,RHBN$STATION_NUMBER))


gages_II<-read.csv("../DATA/1.Spatial_data/regional/North_America/sw_surfacewater_streamflow_runoff_river_network_waterstress/sw1_GAGES-II/basinchar_and_report_sept_2011/conterm_bas_classif.txt")%>%
  
  filter(CLASS=="Ref")

GOF_ls[[18]]<-filter(GOF_ls[[18]],
                     gauge_id %in% gages_II$STAID)


# Figure 3 ############
#seasonal 0.5
stn_comp_ls<-lapply(GOF_ls,function(dat){compareMetricsHighLow_v2(dat,splitVar = "varSeas_fourier",thresh = 0.5)})
stn_all<-bind_rows(stn_comp_ls)
plotComps_v2(stn_all,"3.figures/figure3_SeasonalVar_0.5_v3.svg",groupLabels = c("SeasVar≤0.5",
                                                                                "SeasVar>0.5"))

# Figure 4 #############
# all indicators
stn_comp_ls<-lapply(GOF_ls,function(dat){compareMetricsHighLow(dat,splitVar = "varSeas_fourier",thresh = 0.5)})
stn_all<-bind_rows(stn_comp_ls)
stn_all%>%filter(
  metric == "NSE"
)%>%
  mutate(N = N.above.FALSE+N.above.TRUE)%>%
  mutate(percAbove = round(100*N.above.TRUE/N))%>%
  select(mdl, N,percAbove)

plotComps(stn_all,"3.figures/figure4_SeasonalVar_0.5_v3.svg",
          groupLabels = 
            c("SeasVar≤0.5",
              "SeasVar>0.5"),
          stat= "r")

plotComps(stn_all,"3.figures/figure4_SeasonalVar_0.5_v3.pdf",
          groupLabels = 
            c("SeasVar≤0.5",
              "SeasVar>0.5"),
          stat= "r")

plotComps(stn_all,"3.figures/figure4_SeasonalVar_0.5_v3_rp.pdf",
          groupLabels = 
            c("SeasVar≤0.5",
              "SeasVar>0.5"),
          stat= "rp")

plotComps(stn_all,"3.figures/figure4_SeasonalVar_0.5_v3_nse.svg",
          groupLabels = 
            c("SeasVar≤0.5",
              "SeasVar>0.5"),
          stat= "NSE",
          minXVal = -100,
          breaks = c(-100,-10,0,1),
          minor_breaks = c(-67,-33,-6.7,-3.3,0.33,0.67),
          transform = "reverse_log10_shifted",
          strip.text.y = 4.25)

plotComps(stn_all,"3.figures/figure4_SeasonalVar_0.5_v3_KGE.pdf",
          groupLabels = 
            c("SeasVar≤0.5",
              "SeasVar>0.5"),
          stat= "KGE",
          minXVal = -1,
          # breaks = c(-100,-10,0,1),
          # transform = "reverse_log10_shifted",
          strip.text.y = 4.25)

plotComps(stn_all,"3.figures/figure4_SeasonalVar_0.5_v3_alpha.pdf",
          groupLabels = 
            c("SeasVar≤0.5",
              "SeasVar>0.5"),
          stat= "alpha",
          minXVal = 0,
          strip.text.y = 4.25)
plotComps(stn_all,"3.figures/figure4_SeasonalVar_0.5_v3_bias.pdf",
          groupLabels = 
            c("SeasVar≤0.5",
              "SeasVar>0.5"),
          stat= "bias",
          minXVal = 0,
          strip.text.y = 4.25)

stn_all$metric_prefix<-stn_all$metric%>%
  str_split_fixed("_",2)%>%
  .[,1]

stn_all%>%
  filter(metric_prefix=="r")%>%
  group_by(mdl)%>%
  summarize(mean(above.FALSE>above.TRUE),
            sum(!is.na(above.FALSE)&!is.na(above.TRUE)))


stn_all%>%
  
  filter(metric_prefix=="r")%>%
  group_by(metric)%>%
  summarize(mean(above.FALSE>above.TRUE),
            sum(!is.na(above.FALSE)&!is.na(above.TRUE)))%>%
  print(n = 100)


stn_all%>%
  
  filter(metric_prefix=="r")%>%
  
  # group_by(metric)%>%
  summarize(mean(above.FALSE>above.TRUE),
            mean(above.FALSE>above.TRUE&w.p<0.05),
            mean(above.FALSE<above.TRUE&w.p<0.05),
            sum(!is.na(above.FALSE)&!is.na(above.TRUE)))

stn_all%>%
  
  filter(metric_prefix=="rp")%>%
  # group_by(metric)%>%
  summarize(mean(above.FALSE>above.TRUE),
            mean(above.FALSE>above.TRUE&w.p<0.05),
            mean(above.FALSE<above.TRUE&w.p<0.05),
            sum(!is.na(above.FALSE)&!is.na(above.TRUE)))
stn_all%>%
  
  filter(metric_prefix=="NSE"& !metric =="NSE")%>%
  # group_by(metric)%>%
  summarize(mean(above.FALSE>above.TRUE),
            mean(above.FALSE>above.TRUE&w.p<0.05),
            mean(above.FALSE<above.TRUE&w.p<0.05),
            sum(!is.na(above.FALSE)&!is.na(above.TRUE)))


stn_all%>%
  
  filter(metric_prefix=="NSE"& !metric =="NSE")%>%
  group_by(metric)%>%
  summarize(sum(above.FALSE<0&above.TRUE<0),
            sum(!is.na(above.FALSE)&!is.na(above.TRUE)))%>%
  print(n=100)

stn_all%>%
  
  filter(metric_prefix=="NSE"& !metric =="NSE")%>%
  # group_by(metric)%>%
  summarize(sum(above.FALSE<0&above.TRUE<0),
            sum(!is.na(above.FALSE)&!is.na(above.TRUE)))%>%
  print(n=100)

stn_all%>%
  
  filter(metric_prefix=="KGE"& !metric =="KGE")%>%
  # group_by(metric)%>%
  summarize(mean(above.FALSE>above.TRUE),
            mean(above.FALSE>above.TRUE&w.p<0.05),
            mean(above.FALSE<above.TRUE&w.p<0.05),
            sum(!is.na(above.FALSE)&!is.na(above.TRUE)))


stn_all%>%
  
  filter(metric_prefix=="alpha")%>%
  # group_by(metric)%>%
  summarize(mean(above.FALSE>above.TRUE),
            mean(above.FALSE>above.TRUE&w.p<0.05),
            mean(above.FALSE<above.TRUE&w.p<0.05),
            sum(!is.na(above.FALSE)&!is.na(above.TRUE)))

stn_all%>%
  
  filter(metric_prefix=="bias")%>%
  # group_by(metric)%>%
  summarize(mean(above.FALSE>above.TRUE),
            mean(above.FALSE>above.TRUE&w.p<0.05),
            mean(above.FALSE<above.TRUE&w.p<0.05),
            sum(!is.na(above.FALSE)&!is.na(above.TRUE)))

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
  # filter(!mdl == 'vic-gl')%>%
  # group_by(metric)%>%
  summarize(mean(above.FALSE>above.TRUE),
            mean(above.FALSE>above.TRUE&w.p<0.05),
            mean(above.FALSE<above.TRUE&w.p<0.05),
            sum(!is.na(above.FALSE)&!is.na(above.TRUE)))



stn_comp_ls<-lapply(GOF_ls,function(dat){compareMetricsHighLow_v3(dat,splitVar = "varSeas_fourier",thresh = 0.5)})
stn_all<-bind_rows(stn_comp_ls)

plotComps(stn_all,"3.figures/figure4_SeasonalVar_0.5_v3.pdf",groupLabels = c("SeasVar≤0.5",
                                                                             "SeasVar>0.5"))

## Sensitivity ##########
#seasonal 0.4
stn_comp_ls<-lapply(GOF_ls,function(dat){compareMetricsHighLow_v2(dat,splitVar = "varSeas_fourier",thresh = 0.4)})
stn_all<-bind_rows(stn_comp_ls)
plotComps_v2(stn_all,"3.figures/figure3_SeasonalVar_0.4_v3.png",groupLabels = c("SeasVar≤0.4",
                                                                                "SeasVar>0.4"))

stn_comp_ls<-lapply(GOF_ls,function(dat){compareMetricsHighLow(dat,splitVar = "varSeas_fourier",thresh = 0.4)})
stn_all<-bind_rows(stn_comp_ls)
plotComps(stn_all,"3.figures/figure4_SeasonalVar_0.4_v3.png",groupLabels = c("SeasVar≤0.4",
                                                                             "SeasVar>0.4"))
plotComps(stn_all,"3.figures/figure4_SeasonalVar_0.4_v3.svg",groupLabels = c("SeasVar≤0.4",
                                                                             "SeasVar>0.4"))

#seasonal 0.6
stn_comp_ls<-lapply(GOF_ls,function(dat){compareMetricsHighLow_v2(dat,splitVar = "varSeas_fourier",thresh = 0.6)})
stn_all<-bind_rows(stn_comp_ls)
plotComps_v2(stn_all,"3.figures/figure3_SeasonalVar_0.6_v3.svg",groupLabels = c("SeasVar≤0.6",
                                                                                "SeasVar>0.6"))


stn_comp_ls<-lapply(GOF_ls,function(dat){compareMetricsHighLow(dat,splitVar = "varSeas_fourier",thresh = 0.6)})
stn_all<-bind_rows(stn_comp_ls)
plotComps(stn_all,"3.figures/figure4_SeasonalVar_0.6_v3.svg",groupLabels = c("SeasVar≤0.6",
                                                                             "SeasVar>0.6"))


#interannual
stn_comp_ls<-lapply(GOF_ls,function(dat){compareMetricsHighLow_v2(dat,splitVar = "varInterannual_fourier",thresh = 0.1)})
stn_all<-bind_rows(stn_comp_ls)
plotComps_v2(stn_all,"3.figures/figure3_InterannualVar_0.1_v2.svg",groupLabels = c("InterannualVar≤  0.1",
                                                                                   "InterannualVar> 0.1"))
stn_comp_ls<-lapply(GOF_ls,function(dat){compareMetricsHighLow(dat,splitVar = "varInterannual_fourier",thresh = 0.1)})
stn_all<-bind_rows(stn_comp_ls)
plotComps(stn_all,"3.figures/figure4_InterannualVar_0.1_v2.svg",groupLabels = c("InterannualVar≤  0.1",
                                                                                "InterannualVar> 0.1"))
#interannual
stn_comp_ls<-lapply(GOF_ls,function(dat){compareMetricsHighLow_v2(dat,splitVar = "varInterannual_fourier",thresh = 0.2)})
stn_all<-bind_rows(stn_comp_ls)
plotComps_v2(stn_all,"3.figures/figure3_InterannualVar_0.2_v2.pdf",groupLabels = c("InterannualVar≤  0.2",
                                                                                   "InterannualVar> 0.2"))
#interannual
stn_comp_ls<-lapply(GOF_ls,function(dat){compareMetricsHighLow_v2(dat,splitVar = "varInterannual_fourier",thresh = 0.3)})
stn_all<-bind_rows(stn_comp_ls)
plotComps_v2(stn_all,"3.figures/figure3_InterannualVar_0.3_v2.pdf",groupLabels = c("InterannualVar≤  0.3",
                                                                                   "InterannualVar> 0.3"))




#irregular
stn_comp_ls<-lapply(GOF_ls,function(dat){compareMetricsHighLow_v2(dat,splitVar = "varRem_fourier",thresh = 0.5)})
stn_all<-bind_rows(stn_comp_ls)

plotComps_v2(stn_all,"3.figures/figure3_RemVar_0.5_v2.pdf",groupLabels = c("IrregularVar≤              0.5",
                                                                           "IrregularVar>      0.5"))



# CoV
GOF_ls[[3]]%>%
  
  ggplot(aes(x = CoV))+geom_histogram()

stn_comp_ls<-lapply(GOF_ls,function(dat){compareMetricsHighLow(dat,splitVar = "CoV",thresh = 1)})
stn_all<-bind_rows(stn_comp_ls)
plotComps(stn_all,"3.figures/figure4_CoV_1.svg",groupLabels = c("SeasVar≤0.6",
                                                                "CoV > 1"))

stn_comp_ls<-lapply(GOF_ls,function(dat){compareMetricsHighLow(dat,splitVar = "CoV",thresh = 0.75)})
stn_all<-bind_rows(stn_comp_ls)
plotComps(stn_all,"3.figures/figure4_CoV_0.75.svg",groupLabels = c("SeasVar≤0.6",
                                                                   "CoV > 0.75"))
# I_mr
GOF_ls[[3]]%>%
  
  ggplot(aes(x = I_mr))+geom_histogram()

stn_comp_ls<-lapply(GOF_ls[c(1:15,17:18)],function(dat){compareMetricsHighLow(dat,splitVar = "I_mr",thresh = 1)})
stn_all<-bind_rows(stn_comp_ls)
plotComps(stn_all,"3.figures/figure4_I_mr_1.svg",groupLabels = c("SeasVar≤0.6",
                                                                 "I_mr > 1"))

stn_comp_ls<-lapply(GOF_ls,function(dat){compareMetricsHighLow(dat,splitVar = "I_mr",thresh = 1.5)})
stn_all<-bind_rows(stn_comp_ls)%>%
  filter(!is.na(w.p))%>%
  select(!c(above.NA,N.above.NA))
plotComps(stn_all,"3.figures/figure4_I_mr_1.5.pdf",groupLabels = c("Benchmark NSE≤9.9",
                                                                   "I_mr > 1.5"))
# fs
GOF_ls[[3]]%>%
  
  ggplot(aes(x = fs))+geom_histogram()

stn_comp_ls<-lapply(GOF_ls,function(dat){compareMetricsHighLow(dat,splitVar = "fs",thresh = 0.5)})
stn_all<-bind_rows(stn_comp_ls)%>%
  filter(!is.na(w.p))%>%
  select(!c(above.NA,N.above.NA))
plotComps(stn_all,"3.figures/figure4_fs_0.5.svg",groupLabels = c("SeasVar≤0.6",
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
plotComps(stn_all,"3.figures/figure4_QCI_15.svg",groupLabels = c("SeasVar≤0.6",
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



# the overall NSE is the weighted sum of the component NSEs ############


# reload the data where I subbed in the variance components from the full time series

#Brazil
GOF_ls[[2]]<-left_join(read.csv("2.data/highLowBenchmarkGOF/GOF_camels-br.csv"),
                       read.csv("2.data/worldclim/indices/camels_br_indices.csv"))

GOF_ls[[2]]<-GOF_ls[[2]]%>%
  filter(ID %in% human_int_br$gauge_id[human_int_br$nearNatural])

GOF_ls[[7]]<-left_join(read.csv("2.data/highLowBenchmarkGOF/GOF_kratzert_HBV_ub.csv"),
                       inds,
                       by = c("gauge_id" = "ID"))
GOF_ls[[8]]<-left_join(read.csv("2.data/highLowBenchmarkGOF/GOF_kratzert_mHm_basin.csv"),
                       inds,
                       by = c("gauge_id" = "ID"))
GOF_ls[[9]]<-left_join(read.csv("2.data/highLowBenchmarkGOF/GOF_kratzert_q_sim_fuse_900.csv"),
                       inds,
                       by = c("gauge_id" = "ID"))
GOF_ls[[10]]<-left_join(read.csv("2.data/highLowBenchmarkGOF/GOF_kratzert_SAC_SMA.csv"),
                        inds,
                        by = c("gauge_id" = "ID"))
GOF_ls[[11]]<-left_join(read.csv("2.data/highLowBenchmarkGOF/GOF_kratzert_VIC_basin.csv"),
                        inds,
                        by = c("gauge_id" = "ID"))
GOF_ls[[12]]<-left_join(read.csv("2.data/highLowBenchmarkGOF/GOF_kratzert_lstm.csv"),
                        inds,
                        by = c("gauge_id" = "ID"))

x<-
  lapply(GOF_ls, function(dat){
  dat%>%mutate(
  NSE.sumComps = varSeas_fourier*NSE.Seas+varInterannual_fourier*NSE.Interannual+varRem_fourier*NSE.Rem
)%>%
      select(mdl, varSeas_fourier,NSE.Seas,varInterannual_fourier,NSE.Interannual,varRem_fourier,NSE.Rem,
             NSE.maxRun,NSE.sumComps)
      })%>%
  bind_rows()%>%
  filter(!is.na(NSE.sumComps))%>%
  mutate(NSE.diff = NSE.maxRun-NSE.sumComps)


quantile((x$NSE.diff),c(0,0.95,0.99,0.999,1))
quantile(abs(x$NSE.diff),c(0,0.95,0.99,0.999,1))
quantile(abs(x$NSE.diff/x$NSE.maxRun),c(0.95,0.99,0.999,1))

# Does the NSE Correlate with each of the metrics? ############
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

