### Folder
setwd("/Users/nicholaspanchy/Documents/Work_UTK/DoseTime_FromScratch2/")

### Process All Data
Integration1_AllClusterMarkers_Quart <- readRDS("AllMarkers_pv05_DropFiltered.rds")
dim(Integration1_AllClusterMarkers_Quart)
Integration1_AllClusterMarkers_Quart$rawFC <- 2^Integration1_AllClusterMarkers_Quart$avg_log2FC

RawFC_3rdquart <- quantile(Integration1_AllClusterMarkers_Quart$rawFC,0.75)
RawFC_1stquart <- quantile(Integration1_AllClusterMarkers_Quart$rawFC,0.25)

### 3rd Quarter ###

Integration1_AllClusterMarkers_Quart_FC3rd <- Integration1_AllClusterMarkers_Quart[Integration1_AllClusterMarkers_Quart$rawFC > RawFC_3rdquart,]
Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05 <- Integration1_AllClusterMarkers_Quart_FC3rd[Integration1_AllClusterMarkers_Quart_FC3rd$p_val_adj2 < 0.05,] 
dim(Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05)

Cluster1_Markers <- Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05[Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05$cluster == "1",]
Cluster2_Markers <- Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05[Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05$cluster == "2",]
Cluster3_Markers <- Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05[Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05$cluster == "3",]
Cluster4_Markers <- Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05[Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05$cluster == "4",]
Cluster5_Markers <- Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05[Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05$cluster == "5",]
Cluster6_Markers <- Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05[Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05$cluster == "6",]
Cluster7_Markers <- Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05[Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05$cluster == "7",]
Cluster8_Markers <- Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05[Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05$cluster == "8",]
Cluster9_Markers <- Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05[Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05$cluster == "9",]
Cluster10_Markers <- Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05[Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05$cluster == "10",]
Cluster11_Markers <- Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05[Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05$cluster == "11",]
Cluster12_Markers <- Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05[Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05$cluster == "12",]
Cluster13_Markers <- Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05[Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05$cluster == "13",]
Cluster14_Markers <- Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05[Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05$cluster == "14",]
Cluster15_Markers <- Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05[Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05$cluster == "15",]
Cluster16_Markers <- Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05[Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05$cluster == "16",]

AllGenes <- Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05$gene
length(AllGenes)

Clusters_Integration_Quart_FC3rd <- list(Cluster1_Markers,Cluster2_Markers,
                                         Cluster3_Markers,Cluster4_Markers,Cluster5_Markers,
                                         Cluster6_Markers,Cluster7_Markers,Cluster8_Markers,
                                         Cluster9_Markers,Cluster10_Markers,Cluster11_Markers,
                                         Cluster12_Markers,Cluster13_Markers,Cluster14_Markers,
                                         Cluster15_Markers,Cluster16_Markers)



for (df in Clusters_Integration_Quart_FC3rd){
  print(dim(df))
}

# check
table(Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05$cluster)
sum(table(Integration1_AllClusterMarkers_Quart_FC3rd_adjPV05$cluster))

saveRDS(Clusters_Integration_Quart_FC3rd,'Clusters_IntegrationFirst_Quart_3rdQuartFC_Alt.rds')

###############################################
### All Clusters --- Integration Time Only  ###
###############################################

Integration1_AllClusterMarkers_QuartTime <- readRDS("TimeMarkers_pv05_DropFiltered.rds")
dim(Integration1_AllClusterMarkers_QuartTime)
Integration1_AllClusterMarkers_QuartTime$rawFC <- 2^Integration1_AllClusterMarkers_QuartTime$avg_log2FC

RawFC_3rdquart <- quantile(Integration1_AllClusterMarkers_QuartTime$rawFC,0.75)

### 3rd Quarter ###

Integration1_AllClusterMarkers_QuartTime_FC3rd <- Integration1_AllClusterMarkers_QuartTime[Integration1_AllClusterMarkers_QuartTime$rawFC > RawFC_3rdquart,]
Integration1_AllClusterMarkers_QuartTime_FC3rd_adjPV05 <- Integration1_AllClusterMarkers_QuartTime_FC3rd[Integration1_AllClusterMarkers_QuartTime_FC3rd$p_val_adj2 < 0.05,]
dim(Integration1_AllClusterMarkers_QuartTime_FC3rd_adjPV05)

Cluster1_Markers <- Integration1_AllClusterMarkers_QuartTime_FC3rd_adjPV05[Integration1_AllClusterMarkers_QuartTime_FC3rd_adjPV05$cluster == "1",]
Cluster2_Markers <- Integration1_AllClusterMarkers_QuartTime_FC3rd_adjPV05[Integration1_AllClusterMarkers_QuartTime_FC3rd_adjPV05$cluster == "2",]
Cluster3_Markers <- Integration1_AllClusterMarkers_QuartTime_FC3rd_adjPV05[Integration1_AllClusterMarkers_QuartTime_FC3rd_adjPV05$cluster == "3",]
Cluster4_Markers <- Integration1_AllClusterMarkers_QuartTime_FC3rd_adjPV05[Integration1_AllClusterMarkers_QuartTime_FC3rd_adjPV05$cluster == "4",]
Cluster5_Markers <- Integration1_AllClusterMarkers_QuartTime_FC3rd_adjPV05[Integration1_AllClusterMarkers_QuartTime_FC3rd_adjPV05$cluster == "5",]
Cluster6_Markers <- Integration1_AllClusterMarkers_QuartTime_FC3rd_adjPV05[Integration1_AllClusterMarkers_QuartTime_FC3rd_adjPV05$cluster == "6",]
Cluster7_Markers <- Integration1_AllClusterMarkers_QuartTime_FC3rd_adjPV05[Integration1_AllClusterMarkers_QuartTime_FC3rd_adjPV05$cluster == "7",]
Cluster8_Markers <- Integration1_AllClusterMarkers_QuartTime_FC3rd_adjPV05[Integration1_AllClusterMarkers_QuartTime_FC3rd_adjPV05$cluster == "8",]
Cluster9_Markers <- Integration1_AllClusterMarkers_QuartTime_FC3rd_adjPV05[Integration1_AllClusterMarkers_QuartTime_FC3rd_adjPV05$cluster == "9",]
Cluster10_Markers <- Integration1_AllClusterMarkers_QuartTime_FC3rd_adjPV05[Integration1_AllClusterMarkers_QuartTime_FC3rd_adjPV05$cluster == "10",]
Cluster11_Markers <- Integration1_AllClusterMarkers_QuartTime_FC3rd_adjPV05[Integration1_AllClusterMarkers_QuartTime_FC3rd_adjPV05$cluster == "11",]
Cluster12_Markers <- Integration1_AllClusterMarkers_QuartTime_FC3rd_adjPV05[Integration1_AllClusterMarkers_QuartTime_FC3rd_adjPV05$cluster == "12",]
Cluster13_Markers <- Integration1_AllClusterMarkers_QuartTime_FC3rd_adjPV05[Integration1_AllClusterMarkers_QuartTime_FC3rd_adjPV05$cluster == "13",]
Cluster14_Markers <- Integration1_AllClusterMarkers_QuartTime_FC3rd_adjPV05[Integration1_AllClusterMarkers_QuartTime_FC3rd_adjPV05$cluster == "14",]
Cluster15_Markers <- Integration1_AllClusterMarkers_QuartTime_FC3rd_adjPV05[Integration1_AllClusterMarkers_QuartTime_FC3rd_adjPV05$cluster == "15",]
Cluster16_Markers <- Integration1_AllClusterMarkers_QuartTime_FC3rd_adjPV05[Integration1_AllClusterMarkers_QuartTime_FC3rd_adjPV05$cluster == "16",]

AllGenes <- Integration1_AllClusterMarkers_QuartTime_FC3rd_adjPV05$gene
length(AllGenes)

Clusters_Integration_QuarTime_FC3rd <- list(Cluster1_Markers,Cluster2_Markers,
                                            Cluster3_Markers,Cluster4_Markers,Cluster5_Markers,
                                            Cluster6_Markers,Cluster7_Markers,Cluster8_Markers,
                                            Cluster9_Markers,Cluster10_Markers,Cluster11_Markers,
                                            Cluster12_Markers,Cluster13_Markers,Cluster14_Markers,
                                            Cluster15_Markers,Cluster16_Markers)



for (df in Clusters_Integration_QuarTime_FC3rd){
  print(dim(df))
}

# check
table(Integration1_AllClusterMarkers_QuartTime_FC3rd_adjPV05$cluster)

saveRDS(Clusters_Integration_QuarTime_FC3rd,'Clusters_IntegrationFirst_Quart_3rdQuarFC_TimeOnlyAlt.rds')

####################################
### All Clusters --- Integration ###
####################################

Integration1_AllClusterMarkers_QuarDose <- readRDS("DoseMarkers_pv05_DropFiltered.rds")
dim(Integration1_AllClusterMarkers_QuarDose)
Integration1_AllClusterMarkers_QuarDose$rawFC <- 2^Integration1_AllClusterMarkers_QuarDose$avg_log2FC

RawFC_3rdQuarDose <- quantile(Integration1_AllClusterMarkers_QuarDose$rawFC,0.75)

### 3rd QuarDoseer ###

Integration1_AllClusterMarkers_QuarDose_FC3rd <- Integration1_AllClusterMarkers_QuarDose[Integration1_AllClusterMarkers_QuarDose$rawFC > RawFC_3rdQuarDose,]
Integration1_AllClusterMarkers_QuarDose_FC3rd_adjPV05 <- Integration1_AllClusterMarkers_QuarDose_FC3rd[Integration1_AllClusterMarkers_QuarDose_FC3rd$p_val_adj2  < 0.05,]
dim(Integration1_AllClusterMarkers_QuarDose_FC3rd_adjPV05)

Cluster1_Markers <- Integration1_AllClusterMarkers_QuarDose_FC3rd_adjPV05[Integration1_AllClusterMarkers_QuarDose_FC3rd_adjPV05$cluster == "1",]
Cluster2_Markers <- Integration1_AllClusterMarkers_QuarDose_FC3rd_adjPV05[Integration1_AllClusterMarkers_QuarDose_FC3rd_adjPV05$cluster == "2",]
Cluster3_Markers <- Integration1_AllClusterMarkers_QuarDose_FC3rd_adjPV05[Integration1_AllClusterMarkers_QuarDose_FC3rd_adjPV05$cluster == "3",]
Cluster4_Markers <- Integration1_AllClusterMarkers_QuarDose_FC3rd_adjPV05[Integration1_AllClusterMarkers_QuarDose_FC3rd_adjPV05$cluster == "4",]
Cluster5_Markers <- Integration1_AllClusterMarkers_QuarDose_FC3rd_adjPV05[Integration1_AllClusterMarkers_QuarDose_FC3rd_adjPV05$cluster == "5",]
Cluster6_Markers <- Integration1_AllClusterMarkers_QuarDose_FC3rd_adjPV05[Integration1_AllClusterMarkers_QuarDose_FC3rd_adjPV05$cluster == "6",]
Cluster7_Markers <- Integration1_AllClusterMarkers_QuarDose_FC3rd_adjPV05[Integration1_AllClusterMarkers_QuarDose_FC3rd_adjPV05$cluster == "7",]
Cluster8_Markers <- Integration1_AllClusterMarkers_QuarDose_FC3rd_adjPV05[Integration1_AllClusterMarkers_QuarDose_FC3rd_adjPV05$cluster == "8",]
Cluster9_Markers <- Integration1_AllClusterMarkers_QuarDose_FC3rd_adjPV05[Integration1_AllClusterMarkers_QuarDose_FC3rd_adjPV05$cluster == "9",]
Cluster10_Markers <- Integration1_AllClusterMarkers_QuarDose_FC3rd_adjPV05[Integration1_AllClusterMarkers_QuarDose_FC3rd_adjPV05$cluster == "10",]
Cluster11_Markers <- Integration1_AllClusterMarkers_QuarDose_FC3rd_adjPV05[Integration1_AllClusterMarkers_QuarDose_FC3rd_adjPV05$cluster == "11",]
Cluster12_Markers <- Integration1_AllClusterMarkers_QuarDose_FC3rd_adjPV05[Integration1_AllClusterMarkers_QuarDose_FC3rd_adjPV05$cluster == "12",]
Cluster13_Markers <- Integration1_AllClusterMarkers_QuarDose_FC3rd_adjPV05[Integration1_AllClusterMarkers_QuarDose_FC3rd_adjPV05$cluster == "13",]
Cluster14_Markers <- Integration1_AllClusterMarkers_QuarDose_FC3rd_adjPV05[Integration1_AllClusterMarkers_QuarDose_FC3rd_adjPV05$cluster == "14",]
Cluster15_Markers <- Integration1_AllClusterMarkers_QuarDose_FC3rd_adjPV05[Integration1_AllClusterMarkers_QuarDose_FC3rd_adjPV05$cluster == "15",]
Cluster16_Markers <- Integration1_AllClusterMarkers_QuarDose_FC3rd_adjPV05[Integration1_AllClusterMarkers_QuarDose_FC3rd_adjPV05$cluster == "16",]

AllGenes <- Integration1_AllClusterMarkers_QuarDose_FC3rd_adjPV05$gene
length(AllGenes)

Clusters_Integration_QuarDose_FC3rd <- list(Cluster1_Markers,Cluster2_Markers,
                                            Cluster3_Markers,Cluster4_Markers,Cluster5_Markers,
                                            Cluster6_Markers,Cluster7_Markers,Cluster8_Markers,
                                            Cluster9_Markers,Cluster10_Markers,Cluster11_Markers,
                                            Cluster12_Markers,Cluster13_Markers,Cluster14_Markers,
                                            Cluster15_Markers,Cluster16_Markers)



for (df in Clusters_Integration_QuarDose_FC3rd){
  print(dim(df))
}
table(Integration1_AllClusterMarkers_QuarDose_FC3rd_adjPV05$cluster)

saveRDS(Clusters_Integration_QuarDose_FC3rd,'Clusters_IntegrationFirst_Quart_3rdQuartFC_DoseOnlyAlt.rds')
