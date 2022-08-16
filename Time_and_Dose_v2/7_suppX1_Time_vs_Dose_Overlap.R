########################
### Setup Enviroment ###
########################

### Folder
setwd("/Users/nicholaspanchy/Documents/Work_UTK/DoseTime_FullCorr_FromScratch/")

seuratCon_combo_E_mat_nsprcomp <- readRDS("seuratCon_combo_E_mat_nsprcomp_fullcorr.rds")
seuratCon_combo_M_mat_nsprcomp <- readRDS("seuratCon_combo_M_mat_nsprcomp_fullcorr.rds")

seuratTime_combo_E_mat_nsprcomp <- readRDS("seuratTime_combo_E_mat_nsprcomp_fullcorr.rds")
seuratTime_combo_M_mat_nsprcomp <- readRDS("seuratTime_combo_M_mat_nsprcomp_fullcorr.rds")

integated_meta <- readRDS("Integrated_MetaData.rds")
BatchA_samples <- row.names(integated_meta[integated_meta$Batch=="A",])
BatchB_samples <- row.names(integated_meta[integated_meta$Batch=="B",])
BatchC_samples <- row.names(integated_meta[integated_meta$Batch=="C",])

plot(seuratTime_combo_E_mat_nsprcomp$x[BatchB_samples,1],seuratTime_combo_M_mat_nsprcomp$x[BatchB_samples,1],pch=16,cex=0.5,col="green",xlab="E PC1",ylab="M PC1")
points(seuratTime_combo_E_mat_nsprcomp$x[BatchA_samples,1],seuratTime_combo_M_mat_nsprcomp$x[BatchA_samples,1],pch=16,cex=0.5,col="red")
points(seuratCon_combo_E_mat_nsprcomp$x[BatchC_samples,1],seuratCon_combo_M_mat_nsprcomp$x[BatchC_samples,1],pch=16,cex=0.5,col="blue")
