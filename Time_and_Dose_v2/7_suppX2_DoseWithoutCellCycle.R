########################
### Setup Enviroment ###
########################

### Folder
setwd("/Users/nicholaspanchy/Documents/Work_UTK/DoseTime_FromScratch2/")

# Dependencies
library(Seurat)
# If there is a spatstat conflict from BioConductor packages (this use to be a problem, might not be needed anymore)
# install.packages('https://cran.r-project.org/src/contrib/Archive/spatstat/spatstat_1.64-1.tar.gz', repos=NULL,type="source", INSTALL_opts = "--no-lock") 

#########################
### Process Dose Data ###
#########################

dose_data.combined <- readRDS("SeuratData/Seurat_concentration_DoseOnly_FullCorr_hg38_v2CC.rds")
DimPlot(dose_data.combined, reduction = "umap",group.by = "Phase")


dose_data.combined <- ScaleData(dose_data.combined,vars.to.regress=c("nCount_RNA","percent.mito"),verbose=FALSE)
saveRDS(dose_data.combined,"Dose_NoCellCycle.rds")
dose_data.combined <- readRDS("Dose_NoCellCycle.rds")

dose_data.combined <- RunPCA(dose_data.combined, npcs = 15, verbose = FALSE)
dose_data.combined <- RunUMAP(dose_data.combined, reduction = "pca", dims = 1:15)
dose_data.combined <- FindNeighbors(dose_data.combined, reduction = "pca", dims = 1:15)
dose_data.combined <- FindClusters(dose_data.combined, resolution = 0.5)

DimPlot(dose_data.combined, reduction = "umap",group.by = "Phase")

combined_nointegration.combined <- readRDS("SeuratData/Seurat_NoIntegration_Combined_FullCorreciton_hg38_v2CC.rds")
DimPlot(combined_nointegration.combined, reduction = "umap",group.by = "Batch")

#######################
### Phase Bar Plots ###
#######################

integated_meta <- readRDS("Integrated_MetaData.rds")
dose_samples <- readRDS("DoseSamples.rds")
time_samples <- readRDS("TimeSamples.rds")

dose_phase_table <- table(integated_meta[dose_samples,]$Phase,integated_meta[dose_samples,]$GBC_pM)
time_phase_table <- table(integated_meta[time_samples,]$Phase,integated_meta[time_samples,]$Time)

dose_phase_tablb_freq <- sweep(dose_phase_table,2,colSums(dose_phase_table),`/`)
time_phase_table_freq <- sweep(time_phase_table,2,colSums(time_phase_table),`/`)

dose_phase_tablb_df <- as.data.frame(dose_phase_tablb_freq)
time_phase_table_df <- as.data.frame(time_phase_table_freq)

ggplot(dose_phase_tablb_df,aes(x=Var2,y=Freq,fill=Var1)) + geom_bar(stat="identity", position=position_dodge()) +
  xlab("Dose") + ylab("Count") + labs(fill="Phase") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1)

ggplot(time_phase_table_df,aes(x=Var2,y=Freq,fill=Var1)) + geom_bar(stat="identity", position=position_dodge()) +
  xlab("Time") + ylab("Count") + labs(fill="Phase") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1)
