########################
### Setup Enviroment ###
########################

### Folder
#setwd("/Users/nicholaspanchy/Documents/Work_UTK/DoseTime_FullCorr_FromScratch/")
setwd("/home/panchy/DoseTIme/_fromscratch_server")

# Dependencies
library(Seurat)
# If there is a spatstat conflict from BioConductor packages (this use to be a problem, might not be needed anymore)
# install.packages('https://cran.r-project.org/src/contrib/Archive/spatstat/spatstat_1.64-1.tar.gz', repos=NULL,type="source", INSTALL_opts = "--no-lock") 

#########################
### Process Dose Data ###
#########################

dose_data.combined <- readRDS("SeuratData/Seurat_concentration_DoseOnly_FullCorr_hg38.rds")
dose_data.combined <- ScaleData(dose_data.combined,vars.to.regress=c("nCount_RNA","percent.mito"),verbose=FALSE)
saveRDS(dose_data.combined,"Dose_NoCellCycle.rds")

dose_data.combined <- RunPCA(dose_data.combined, npcs = 15, verbose = FALSE)
dose_data.combined <- RunUMAP(dose_data.combined, reduction = "pca", dims = 1:15)
dose_data.combined <- FindNeighbors(dose_data.combined, reduction = "pca", dims = 1:15)
dose_data.combined <- FindClusters(dose_data.combined, resolution = 0.5)

DimPlot(ConCC_integration.combined, reduction = "umap",group.by = "Phase")
