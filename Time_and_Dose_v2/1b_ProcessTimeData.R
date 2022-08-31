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
### Process Time Data ###
#########################

time_data.combined <- readRDS("SeuratData/Seurat_integrated_TimeOnly_FullCorr_hg38_v2CC.rds")

### Save Processed Data ###

# Subset and save parts of the Seurat Obects so we don't have to use the full
# object all of the time. This is useful if you are working on a latop

Scaled_data <- time_data.combined@assays$integrated@scale.data
dose_meta <- time_data.combined@meta.data
UMAP_values <- time_data.combined[["umap"]]@cell.embeddings

saveRDS(Scaled_data,"Time_ScaledData.rds")
saveRDS(dose_meta,"Time_MetaData.rds")
saveRDS(UMAP_values,"Time_UMAP_values.rds")

### Generate and Save nnPCA Components ###

# Load Genes Sets #
EMT_genes <- read.table("EMTGeneUpdate3.txt",header=TRUE)
E_genes <- subset(EMT_genes, EMT_genes$Annotation == 'E')
M_genes <- subset(EMT_genes, EMT_genes$Annotation == 'M')
Gene_Sets <- list('Egenes'=E_genes$Gene,'Mgenes'=M_genes$Gene)

### Run NN-PCA ###
library(nsprcomp)
set.seed(5) # Controls randomness for nnPC, minor variations , mainly to keep enriched gene lists consistent (+/-1 genes)

# Subset data
Scaled_data_E <- Scaled_data[row.names(Scaled_data) %in% E_genes$Gene,]
Scaled_data_M <- Scaled_data[row.names(Scaled_data) %in% M_genes$Gene,]
dim(Scaled_data_E)
dim(Scaled_data_M)

# Run nsprcomp
seuratTime_combo_E_mat_nsprcomp <- nsprcomp(t(Scaled_data_E),nneg=TRUE,ncomp=25,em_maxiter = 10000,em_tol = 0.00001)
seuratTime_combo_M_mat_nsprcomp <- nsprcomp(t(Scaled_data_M),nneg=TRUE,ncomp=25,em_maxiter = 10000,em_tol = 0.00001)

saveRDS(seuratTime_combo_E_mat_nsprcomp,"seuratTime_combo_E_mat_nsprcomp_fullcorr.rds")
saveRDS(seuratTime_combo_M_mat_nsprcomp,"seuratTime_combo_M_mat_nsprcomp_fullcorr.rds")
seuratTime_combo_E_mat_nsprcomp <- readRDS("seuratTime_combo_E_mat_nsprcomp_fullcorr.rds")
seuratTime_combo_M_mat_nsprcomp <- readRDS("seuratTime_combo_M_mat_nsprcomp_fullcorr.rds")