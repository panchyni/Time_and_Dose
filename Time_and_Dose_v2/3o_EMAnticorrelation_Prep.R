########################
### Setup Enviroment ###
########################

### Folder
setwd("/Users/nicholaspanchy/Documents/Work_UTK/DoseTime_FromScratch2/")

### IMPORTANT NOTE ###
# Some wonkiness is expected with UMAP across systems, even with set seed
# https://stackoverflow.com/questions/67101829/seurat-umap-visualization-result-is-mirrored-after-running-in-two-identical-envi

# Dependencies
# install.packages('https://cran.r-project.org/src/contrib/Archive/spatstat/spatstat_1.64-1.tar.gz', repos=NULL,type="source", INSTALL_opts = "--no-lock") 
library(Seurat) #v3.x
library(viridis)
library(ggplot2)
library(scales)
library(gridExtra)

set.seed(5) # Controls randomness for nnPC, minor variations , mainly to keep enriched gene lists consistent (+/-1 genes)

########################
### Generate Figures ###
########################

### Load Processed Data ###

integated_meta <- readRDS("Integrated_MetaData.rds")
UMAP_values <- readRDS("Integrated_UMAP_values.rds")

dose_samples <- readRDS("DoseSamples.rds")
time_samples <- readRDS("TimeSamples.rds")

Integrated_E_mat_nsprcomp <- readRDS("Integrated_E_mat_nsprcomp.rds")
Integrated_M_mat_nsprcomp <- readRDS("Integrated_M_mat_nsprcomp.rds")

#######################################################
### Correlation of NONEMT Genes with E and M-scores ###
#######################################################

### Load Genes Sets ###
EMT_genes <- read.table("EMTGeneUpdate3.txt",header=TRUE)
E_genes <- subset(EMT_genes, EMT_genes$Annotation == 'E')
M_genes <- subset(EMT_genes, EMT_genes$Annotation == 'M')
Gene_Sets <- list('Egenes'=E_genes$Gene,'Mgenes'=M_genes$Gene)

### Subset Integrated Data ###
integated_scale <- readRDS("Integrated_ScaledData.rds")
seurat_combo_All_mat_NOEMT <- integated_scale[-which(row.names(integated_scale) %in% EMT_genes$Gene),]
dim(seurat_combo_All_mat_NOEMT)
saveRDS(seurat_combo_All_mat_NOEMT,"seurat_combo_All_mat_NOEMT.rds")

### Split Dose and Time Data ###
seurat_combo_All_mat_NOEMT <- readRDS("seurat_combo_All_mat_NOEMT.rds")
seurat_combo_All_mat_NOEMT <- t(seurat_combo_All_mat_NOEMT)

seurat_combo_All_mat_NOEMT_TimeT <- seurat_combo_All_mat_NOEMT[time_samples,]
seurat_combo_All_mat_NOEMT_ConT <- seurat_combo_All_mat_NOEMT[dose_samples,]
dim(seurat_combo_All_mat_NOEMT_TimeT)
dim(seurat_combo_All_mat_NOEMT_ConT)
rm(seurat_combo_All_mat_NOEMT)
saveRDS(seurat_combo_All_mat_NOEMT_TimeT,"seurat_combo_All_mat_NOEMT_TimeT.rds")
saveRDS(seurat_combo_All_mat_NOEMT_ConT,"seurat_combo_All_mat_NOEMT_ConT.rds")


### Update NonEMT Gene Symbols with whatever is current ###

seurat_combo_All_mat_NOEMT_TimeT <- readRDS("seurat_combo_All_mat_NOEMT_TimeT.rds")
seurat_combo_All_mat_NOEMT_ConT <- readRDS("seurat_combo_All_mat_NOEMT_ConT.rds")
seurat_combo_All_mat_NOEMT_TimeT <- seurat_combo_All_mat_NOEMT_TimeT[time_samples,]

# Update gene map
library(HGNChelper)
source("NonOverlappingReplacement.R")
currentmamp <- getCurrentHumanMap()

check_time <- checkGeneSymbols(colnames(seurat_combo_All_mat_NOEMT_TimeT),unmapped.as.na = FALSE,map=currentmamp)
check_dose <- checkGeneSymbols(colnames(seurat_combo_All_mat_NOEMT_ConT),unmapped.as.na = FALSE,map=currentmamp)

check_time_NoOverlap <- NonOverlappingReplacement(check_time)
check_dose_NoOverlap <- NonOverlappingReplacement(check_dose)

colnames(seurat_combo_All_mat_NOEMT_TimeT) <- check_time_NoOverlap$Suggested.Symbol
colnames(seurat_combo_All_mat_NOEMT_ConT) <- check_dose_NoOverlap$Suggested.Symbol

# Calcluate Correlation with Time Samples
correl_TimeE <- cor(seurat_combo_All_mat_NOEMT_TimeT,Integrated_E_mat_nsprcomp$x[time_samples,1])
hist(correl_TimeE)
correl_TimeM <- cor(seurat_combo_All_mat_NOEMT_TimeT,Integrated_M_mat_nsprcomp$x[time_samples,1])
hist(correl_TimeM)
correl_TimeM2 <- cor(seurat_combo_All_mat_NOEMT_TimeT,Integrated_M_mat_nsprcomp$x[time_samples,2])
hist(correl_TimeM2)

saveRDS(correl_TimeE,"correl_TimeE.rds")
saveRDS(correl_TimeM,"correl_TimeM.rds")
saveRDS(correl_TimeM2,"correl_TimeM2.rds")

# Calcualte Correlation with Dose Samples
correl_ConE <- cor(seurat_combo_All_mat_NOEMT_ConT,Integrated_E_mat_nsprcomp$x[dose_samples,1])
hist(correl_ConE)
correl_ConM <- cor(seurat_combo_All_mat_NOEMT_ConT,Integrated_M_mat_nsprcomp$x[dose_samples,1])
hist(correl_ConM)
correl_ConM2 <- cor(seurat_combo_All_mat_NOEMT_ConT,Integrated_M_mat_nsprcomp$x[dose_samples,2])
hist(correl_ConM2)

saveRDS(correl_ConE,"correl_ConE.rds")
saveRDS(correl_ConM,"correl_ConM.rds")
saveRDS(correl_ConM2,"correl_ConM2.rds")

dim(seurat_combo_All_mat_NOEMT_ConT)
all_genes <- colnames(seurat_combo_All_mat_NOEMT_ConT)
saveRDS(all_genes,"all_genes.rds")