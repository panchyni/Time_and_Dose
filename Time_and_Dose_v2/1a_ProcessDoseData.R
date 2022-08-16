########################
### Setup Enviroment ###
########################

### Folder
setwd("/Users/nicholaspanchy/Documents/Work_UTK/DoseTime_FullCorr_FromScratch/")

# Dependencies
library(Seurat)
# If there is a spatstat conflict from BioConductor packages (this use to be a problem, might not be needed anymore)
# install.packages('https://cran.r-project.org/src/contrib/Archive/spatstat/spatstat_1.64-1.tar.gz', repos=NULL,type="source", INSTALL_opts = "--no-lock") 

#########################
### Process Dose Data ###
#########################

dose_data.combined <- readRDS("SeuratData/Seurat_concentration_DoseOnly_FullCorr_hg38.rds")

### Fix annotation ##

# Check against original list
check_table <- read.table("GBCannot/CBC_GBC_summary.txt.alt_names",head=TRUE)
check_table$CBC_alt <- cbind(lapply(check_table$CBC.1, function(x) {sub(".1", "-1", x)}))
row.names(check_table) <- check_table$CBC_alt
overlap <- intersect(check_table$CBC_alt,row.names(dose_data.combined@meta.data))
overlap <- unlist(overlap)

check_table <- as.data.frame(check_table)
check_table$meta_annot <- 0
check_table[overlap,]$meta_annot <- dose_data.combined@meta.data[overlap ,]$GBC

table(check_table$meta_annot,check_table$GBC)
# Note: 0s are samples present in the annotation, but not post filtering
# 149 samples in check table but not in processed data == 149 0s
# 4 samples that are GBC0 in the processed, but not the annotation (these were dropped in a previous filtering)

# Make replacements
#check_table[check_table$GBC == "GBC1" & check_table$meta_annot == "GBC0",]
dose_data.combined@meta.data["CTGAAACAGTTGTCGT-1",]$GBC = "GBC1"
dose_data.combined@meta.data["CTGAAACAGTTGTCGT-1",]$GBC_values = "1"
dose_data.combined@meta.data["CTGAAACAGTTGTCGT-1",]$GBC_pM = 0
dose_data.combined@meta.data["TGTATTCCATTTGCTT-1",]$GBC = "GBC1"
dose_data.combined@meta.data["TGTATTCCATTTGCTT-1",]$GBC_values = "1"
dose_data.combined@meta.data["TGTATTCCATTTGCTT-1",]$GBC_pM = 0
#check_table[check_table$GBC == "GBC4" & check_table$meta_annot == "GBC0",]
dose_data.combined@meta.data["ACGCCGAGTGGTACAG-1",]$GBC = "GBC4"
dose_data.combined@meta.data["ACGCCGAGTGGTACAG-1",]$GBC_values = "4"
dose_data.combined@meta.data["ACGCCGAGTGGTACAG-1",]$GBC_pM = 50
#check_table[check_table$GBC == "GBC9" & check_table$meta_annot == "GBC0",]
dose_data.combined@meta.data["CCGGGATCAGGCGATA-1",]$GBC = "GBC9"
dose_data.combined@meta.data["CCGGGATCAGGCGATA-1",]$GBC_values = "9"
dose_data.combined@meta.data["CCGGGATCAGGCGATA-1",]$GBC_pM = 200

check_table[overlap,]$meta_annot <- dose_data.combined@meta.data[overlap ,]$GBC
table(check_table$meta_annot,check_table$GBC)
# The are 4 are fixed and annotations line up except for 0s, as expected

# Check correspondence post correction
table(dose_data.combined@meta.data$GBC,dose_data.combined@meta.data$GBC_values)
table(dose_data.combined@meta.data$GBC,dose_data.combined@meta.data$GBC_pM)

# Check missing samples
check_table_0 <- check_table[check_table$meta_annot ==0,]
dim(check_table_0)
kazu_barcodes <- read.table("data_KazuConcentration/barcodes.tsv")
length(intersect(kazu_barcodes$V1,check_table_0$CBC_alt))
# All but ten are missing from the filtered barcode set
# The missing ten are filtered (these are the 10 NA sample post stringent filtering)

### Save Processed Data ###

# Subset and save parts of the Seurat Obects so we don't have to use the full
# object all of the time. This is useful if you are working on a latop

Scaled_data <- dose_data.combined@assays$RNA@scale.data
dose_meta <- dose_data.combined@meta.data
UMAP_values <- dose_data.combined[["umap"]]@cell.embeddings

saveRDS(Scaled_data,"Dose_ScaledData.rds")
saveRDS(dose_meta,"Dose_MetaData.rds")
saveRDS(UMAP_values,"Dose_UMAP_values.rds")

# Save the Processed Seurat Obejct
saveRDS(dose_data.combined,"SeuratData/Seurat_concentration_DoseOnly_FullCorr_hg38_Processed.rds")

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
seuratCon_combo_E_mat_nsprcomp <- nsprcomp(t(Scaled_data_E),nneg=TRUE,ncomp=25,em_maxiter = 10000,em_tol = 0.00001)
seuratCon_combo_M_mat_nsprcomp <- nsprcomp(t(Scaled_data_M),nneg=TRUE,ncomp=25,em_maxiter = 10000,em_tol = 0.00001)

saveRDS(seuratCon_combo_E_mat_nsprcomp,"seuratCon_combo_E_mat_nsprcomp_fullcorr.rds")
saveRDS(seuratCon_combo_M_mat_nsprcomp,"seuratCon_combo_M_mat_nsprcomp_fullcorr.rds")
seuratCon_combo_E_mat_nsprcomp <- readRDS("seuratCon_combo_E_mat_nsprcomp_fullcorr.rds")
seuratCon_combo_M_mat_nsprcomp <- readRDS("seuratCon_combo_M_mat_nsprcomp_fullcorr.rds")

# Without UnLabeled Samples (PC > 0.999)
dose_samples <- row.names(dose_meta[dose_meta$GBC_values > 0,])
Scaled_data_Labeled_E <- Scaled_data[row.names(Scaled_data) %in% E_genes$Gene,dose_samples]
Scaled_data_Labeled_M <- Scaled_data[row.names(Scaled_data) %in% M_genes$Gene,dose_samples]

seuratCon_combo_E_mat_Labeled_nsprcomp <- nsprcomp(t(Scaled_data_Labeled_E),nneg=TRUE,ncomp=25,em_maxiter = 10000,em_tol = 0.00001)
seuratCon_combo_M_mat_Labeled_nsprcomp <- nsprcomp(t(Scaled_data_Labeled_M),nneg=TRUE,ncomp=25,em_maxiter = 10000,em_tol = 0.00001)

saveRDS(seuratCon_combo_E_mat_Labeled_nsprcomp,"seuratCon_combo_E_mat_Labeled_nsprcomp.rds")
saveRDS(seuratCon_combo_M_mat_Labeled_nsprcomp,"seuratCon_combo_M_mat_Labeled_nsprcomp.rds")