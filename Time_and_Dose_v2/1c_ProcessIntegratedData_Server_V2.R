########################
### Setup Enviroment ###
########################

### Folder
setwd("/Users/nicholaspanchy/Documents/Work_UTK/DoseTime_FromScratch2/")
#setwd("/home/panchy/DoseTIme/_fromscratch_2")

# Dependencies
library(Seurat)

###############################
### Process Integrated Data ###
###############################

integration.combined.All <- readRDS("SeuratData/Seurat_integrated_TimeAndCon_FullCorrection_hg38_v2CC.rds")

# Check against original list
check_table <- read.table("GBCannot/CBC_GBC_summary.txt.alt_names",head=TRUE)
check_table$CBC_alt <- cbind(lapply(check_table$CBC.1, function(x) {sub(".1", "-1", x)}))
row.names(check_table) <- check_table$CBC_alt
overlap <- intersect(check_table$CBC_alt,row.names(integration.combined.All@meta.data))
overlap <- unlist(overlap)

check_table <- as.data.frame(check_table)
check_table$meta_annot <- 0
check_table[overlap,]$meta_annot <- integration.combined.All@meta.data[overlap ,]$GBC

table(check_table$meta_annot,check_table$GBC)
# Note: 0s are samples present in the annotation, but not post filtering
# 149 samples in check table but not in processed data == 149 0s
# 4 samples that are GBC0 in the processed, but not the annotation (these were dropped in a previous filtering)

# Make replacements
#check_table[check_table$GBC == "GBC1" & check_table$meta_annot == "GBC0",]
integration.combined.All@meta.data["CTGAAACAGTTGTCGT-1",]$GBC = "GBC1"
integration.combined.All@meta.data["CTGAAACAGTTGTCGT-1",]$GBC_values = "1"
integration.combined.All@meta.data["CTGAAACAGTTGTCGT-1",]$GBC_pM = 0
integration.combined.All@meta.data["TGTATTCCATTTGCTT-1",]$GBC = "GBC1"
integration.combined.All@meta.data["TGTATTCCATTTGCTT-1",]$GBC_values = "1"
integration.combined.All@meta.data["TGTATTCCATTTGCTT-1",]$GBC_pM = 0
#check_table[check_table$GBC == "GBC4" & check_table$meta_annot == "GBC0",]
integration.combined.All@meta.data["ACGCCGAGTGGTACAG-1",]$GBC = "GBC4"
integration.combined.All@meta.data["ACGCCGAGTGGTACAG-1",]$GBC_values = "4"
integration.combined.All@meta.data["ACGCCGAGTGGTACAG-1",]$GBC_pM = 50
#check_table[check_table$GBC == "GBC9" & check_table$meta_annot == "GBC0",]
integration.combined.All@meta.data["CCGGGATCAGGCGATA-1",]$GBC = "GBC9"
integration.combined.All@meta.data["CCGGGATCAGGCGATA-1",]$GBC_values = "9"
integration.combined.All@meta.data["CCGGGATCAGGCGATA-1",]$GBC_pM = 200

check_table[overlap,]$meta_annot <- integration.combined.All@meta.data[overlap ,]$GBC
table(check_table$meta_annot,check_table$GBC)
# The are 4 are fixed and annotations line up except for 0s, as expected

# Check correspondence post correction
table(integration.combined.All@meta.data$GBC,integration.combined.All@meta.data$GBC_values)
table(integration.combined.All@meta.data$GBC,integration.combined.All@meta.data$GBC_pM)

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

integated_meta <- integration.combined.All@meta.data
integated_scale <- integration.combined.All@assays$integrated@scale.data
UMAP_values <- integration.combined.All[["umap"]]@cell.embeddings

saveRDS(integated_scale,"Integrated_ScaledData.rds")
saveRDS(integated_meta,"Integrated_MetaData.rds")
saveRDS(UMAP_values,"Integrated_UMAP_values.rds")

# Save the Processed Seurat Obejct
saveRDS(integration.combined.All,"SeuratData/Seurat_integrated_TimeAndCon_FullCorrection_hg38_Processed.rds")
dose_samples <- rownames(integration.combined.All@meta.data[integration.combined.All$GBC_pM > -1,])
time_samples <- rownames(integration.combined.All@meta.data[integration.combined.All$Time < 100,])
saveRDS(dose_samples,"DoseSamples.rds")
saveRDS(time_samples,"TimeSamples.rds")
#rm(integration.combined.All)

### Generate and Save nnPCA Components ###

# Load Genes Sets #
EMT_genes <- read.table("EMTGeneUpdate3.txt",header=TRUE)
E_genes <- subset(EMT_genes, EMT_genes$Annotation == 'E')
M_genes <- subset(EMT_genes, EMT_genes$Annotation == 'M')
Gene_Sets <- list('Egenes'=E_genes$Gene,'Mgenes'=M_genes$Gene)

### Run NN-PCA ###
library(nsprcomp)
set.seed(5) # Controls randomness for nnPC, minor variations , mainly to keep enriched gene lists consistent (+/-1 genes)

integated_scale <- readRDS("Integrated_ScaledData.rds")
integated_meta <- readRDS("Integrated_MetaData.rds")
dose_samples <- readRDS("DoseSamples.rds")
time_samples <- readRDS("TimeSamples.rds")

# Subset Data
Scaled_Integrated_data_E <- integated_scale[row.names(integated_scale) %in% E_genes$Gene,]
Scaled_Integrated_data_M <- integated_scale[row.names(integated_scale) %in% M_genes$Gene,]
dim(Scaled_Integrated_data_E)
dim(Scaled_Integrated_data_M)

# Run nsprcomp
Integrated_E_mat_nsprcomp <- nsprcomp(t(Scaled_Integrated_data_E),nneg=TRUE,ncomp=25,em_maxiter = 10000,em_tol = 0.00001)
Integrated_M_mat_nsprcomp <- nsprcomp(t(Scaled_Integrated_data_M),nneg=TRUE,ncomp=25,em_maxiter = 10000,em_tol = 0.00001)

saveRDS(Integrated_E_mat_nsprcomp,"Integrated_E_mat_nsprcomp.rds")
saveRDS(Integrated_M_mat_nsprcomp,"Integrated_M_mat_nsprcomp.rds")
Integrated_E_mat_nsprcomp <- readRDS("Integrated_E_mat_nsprcomp.rds")
Integrated_M_mat_nsprcomp <- readRDS("Integrated_M_mat_nsprcomp.rds")

### Additional Files ###

saveRDS(Scaled_Integrated_data_E,"Scaled_Integrated_data_E.rds")
saveRDS(Scaled_Integrated_data_M,"Scaled_Integrated_data_M.rds")

#integration.combined.All <- readRDS("SeuratData/Seurat_integrated_TimeAndCon_FullCorrection_hg38_Processed.rds")
integration.combined.EMT <- subset(integration.combined.All,features = EMT_genes$Gene)
print(dim(integration.combined.EMT))
saveRDS(integration.combined.EMT,"integration.combined.EMT.rds")

# Make and save the more complete UMAP object because its slow
UMAP_annot <- merge(UMAP_values,integated_meta,by=0)
row.names(UMAP_annot) <- UMAP_annot$Row.names
UMAP_samples <- row.names(UMAP_annot)

UMAP_annot$VIM <- integated_scale["VIM",UMAP_annot$Row.names]
UMAP_annot$EPCAM <- integated_scale["EPCAM",UMAP_annot$Row.names]
UMAP_annot$FN1 <- integated_scale["FN1",UMAP_annot$Row.names]
rm(integated_scale)

UMAP_annot$E_PC1 <- Integrated_E_mat_nsprcomp$x[UMAP_samples,1]
UMAP_annot$M_PC1 <- Integrated_M_mat_nsprcomp$x[UMAP_samples,1]
UMAP_annot$M_PC2 <- Integrated_M_mat_nsprcomp$x[UMAP_samples,2]

saveRDS(UMAP_annot,"UMAP_annot.rds")

############################################################
### Supplemental --> Affect of droping samples on scores ###
############################################################

integated_scale <- readRDS("Integrated_ScaledData.rds")
dose_samples <- readRDS("DoseSamples.rds")
time_samples <- readRDS("TimeSamples.rds")

# Load Genes Sets #
EMT_genes <- read.table("EMTGeneUpdate3.txt",header=TRUE)
E_genes <- subset(EMT_genes, EMT_genes$Annotation == 'E')
M_genes <- subset(EMT_genes, EMT_genes$Annotation == 'M')
Gene_Sets <- list('Egenes'=E_genes$Gene,'Mgenes'=M_genes$Gene)

### Check unlabeled data
Scaled_data_Labeled_E <- integated_scale[row.names(integated_scale) %in% E_genes$Gene,union(dose_samples,time_samples)]
Scaled_data_Labeled_M <- integated_scale[row.names(integated_scale) %in% M_genes$Gene,union(dose_samples,time_samples)]

# Check against scores without unlabeld (little diff; cor > 0.999)
library(nsprcomp)
set.seed(5) # Controls randomness for nnPC, minor variations , mainly to keep enriched gene lists consistent (+/-1 genes)

Integrated_E_mat_Labeled_nsprcomp <- nsprcomp(t(Scaled_data_Labeled_E),nneg=TRUE,ncomp=5,em_maxiter = 10000,em_tol = 0.00001)
Integrated_M_mat_Labeled_nsprcomp <- nsprcomp(t(Scaled_data_Labeled_M),nneg=TRUE,ncomp=5,em_maxiter = 10000,em_tol = 0.00001)

saveRDS(Integrated_E_mat_Labeled_nsprcomp,"Integrated_E_mat_Labeled_nsprcomp.rds")
saveRDS(Integrated_M_mat_Labeled_nsprcomp,"Integrated_M_mat_Labeled_nsprcomp.rds")
Integrated_E_mat_Labeled_nsprcomp <- readRDS("Integrated_E_mat_Labeled_nsprcomp.rds")
Integrated_M_mat_Labeled_nsprcomp <- readRDS("Integrated_M_mat_Labeled_nsprcomp.rds")

# Check against scores without unlabeld (little diff cor > 0.999)
labeled_sampels <- union(dose_samples,time_samples)
cor(Integrated_E_mat_Labeled_nsprcomp$x[labeled_sampels,1],Integrated_E_mat_nsprcomp$x[labeled_sampels,1])
cor(Integrated_M_mat_Labeled_nsprcomp$x[labeled_sampels,1],Integrated_M_mat_nsprcomp$x[labeled_sampels,1])

### Check unlabeled data + dose samples
# Uncomment and run after running 7d

#integated_meta <- readRDS("Integrated_MetaData.rds")
#dose_samples <- readRDS("DoseSamples.rds")
#time_samples <- readRDS("TimeSamples.rds")

#time_doublets <- readRDS("time_doublets.rds")
#time_wo_doublets <- setdiff(time_samples,time_doublets)
#time_samples <- time_wo_doublets

#Scaled_data_Labeled_E <- readRDS("Scaled_Integrated_data_E.rds")
#Scaled_data_Labeled_M <- readRDS("Scaled_Integrated_data_M.rds")

#Scaled_data_Labeled_E_fliter <- Scaled_data_Labeled_E[,union(dose_samples,time_samples)]
#Scaled_data_Labeled_M_filter <- Scaled_data_Labeled_M[,union(dose_samples,time_samples)]

#dim(Scaled_data_Labeled_E_fliter)
#dim(Scaled_data_Labeled_M_filter)

# Quick Check so use fewer npcs
#library(nsprcomp)
#set.seed(5) # Controls randomness for nnPC, minor variations , mainly to keep enriched gene lists consistent (+/-1 genes)

#Integrated_E_mat_LabeledFilter_nsprcomp <- nsprcomp(t(Scaled_data_Labeled_E_fliter),nneg=TRUE,ncomp=5,em_maxiter = 10000,em_tol = 0.00001)
#Integrated_M_mat_LabeledFilter_nsprcomp <- nsprcomp(t(Scaled_data_Labeled_M_filter),nneg=TRUE,ncomp=5,em_maxiter = 10000,em_tol = 0.00001)

#saveRDS(Integrated_E_mat_LabeledFilter_nsprcomp,"Integrated_E_mat_LabeledFilter_nsprcomp.rds")
#saveRDS(Integrated_M_mat_LabeledFilter_nsprcomp,"Integrated_M_mat_LabeledFilter_nsprcomp.rds")

# Check against scores without unlabeld (little diff cor > 0.9999)
#labeledfiltered_samples <- union(dose_samples,time_samples)
#cor(Integrated_E_mat_LabeledFilter_nsprcomp$x[labeledfiltered_samples,1],Integrated_E_mat_nsprcomp$x[labeledfiltered_samples,1])
#cor(Integrated_M_mat_LabeledFilter_nsprcomp$x[labeledfiltered_samples,1],Integrated_M_mat_nsprcomp$x[labeledfiltered_samples,1])
