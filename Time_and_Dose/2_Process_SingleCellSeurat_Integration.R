### Folder
setwd("<path_to_folder>/Time_and_Dose")

### IMPORTANT NOTE ###
# Some wonkiness is expected with UMAP across systems, even with set seed
# https://stackoverflow.com/questions/67101829/seurat-umap-visualization-result-is-mirrored-after-running-in-two-identical-envi

# Dependencies
# install.packages('https://cran.r-project.org/src/contrib/Archive/spatstat/spatstat_1.64-1.tar.gz', repos=NULL,type="source", INSTALL_opts = "--no-lock") 
library(Seurat) #v3.x
library(viridis)
library(dplyr)
library(cowplot)
library(HGNChelper)

# Start the clock!
ptm <- proc.time()
#   user     system elapsed

# Load in the Cellranger output
SRR13606301_filtered_mat <- Read10X(data.dir = "_10XData/SRR13606301_Run_hg38_NoExpect/outs/filtered_feature_bc_matrix")
SRR13606302_filtered_mat <- Read10X(data.dir = "_10XData/SRR13606302_Run_hg38_NoExpect/outs/filtered_feature_bc_matrix")
SRR13606303_filtered_mat <- Read10X(data.dir = "_10XData/SRR13606303_Run_hg38_NoExpect/outs/filtered_feature_bc_matrix")
SRR13606304_filtered_mat <- Read10X(data.dir = "_10XData/SRR13606304_Run_hg38_NoExpect/outs/filtered_feature_bc_matrix")
SRR13606305_filtered_mat <- Read10X(data.dir = "_10XData/SRR13606305_Run_hg38_NoExpect/outs/filtered_feature_bc_matrix")
SRR13606306_filtered_mat <- Read10X(data.dir = "_10XData/SRR13606306_Run_hg38_NoExpect/outs/filtered_feature_bc_matrix")
SRR13606307_filtered_mat <- Read10X(data.dir = "_10XData/SRR13606307_Run_hg38_NoExpect/outs/filtered_feature_bc_matrix")
KazuData_filtered_mat <- Read10X(data.dir = "data_KazuConcentration/")
length(intersect(SRR13606301_filtered_mat@Dimnames[[1]],KazuData_filtered_mat@Dimnames[[1]]))

# Update Gene Symbols
getCurrentHumanMap()
check1 <- checkGeneSymbols(SRR13606301_filtered_mat@Dimnames[[1]],unmapped.as.na = FALSE,species='human')
check2 <- checkGeneSymbols(SRR13606302_filtered_mat@Dimnames[[1]],unmapped.as.na = FALSE,species='human')
check3 <- checkGeneSymbols(SRR13606303_filtered_mat@Dimnames[[1]],unmapped.as.na = FALSE,species='human')
check4 <- checkGeneSymbols(SRR13606304_filtered_mat@Dimnames[[1]],unmapped.as.na = FALSE,species='human')
check5 <- checkGeneSymbols(SRR13606305_filtered_mat@Dimnames[[1]],unmapped.as.na = FALSE,species='human')
check6 <- checkGeneSymbols(SRR13606306_filtered_mat@Dimnames[[1]],unmapped.as.na = FALSE,species='human')
check7 <- checkGeneSymbols(SRR13606307_filtered_mat@Dimnames[[1]],unmapped.as.na = FALSE,species='human')
check_kazu <- checkGeneSymbols(KazuData_filtered_mat@Dimnames[[1]],unmapped.as.na = FALSE,species='human')

NonOverlappingReplacement <- function(check){
  # Go through suggested symbols and replaced the suggested
  # symbol with the exist symbol if the suggested symbol already exists
  # because duplicate name are not allowed and visiaul inspection suggests
  # at least some of these same symbol genes have different expression patterns
  
  check_suggested_table <- table(check$Suggested.Symbol)
  
  for (i in 1:dim(check)[1]){
    #print(i)
    
    # If this is a replacement row
    if(check[i,]$Approved == 'FALSE'){
      # If its actually an new symbol
      if (check[i,]$Suggested.Symbol != check[i,]$x){
        if (check_suggested_table[check[i,]$Suggested.Symbol] > 1){
          check[i,]$Suggested.Symbol <- check[i,]$x
          #print(check[i,]$x)
        }
      }
    }
  }
  return(check)
}

# Correct for repeated gene symbols
check1_NoOverlap <- NonOverlappingReplacement(check1)
check2_NoOverlap <- NonOverlappingReplacement(check2)
check3_NoOverlap <- NonOverlappingReplacement(check3)
check4_NoOverlap <- NonOverlappingReplacement(check4)
check5_NoOverlap <- NonOverlappingReplacement(check5)
check6_NoOverlap <- NonOverlappingReplacement(check6)
check7_NoOverlap <- NonOverlappingReplacement(check7)

check_kazu_NoOverlap <- NonOverlappingReplacement(check_kazu)

# Update Gene Names
SRR13606301_filtered_mat@Dimnames[[1]] <- check1_NoOverlap$Suggested.Symbol
SRR13606302_filtered_mat@Dimnames[[1]] <- check2_NoOverlap$Suggested.Symbol
SRR13606303_filtered_mat@Dimnames[[1]] <- check3_NoOverlap$Suggested.Symbol
SRR13606304_filtered_mat@Dimnames[[1]] <- check4_NoOverlap$Suggested.Symbol
SRR13606305_filtered_mat@Dimnames[[1]] <- check5_NoOverlap$Suggested.Symbol
SRR13606306_filtered_mat@Dimnames[[1]] <- check6_NoOverlap$Suggested.Symbol
SRR13606307_filtered_mat@Dimnames[[1]] <- check7_NoOverlap$Suggested.Symbol
KazuData_filtered_mat@Dimnames[[1]] <- check_kazu_NoOverlap$Suggested.Symbol

# Check Intersection
length(intersect(SRR13606301_filtered_mat@Dimnames[[1]],KazuData_filtered_mat@Dimnames[[1]]))

# Check dimensions
dim(SRR13606301_filtered_mat)
dim(SRR13606302_filtered_mat)
dim(SRR13606303_filtered_mat)
dim(SRR13606304_filtered_mat)
dim(SRR13606305_filtered_mat)
dim(SRR13606306_filtered_mat)
dim(SRR13606307_filtered_mat)
dim(KazuData_filtered_mat)

# Make Seurat Objects
seurat_01 <- CreateSeuratObject(SRR13606301_filtered_mat, project = "TGFB_Timecourse", min.cells = 3, min.features = 200)
seurat_02 <- CreateSeuratObject(SRR13606302_filtered_mat, project = "TGFB_Timecourse", min.cells = 3, min.features = 200)
seurat_03 <- CreateSeuratObject(SRR13606303_filtered_mat, project = "TGFB_Timecourse", min.cells = 3, min.features = 200)
seurat_04 <- CreateSeuratObject(SRR13606304_filtered_mat, project = "TGFB_Timecourse", min.cells = 3, min.features = 200)
seurat_05 <- CreateSeuratObject(SRR13606305_filtered_mat, project = "TGFB_Timecourse", min.cells = 3, min.features = 200)
seurat_06 <- CreateSeuratObject(SRR13606306_filtered_mat, project = "TGFB_Timecourse", min.cells = 3, min.features = 200)
seurat_07 <- CreateSeuratObject(SRR13606307_filtered_mat, project = "TGFB_Timecourse", min.cells = 3, min.features = 200)
seurat_Kazu <- CreateSeuratObject(KazuData_filtered_mat, project = "Kazu_Concentration", min.cells = 3, min.features = 200)

# Annotate
# see: https://www.ncbi.nlm.nih.gov/sra?LinkName=bioproject_sra_all&from_uid=698642
seurat_01[["Time"]] <- 8
seurat_01[["Batch"]] <- "A"

seurat_02[["Time"]] <- 4
seurat_02[["Batch"]] <- "A"

seurat_03[["Time"]] <- 3
seurat_03[["Batch"]] <- "B"

seurat_04[["Time"]] <- 2
seurat_04[["Batch"]] <- "B"

seurat_05[["Time"]] <- 1
seurat_05[["Batch"]] <- "B"

seurat_06[["Time"]] <- 0
seurat_06[["Batch"]] <- "A"

seurat_07[["Time"]] <- 0
seurat_07[["Batch"]] <- "B"

seurat_Kazu[["Time"]] <- 100 # Placerholder
seurat_Kazu[["Batch"]] <- "C"


# Merge Common Batches
seurat_A <- merge(seurat_01, y = c(seurat_02 ,seurat_06), add.cell.ids = c("SRR13606301", "SRR13606302", "SRR13606306"), project = "TGFB_timecourse_A")
seurat_B <- merge(seurat_03, y = c(seurat_04,seurat_05,seurat_07), add.cell.ids = c("SRR13606303", "SRR13606304", "SRR13606305", "SRR13606307"), project = "TGFB_timecourse_B")
seurat_C <- seurat_Kazu

dim(seurat_A)
dim(seurat_B)
dim(seurat_C)

# Check Feature Distribution
FeatureScatter(seurat_A, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
summary(seurat_A@meta.data$nFeature_RNA)
# Filter @ 1.5*IRQ
iqr_A <- 6019 + (6019-4880)*1.5
# Check Percentile
sum(seurat_A@meta.data$nFeature_RNA < iqr_A)/length(seurat_A@meta.data$nFeature_RNA)
# Try MAD
mad_A <- stats::mad(seurat_A@meta.data$nFeature_RNA)
mad_thresh_A <- median(seurat_A@meta.data$nFeature_RNA) + 3*mad_A
sum(seurat_A@meta.data$nFeature_RNA < mad_thresh_A)/length(seurat_A@meta.data$nFeature_RNA)

FeatureScatter(seurat_B, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
summary(seurat_B@meta.data$nFeature_RNA)
# Filter @ 1.5*IRQ
iqr_B <- 4667 + (4667-2525)*1.5
# Check percentile
sum(seurat_B@meta.data$nFeature_RNA < iqr_B)/length(seurat_B@meta.data$nFeature_RNA)
# Try MAD
mad_B <- stats::mad(seurat_B@meta.data$nFeature_RNA)
mad_thresh_B <- median(seurat_B@meta.data$nFeature_RNA) + 3*mad_B
sum(seurat_B@meta.data$nFeature_RNA < mad_thresh_B)/length(seurat_B@meta.data$nFeature_RNA)

FeatureScatter(seurat_C, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
summary(seurat_C@meta.data$nFeature_RNA)
# Filter @ 1.5*IRQ
iqr_C <- 4105 + (4105-3027)*1.5
# Check percentile
sum(seurat_C@meta.data$nFeature_RNA < iqr_C)/length(seurat_C@meta.data$nFeature_RNA)
# Try MAD
mad_C <- stats::mad(seurat_C@meta.data$nFeature_RNA)
mad_thresh_C <- median(seurat_C@meta.data$nFeature_RNA) + 3*mad_C
sum(seurat_C@meta.data$nFeature_RNA < mad_thresh_C)/length(seurat_C@meta.data$nFeature_RNA)

# Going to Use MAD threshold. See:
# https://www.biostars.org/p/446554/
# https://www.biorxiv.org/content/10.1101/2020.02.02.930578v3.full

# Filter
# These features are the same across all data sets as we would hope
mito.features <- grep(pattern = "^MT-", x = rownames(seurat_A), value = TRUE)

# A
percent.mito <- Matrix::colSums(GetAssayData(seurat_A, slot = 'counts')[mito.features, ]) / Matrix::colSums(GetAssayData(seurat_A, slot = 'counts'))
seurat_A[['percent.mito']] <- percent.mito
summary(seurat_A[['percent.mito']])

seurat_A <- subset(seurat_A, subset = nFeature_RNA > 500 & nFeature_RNA < mad_thresh_A & percent.mito < 0.2)
dim(seurat_A)

# B
percent.mito <- Matrix::colSums(GetAssayData(seurat_B, slot = 'counts')[mito.features, ]) / Matrix::colSums(GetAssayData(seurat_B, slot = 'counts'))
seurat_B[['percent.mito']] <- percent.mito
summary(seurat_B[['percent.mito']])

seurat_B <- subset(seurat_B, subset = nFeature_RNA > 500 & nFeature_RNA < mad_thresh_B & percent.mito < 0.2)
dim(seurat_B)

# C
#mito.features <- grep(pattern = "^MT-", x = rownames(seurat_C), value = TRUE)
percent.mito <- Matrix::colSums(GetAssayData(seurat_C, slot = 'counts')[mito.features, ]) / Matrix::colSums(GetAssayData(seurat_C, slot = 'counts'))
seurat_C[['percent.mito']] <- percent.mito
summary(seurat_C[['percent.mito']])

seurat_C <- subset(seurat_C, subset = nFeature_RNA > 500 & nFeature_RNA < mad_thresh_C & percent.mito < 0.2)
dim(seurat_C)

# Merge All Data
seurat_list <- list(seurat_A,seurat_B,seurat_C)

# Cell Cycle Prep
cc.genes <- readLines(con = "regev_lab_cell_cycle_genes.txt")
cc.genes_check <- checkGeneSymbols(cc.genes)
cc.genes <- cc.genes_check$Suggested.Symbol
#Split these genes into S markers and G2M markers
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

# Integrated Time data as in previous paper 
# seurat_list and identify variable features for each dataset independently
seurat_list <- lapply(X = seurat_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 15000)
  x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
})

# select features that are repeatedly variable across datasets for integration
vst_features <- SelectIntegrationFeatures(object.list = seurat_list,nfeatures = 15000)

integration.anchors <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = vst_features,normalization.method="LogNormalize")
integration.combined <- IntegrateData(anchorset = integration.anchors)

integration.combined <- ScaleData(integration.combined, vars.to.regress=c("S.Score","G2M.Score") ,verbose = FALSE)
integration.combined <- RunPCA(integration.combined, npcs = 15, verbose = FALSE)
integration.combined <- RunUMAP(integration.combined, reduction = "pca", dims = 1:15)
integration.combined <- FindNeighbors(integration.combined, reduction = "pca", dims = 1:15)
integration.combined <- FindClusters(integration.combined, resolution = 0.5)

##########################
### Add GBC Annotation ###

GBC_annot_columns <- read.csv("GBCannot/KazuConcentration_GBC_annotation_Samples.csv")
table(GBC_annot_columns$GBC_annot)
table(GBC_annot_columns$Values)
# 1    2    3    4    8    9   12   13 
# 1337 1422 1000 1020 1172 1020  862 1049

metadata <- integration.combined@meta.data
metadata$GBC <- "GBC0"
metadata$GBC_values <- 0

for (index in 1: dim(GBC_annot_columns)[1]){
  sample <- GBC_annot_columns[index,1]
  GBC <- GBC_annot_columns[index,2]
  GBC_value <- GBC_annot_columns[index,3]
  
  metadata[sample,]$GBC <- GBC
  metadata[sample,]$GBC_values <- GBC_value
}

table(metadata$GBC)
table(metadata$GBC_values)

# Check
check = matrix(, nrow = 8882, ncol = 4)
for (row in 1:8882){
  name = GBC_annot_columns[row,1]
  value = GBC_annot_columns[row,3]
  check[row,1] = name
  check[row,2] = value
  check[row,3] = metadata[name,]$GBC_values 
}
sum(check[,2] == check[,3])

# Concentration
#GBC1: 0
#GBC2: 12.5
#GBC3: 25
#GBC4: 50
#GBC8: 100
#GBC9: 200
#GBC12: 400
#GBC13: 800

metadata$GBC_pM <- metadata$GBC_values
metadata[metadata$GBC_pM == 0,]$GBC_pM <- -100
metadata[metadata$GBC_pM == 1,]$GBC_pM <- 0
metadata[metadata$GBC_pM == 2,]$GBC_pM <- 12.5
metadata[metadata$GBC_pM == 3,]$GBC_pM <- 25
metadata[metadata$GBC_pM == 4,]$GBC_pM <- 50
metadata[metadata$GBC_pM == 8,]$GBC_pM <- 100
metadata[metadata$GBC_pM == 9,]$GBC_pM <- 200
metadata[metadata$GBC_pM == 12,]$GBC_pM <- 400
metadata[metadata$GBC_pM == 13,]$GBC_pM <- 800

integration.combined@meta.data <- metadata

##########################

saveRDS(integration.combined,'SeuratData/Seurat_integrated_TimeAndCon_CellCycle_hg38.rds')
