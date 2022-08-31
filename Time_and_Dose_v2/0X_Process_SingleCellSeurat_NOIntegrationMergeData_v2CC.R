### Folder
setwd("/home/panchy/DoseTIme/_fullcorrection")

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
currentmap <- getCurrentHumanMap()
check1 <- checkGeneSymbols(SRR13606301_filtered_mat@Dimnames[[1]],unmapped.as.na = FALSE,species='human',map=currentmap)
check2 <- checkGeneSymbols(SRR13606302_filtered_mat@Dimnames[[1]],unmapped.as.na = FALSE,species='human',map=currentmap)
check3 <- checkGeneSymbols(SRR13606303_filtered_mat@Dimnames[[1]],unmapped.as.na = FALSE,species='human',map=currentmap)
check4 <- checkGeneSymbols(SRR13606304_filtered_mat@Dimnames[[1]],unmapped.as.na = FALSE,species='human',map=currentmap)
check5 <- checkGeneSymbols(SRR13606305_filtered_mat@Dimnames[[1]],unmapped.as.na = FALSE,species='human',map=currentmap)
check6 <- checkGeneSymbols(SRR13606306_filtered_mat@Dimnames[[1]],unmapped.as.na = FALSE,species='human',map=currentmap)
check7 <- checkGeneSymbols(SRR13606307_filtered_mat@Dimnames[[1]],unmapped.as.na = FALSE,species='human',map=currentmap)
check_kazu <- checkGeneSymbols(KazuData_filtered_mat@Dimnames[[1]],unmapped.as.na = FALSE,species='human',map=currentmap)

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

print(dim(seurat_01))
print(dim(seurat_02))
print(dim(seurat_03))
print(dim(seurat_04))
print(dim(seurat_05))
print(dim(seurat_06))
print(dim(seurat_07))
print(dim(seurat_Kazu))

individual_samples <- list(seurat_01,seurat_02,seurat_03,seurat_04,seurat_05,seurat_06,seurat_07,seurat_Kazu)

# Filter

# Going to Use MAD threshold. See:
# https://www.biostars.org/p/446554/
# https://www.biorxiv.org/content/10.1101/2020.02.02.930578v3.full

for (i in 1:length(individual_samples)){
  seurat_sample <- individual_samples[[i]]
  
  # MAD thresholding for features 
  mad_A <- stats::mad(seurat_sample@meta.data$nFeature_RNA)
  mad_thresh_A <- median(seurat_sample@meta.data$nFeature_RNA) + 3*mad_A
  print(sum(seurat_sample@meta.data$nFeature_RNA < mad_thresh_A)/length(seurat_sample@meta.data$nFeature_RNA))
  
  # Mito genes
  # These features are the same across all data sets as we would hope
  mito.features <- grep(pattern = "^MT-", x = rownames(seurat_sample), value = TRUE)
  percent.mito <- Matrix::colSums(GetAssayData(seurat_sample, slot = 'counts')[mito.features, ]) / Matrix::colSums(GetAssayData(seurat_sample, slot = 'counts'))
  seurat_sample[['percent.mito']] <- percent.mito
  summary(seurat_sample[['percent.mito']])
  
  seurat_sample <- subset(seurat_sample, subset = nFeature_RNA > 500 & nFeature_RNA < mad_thresh_A & percent.mito < 0.2)
  print(dim(seurat_sample))
  individual_samples[[i]] <- seurat_sample
}

seurat_01 <- individual_samples[[1]]
seurat_02 <- individual_samples[[2]] 
seurat_03 <- individual_samples[[3]] 
seurat_04 <- individual_samples[[4]] 
seurat_05 <- individual_samples[[5]] 
seurat_06 <- individual_samples[[6]] 
seurat_07 <- individual_samples[[7]] 
seurat_Kazu <- individual_samples[[8]] 

print(dim(seurat_01))
print(dim(seurat_02))
print(dim(seurat_03))
print(dim(seurat_04))
print(dim(seurat_05))
print(dim(seurat_06))
print(dim(seurat_07))
print(dim(seurat_Kazu))

# Merge Common Batches
seurat_A <- merge(seurat_01, y = c(seurat_02 ,seurat_06), add.cell.ids = c("SRR13606301", "SRR13606302", "SRR13606306"), project = "TGFB_timecourse_A")
seurat_B <- merge(seurat_03, y = c(seurat_04,seurat_05,seurat_07), add.cell.ids = c("SRR13606303", "SRR13606304", "SRR13606305", "SRR13606307"), project = "TGFB_timecourse_B")
seurat_C <- seurat_Kazu

dim(seurat_A)
dim(seurat_B)
dim(seurat_C)

############################
### Only Merge Time Data ###
############################

#seurat_time_list <- list(seurat_A,seurat_B)

# Integrated Time data as in previous paper 
# seurat_list and identify variable features for each dataset independently
#seurat_time_list <- lapply(X = seurat_time_list, FUN = function(x) {
#  x <- NormalizeData(x)
#  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 15000)
#  x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
#})

# select features that are repeatedly variable across datasets for integration
#vst_features <- SelectIntegrationFeatures(object.list = seurat_time_list,nfeatures = 15000)

#integration_time.anchors <- FindIntegrationAnchors(object.list = seurat_time_list, anchor.features = vst_features,normalization.method="LogNormalize")
#integration_time.combined <- IntegrateData(anchorset = integration_time.anchors)

#integration_time.combined <- ScaleData(integration_time.combined, vars.to.regress=c("nCount_RNA","percent.mito","S.Score","G2M.Score") ,verbose = FALSE)
#integration_time.combined <- RunPCA(integration_time.combined, npcs = 15, verbose = FALSE)
#integration_time.combined <- RunUMAP(integration_time.combined, reduction = "pca", dims = 1:15)
#integration_time.combined <- FindNeighbors(integration_time.combined, reduction = "pca", dims = 1:15)
#integration_time.combined <- FindClusters(integration_time.combined, resolution = 0.5)

### Save Integrated Time Data ###

#saveRDS(integration_time.combined,'SeuratData/Seurat_integrated_TimeOnly_FullCorr_hg38.rds')

##################################
### Correct Concentration Data ###
##################################

#seurat_C <- NormalizeData(seurat_C)
#seurat_C <- FindVariableFeatures(seurat_C, selection.method = "vst", nfeatures = 15000)
#seurat_C <- CellCycleScoring(seurat_C, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
#seurat_concentration_scale <- ScaleData(seurat_C,vars.to.regress=c("nCount_RNA","percent.mito","S.Score","G2M.Score"),verbose=FALSE)

#seurat_concentration_scale <- RunPCA(seurat_concentration_scale, npcs = 15, verbose = FALSE)
#seurat_concentration_scale <- RunUMAP(seurat_concentration_scale, reduction = "pca", dims = 1:15)
#seurat_concentration_scale <- FindNeighbors(seurat_concentration_scale, reduction = "pca", dims = 1:15)
#seurat_concentration_scale <- FindClusters(seurat_concentration_scale, resolution = 0.5)

### Add GBC Annotation ###

#GBC_annot_columns <- read.csv("GBCannot/KazuConcentration_GBC_annotation_Samples.csv")
#table(GBC_annot_columns$GBC_annot)
#table(GBC_annot_columns$Values)
# 1    2    3    4    8    9   12   13 
# 1337 1422 1000 1020 1172 1020  862 1049


#metadata <- seurat_concentration_scale@meta.data
#metadata$GBC <- "GBC0"
#metadata$GBC_values <- 0

#for (index in 1: dim(GBC_annot_columns)[1]){
#  sample <- GBC_annot_columns[index,1]
#  GBC <- GBC_annot_columns[index,2]
#  GBC_value <- GBC_annot_columns[index,3]
  
#  if (sample %in% row.names(metadata)){
#    metadata[sample,]$GBC <- GBC
#    metadata[sample,]$GBC_values <- GBC_value
#  }
#}

#table(metadata$GBC)
#table(metadata$GBC_values)

# Check
#check = matrix(, nrow = 8882, ncol = 4)
#for (row in 1:8882){
#  name = GBC_annot_columns[row,1]
#  value = GBC_annot_columns[row,3]
#  check[row,1] = name
#  check[row,2] = value
#  check[row,3] = metadata[name,]$GBC_values 
#}
#sum(check[,2] == check[,3])

# Concentration
#GBC1: 0
#GBC2: 12.5
#GBC3: 25
#GBC4: 50
#GBC8: 100
#GBC9: 200
#GBC12: 400
#GBC13: 800

#metadata$GBC_pM <- metadata$GBC_values
#metadata[metadata$GBC_pM == 0,]$GBC_pM <- -100
#metadata[metadata$GBC_pM == 1,]$GBC_pM <- 0
#metadata[metadata$GBC_pM == 2,]$GBC_pM <- 12.5
#metadata[metadata$GBC_pM == 3,]$GBC_pM <- 25
#metadata[metadata$GBC_pM == 4,]$GBC_pM <- 50
#metadata[metadata$GBC_pM == 8,]$GBC_pM <- 100
#metadata[metadata$GBC_pM == 9,]$GBC_pM <- 200
#metadata[metadata$GBC_pM == 12,]$GBC_pM <- 400
#metadata[metadata$GBC_pM == 13,]$GBC_pM <- 800

#seurat_concentration_scale@meta.data <- metadata

### Save Dose Data ###

#saveRDS(seurat_concentration_scale,"SeuratData/Seurat_concentration_DoseOnly_FullCorr_hg38.rds")

############################
##### FULL INTEGRATION #####
############################

# Cell Cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

s_check <- checkGeneSymbols(s.genes,unmapped.as.na = FALSE,species='human',map=currentmap)
g2m_check <- checkGeneSymbols(g2m.genes,unmapped.as.na = FALSE,species='human',map=currentmap)

s.genes <- s_check$Suggested.Symbol
g2m.genes <- g2m_check$Suggested.Symbol

print(length(intersect(s.genes,row.names(seurat_A@assays$RNA@data))))
print(length(intersect(g2m.genes,row.names(seurat_A@assays$RNA@data))))

# Merge All Data
seurat_all <- merge(seurat_A, y=c(seurat_B,seurat_C))

# Integrated Time data as in previous paper 
# seurat_list and identify variable features for each dataset independently
seurat_all <- NormalizeData(seurat_all)
seurat_all <- FindVariableFeatures(seurat_all, selection.method = "vst", nfeatures = 15000)
seurat_all <- CellCycleScoring(seurat_all, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

integration.combined <- ScaleData(seurat_all, vars.to.regress=c("nCount_RNA","percent.mito","S.Score","G2M.Score") ,verbose = FALSE)
integration.combined <- RunPCA(integration.combined, npcs = 15, verbose = FALSE)
integration.combined <- RunUMAP(integration.combined, reduction = "pca", dims = 1:15)
integration.combined <- FindNeighbors(integration.combined, reduction = "pca", dims = 1:15)
integration.combined <- FindClusters(integration.combined, resolution = 0.5)

### Save Fully Integrated Data ###

saveRDS(integration.combined,'SeuratData/Seurat_NoIntegration_Combined_FullCorreciton_hg38_v2CC.rds')

print("Done")
