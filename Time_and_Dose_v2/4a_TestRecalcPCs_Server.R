########################
### Setup Enviroment ###
########################

### Folder
#setwd("/Users/nicholaspanchy/Documents/Work_UTK/DoseTime_FullCorr_FromScratch//")
setwd("/home/panchy/DoseTIme/_fromscratch_server")

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

###########################
### Load Processed Data ###
###########################

integated_meta <- readRDS("Integrated_MetaData.rds")
UMAP_values <- readRDS("Integrated_UMAP_values.rds")

dose_samples <- readRDS("DoseSamples.rds")
time_samples <- readRDS("TimeSamples.rds")

Integrated_E_mat_nsprcomp <- readRDS("Integrated_E_mat_nsprcomp.rds")
Integrated_M_mat_nsprcomp <- readRDS("Integrated_M_mat_nsprcomp.rds")

### Load EMT Data only ###
integration.combined.EMT.rds <- readRDS("integration.combined.EMT.rds")

# Load Genes Sets #
EMT_genes <- read.table("EMTGeneUpdate3.txt",header=TRUE)
E_genes <- subset(EMT_genes, EMT_genes$Annotation == 'E')
M_genes <- subset(EMT_genes, EMT_genes$Annotation == 'M')
Gene_Sets <- list('Egenes'=E_genes$Gene,'Mgenes'=M_genes$Gene)

### Try recalculating PC ###

recalcM_results <- matrix(nrow=length(M_genes$Gene),ncol=3)
index = 1
for(gene in M_genes$Gene){
  gene_order <- names(Integrated_M_mat_nsprcomp$rotation[,1])
  if (gene %in% gene_order){
    recalcM_results[index,1] <- gene
    recalcM_results[index,2] <- Integrated_M_mat_nsprcomp$rotation[gene,1]
    gene_order_drop <- setdiff(gene_order,gene)
    print(length(gene_order_drop))
    recalc_PC <- t(integration.combined.EMT.rds@assays$integrated@scale.data[gene_order_drop,]) %*% Integrated_M_mat_nsprcomp$rotation[gene_order_drop,1]
    recalc_PC = recalc_PC - mean(recalc_PC)
    recalcM_results[index,3] <- cor(Integrated_M_mat_nsprcomp$x[,1],recalc_PC)
  }
  else{
    recalcM_results[index,1] <- gene
    recalcM_results[index,2] <- NA
    recalcM_results[index,3] <- NA
  }
  index=index+1
}

recalcE_results <- matrix(nrow=length(E_genes$Gene),ncol=3)
index = 1
for(gene in E_genes$Gene){
  gene_order <- names(Integrated_E_mat_nsprcomp$rotation[,1])
  if (gene %in% gene_order){
    recalcE_results[index,1] <- gene
    recalcE_results[index,2] <- Integrated_E_mat_nsprcomp$rotation[gene,1]
    gene_order_drop <- setdiff(gene_order,gene)
    print(length(gene_order_drop))
    recalc_PC <- t(integration.combined.EMT.rds@assays$integrated@scale.data[gene_order_drop,]) %*% Integrated_E_mat_nsprcomp$rotation[gene_order_drop,1]
    recalc_PC = recalc_PC - mean(recalc_PC)
    recalcE_results[index,3] <- cor(Integrated_E_mat_nsprcomp$x[,1],recalc_PC)
  }
  else{
    recalcE_results[index,1] <- gene
    recalcE_results[index,2] <- NA
    recalcE_results[index,3] <- NA
  }
  index=index+1
}

### Try recreating PCs without individual genes ###
GetBoundarySamples <- function(Scores,row,column,bonds_df){
  
  # Select samples based on score
  if (row == 1){ # If first values, get less than
    samples <- names(Scores[Scores < bonds_df[row,column]])
  }
  else if (row == dim(bonds_df)[1]+1){ # if beyond the last value, > last value
    samples <- names(Scores[Scores >= bonds_df[row-1,column]])
  }
  else { # else between current value and previous value
    tmp <- Scores[Scores >= bonds_df[row-1,column]]
    samples <- names(tmp[tmp< bonds_df[row,column]])
  }
  return(samples)
}

Con_Samples <- dose_samples
Time_Samples <- time_samples
ConAll_samples <- rownames(integated_meta[integated_meta$Batch == "C",])
length(Time_Samples)
length(Con_Samples)
length(ConAll_samples)
Labeled_Samples <- c(Time_Samples,Con_Samples)

library(nsprcomp)
integated_scale <- integration.combined.EMT.rds@assays$integrated@scale.data

present_M_genes <- names(Integrated_M_mat_nsprcomp$rotation[,1])
rerun_nnPCA_results_M <- matrix(,nrow=length(present_M_genes),ncol=2)
indexM = 1
for (gene in present_M_genes){

  Mgenes_drop_gene <- setdiff(M_genes$Gene,gene)

  Scaled_Integrated_data_M <- integated_scale[row.names(integated_scale) %in% Mgenes_drop_gene,]
  print(dim(Scaled_Integrated_data_M))
  set.seed(5)
  Integrated_M_mat_nsprcomp_drop <- nsprcomp(t(Scaled_Integrated_data_M),nneg=TRUE,ncomp=5,em_maxiter = 10000,em_tol = 0.00001)
  rerun_nnPCA_results_M[indexM,1] <- gene
  rerun_nnPCA_results_M[indexM,2] <- cor(Integrated_M_mat_nsprcomp_drop$x[,1],Integrated_M_mat_nsprcomp$x[,1])
  indexM=indexM+1
  
  print(gene)
  E_scores <- Integrated_E_mat_nsprcomp$x[,1]
  M_scores <- Integrated_M_mat_nsprcomp_drop$x[,1]
  
  Escore_IQR_med <- quantile(E_scores[Labeled_Samples],probs=c(0.25,0.50,0.75))
  Mscore_IQR_med <- quantile(M_scores[Labeled_Samples],probs=c(0.25,0.50,0.75))
  Bounds_IQR_med <- as.data.frame(cbind(Escore_IQR_med,Mscore_IQR_med))
  
  quart_clusters <- rep(NA, 16)
  index = 1  
  for(i in 1:(dim(Bounds_IQR_med)[1]+1)){
    for(j in 1:(dim(Bounds_IQR_med)[1]+1)){
      E_boundary_samples <- GetBoundarySamples(E_scores,i,1,Bounds_IQR_med)
      M_boundary_sampels <- GetBoundarySamples(M_scores,j,2,Bounds_IQR_med)
      intersect_samples = intersect(E_boundary_samples, M_boundary_sampels)
      filter_samples <- intersect(intersect_samples,union(Time_Samples,Con_Samples))
      quart_clusters[index] <- list(filter_samples)
      index = index + 1
    }
  }
  integration.combined.EMT.labeled <- subset(integration.combined.EMT.rds,cells=Labeled_Samples)
  
  integration.combined.EMT.labeled@meta.data$label <- 0
  integration.combined.EMT.labeled@meta.data[quart_clusters[[1]],]$label <- 1
  integration.combined.EMT.labeled@meta.data[quart_clusters[[2]],]$label <- 2
  integration.combined.EMT.labeled@meta.data[quart_clusters[[3]],]$label <- 3
  integration.combined.EMT.labeled@meta.data[quart_clusters[[4]],]$label <- 4
  integration.combined.EMT.labeled@meta.data[quart_clusters[[5]],]$label <- 5
  integration.combined.EMT.labeled@meta.data[quart_clusters[[6]],]$label <- 6
  integration.combined.EMT.labeled@meta.data[quart_clusters[[7]],]$label <- 7
  integration.combined.EMT.labeled@meta.data[quart_clusters[[8]],]$label <- 8
  integration.combined.EMT.labeled@meta.data[quart_clusters[[9]],]$label <- 9
  integration.combined.EMT.labeled@meta.data[quart_clusters[[10]],]$label <- 10
  integration.combined.EMT.labeled@meta.data[quart_clusters[[11]],]$label <- 11
  integration.combined.EMT.labeled@meta.data[quart_clusters[[12]],]$label <- 12
  integration.combined.EMT.labeled@meta.data[quart_clusters[[13]],]$label <- 13
  integration.combined.EMT.labeled@meta.data[quart_clusters[[14]],]$label <- 14
  integration.combined.EMT.labeled@meta.data[quart_clusters[[15]],]$label <- 15
  integration.combined.EMT.labeled@meta.data[quart_clusters[[16]],]$label <- 16
  DefaultAssay(integration.combined.EMT.labeled) <- "RNA"
  #table(integration.combined.EMT.labeled$label)
  
  integration.combined.EMT.labeled@active.ident <- as.factor(integration.combined.EMT.labeled$label)
  AllClusterMarkers_LR_QuartClusters <- FindAllMarkers(integration.combined.EMT.labeled,features=gene,test.use="LR",latent.vars=c("nCount_RNA","percent.mito","S.Score","G2M.Score","Batch"),logfc.threshold=0,min.pct=0,return.thresh=1)
  saveRDS(AllClusterMarkers_LR_QuartClusters,paste0("GeneDropoutM_",gene,"_LRTest.rds"))
  
  ### Alternate Quarts ###
  integration.combined.EMT.labeled.Time <- subset(integration.combined.EMT.labeled, subset = Batch %in% c("A","B"))
  integration.combined.EMT.labeled.Con <- subset(integration.combined.EMT.labeled, subset = Batch %in% c("C"))
  
  AllClusterMarkers_TimeOnly_Quart <- FindAllMarkers(integration.combined.EMT.labeled.Time,features=gene,test.use="LR",latent.vars=c("nCount_RNA","percent.mito","S.Score","G2M.Score","Batch"),logfc.threshold=0,min.pct=0,return.thresh=1)
  saveRDS(AllClusterMarkers_TimeOnly_Quart,paste0("GeneDropoutM_",gene,"_Time_LRTest.rds"))
  AllClusterMarkers_ConOnly_Quart <- FindAllMarkers(integration.combined.EMT.labeled.Con,features=gene,test.use="LR",latent.vars=c("nCount_RNA","percent.mito","S.Score","G2M.Score"),logfc.threshold=0,min.pct=0,return.thresh=1)
  saveRDS(AllClusterMarkers_ConOnly_Quart,paste0("GeneDropoutM_",gene,"_Dose_LRTest.rds"))
}
saveRDS(rerun_nnPCA_results_M,"rerun_nnPCA_results_M.rds")

present_E_genes <- names(Integrated_E_mat_nsprcomp$rotation[,1])
rerun_nnPCA_results_E <- matrix(,nrow=length(present_E_genes),ncol=2)
indexE = 1
for (gene in present_E_genes){
  
  Egenes_drop_gene <- setdiff(E_genes$Gene,gene)
  
  Scaled_Integrated_data_E <- integated_scale[row.names(integated_scale) %in% Egenes_drop_gene,]
  print(dim(Scaled_Integrated_data_E))
  set.seed(5)
  Integrated_E_mat_nsprcomp_drop <- nsprcomp(t(Scaled_Integrated_data_E),nneg=TRUE,ncomp=5,em_maxiter = 10000,em_tol = 0.00001)
  rerun_nnPCA_results_E[indexE,1] <- gene
  rerun_nnPCA_results_E[indexE,2] <- cor(Integrated_E_mat_nsprcomp_drop$x[,1],Integrated_E_mat_nsprcomp$x[,1])
  indexE=indexE+1
  
  print(gene)
  E_scores <- Integrated_E_mat_nsprcomp_drop$x[,1]
  M_scores <- Integrated_M_mat_nsprcomp$x[,1]
  
  Escore_IQR_med <- quantile(E_scores[Labeled_Samples],probs=c(0.25,0.50,0.75))
  Mscore_IQR_med <- quantile(M_scores[Labeled_Samples],probs=c(0.25,0.50,0.75))
  Bounds_IQR_med <- as.data.frame(cbind(Escore_IQR_med,Mscore_IQR_med))
  
  quart_clusters <- rep(NA, 16)
  index = 1  
  for(i in 1:(dim(Bounds_IQR_med)[1]+1)){
    for(j in 1:(dim(Bounds_IQR_med)[1]+1)){
      E_boundary_samples <- GetBoundarySamples(E_scores,i,1,Bounds_IQR_med)
      M_boundary_sampels <- GetBoundarySamples(M_scores,j,2,Bounds_IQR_med)
      intersect_samples = intersect(E_boundary_samples, M_boundary_sampels)
      filter_samples <- intersect(intersect_samples,union(Time_Samples,Con_Samples))
      quart_clusters[index] <- list(filter_samples)
      index = index + 1
    }
  }
  integration.combined.EMT.labeled <- subset(integration.combined.EMT.rds,cells=Labeled_Samples)
  
  integration.combined.EMT.labeled@meta.data$label <- 0
  integration.combined.EMT.labeled@meta.data[quart_clusters[[1]],]$label <- 1
  integration.combined.EMT.labeled@meta.data[quart_clusters[[2]],]$label <- 2
  integration.combined.EMT.labeled@meta.data[quart_clusters[[3]],]$label <- 3
  integration.combined.EMT.labeled@meta.data[quart_clusters[[4]],]$label <- 4
  integration.combined.EMT.labeled@meta.data[quart_clusters[[5]],]$label <- 5
  integration.combined.EMT.labeled@meta.data[quart_clusters[[6]],]$label <- 6
  integration.combined.EMT.labeled@meta.data[quart_clusters[[7]],]$label <- 7
  integration.combined.EMT.labeled@meta.data[quart_clusters[[8]],]$label <- 8
  integration.combined.EMT.labeled@meta.data[quart_clusters[[9]],]$label <- 9
  integration.combined.EMT.labeled@meta.data[quart_clusters[[10]],]$label <- 10
  integration.combined.EMT.labeled@meta.data[quart_clusters[[11]],]$label <- 11
  integration.combined.EMT.labeled@meta.data[quart_clusters[[12]],]$label <- 12
  integration.combined.EMT.labeled@meta.data[quart_clusters[[13]],]$label <- 13
  integration.combined.EMT.labeled@meta.data[quart_clusters[[14]],]$label <- 14
  integration.combined.EMT.labeled@meta.data[quart_clusters[[15]],]$label <- 15
  integration.combined.EMT.labeled@meta.data[quart_clusters[[16]],]$label <- 16
  DefaultAssay(integration.combined.EMT.labeled) <- "RNA"
  #table(integration.combined.EMT.labeled$label)
  
  integration.combined.EMT.labeled@active.ident <- as.factor(integration.combined.EMT.labeled$label)
  AllClusterMarkers_LR_QuartClusters <- FindAllMarkers(integration.combined.EMT.labeled,features=gene,test.use="LR",latent.vars=c("nCount_RNA","percent.mito","S.Score","G2M.Score","Batch"),logfc.threshold=0,min.pct=0,return.thresh=1)
  saveRDS(AllClusterMarkers_LR_QuartClusters,paste0("GeneDropoutE_",gene,"_LRTest.rds"))
  
  ### Alternate Quarts ###
  integration.combined.EMT.labeled.Time <- subset(integration.combined.EMT.labeled, subset = Batch %in% c("A","B"))
  integration.combined.EMT.labeled.Con <- subset(integration.combined.EMT.labeled, subset = Batch %in% c("C"))
  
  AllClusterMarkers_TimeOnly_Quart <- FindAllMarkers(integration.combined.EMT.labeled.Time,features=gene,test.use="LR",latent.vars=c("nCount_RNA","percent.mito","S.Score","G2M.Score","Batch"),logfc.threshold=0,min.pct=0,return.thresh=1)
  saveRDS(AllClusterMarkers_TimeOnly_Quart,paste0("GeneDropoutE_",gene,"_Time_LRTest.rds"))
  AllClusterMarkers_ConOnly_Quart <- FindAllMarkers(integration.combined.EMT.labeled.Con,features=gene,test.use="LR",latent.vars=c("nCount_RNA","percent.mito","S.Score","G2M.Score"),logfc.threshold=0,min.pct=0,return.thresh=1)  
  saveRDS(AllClusterMarkers_ConOnly_Quart,paste0("GeneDropoutE_",gene,"_Dose_LRTest.rds"))
}
saveRDS(rerun_nnPCA_results_E,"rerun_nnPCA_results_E.rds")
