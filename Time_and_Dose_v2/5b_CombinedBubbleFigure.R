########################
### Setup Enviroment ###
########################

### Folder
setwd("/Users/nicholaspanchy/Documents/Work_UTK/DoseTime_FullCorr_FromScratch//")

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

# Bubble Data
Con_Samples <- dose_samples
Time_Samples <- time_samples
Labeled_Samples <- c(Time_Samples,Con_Samples)

Bounds_IQR_med <- readRDS("Bounds_IQR_med.rds")
odd_mat_long <- readRDS("BubbleOddsMat.rds")
quart_clusters <- readRDS("EMQuartileClusters.rds")

####################
### Bubble Plots ###
####################

integated_meta_labeled<- integated_meta[Labeled_Samples,]
integated_meta_labeled$label <- 0
integated_meta_labeled[quart_clusters[[1]],]$label <- 1
integated_meta_labeled[quart_clusters[[2]],]$label <- 2
integated_meta_labeled[quart_clusters[[3]],]$label <- 3
integated_meta_labeled[quart_clusters[[4]],]$label <- 4
integated_meta_labeled[quart_clusters[[5]],]$label <- 5
integated_meta_labeled[quart_clusters[[6]],]$label <- 6
integated_meta_labeled[quart_clusters[[7]],]$label <- 7
integated_meta_labeled[quart_clusters[[8]],]$label <- 8
integated_meta_labeled[quart_clusters[[9]],]$label <- 9
integated_meta_labeled[quart_clusters[[10]],]$label <- 10
integated_meta_labeled[quart_clusters[[11]],]$label <- 11
integated_meta_labeled[quart_clusters[[12]],]$label <- 12
integated_meta_labeled[quart_clusters[[13]],]$label <- 13
integated_meta_labeled[quart_clusters[[14]],]$label <- 14
integated_meta_labeled[quart_clusters[[15]],]$label <- 15
integated_meta_labeled[quart_clusters[[16]],]$label <- 16

table(integated_meta_labeled$label)

### Contour Plot ### 
UMAP_annot <- merge(UMAP_values,integated_meta,by=0)
row.names(UMAP_annot) <- UMAP_annot$Row.names
UMAP_annot[UMAP_annot$Batch=="A",]$Batch <- "Time"
UMAP_annot[UMAP_annot$Batch=="B",]$Batch  <- "Time"
UMAP_annot[UMAP_annot$Batch=="C",]$Batch  <- "Dose"

UMAP_samples <- row.names(UMAP_annot)
UMAP_annot$E_PC1 <- Integrated_E_mat_nsprcomp$x[UMAP_samples,1]
UMAP_annot$M_PC1 <- Integrated_M_mat_nsprcomp$x[UMAP_samples,1]

ggplot(UMAP_annot, aes(x = E_PC1, y = M_PC1, colour=Batch)) +
  geom_density_2d(bins=10,alpha=0.5,size=1) + 
  xlab("E score") + ylab("M score") + labs(fill="Log2FC") +
  theme_bw() + theme(axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1)


### Base Bubble Plot ###
bubble_p <- ggplot(odd_mat_long, aes(x = X, y = Y,size=Count,fill=Odds)) + 
  geom_point(alpha=1.0, shape=21, color="black") + ggtitle("Time vs. Dose") +
  scale_size(range = c(5, 15)) + xlab("E Quartile") + ylab("M Quartile") + labs(fill="Log2Odds") +
  scale_fill_gradient2(low='blue',high='red',mid='white',midpoint=0,limits=c(-3,3),oob=squish) +
  theme_bw() + coord_fixed(ratio = 1) + ylim(c(0.75,4.25)) + xlim(c(0.75,4.25)) +
  theme(axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1)
bubble_p

# Read Marker Data

AllClusterMarkers <- readRDS("AllClusterMarkers_EMTGenes_AllResults.rds")


odd_mat_long$VIM <- NA
odd_mat_long$FN1 <- NA
odd_mat_long$EPCAM <- NA
odd_mat_long$CDH1 <- NA

for(i in 1:16){
  AllClusterMarkers_CurrentCluster <- AllClusterMarkers[AllClusterMarkers$cluster == i,]
  if("VIM" %in% AllClusterMarkers_CurrentCluster$gene){
    VIM_results <- AllClusterMarkers_CurrentCluster[AllClusterMarkers_CurrentCluster$gene=="VIM",]
    odd_mat_long[i,]$VIM <- VIM_results$avg_log2FC
  }
  if("FN1" %in% AllClusterMarkers_CurrentCluster$gene){
    FN1_results <- AllClusterMarkers_CurrentCluster[AllClusterMarkers_CurrentCluster$gene=="FN1",]
    odd_mat_long[i,]$FN1 <- FN1_results$avg_log2FC
  }
  if("CDH1" %in% AllClusterMarkers_CurrentCluster$gene){
    CDH1_results <- AllClusterMarkers_CurrentCluster[AllClusterMarkers_CurrentCluster$gene=="CDH1",]
    odd_mat_long[i,]$CDH1 <- CDH1_results$avg_log2FC
  }
  if("EPCAM" %in% AllClusterMarkers_CurrentCluster$gene){
    EPCAM_results <- AllClusterMarkers_CurrentCluster[AllClusterMarkers_CurrentCluster$gene=="EPCAM",]
    odd_mat_long[i,]$EPCAM <- EPCAM_results$avg_log2FC
  }
}

bubble_VIM <- ggplot(odd_mat_long, aes(x = X, y = Y,size=Count,fill=VIM)) + 
  geom_point(alpha=1.0, shape=21, color="black") + ggtitle("VIM") +
  scale_size(range = c(5, 15)) + xlab("E Quartile") + ylab("M Quartile") + labs(fill="Log2FC") +
  scale_fill_gradient2(low='blue',high='red',mid='white',midpoint=0,limits=c(-1,1),oob=squish) +
  theme_bw() + coord_fixed(ratio = 1) + ylim(c(0.75,4.25)) + xlim(c(0.75,4.25)) +
  theme(axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1)
bubble_VIM

bubble_FN1<- ggplot(odd_mat_long, aes(x = X, y = Y,size=Count,fill=FN1)) + 
  geom_point(alpha=1.0, shape=21, color="black") + ggtitle("FN1") +
  scale_size(range = c(5, 15)) + xlab("E Quartile") + ylab("M Quartile") + labs(color="Log2FC") +
  scale_fill_gradient2(low='blue',high='red',mid='white',midpoint=0,limits=c(-3,2),oob=squish) +
  theme_bw() + coord_fixed(ratio = 1) + ylim(c(0.75,4.25)) + xlim(c(0.75,4.25)) +
  theme(axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1)
bubble_FN1

bubble_EPCAM <- ggplot(odd_mat_long, aes(x = X, y = Y,size=Count,fill=EPCAM)) + 
  geom_point(alpha=1.0, shape=21, color="black") + ggtitle("EPCAM") +
  scale_size(range = c(5, 15)) + xlab("E Quartile") + ylab("M Quartile") + labs(color="Log2FC") +
  scale_fill_gradient2(low='blue',high='red',mid='white',midpoint=0,limits=c(-1,1),oob=squish) +
  theme_bw() + coord_fixed(ratio = 1) + ylim(c(0.75,4.25)) + xlim(c(0.75,4.25)) +
  theme(axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1)
bubble_EPCAM

bubble_CDH1 <- ggplot(odd_mat_long, aes(x = X, y = Y,size=Count,fill=CDH1)) + 
  geom_point(alpha=1.0, shape=21, color="black") + ggtitle("CDH1") +
  scale_size(range = c(5, 15)) + xlab("E Quartile") + ylab("M Quartile") + labs(color="Log2FC") +
  scale_fill_gradient2(low='blue',high='red',mid='white',midpoint=0,limits=c(-1,1),oob=squish) +
  theme_bw() + coord_fixed(ratio = 1) + ylim(c(0.75,4.25)) + xlim(c(0.75,4.25)) +
  theme(axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1)
bubble_CDH1

### Nearest Neighbor Clusters by Quart Clusters ###
table(integated_meta_labeled$label,integated_meta_labeled$seurat_clusters)
table(integated_meta_labeled$Time,integated_meta_labeled$seurat_clusters)
table(integated_meta_labeled$GBC_pM,integated_meta_labeled$seurat_clusters)

Quart_by_NN <- table(integated_meta_labeled$label,integated_meta_labeled$seurat_clusters)
Quart_by_NN_Per <- sweep(Quart_by_NN,2,colSums(Quart_by_NN),`/`)
Quart_by_NN_PerT <- sweep(t(Quart_by_NN),2,colSums(t(Quart_by_NN)),`/`)

odd_mat_long$Cluster_1 <- -1
for(i in 1:16){
  odd_mat_long[i,]$Cluster_1 <- Quart_by_NN_Per[i,2] # Add 1 because starts at cluster 0
}

bubble_Cluster1 <- ggplot(odd_mat_long, aes(x = X, y = Y,size=Count,fill=Cluster_1)) + 
  geom_point(alpha=1.0, shape=21, color="black") +
  scale_size(range = c(5, 15)) +
  scale_fill_gradient2(low='black',high='red',mid='yellow',midpoint=0.1250,limits=c(0.0,0.5),oob=squish) +
  theme_bw() + coord_fixed(ratio = 1) + ylim(c(0.75,4.25)) + xlim(c(0.75,4.25))
bubble_Cluster1

odd_mat_long$Cluster_8 <- -1
for(i in 1:16){
  odd_mat_long[i,]$Cluster_8 <- Quart_by_NN_Per[i,9] # Add 1 because starts at cluster 0
}

bubble_Cluster8 <- ggplot(odd_mat_long, aes(x = X, y = Y,size=Count,fill=Cluster_8)) + 
  geom_point(alpha=1.0, shape=21, color="black") +
  scale_size(range = c(5, 15)) +
  scale_fill_gradient2(low='black',high='red',mid='yellow',midpoint=0.1250,limits=c(0.0,0.5),oob=squish) +
  theme_bw() + coord_fixed(ratio = 1) + ylim(c(0.75,4.25)) + xlim(c(0.75,4.25))
bubble_Cluster8

odd_mat_long$Cluster_3 <- -1
for(i in 1:16){
  odd_mat_long[i,]$Cluster_3 <- Quart_by_NN_Per[i,4] # Add 1 because starts at cluster 0
}

### Middle 50% of E and M Scores ###

MiddleMarkers <- readRDS("MiddleMarkers_EMTGenes_AllResults.rds")
# Haven't re-corrected these, but P-values are well in affected of the change which
# is roughly dividing by 16