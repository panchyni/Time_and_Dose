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
library(ggplot2)
library(scales)

set.seed(5) # Controls randomness for nnPC, minor variations , mainly to keep enriched gene lists consistent (+/-1 genes)

########################################################
### Integrated Data --- Seurat/Quad Clusters Overlap ###
########################################################

integration.combined.All <- readRDS('SeuratData/Seurat_integrated_TimeAndCon_CellCycle_hg38.rds')

### Fix annotation ###

integration.combined.All@meta.data <- integration.combined.All@meta.data[1:25109,]

# Check against original list
check_table <- read.table("GBCannot/CBC_GBC_summary.txt.alt_names",head=TRUE)
check_table$CBC_alt <- cbind(lapply(check_table$CBC.1, function(x) {sub(".1", "-1", x)}))
row.names(check_table) <- check_table$CBC_alt
overlap <- intersect(check_table$CBC_alt,row.names(integration.combined.All@meta.data))

check_table <- as.data.frame(check_table)
check_table$meta_annot <- 0
check_table[overlap,]$meta_annot <- integration.combined.All@meta.data[overlap ,]$GBC

table(check_table$meta_annot,check_table$GBC)

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

# Check correspondence post correction
table(integration.combined.All@meta.data$GBC,integration.combined.All@meta.data$GBC_values)
table(integration.combined.All@meta.data$GBC,integration.combined.All@meta.data$GBC_pM)

# Check missing samples
check_table_0 <- check_table[check_table$meta_annot ==0,]
dim(check_table_0)
kazu_barcodes <- read.table("data_KazuConcentration/barcodes.tsv")
length(intersect(kazu_barcodes$V1,check_table_0$CBC_alt))
# All but seven are missing from the filtered barcode set
# The missing seven are filtered (these are the 7 NA sample post stringent filtering)

p_clusters <- DimPlot(integration.combined.All, reduction = "umap", label = TRUE, repel = TRUE)

# Genes Sets #
EMT_genes <- read.table("EMTGenesUpdateV2.txt",header=TRUE)
E_genes <- subset(EMT_genes, EMT_genes$Annotation == 'E')
M_genes <- subset(EMT_genes, EMT_genes$Annotation == 'M')
Gene_Sets <- list('Egenes'=E_genes$Gene,'Mgenes'=M_genes$Gene)

seurat_combo_E <- ScaleData(integration.combined.All, ,vars.to.regress=c("S.Score","G2M.Score"),verbose = FALSE, features=E_genes$Gene)
seurat_combo_M <- ScaleData(integration.combined.All, ,vars.to.regress=c("S.Score","G2M.Score"), verbose = FALSE, features=M_genes$Gene)

seurat_combo_E_mat <- GetAssayData(object = seurat_combo_E, slot = "scale.data")
seurat_combo_M_mat <- GetAssayData(object = seurat_combo_M, slot = "scale.data")

dim(seurat_combo_E_mat)
dim(seurat_combo_M_mat)

### NN-PCA ###

library(nsprcomp)
set.seed(5)

# E
seurat_combo_E_mat_nsprcomp <- nsprcomp(t(seurat_combo_E_mat),nneg=TRUE,ncomp=5)
seurat_combo_E_mat_nsprcomp_annot <- cbind(seurat_combo_E_mat_nsprcomp$x,seurat_combo_E@meta.data)

Time_Samples <- row.names(seurat_combo_E_mat_nsprcomp_annot[seurat_combo_E_mat_nsprcomp_annot$Time < 100,])
ConGood_Samples <- row.names(seurat_combo_E_mat_nsprcomp_annot[seurat_combo_E_mat_nsprcomp_annot$GBC_pM >= 0,])


cor(seurat_combo_E_mat_nsprcomp_annot[Time_Samples,]$PC1,seurat_combo_E_mat_nsprcomp_annot[Time_Samples,]$Time)
cor(seurat_combo_E_mat_nsprcomp_annot[ConGood_Samples,]$PC1,seurat_combo_E_mat_nsprcomp_annot[ConGood_Samples,]$GBC_pM)

# M 
seurat_combo_M_mat_nsprcomp <- nsprcomp(t(seurat_combo_M_mat),nneg=TRUE,ncomp=5)
seurat_combo_M_mat_nsprcomp_annot <- cbind(seurat_combo_M_mat_nsprcomp$x,seurat_combo_M@meta.data)
cor(seurat_combo_M_mat_nsprcomp_annot[Time_Samples,]$PC1,seurat_combo_M_mat_nsprcomp_annot[Time_Samples,]$Time)
cor(seurat_combo_M_mat_nsprcomp_annot[ConGood_Samples,]$PC1,seurat_combo_M_mat_nsprcomp_annot[ConGood_Samples,]$GBC_pM)

model <- lm(seurat_combo_E_mat_nsprcomp_annot[Time_Samples,]$Time ~ seurat_combo_E_mat_nsprcomp_annot[Time_Samples,]$PC1 + seurat_combo_M_mat_nsprcomp_annot[Time_Samples,]$PC1)
summary(model)

model <- lm(seurat_combo_E_mat_nsprcomp_annot[ConGood_Samples,]$GBC_pM ~ seurat_combo_E_mat_nsprcomp_annot[ConGood_Samples,]$PC1 + seurat_combo_M_mat_nsprcomp_annot[ConGood_Samples,]$PC1)
summary(model)

### Overlapping Scores ###

overlapping_Escores <- seurat_combo_E_mat_nsprcomp_annot
overlapping_Mscores <- seurat_combo_M_mat_nsprcomp_annot

BatchA <- row.names(overlapping_Escores[overlapping_Escores$Batch == "A",])
BatchB <- row.names(overlapping_Escores[overlapping_Escores$Batch == "B",])
BatchC <- row.names(overlapping_Escores[overlapping_Escores$Batch == "C",])


Time_Samples <- rownames(integration.combined.All@meta.data[integration.combined.All@meta.data$Time < 100,])
Con_Samples <- rownames(integration.combined.All@meta.data[integration.combined.All@meta.data$GBC_values > 0,])
ConAll_samples <- rownames(integration.combined.All@meta.data[integration.combined.All@meta.data$Batch == "C",])
Labeled_Samples <- c(Time_Samples,Con_Samples)

Escore_IQR_med <- quantile(overlapping_Escores[Labeled_Samples,]$PC1,probs=c(0.25,0.50,0.75))
Mscore_IQR_med <- quantile(overlapping_Mscores[Labeled_Samples,]$PC1,probs=c(0.25,0.50,0.75))
Bounds_IQR_med <- as.data.frame(cbind(Escore_IQR_med,Mscore_IQR_med))

GetBoundarySamples <- function(row,column,bonds_df){
  
  # Set E or M score based on column
  if(column == 1){
    score_annot <- overlapping_Escores
  }
  else if (column == 2){
    score_annot <- overlapping_Mscores
  }
  
  # Select samples based on score
  if (row == 1){ # If first values, get less than
    samples <- row.names(score_annot[score_annot$PC1 < bonds_df[row,column],])
  }
  else if (row == dim(bonds_df)[1]+1){ # if beyond the last value, > last value
    samples <- row.names(score_annot[score_annot$PC1 > bonds_df[row-1,column],])
  }
  else { # else between current value and previous value
    tmp <- score_annot[score_annot$PC1 > bonds_df[row-1,column],]
    samples <- row.names(tmp[tmp$PC1 < bonds_df[row,column],])
  }
  return(samples)
}

### Quartile Bounary Clusters ###

quart_clusters <- rep(NA, 16)
index = 1  
for(i in 1:(dim(Bounds_IQR_med)[1]+1)){
  for(j in 1:(dim(Bounds_IQR_med)[1]+1)){
    E_boundary_samples <- GetBoundarySamples(i,1,Bounds_IQR_med)
    M_boundary_sampels <- GetBoundarySamples(j,2,Bounds_IQR_med)
    intersect_samples = intersect(E_boundary_samples, M_boundary_sampels)
    filter_samples <- intersect(intersect_samples,union(Time_Samples,Con_Samples))
    quart_clusters[index] <- list(filter_samples)
    index = index + 1
  }
}

### Seruat/Quad Overlap ###
integration.combined.All.QuartOverlap <- integration.combined.All
integration.combined.All.QuartOverlap <- subset(integration.combined.All.QuartOverlap, subset = GBC != 'GBC0' | Batch == 'A' | Batch == 'B')
integration.combined.All.QuartOverlap@meta.data$QuadCluster <- 0

table(integration.combined.All.QuartOverlap@meta.data$QuadCluster)
table(integration.combined.All.QuartOverlap@meta.data$seurat_clusters)

integration.combined.All.QuartOverlap@meta.data[quart_clusters[[1]],]$QuadCluster <- 1
integration.combined.All.QuartOverlap@meta.data[quart_clusters[[2]],]$QuadCluster <- 2
integration.combined.All.QuartOverlap@meta.data[quart_clusters[[3]],]$QuadCluster <- 3
integration.combined.All.QuartOverlap@meta.data[quart_clusters[[4]],]$QuadCluster <- 4
integration.combined.All.QuartOverlap@meta.data[quart_clusters[[5]],]$QuadCluster <- 5
integration.combined.All.QuartOverlap@meta.data[quart_clusters[[6]],]$QuadCluster <- 6
integration.combined.All.QuartOverlap@meta.data[quart_clusters[[7]],]$QuadCluster <- 7
integration.combined.All.QuartOverlap@meta.data[quart_clusters[[8]],]$QuadCluster <- 8
integration.combined.All.QuartOverlap@meta.data[quart_clusters[[9]],]$QuadCluster <- 9
integration.combined.All.QuartOverlap@meta.data[quart_clusters[[10]],]$QuadCluster <- 10
integration.combined.All.QuartOverlap@meta.data[quart_clusters[[11]],]$QuadCluster <- 11
integration.combined.All.QuartOverlap@meta.data[quart_clusters[[12]],]$QuadCluster <- 12
integration.combined.All.QuartOverlap@meta.data[quart_clusters[[13]],]$QuadCluster <- 13
integration.combined.All.QuartOverlap@meta.data[quart_clusters[[14]],]$QuadCluster <- 14
integration.combined.All.QuartOverlap@meta.data[quart_clusters[[15]],]$QuadCluster <- 15
integration.combined.All.QuartOverlap@meta.data[quart_clusters[[16]],]$QuadCluster <- 16

table(integration.combined.All.QuartOverlap@meta.data$QuadCluster)

# Remove borderline samples (6 of 22765)
integration.combined.All.QuartOverlap <- subset(integration.combined.All.QuartOverlap, subset = QuadCluster != 0)
table(integration.combined.All.QuartOverlap@meta.data$QuadCluster)

overlap_table <- table(integration.combined.All.QuartOverlap@meta.data$QuadCluster,
                       integration.combined.All.QuartOverlap@meta.data$seurat_clusters)

overlap_table_batch <- table(integration.combined.All.QuartOverlap@meta.data$Batch,
                             integration.combined.All.QuartOverlap@meta.data$seurat_clusters)

overlap_table_phase <- table(integration.combined.All.QuartOverlap@meta.data$Phase,
                             integration.combined.All.QuartOverlap@meta.data$seurat_clusters)

overlap_table_dose <- table(integration.combined.All.QuartOverlap@meta.data$GBC_pM,
                            integration.combined.All.QuartOverlap@meta.data$seurat_clusters)

overlap_table_time <- table(integration.combined.All.QuartOverlap@meta.data$Time,
                            integration.combined.All.QuartOverlap@meta.data$seurat_clusters)

# 4x4 --- Bubble Plot
overlap_mat_long <- matrix(, nrow = 16, ncol =17)
index = 1 
for (i in 1:4){
  for (j in 1:4){
    overlap_mat_long[index,1] = i
    overlap_mat_long[index,2] = j
    overlap_mat_long[index,3] <- length(quart_clusters[[index]])
    
    overlap_mat_long[index,4] = overlap_table[index,1]/sum(overlap_table[,1])
    overlap_mat_long[index,5] = overlap_table[index,2]/sum(overlap_table[,2])
    overlap_mat_long[index,6] = overlap_table[index,3]/sum(overlap_table[,3])
    overlap_mat_long[index,7] = overlap_table[index,4]/sum(overlap_table[,4])
    overlap_mat_long[index,8] = overlap_table[index,5]/sum(overlap_table[,5])
    overlap_mat_long[index,9] = overlap_table[index,6]/sum(overlap_table[,6])
    overlap_mat_long[index,10] = overlap_table[index,7]/sum(overlap_table[,7])
    overlap_mat_long[index,11] = overlap_table[index,8]/sum(overlap_table[,8])
    overlap_mat_long[index,12] = overlap_table[index,9]/sum(overlap_table[,9])
    overlap_mat_long[index,13] = overlap_table[index,10]/sum(overlap_table[,10])
    overlap_mat_long[index,14] = overlap_table[index,11]/sum(overlap_table[,11])
    overlap_mat_long[index,15] = overlap_table[index,12]/sum(overlap_table[,12])
    overlap_mat_long[index,16] = overlap_table[index,13]/sum(overlap_table[,13])
    overlap_mat_long[index,17] = overlap_table[index,14]/sum(overlap_table[,14])
    
    index = index + 1
  }
}
overlap_mat_long_df <- as.data.frame(overlap_mat_long)
colnames(overlap_mat_long_df) <- c("X","Y","Count","S0","S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12","S13")

# Dose
ggplot(overlap_mat_long_df, aes(x = X, y = Y,size=Count,fill=S0)) + 
  geom_point(alpha=1.0, shape=21, color="black") +
  scale_size(range = c(5, 15)) +
  scale_fill_gradient2(low='black',high='red',mid='yellow',midpoint=0.25,limits=c(0,0.5),oob=squish) +
  theme_bw() + coord_fixed(ratio = 1) + ylim(c(0.75,4.25)) + xlim(c(0.75,4.25))

ggplot(overlap_mat_long_df, aes(x = X, y = Y,size=Count,fill=S3)) + 
  geom_point(alpha=1.0, shape=21, color="black") +
  scale_size(range = c(5, 15)) +
  scale_fill_gradient2(low='black',high='red',mid='yellow',midpoint=0.25,limits=c(0,0.5),oob=squish) +
  theme_bw() + coord_fixed(ratio = 1) + ylim(c(0.75,4.25)) + xlim(c(0.75,4.25))

ggplot(overlap_mat_long_df, aes(x = X, y = Y,size=Count,fill=S6)) + 
  geom_point(alpha=1.0, shape=21, color="black") +
  scale_size(range = c(5, 15)) +
  scale_fill_gradient2(low='black',high='red',mid='yellow',midpoint=0.25,limits=c(0,0.5),oob=squish) +
  theme_bw() + coord_fixed(ratio = 1) + ylim(c(0.75,4.25)) + xlim(c(0.75,4.25))

# Time
ggplot(overlap_mat_long_df, aes(x = X, y = Y,size=Count,fill=S4)) + 
  geom_point(alpha=1.0, shape=21, color="black") +
  scale_size(range = c(5, 15)) +
  scale_fill_gradient2(low='black',high='red',mid='yellow',midpoint=0.25,limits=c(0,0.5),oob=squish) +
  theme_bw() + coord_fixed(ratio = 1) + ylim(c(0.75,4.25)) + xlim(c(0.75,4.25))

ggplot(overlap_mat_long_df, aes(x = X, y = Y,size=Count,fill=S7)) + 
  geom_point(alpha=1.0, shape=21, color="black") +
  scale_size(range = c(5, 15)) +
  scale_fill_gradient2(low='black',high='red',mid='yellow',midpoint=0.25,limits=c(0,0.5),oob=squish) +
  theme_bw() + coord_fixed(ratio = 1) + ylim(c(0.75,4.25)) + xlim(c(0.75,4.25))

ggplot(overlap_mat_long_df, aes(x = X, y = Y,size=Count,fill=S12)) + 
  geom_point(alpha=1.0, shape=21, color="black") +
  scale_size(range = c(5, 15)) +
  scale_fill_gradient2(low='black',high='red',mid='yellow',midpoint=0.25,limits=c(0,0.5),oob=squish) +
  theme_bw() + coord_fixed(ratio = 1) + ylim(c(0.75,4.25)) + xlim(c(0.75,4.25))

# Zeros
ggplot(overlap_mat_long_df, aes(x = X, y = Y,size=Count,fill=S5)) + 
  geom_point(alpha=1.0, shape=21, color="black") +
  scale_size(range = c(5, 15)) +
  scale_fill_gradient2(low='black',high='red',mid='yellow',midpoint=0.25,limits=c(0,0.5),oob=squish) +
  theme_bw() + coord_fixed(ratio = 1) + ylim(c(0.75,4.25)) + xlim(c(0.75,4.25))

ggplot(overlap_mat_long_df, aes(x = X, y = Y,size=Count,fill=S9)) + 
  geom_point(alpha=1.0, shape=21, color="black") +
  scale_size(range = c(5, 15)) +
  scale_fill_gradient2(low='black',high='red',mid='yellow',midpoint=0.25,limits=c(0,0.5),oob=squish) +
  theme_bw() + coord_fixed(ratio = 1) + ylim(c(0.75,4.25)) + xlim(c(0.75,4.25))

# Barplots #

overlap_table_batch <- as.data.frame(table(integration.combined.All.QuartOverlap@meta.data$Batch,
                                           integration.combined.All.QuartOverlap@meta.data$seurat_clusters))
colnames(overlap_table_batch) <- c("Batch","Cluster","Count")
overlap_table_phase <- as.data.frame(table(integration.combined.All.QuartOverlap@meta.data$Phase,
                                           integration.combined.All.QuartOverlap@meta.data$seurat_clusters))
colnames(overlap_table_phase) <- c("Phase","Cluster","Count")
overlap_table_dose <- as.data.frame(table(integration.combined.All.QuartOverlap@meta.data$GBC_pM,
                                          integration.combined.All.QuartOverlap@meta.data$seurat_clusters))
colnames(overlap_table_dose) <- c("Dose","Cluster","Count")
overlap_table_time <- as.data.frame(table(integration.combined.All.QuartOverlap@meta.data$Time,
                                          integration.combined.All.QuartOverlap@meta.data$seurat_clusters))
colnames(overlap_table_time) <- c("Time","Cluster","Count")

overlap_table_batch_subset <- overlap_table_batch[overlap_table_batch$Cluster %in% c(0,3,6,4,7,12,5,9),]
ggplot(data=overlap_table_batch_subset, aes(x=Cluster, y=Count, fill=Batch)) +
  geom_bar(stat="identity", position=position_dodge()) + ggtitle("Batch")

overlap_table_phase_subset <- overlap_table_phase[overlap_table_phase$Cluster %in% c(0,3,6,4,7,12,5,9),]
ggplot(data=overlap_table_phase_subset, aes(x=Cluster, y=Count, fill=Phase)) +
  geom_bar(stat="identity", position=position_dodge()) + ggtitle("Phase")

overlap_table_time_subset <- overlap_table_time[as.numeric(as.character(overlap_table_time$Time)) < 100,]
overlap_table_time_subset <- overlap_table_time_subset[overlap_table_time_subset$Cluster %in% c(4,5,7,12),]
ggplot(data=overlap_table_time_subset, aes(x=Cluster, y=Count, fill=Time)) +
  geom_bar(stat="identity", position=position_dodge()) + ggtitle("Time")

overlap_table_dose_subset <- overlap_table_dose[as.numeric(as.character(overlap_table_dose$Dose)) > -1,]
overlap_table_dose_subset <- overlap_table_dose_subset[overlap_table_dose_subset$Cluster %in% c(0,3,9),]
ggplot(data=overlap_table_dose_subset, aes(x=Cluster, y=Count, fill=Dose)) +
  geom_bar(stat="identity", position=position_dodge()) + ggtitle("Dose")

batch_vs_cluster <- table(integration.combined.All.QuartOverlap@meta.data$Batch, integration.combined.All.QuartOverlap@meta.data$seurat_clusters)
batch_vs_cluster[3,1]/sum(batch_vs_cluster[,1]) # Cluster 0
batch_vs_cluster[3,4]/sum(batch_vs_cluster[,4]) # Cluster 3
batch_vs_cluster[3,10]/sum(batch_vs_cluster[,10]) # Cluster 9

sum(batch_vs_cluster[1:2,5])/sum(batch_vs_cluster[,5]) # Cluster 4
sum(batch_vs_cluster[1:2,6])/sum(batch_vs_cluster[,6]) # Cluster 5
sum(batch_vs_cluster[1:2,8])/sum(batch_vs_cluster[,8]) # Cluster 7

#############################################
### Dose Data Only --- Multiple Endpoints ###
#############################################

gg_means_plot <- function(plot_data,color_n,x_down,x_up,y_down,y_up){
  manual_colors <- hue_pal()(color_n)
  
  base_plot <- ggplot(plot_data, aes(x=Egenes, y=Mgenes)) + geom_point(size=0.1, alpha=0.75, color='lightgrey') + stat_density_2d(aes(alpha=(..level..)^1, fill=key_label),geom="polygon",show.legend=F) + scale_fill_manual(values=manual_colors) + scale_alpha_continuous(range=c(0,0.4)) + theme_bw()
  #base_plot
  
  data_Emeans <- r1<-with(plot_data, tapply(Egenes, key, mean))
  data_Mmeans <- r1<-with(plot_data, tapply(Mgenes, key, mean))
  data_Esds <- r1<-with(plot_data, tapply(Egenes, key, sd))
  data_Msds <- r1<-with(plot_data, tapply(Mgenes, key, sd))
  
  data_ceneters <- as.data.frame(t(rbind(data_Emeans,data_Mmeans,data_Esds,data_Msds)))
  names(data_ceneters) <- c("Emean","Mmean","Esd","Msd")
  
  #x_down = -5
  #x_up = 10
  #y_down = -5
  #y_up = 7.5
  
  E_20th = (x_up - x_down)/20
  M_20th  = (y_up - y_down)/20
  
  base_plot_centers <- base_plot + geom_point(data=data_ceneters, aes(x=Emean, y=Mmean),fill=manual_colors,pch=21,cex=4.0) + geom_errorbar(data=data_ceneters, aes(x=Emean, y=Mmean, ymin=Mmean-Msd, ymax=Mmean+Msd), colour="black", width=E_20th) + geom_errorbarh(data=data_ceneters, aes(x=Emean, y=Mmean, xmin=Emean-Esd, xmax=Emean+Esd), colour="black", height=M_20th) + xlim(x_down,x_up) + ylim(y_down,y_up) + xlab("Egenes") + ylab("Mgenes")
  #base_plot_centers
  return(base_plot_centers)
}


### Dose Only ###

ConCC_integration.combined <- readRDS("SeuratData/Seurat_concentration_DoseOnly_cellcycle_correct_hg38.rds")

### Fix annotation ###

ConCC_integration.combined@meta.data <- ConCC_integration.combined@meta.data[1:11220,]

# Check against original list
check_table <- read.table("GBCannot/CBC_GBC_summary.txt.alt_names",head=TRUE)
check_table$CBC_alt <- cbind(lapply(check_table$CBC.1, function(x) {sub(".1", "-1", x)}))
row.names(check_table) <- check_table$CBC_alt
overlap <- intersect(check_table$CBC_alt,row.names(ConCC_integration.combined@meta.data))

check_table <- as.data.frame(check_table)
check_table$meta_annot <- 0
check_table[overlap,]$meta_annot <- ConCC_integration.combined@meta.data[overlap ,]$GBC

table(check_table$meta_annot,check_table$GBC)

# Make replacements
#check_table[check_table$GBC == "GBC1" & check_table$meta_annot == "GBC0",]
ConCC_integration.combined@meta.data["CTGAAACAGTTGTCGT-1",]$GBC = "GBC1"
ConCC_integration.combined@meta.data["CTGAAACAGTTGTCGT-1",]$GBC_values = "1"
ConCC_integration.combined@meta.data["CTGAAACAGTTGTCGT-1",]$GBC_pM = 0
ConCC_integration.combined@meta.data["TGTATTCCATTTGCTT-1",]$GBC = "GBC1"
ConCC_integration.combined@meta.data["TGTATTCCATTTGCTT-1",]$GBC_values = "1"
ConCC_integration.combined@meta.data["TGTATTCCATTTGCTT-1",]$GBC_pM = 0
#check_table[check_table$GBC == "GBC4" & check_table$meta_annot == "GBC0",]
ConCC_integration.combined@meta.data["ACGCCGAGTGGTACAG-1",]$GBC = "GBC4"
ConCC_integration.combined@meta.data["ACGCCGAGTGGTACAG-1",]$GBC_values = "4"
ConCC_integration.combined@meta.data["ACGCCGAGTGGTACAG-1",]$GBC_pM = 50
#check_table[check_table$GBC == "GBC9" & check_table$meta_annot == "GBC0",]
ConCC_integration.combined@meta.data["CCGGGATCAGGCGATA-1",]$GBC = "GBC9"
ConCC_integration.combined@meta.data["CCGGGATCAGGCGATA-1",]$GBC_values = "9"
ConCC_integration.combined@meta.data["CCGGGATCAGGCGATA-1",]$GBC_pM = 200

check_table[overlap,]$meta_annot <- ConCC_integration.combined@meta.data[overlap ,]$GBC
table(check_table$meta_annot,check_table$GBC)

# Check correspondence post correction
table(ConCC_integration.combined@meta.data$GBC,ConCC_integration.combined@meta.data$GBC_values)
table(ConCC_integration.combined@meta.data$GBC,ConCC_integration.combined@meta.data$GBC_pM)

# Check missing samples
check_table_0 <- check_table[check_table$meta_annot ==0,]
dim(check_table_0)
kazu_barcodes <- read.table("data_KazuConcentration/barcodes.tsv")
length(intersect(kazu_barcodes$V1,check_table_0$CBC_alt))
# All but seven are missing from the filtered barcode set
# The missing seven are filtered (these are the 7 NA sample post stringent filtering)

p_dose_clusters <- DimPlot(ConCC_integration.combined, reduction = "umap", label = TRUE, repel = TRUE)

seurat_combo_Con_mat <- GetAssayData(object = ConCC_integration.combined, slot = "scale.data")

# Genes Sets #
EMT_genes <- read.table("EMTGenesUpdateV2.txt",header=TRUE)
E_genes <- subset(EMT_genes, EMT_genes$Annotation == 'E')
M_genes <- subset(EMT_genes, EMT_genes$Annotation == 'M')
Gene_Sets <- list('Egenes'=E_genes$Gene,'Mgenes'=M_genes$Gene)


### NN-PCA ###
library(nsprcomp)
set.seed(5) # Controls randomness for nnPC, minor variations , mainly to keep enriched gene lists consistent (+/-1 genes)

# For annotation data
seuratCon_combo_E <- ScaleData(ConCC_integration.combined, verbose = FALSE, vars.to.regress=c("S.Score","G2M.Score"), features=E_genes$Gene)
seuratCon_combo_M <- ScaleData(ConCC_integration.combined, verbose = FALSE, vars.to.regress=c("S.Score","G2M.Score"), features=M_genes$Gene)

# For consistent gene coverage
seuratCon_combo_E_mat <- seurat_combo_Con_mat[E_genes$Gene[E_genes$Gene %in% row.names(seurat_combo_Con_mat)],]
seuratCon_combo_M_mat <- seurat_combo_Con_mat[M_genes$Gene[M_genes$Gene %in% row.names(seurat_combo_Con_mat)],]

dim(seuratCon_combo_E_mat)
dim(seuratCon_combo_M_mat)

# E
seuratCon_combo_E_mat_nsprcomp <- nsprcomp(t(seuratCon_combo_E_mat),nneg=TRUE,ncomp=5)
seuratCon_combo_E_mat_mat_nsprcomp_annot <- cbind(seuratCon_combo_E_mat_nsprcomp$x,seuratCon_combo_E@meta.data)

# M 
seuratCon_combo_M_mat_nsprcomp <- nsprcomp(t(seuratCon_combo_M_mat),nneg=TRUE,ncomp=5)
seuratCon_combo_M_mat_mat_nsprcomp_annot <- cbind(seuratCon_combo_M_mat_nsprcomp$x,seuratCon_combo_M@meta.data)

seuratCon_combo_E_mat_mat_nsprcomp_annot_GBCGood <- seuratCon_combo_E_mat_mat_nsprcomp_annot[seuratCon_combo_E_mat_mat_nsprcomp_annot$GBC_pM >= 0,]
seuratCon_combo_M_mat_mat_nsprcomp_annot_GBCGood <- seuratCon_combo_M_mat_mat_nsprcomp_annot[seuratCon_combo_M_mat_mat_nsprcomp_annot$GBC_pM >= 0,]

### Supplemental --- Two Endpoints ###
p_dose <- DimPlot(ConCC_integration.combined, reduction = "umap",group.by = "GBC_pM")
p_clusters <- DimPlot(ConCC_integration.combined, reduction = "umap", label = TRUE, repel = TRUE)
p_dose+p_clusters

plot_data <- as.data.frame(cbind(seuratCon_combo_E_mat_mat_nsprcomp_annot_GBCGood$PC1,seuratCon_combo_M_mat_mat_nsprcomp_annot_GBCGood$PC1,seuratCon_combo_E_mat_mat_nsprcomp_annot_GBCGood$seurat_clusters))
colnames(plot_data) <- c("Egenes","Mgenes","key")
plot_data$key_label <- factor(plot_data$key)

manual_colors <- hue_pal()(7)
dose_means_plot <- gg_means_plot(plot_data,7,-7.5,10.5,-7.5,10)
dose_means_plot

cluster_0 <- row.names(ConCC_integration.combined@meta.data[ConCC_integration.combined@meta.data$seurat_clusters==0,])
cluster_1 <- row.names(ConCC_integration.combined@meta.data[ConCC_integration.combined@meta.data$seurat_clusters==1,])

cluster_0_data <- seuratCon_combo_E_mat_mat_nsprcomp_annot_GBCGood[row.names(seuratCon_combo_E_mat_mat_nsprcomp_annot_GBCGood) %in% cluster_0,]
cluster_1_data <- seuratCon_combo_E_mat_mat_nsprcomp_annot_GBCGood[row.names(seuratCon_combo_E_mat_mat_nsprcomp_annot_GBCGood) %in% cluster_1,]

cluster_0_DoseByPhase <- as.data.frame(table(cluster_0_data$Phase,cluster_0_data$GBC_pM))
cluster_0_DoseByPhase$Freq <- cluster_0_DoseByPhase$Freq/sum(cluster_0_DoseByPhase$Freq)
colnames(cluster_0_DoseByPhase) <- c("Phase","Dose","Freq")
cluster_1_DoseByPhase <- as.data.frame(table(cluster_1_data$Phase,cluster_1_data$GBC_pM))
cluster_1_DoseByPhase$Freq <- cluster_1_DoseByPhase$Freq/sum(cluster_1_DoseByPhase$Freq)
colnames(cluster_1_DoseByPhase) <- c("Phase","Dose","Freq")

ggplot(data=cluster_0_DoseByPhase, aes(x=Dose, y=Freq, fill=Phase)) +
  geom_bar(stat="identity", position=position_dodge()) + ggtitle("Cluster 0")

sum(cluster_0_DoseByPhase[as.numeric(as.character(cluster_0_DoseByPhase$Dose)) < 100,]$Freq)
# 0.1526055
sum(cluster_0_DoseByPhase[cluster_0_DoseByPhase$Phase=="G1",]$Freq)
# 0.3846154
sum(cluster_0_DoseByPhase[as.numeric(as.character(cluster_0_DoseByPhase$Dose)) < 100 & cluster_0_DoseByPhase$Phase=="G1",]$Freq)
# 0.05334988
sum(cluster_0_DoseByPhase[as.numeric(as.character(cluster_0_DoseByPhase$Dose)) < 100 & cluster_0_DoseByPhase$Phase=="G1",]$Freq)/sum(cluster_0_DoseByPhase[as.numeric(as.character(cluster_0_DoseByPhase$Dose)) < 100,]$Freq)
# 0.3495935

ggplot(data=cluster_1_DoseByPhase, aes(x=Dose, y=Freq, fill=Phase)) +
  geom_bar(stat="identity", position=position_dodge()) + ggtitle("Cluster 1")

sum(cluster_1_DoseByPhase[as.numeric(as.character(cluster_1_DoseByPhase$Dose)) < 100,]$Freq)
# 0.4682484
sum(cluster_1_DoseByPhase[cluster_1_DoseByPhase$Phase=="G1",]$Freq)
# 0.6608514
sum(cluster_1_DoseByPhase[as.numeric(as.character(cluster_1_DoseByPhase$Dose)) < 100 & cluster_1_DoseByPhase$Phase=="G1",]$Freq)
# 0.3251919
