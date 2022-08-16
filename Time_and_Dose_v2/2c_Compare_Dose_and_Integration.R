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

########################
### Generate Figures ###
########################

### Load Processed Data ###

dose_meta <- readRDS("Dose_MetaData.rds")
UMAP_dose_values <- readRDS("Dose_UMAP_values.rds")

dose_samples <- readRDS("DoseSamples.rds")
time_samples <- readRDS("TimeSamples.rds")

Dose_E_mat_nsprcomp <- readRDS("seuratCon_combo_E_mat_nsprcomp_fullcorr.rds")
Dose_M_mat_nsprcomp <- readRDS("seuratCon_combo_M_mat_nsprcomp_fullcorr.rds")

integated_meta <- readRDS("Integrated_MetaData.rds")
UMAP_integated_values <- readRDS("Integrated_UMAP_values.rds")

dose_samples <- readRDS("DoseSamples.rds")
time_samples <- readRDS("TimeSamples.rds")

Integrated_E_mat_nsprcomp <- readRDS("Integrated_E_mat_nsprcomp.rds")
Integrated_M_mat_nsprcomp <- readRDS("Integrated_M_mat_nsprcomp.rds")

### Correlation between Dose and Integrated PCs ###
cor(Dose_E_mat_nsprcomp$x[dose_samples,1],Integrated_E_mat_nsprcomp$x[dose_samples,1]) # 0.8368308

cor(Dose_M_mat_nsprcomp$x[dose_samples,1],Integrated_M_mat_nsprcomp$x[dose_samples,1]) # 0.8104603
cor(Dose_M_mat_nsprcomp$x[dose_samples,1],Integrated_M_mat_nsprcomp$x[dose_samples,2]) # 0.7823477
cor(Dose_M_mat_nsprcomp$x[dose_samples,2],Integrated_M_mat_nsprcomp$x[dose_samples,1]) # 0.7872044
cor(Dose_M_mat_nsprcomp$x[dose_samples,2],Integrated_M_mat_nsprcomp$x[dose_samples,2]) # 0.2135162

######################################
### Integration Data with Dose PCs ###
######################################

# Add annotation data to UMAP valuse
UMAP_integrated_annot <- merge(UMAP_integated_values,integated_meta,by=0)
row.names(UMAP_integrated_annot) <- UMAP_integrated_annot$Row.names
UMAP_integrated_annot$GBC_pM_factor <- factor(UMAP_integrated_annot$GBC_pM)
UMAP_integrated_annot <- UMAP_integrated_annot[dose_samples,]

# Add PCs
UMAP_integrated_annot$E_PC1_Dose <- Dose_E_mat_nsprcomp$x[dose_samples,1]
UMAP_integrated_annot$M_PC1_Dose <- Dose_M_mat_nsprcomp$x[dose_samples,1]
UMAP_integrated_annot$M_PC2_Dose <- Dose_M_mat_nsprcomp$x[dose_samples,2]

UMAP_integrated_annot$E_PC1 <- Integrated_E_mat_nsprcomp$x[dose_samples,1]
UMAP_integrated_annot$M_PC1 <- Integrated_M_mat_nsprcomp$x[dose_samples,1]
UMAP_integrated_annot$M_PC2 <- Integrated_M_mat_nsprcomp$x[dose_samples,2]

# E PCs
UMAP_E1 <- ggplot(UMAP_integrated_annot, aes(x = UMAP_1, y = UMAP_2, colour=E_PC1)) +
  geom_point(size=0.5) + theme_bw() + scale_color_gradient2(midpoint=0, low="darkblue", mid="white",
                                                            high="red", space ="Lab") + labs(color="E_PC1") + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1)
UMAP_E1

UMAP_E1_Dose <- ggplot(UMAP_integrated_annot, aes(x = UMAP_1, y = UMAP_2, colour=E_PC1_Dose)) +
  geom_point(size=0.5) + theme_bw() + scale_color_gradient2(midpoint=0, low="darkblue", mid="white",
                                                            high="red", space ="Lab") + labs(color="E_PC1") + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1)
UMAP_E1_Dose


# M PCs
UMAP_M1 <- ggplot(UMAP_integrated_annot, aes(x = UMAP_1, y = UMAP_2, colour=M_PC1)) +
  geom_point(size=0.5) + theme_bw() + scale_color_gradient2(midpoint=0, low="darkblue", mid="white",
                                                            high="red", space ="Lab") + labs(color="M PC1") + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1)
UMAP_M1
UMAP_M2_Dose <- ggplot(UMAP_integrated_annot, aes(x = UMAP_1, y = UMAP_2, colour=M_PC2_Dose)) +
  geom_point(size=0.5) + theme_bw() + scale_color_gradient2(midpoint=0, low="darkblue", mid="white",
                                                            high="red", space ="Lab") + labs(color="M PC2") + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1)
UMAP_M2_Dose


UMAP_M2 <- ggplot(UMAP_integrated_annot, aes(x = UMAP_1, y = UMAP_2, colour=M_PC2)) +
  geom_point(size=0.5) + theme_bw() + scale_color_gradient2(midpoint=0, low="darkblue", mid="white",
                                                            high="red", space ="Lab") + labs(color="M PC2") + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1)
UMAP_M2
UMAP_M1_Dose <- ggplot(UMAP_integrated_annot, aes(x = UMAP_1, y = UMAP_2, colour=M_PC1_Dose)) +
  geom_point(size=0.5) + theme_bw() + scale_color_gradient2(midpoint=0, low="darkblue", mid="white",
                                                            high="red", space ="Lab") + labs(color="M PC1") + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1)
UMAP_M1_Dose

par(mfrow = c(3, 2))
grid.arrange(UMAP_E1,UMAP_E1_Dose,UMAP_M1,UMAP_M2_Dose,UMAP_M2,UMAP_M1_Dose, nrow=3)
# Fig_M1vM2_Int, save as 8x8

#####################################
### Dose Data with Integrated PCs ###
#####################################

# Add annotation to UMAP values
UMAP_dose_annot <- merge(UMAP_dose_values,dose_meta,by=0)
row.names(UMAP_dose_annot) <- UMAP_dose_annot$Row.names
UMAP_dose_annot$GBC_pM_factor <- factor(UMAP_dose_annot$GBC_pM)
UMAP_dose_annot <- UMAP_dose_annot[dose_samples,]

# Add PCs
UMAP_dose_annot$E_PC1 <- Dose_E_mat_nsprcomp$x[dose_samples,1]
UMAP_dose_annot$M_PC1 <- Dose_M_mat_nsprcomp$x[dose_samples,1]
UMAP_dose_annot$M_PC2 <- Dose_M_mat_nsprcomp$x[dose_samples,2]

UMAP_dose_annot$E_PC1_Int <- Integrated_E_mat_nsprcomp$x[dose_samples,1]
UMAP_dose_annot$M_PC1_Int <- Integrated_M_mat_nsprcomp$x[dose_samples,1]
UMAP_dose_annot$M_PC2_Int <- Integrated_M_mat_nsprcomp$x[dose_samples,2]

# E PCs
UMAP_E1 <- ggplot(UMAP_dose_annot, aes(x = UMAP_1, y = UMAP_2, colour=E_PC1)) +
  geom_point(size=0.5) + theme_bw() + scale_color_gradient2(midpoint=0, low="darkblue", mid="white",
                                                            high="red", space ="Lab") + labs(color="E PC1") + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1)
UMAP_E1

UMAP_E1_Int <- ggplot(UMAP_dose_annot, aes(x = UMAP_1, y = UMAP_2, colour=E_PC1_Int)) +
  geom_point(size=0.5) + theme_bw() + scale_color_gradient2(midpoint=0, low="darkblue", mid="white",
                                                            high="red", space ="Lab") + labs(color="E PC1") + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1)
UMAP_E1_Int


# M PCs
UMAP_M1 <- ggplot(UMAP_dose_annot, aes(x = UMAP_1, y = UMAP_2, colour=M_PC1)) +
  geom_point(size=0.5) + theme_bw() + scale_color_gradient2(midpoint=0, low="darkblue", mid="white",
                                                            high="red", space ="Lab") + labs(color="M PC1") + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1)
UMAP_M1
UMAP_M2_Int <- ggplot(UMAP_dose_annot, aes(x = UMAP_1, y = UMAP_2, colour=M_PC2_Int)) +
  geom_point(size=0.5) + theme_bw() + scale_color_gradient2(midpoint=0, low="darkblue", mid="white",
                                                            high="red", space ="Lab") + labs(color="M PC2") + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1)
UMAP_M2_Int


UMAP_M2 <- ggplot(UMAP_dose_annot, aes(x = UMAP_1, y = UMAP_2, colour=M_PC2)) +
  geom_point(size=0.5) + theme_bw() + scale_color_gradient2(midpoint=0, low="darkblue", mid="white",
                                                            high="red", space ="Lab") + labs(color="M PC2") + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1)
UMAP_M2
UMAP_M1_Int <- ggplot(UMAP_dose_annot, aes(x = UMAP_1, y = UMAP_2, colour=M_PC1_Int)) +
  geom_point(size=0.5) + theme_bw() + scale_color_gradient2(midpoint=0, low="darkblue", mid="white",
                                                            high="red", space ="Lab") + labs(color="M PC1") + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1)
UMAP_M1_Int

par(mfrow = c(3, 2))
grid.arrange(UMAP_E1,UMAP_E1_Int,UMAP_M1,UMAP_M2_Int,UMAP_M2,UMAP_M1_Int, nrow=3)
# Fig_M1vM2_Dose, save as 8x8

###############################
### Overlapping NN Clusters ###
###############################

# Add dose clusters to integrated data
UMAP_integrated_annot$dose_clusters <- -1
UMAP_integrated_annot[dose_samples,]$dose_clusters <- UMAP_dose_annot[dose_samples,]$seurat_clusters
UMAP_integrated_annot$dose_clusters <- as.factor(UMAP_integrated_annot$dose_clusters-1)
table(UMAP_integrated_annot$dose_clusters)

# Add integrated clusters to dose data
UMAP_dose_annot$int_clusters <- -1
UMAP_dose_annot[dose_samples,]$int_clusters <- UMAP_integrated_annot[dose_samples,]$seurat_clusters
UMAP_dose_annot$int_clusters <- as.factor(UMAP_dose_annot$int_clusters-1)
table(UMAP_dose_annot$int_clusters)

# Clusters on Dose Data
UMAP_Dose_by_DoseClusters <- ggplot(UMAP_dose_annot, aes(x = UMAP_1, y = UMAP_2, colour=seurat_clusters)) +
  geom_point(size=0.5) + theme_bw() + labs(color="clusters") + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1,legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 3)))
UMAP_Dose_by_DoseClusters
UMAP_Dose_by_IntClusters <- ggplot(UMAP_dose_annot, aes(x = UMAP_1, y = UMAP_2, colour=int_clusters)) +
  geom_point(size=0.5) + theme_bw() + labs(color="clusters") + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1,legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 3)))
UMAP_Dose_by_IntClusters

# Clusters on Integrated Data
UMAP_dose_annot_subset <- UMAP_dose_annot[UMAP_dose_annot$int_clusters %in% c(0,1,4,8),]
UMAP_Dose_by_IntClusters_subset <- ggplot(UMAP_dose_annot_subset, aes(x = UMAP_1, y = UMAP_2, colour=int_clusters)) +
  geom_point(size=0.5) + theme_bw() + labs(color="clusters") + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1,legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 3)))
UMAP_Dose_by_IntClusters_subset

UMAP_Int_by_IntClusters <- ggplot(UMAP_integrated_annot, aes(x = UMAP_1, y = UMAP_2, colour=seurat_clusters)) +
  geom_point(size=0.5) + theme_bw() + labs(color="clusters") + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1,legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 3)))
UMAP_Int_by_IntClusters
UMAP_Int_by_DoseClusters <- ggplot(UMAP_integrated_annot, aes(x = UMAP_1, y = UMAP_2, colour=dose_clusters)) +
  geom_point(size=0.5) + theme_bw() + labs(color="clusters") + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1,legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 3)))
UMAP_Int_by_DoseClusters

UMAP_integrated_annot_subset <- UMAP_integrated_annot[UMAP_integrated_annot$dose_clusters %in% c(0,2),]
UMAP_Int_by_DoseClusters_subset <- ggplot(UMAP_integrated_annot_subset, aes(x = UMAP_1, y = UMAP_2, colour=dose_clusters)) +
  geom_point(size=0.5) + theme_bw() + labs(color="clusters") + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1,legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 3)))
UMAP_Int_by_DoseClusters_subset

par(mfrow = c(2, 2))
grid.arrange(UMAP_Dose_by_DoseClusters,UMAP_Int_by_IntClusters,
             UMAP_Dose_by_IntClusters,UMAP_Int_by_DoseClusters, nrow=2)

# Overlap in Clusters 
int_clusterA <- UMAP_integrated_annot[UMAP_integrated_annot$seurat_clusters %in% c(0,4),]
int_clusterB <- UMAP_integrated_annot[UMAP_integrated_annot$seurat_clusters %in% c(1,8),]

dose_cluster0 <- UMAP_integrated_annot[UMAP_integrated_annot$dose_clusters == 0,]
dose_cluster2 <- UMAP_integrated_annot[UMAP_integrated_annot$dose_clusters == 2,]

# Cacluate Jacard Index
inter_A <- length(intersect(row.names(int_clusterA),row.names(dose_cluster0)))
union_A <- length(union(row.names(int_clusterA),row.names(dose_cluster0)))
JI_A <- inter_A/union_A
Jsim_A = inter_A/(dim(int_clusterA)[1] + dim(dose_cluster0)[1] - inter_A)
  
inter_B <- length(intersect(row.names(int_clusterB),row.names(dose_cluster2)))
union_B <- length(union(row.names(int_clusterB),row.names(dose_cluster2)))
JI_B <- inter_B/union_B
Jsim_B = inter_B/(dim(int_clusterB)[1] + dim(dose_cluster2)[1] - inter_B)

# Test Jaccard Index
library(jaccard)
int_cluster_A_binary <- dose_samples %in% row.names(int_clusterA)
int_cluster_B_binary <- dose_samples %in% row.names(int_clusterB)

dose_cluster_0_binary <- dose_samples %in% row.names(dose_cluster0)
dose_cluster_2_binary <- dose_samples %in% row.names(dose_cluster2)

# For details on test, see https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-3118-5#Sec2
JTest_A <- jaccard.test(int_cluster_A_binary,dose_cluster_0_binary,method='exact')
JTest_B <- jaccard.test(int_cluster_B_binary,dose_cluster_2_binary,method='exact')

JTest_Aboot <- jaccard.test(int_cluster_A_binary,dose_cluster_0_binary,method='bootstrap',B=100000)
JTest_Bboot <- jaccard.test(int_cluster_B_binary,dose_cluster_2_binary,method='bootstrap',B=100000)
