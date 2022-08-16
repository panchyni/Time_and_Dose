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

#################
### Functions ###
#################

gg_means_plot <- function(plot_data,color_n,x_down,x_up,y_down,y_up){
  manual_colors <- hue_pal()(color_n)
  
  base_plot <- ggplot(plot_data, aes(x=Egenes, y=Mgenes)) + geom_point(size=0.1, alpha=0.75, color='lightgrey') + stat_density_2d(aes(alpha=(..level..)^1, fill=key_label),geom="polygon",show.legend=F) + scale_fill_manual(values=manual_colors) + scale_alpha_continuous(range=c(0,0.4)) + theme_bw() + theme(axis.text = element_text(size = 10), axis.title=element_text(size=11), aspect.ratio=1)
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
  
  base_plot_centers <- base_plot + geom_point(data=data_ceneters, aes(x=Emean, y=Mmean),fill=manual_colors,pch=21,cex=4.0) + geom_errorbar(data=data_ceneters, aes(x=Emean, y=Mmean, ymin=Mmean-Msd, ymax=Mmean+Msd), colour="black", width=E_20th) + geom_errorbarh(data=data_ceneters, aes(x=Emean, y=Mmean, xmin=Emean-Esd, xmax=Emean+Esd), colour="black", height=M_20th) + xlim(x_down,x_up) + ylim(y_down,y_up) + xlab("E PC") + ylab("M PC")
  #base_plot_centers
  return(base_plot_centers)
}

gg_means <- function(plot_data,color_n,x_down,x_up,y_down,y_up){
  
  data_Emeans <- r1<-with(plot_data, tapply(Egenes, key, mean))
  data_Mmeans <- r1<-with(plot_data, tapply(Mgenes, key, mean))
  data_Esds <- r1<-with(plot_data, tapply(Egenes, key, sd))
  data_Msds <- r1<-with(plot_data, tapply(Mgenes, key, sd))
  
  data_ceneters <- as.data.frame(t(rbind(data_Emeans,data_Mmeans,data_Esds,data_Msds)))
  names(data_ceneters) <- c("Emean","Mmean","Esd","Msd")
  return(data_ceneters)
}

########################
### Generate Figures ###
########################

### Load Processed Data ###

dose_meta <- readRDS("Dose_MetaData.rds")
dose_scale <- readRDS("Dose_ScaledData.rds")
UMAP_values <- readRDS("Dose_UMAP_values.rds")

dose_samples <- readRDS("DoseSamples.rds")
time_samples <- readRDS("TimeSamples.rds")

Dose_E_mat_nsprcomp <- readRDS("seuratCon_combo_E_mat_nsprcomp_fullcorr.rds")
Dose_M_mat_nsprcomp <- readRDS("seuratCon_combo_M_mat_nsprcomp_fullcorr.rds")

### PC UMAP Plots ###

# Add annotation data to UMAP
UMAP_annot <- merge(UMAP_values,dose_meta,by=0)
row.names(UMAP_annot) <- UMAP_annot$Row.names
UMAP_annot$GBC_pM_factor <- factor(UMAP_annot$GBC_pM)

# Add scaled genes values
UMAP_annot$VIM <- dose_scale["VIM",UMAP_annot$Row.names]
UMAP_annot$CDH1 <- dose_scale["CDH1",UMAP_annot$Row.names]
UMAP_annot$FN1 <- dose_scale["FN1",UMAP_annot$Row.names]

# Add PCs
UMAP_samples <- row.names(UMAP_annot)
UMAP_annot$E_PC1 <- Dose_E_mat_nsprcomp$x[UMAP_samples,1]

UMAP_annot$M_PC1 <- Dose_M_mat_nsprcomp$x[UMAP_samples,1]
UMAP_annot$M_PC2 <- Dose_M_mat_nsprcomp$x[UMAP_samples,2]

# Subset only labeled samples
UMAP_annot <- UMAP_annot[dose_samples,]

# E PCs
UMAP_E1 <- ggplot(UMAP_annot, aes(x = UMAP_1, y = UMAP_2, colour=E_PC1)) +
  geom_point(size=0.5) + theme_bw() + scale_color_gradient2(midpoint=0, low="darkblue", mid="white",
                                                            high="red", space ="Lab",guide=guide_colorbar(direction="horizontal")) + labs(color="E PC1") + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1,legend.position = "bottom")
UMAP_E1

# M PCs
UMAP_M1 <- ggplot(UMAP_annot, aes(x = UMAP_1, y = UMAP_2, colour=M_PC1)) +
  geom_point(size=0.5) + theme_bw() + scale_color_gradient2(midpoint=0, low="darkblue", mid="white",
                                                            high="red", space ="Lab",guide=guide_colorbar(direction="horizontal")) + labs(color="M PC1") + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1,legend.position = "bottom")
UMAP_M1

UMAP_M2 <- ggplot(UMAP_annot, aes(x = UMAP_1, y = UMAP_2, colour=M_PC2)) +
  geom_point(size=0.5) + theme_bw() + scale_color_gradient2(midpoint=0, low="darkblue", mid="white",
                                                            high="red", space ="Lab",guide=guide_colorbar(direction="horizontal")) + labs(color="M PC2") + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1,legend.position = "bottom")
UMAP_M2

# Genes
# M PCs
UMAP_VIM <- ggplot(UMAP_annot, aes(x = UMAP_1, y = UMAP_2, colour=VIM)) +
  geom_point(size=0.5) + theme_bw() + scale_color_gradient2(midpoint=0, low="darkblue", mid="white",
                                                            high="red", space ="Lab",limits=c(-3,3),guide=guide_colorbar(direction="horizontal")) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1,legend.position = "bottom")
UMAP_VIM

UMAP_FN1 <- ggplot(UMAP_annot, aes(x = UMAP_1, y = UMAP_2, colour=FN1)) +
  geom_point(size=0.5) + theme_bw() + scale_color_gradient2(midpoint=0, low="darkblue", mid="white",
                                                            high="red", space ="Lab",limits=c(-3,3),guide=guide_colorbar(direction="horizontal")) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1,legend.position = "bottom")
UMAP_FN1

UMAP_CDH1 <- ggplot(UMAP_annot, aes(x = UMAP_1, y = UMAP_2, colour=CDH1)) +
  geom_point(size=0.5) + theme_bw() + scale_color_gradient2(midpoint=0, low="darkblue", mid="white",
                                                            high="red", space ="Lab",limits=c(-3,3),guide=guide_colorbar(direction="horizontal")) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1,legend.position = "bottom")
UMAP_CDH1 

### Progression Plots ### 

# Dosage UAMP
UMAP_dose <- ggplot(UMAP_annot, aes(x = UMAP_1, y = UMAP_2, colour=GBC_pM_factor)) +
  geom_point(size=0.5) + theme_bw() + labs(color="pM") + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1,legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 3)))

# E1 vs M1
plot_data <- as.data.frame(cbind(UMAP_annot$E_PC1,UMAP_annot$M_PC1,UMAP_annot$GBC_pM))
colnames(plot_data) <- c("Egenes","Mgenes","key")
plot_data$key_label <- factor(plot_data$key)

manual_colors <- hue_pal()(8)
base_plot_centers_M1 <- gg_means_plot(plot_data,8,-7.5,10,-7.5,10)
base_plot_centers_M1

# E1 vs M2
plot_data <- as.data.frame(cbind(UMAP_annot[dose_samples,]$E_PC1,UMAP_annot[dose_samples,]$M_PC2,UMAP_annot[dose_samples,]$GBC_pM))
colnames(plot_data) <- c("Egenes","Mgenes","key")
plot_data$key_label <- factor(plot_data$key)

manual_colors <- hue_pal()(8)
base_plot_centers_M2 <- gg_means_plot(plot_data,8,-7.5,10,-7.5,10)
base_plot_centers_M2

# Condense plot to single plots

par(mfrow = c(3, 3))
grid.arrange(UMAP_E1,UMAP_M1,UMAP_M2,UMAP_CDH1,UMAP_VIM,UMAP_FN1,
             UMAP_dose, base_plot_centers_M1, base_plot_centers_M2, nrow=3)
# Figure 1

### Seperation Table ### 
mwu_greater <- function(A, B){
  mwu_test <- wilcox.test(A, B,alternative='greater')
  auc_roc <- mwu_test$statistic/(length(A) * length(B))
  #print(mwu_test$p.value)
  #print(auc_roc)
  return (c(mwu_test$p.value,auc_roc))
}

# Go through pairs of dose samples
results <- matrix(,nrow=7,ncol=8)
dose_values <- c(1,2,3,4,8,9,12,13)
for (index in 1:7){
  
  # Get each value and the next
  v1 = dose_values[index]
  v2 = dose_values[index+1]
  
  # Get samples
  lesser_treatmeant = row.names(dose_meta[dose_meta$GBC_values == v1,])
  greater_treatment = row.names(dose_meta[dose_meta$GBC_values == v2,])
  
  # Check
  #print(paste0(v1," ",v2))
  #print(table(dose_meta[lesser_treatmeant,]$GBC_pM))
  #print(table(dose_meta[greater_treatment,]$GBC_pM))
  
  results[index,1] <- names(table(dose_meta[lesser_treatmeant,]$GBC_pM))
  results[index,2] <- names(table(dose_meta[greater_treatment,]$GBC_pM))
  
  # Run Tests for E1, M1, and M2
  results[index,3:4] <- mwu_greater(Dose_E_mat_nsprcomp$x[lesser_treatmeant,1],Dose_E_mat_nsprcomp$x[greater_treatment,1]) # use lesser for E
  results[index,5:6] <- mwu_greater(Dose_M_mat_nsprcomp$x[greater_treatment,1],Dose_M_mat_nsprcomp$x[lesser_treatmeant,1])
  results[index,7:8] <- mwu_greater(Dose_M_mat_nsprcomp$x[greater_treatment,2],Dose_M_mat_nsprcomp$x[lesser_treatmeant,2])
}
results

####################
### Supplemental ###
####################

### NN Clusters ###

# UMAP
UMAP_NN <- ggplot(UMAP_annot, aes(x = UMAP_1, y = UMAP_2, colour=seurat_clusters)) +
  geom_point(size=0.5) + theme_bw() + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1,legend.position = "bottom")
UMAP_NN 

# Progression
plot_data <- as.data.frame(cbind(UMAP_annot[dose_samples,]$E_PC1,UMAP_annot[dose_samples,]$M_PC1,UMAP_annot[dose_samples,]$seurat_clusters))
colnames(plot_data) <- c("Egenes","Mgenes","key")
plot_data$key_label <- factor(plot_data$key)

manual_colors <- hue_pal()(7)
base_plot_centers_NN <- gg_means_plot(plot_data,7,-7.5,10.5,-7.5,10.5)
base_plot_centers_NN

# Barplot
NN_Cluster0 <- UMAP_annot[UMAP_annot$seurat_clusters ==0,]
NN_Cluster2 <- UMAP_annot[UMAP_annot$seurat_clusters ==2,]

Cluster0_table <- as.data.frame(table(NN_Cluster0$GBC_pM,NN_Cluster0$Phase))
Cluster2_table <- as.data.frame(table(NN_Cluster2$GBC_pM,NN_Cluster2$Phase))

Cluster0_table$Var1 <- as.factor(Cluster0_table$Var1)
Cluster2_table$Var1 <- as.factor(Cluster2_table$Var1)

ggplot(Cluster0_table,aes(x=Var1,y=Freq,fill=Var2)) + geom_bar(stat="identity", position=position_dodge()) +
  xlab("Dose") + ylab("Count") + labs(fill="Phase") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1)

ggplot(Cluster2_table,aes(x=Var1,y=Freq,fill=Var2)) + geom_bar(stat="identity", position=position_dodge()) +
  xlab("Dose") + ylab("Count") + labs(fill="Phase") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1)

### Unlabeled Points ###

# Label Data
plot_data <- as.data.frame(cbind(UMAP_annot$E_PC1,UMAP_annot$M_PC1,UMAP_annot$seurat_clusters))
colnames(plot_data) <- c("Egenes","Mgenes","key")
plot_data$key_label <- factor(plot_data$key)

#manual_colors <- hue_pal()(7)
#dose_means_plot <- gg_means_plot(plot_data,7,-7.5,10.5,-7.5,10)
#dose_means_plot
dose_means <- gg_means(plot_data)

# Define and annot UMAP values for unalebed data
UMAP_annot_unlabeled <- merge(UMAP_values,dose_meta,by=0)
UMAP_annot_unlabeled <- UMAP_annot_unlabeled[UMAP_annot_unlabeled$GBC_values==0,]
row.names(UMAP_annot_unlabeled) <- UMAP_annot_unlabeled$Row.names
UMAP_annot_unlabeled$GBC_pM_factor <- factor(UMAP_annot_unlabeled$GBC_pM)

# Add PCs
UMAP_samples <- row.names(UMAP_annot_unlabeled)
UMAP_annot_unlabeled$E_PC1 <- Dose_E_mat_nsprcomp$x[UMAP_samples,1]

UMAP_annot_unlabeled$M_PC1 <- Dose_M_mat_nsprcomp$x[UMAP_samples,1]
UMAP_annot_unlabeled$M_PC2 <- Dose_M_mat_nsprcomp$x[UMAP_samples,2]

plot_data_unlabeled <- as.data.frame(cbind(UMAP_annot_unlabeled$E_PC1,UMAP_annot_unlabeled$M_PC1,UMAP_annot_unlabeled$seurat_clusters))
colnames(plot_data_unlabeled) <- c("Egenes","Mgenes","key")
plot_data_unlabeled$key_label <- factor(plot_data_unlabeled$key)

#manual_colors <- hue_pal()(7)
#dose_means_unlabeled_plot <- gg_means_plot(plot_data_unlabeled,7,-7.5,10.5,-7.5,10)
#dose_means_unlabeled_plot
dose_means_unlabeled <- gg_means(plot_data_unlabeled)

plot_data <- as.data.frame(cbind(UMAP_annot$E_PC1,UMAP_annot$M_PC1))
colnames(plot_data) <- c("Egenes","Mgenes") 
base_plot <- ggplot(plot_data, aes(x=Egenes, y=Mgenes)) + geom_point(size=0.1, alpha=0.75, color='lightgrey') + theme_bw()
base_plot_centers <- base_plot + geom_point(data=dose_means, aes(x=Emean, y=Mmean),fill='black',pch=21,cex=4.0) + geom_errorbar(data=dose_means, aes(x=Emean, y=Mmean, ymin=Mmean-Msd, ymax=Mmean+Msd), colour="black") + geom_errorbarh(data=dose_means, aes(x=Emean, y=Mmean, xmin=Emean-Esd, xmax=Emean+Esd), colour="black")
base_plot_centers <- base_plot_centers + geom_point(data=dose_means_unlabeled, aes(x=Emean, y=Mmean),fill='red',pch=21,cex=4.0) + geom_errorbar(data=dose_means_unlabeled, aes(x=Emean, y=Mmean, ymin=Mmean-Msd, ymax=Mmean+Msd), colour="black") + geom_errorbarh(data=dose_means_unlabeled, aes(x=Emean, y=Mmean, xmin=Emean-Esd, xmax=Emean+Esd), colour="black")
base_plot_centers

### Plot Labeled and Unlabeled points ###
UMAP_annot_both <- merge(UMAP_values,dose_meta,by=0)
row.names(UMAP_annot_both) <- UMAP_annot_both$Row.names
UMAP_samples <- row.names(UMAP_annot_both)
UMAP_annot_both$Label <- "Unlabeled"
UMAP_annot_both[dose_samples,]$Label <- "Labeled"

ggplot(UMAP_annot_both, aes(x=UMAP_1, y=UMAP_2,colour=Label)) + geom_point(size=0.25) + 
  scale_colour_manual(values = c("black", "red")) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1) 
