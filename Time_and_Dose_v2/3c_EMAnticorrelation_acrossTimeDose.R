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
integated_scale <- readRDS("Integrated_ScaledData.rds")
UMAP_values <- readRDS("Integrated_UMAP_values.rds")

dose_samples <- readRDS("DoseSamples.rds")
time_samples <- readRDS("TimeSamples.rds")

Integrated_E_mat_nsprcomp <- readRDS("Integrated_E_mat_nsprcomp.rds")
Integrated_M_mat_nsprcomp <- readRDS("Integrated_M_mat_nsprcomp.rds")

# Add annotation data to UMAP
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

cor(UMAP_annot[time_samples,]$E_PC1,UMAP_annot[time_samples,]$M_PC1)
cor(UMAP_annot[dose_samples,]$E_PC1,UMAP_annot[dose_samples,]$M_PC1)

cor(UMAP_annot[time_samples,]$E_PC1,UMAP_annot[time_samples,]$M_PC2)
cor(UMAP_annot[dose_samples,]$E_PC1,UMAP_annot[dose_samples,]$M_PC2)

############################################
### E and M Correlation across treatment ###
############################################

### Time Samples ###
time_E1_vs_M1 <- matrix(,nrow=6,ncol=1)
time_E1_vs_M2 <- matrix(,nrow=6,ncol=1)
index = 1
for (i in c(0,1,2,3,4,8)){
  # Get each time treatment
  time_subset <- row.names(integated_meta[integated_meta$Time==i,])
  # E1 vs M1 and M2 correlation
  e1m1 <- cor(Integrated_E_mat_nsprcomp$x[time_subset,1],Integrated_M_mat_nsprcomp$x[time_subset,1])
  e1m2 <- cor(Integrated_E_mat_nsprcomp$x[time_subset,1],Integrated_M_mat_nsprcomp$x[time_subset,2])
  time_E1_vs_M1[index,1] <- e1m1
  time_E1_vs_M2[index,1] <- e1m2
  index = index + 1
}

### Dose Sampels ###
dose_E1_vs_M1 <- matrix(,nrow=8,ncol=1)
dose_E1_vs_M2 <- matrix(,nrow=8,ncol=1)
index = 1
for (i in c(1,2,3,4,8,9,12,13)){
  # Get each dose sample
  dose_subset <- row.names(integated_meta[integated_meta$GBC_values==i,])
  # E1 vs M1 and M2 correlation
  e1m1 <- cor(Integrated_E_mat_nsprcomp$x[dose_subset,1],Integrated_M_mat_nsprcomp$x[dose_subset,1])
  e1m2 <- cor(Integrated_E_mat_nsprcomp$x[dose_subset,1],Integrated_M_mat_nsprcomp$x[dose_subset,2])
  dose_E1_vs_M1[index,1] <- e1m1
  dose_E1_vs_M2[index,1] <- e1m2
  index = index + 1
}

# Make data frames
time_heat_data <- as.data.frame(cbind(time_E1_vs_M1,time_E1_vs_M2))
colnames(time_heat_data) <- c("E1_vs_M1","E1_vs_M2")
time_heat_data$E1_vs_M1 <- as.numeric(time_heat_data$E1_vs_M1)
time_heat_data$E1_vs_M2 <- as.numeric(time_heat_data$E1_vs_M2)

dose_heat_data <- as.data.frame(cbind(dose_E1_vs_M1,dose_E1_vs_M2))
colnames(dose_heat_data) <- c("E1_vs_M1","E1_vs_M2")
dose_heat_data$E1_vs_M1 <- as.numeric(dose_heat_data$E1_vs_M1)
dose_heat_data$E1_vs_M2 <- as.numeric(dose_heat_data$E1_vs_M2)

### Make heatmaps ###
library(reshape2)

# Make data long format
time_heat_data_long <- melt(time_heat_data)

# Lable data
time_heat_data_long$Y <- "M2"
time_heat_data_long[time_heat_data_long$variable=="E1_vs_M1",]$Y <- "M1"
time_heat_data_long$X <- factor(c("0","1","2","3","4","8","0","1","2","3","4","8"), levels = c("0","1","2","3","4","8"))

# Make Heatmap
ggplot(time_heat_data_long, aes(X, Y, fill= value)) + geom_tile() + 
  scale_fill_gradient2(midpoint=-0.3, low="darkblue", mid="white", high="red", space ="Lab",limits=c(-0.85,0)) + theme_bw() + labs(fill="Corr.") + xlab("Time (days)") + ylab("PCs") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 11),aspect.ratio=1,axis.text.x = element_text(angle = 90))

# Make data long format
dose_heat_data_long <- melt(dose_heat_data)

# Label data
dose_heat_data_long$Y <- "M2"
dose_heat_data_long[dose_heat_data_long$variable=="E1_vs_M1",]$Y <- "M1"
dose_heat_data_long$X <-factor(c("0","12.5","25","50","100","200","400","800","0","12.5","25","50","100","200","400","800"), levels = c("0","12.5","25","50","100","200","400","800"))

# Make heatmap
ggplot(dose_heat_data_long, aes(X, Y, fill= value)) + geom_tile() + 
  scale_fill_gradient2(midpoint=-0.3, low="darkblue", mid="white", high="red", space ="Lab",limits=c(-0.85,0)) + theme_bw() + labs(fill="Corr.") + xlab("Dose (pm)") + ylab("PCs") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 11),aspect.ratio=1,axis.text.x = element_text(angle = 90))

#############################
### Without Doublets Time ###
#############################
# Unomment and run after 7X4

#time_doublets <- readRDS("time_doublets.rds")
#time_wo_doublets <- setdiff(time_samples,time_doublets)

#integated_meta_drop_doublet <- integated_meta[time_wo_doublets,]

### Time Samples ###
#time_E1_vs_M1_drop <- matrix(,nrow=6,ncol=1)
#time_E1_vs_M2_drop <- matrix(,nrow=6,ncol=1)
#index = 1
#for (i in c(0,1,2,3,4,8)){
  # Get time samples
#  time_subset <- row.names(integated_meta_drop_doublet[integated_meta_drop_doublet$Time==i,])
  
  # E1 vs M1 and M2 correlation
#  e1m1 <- cor(Integrated_E_mat_nsprcomp$x[time_subset,1],Integrated_M_mat_nsprcomp$x[time_subset,1])
#  e1m2 <- cor(Integrated_E_mat_nsprcomp$x[time_subset,1],Integrated_M_mat_nsprcomp$x[time_subset,2])
#  time_E1_vs_M1_drop[index,1] <- e1m1
#  time_E1_vs_M2_drop[index,1] <- e1m2
#  index = index + 1
#}
#cor(time_E1_vs_M1_drop,time_E1_vs_M1)
#cor(time_E1_vs_M2_drop,time_E1_vs_M2)
