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

### Time Samples vs Genes ###
time_genes_results <- matrix(,nrow=6,ncol=9)
index = 1
for (i in c(0,1,2,3,4,8)){
  # Get each time subset
  time_subset <- row.names(integated_meta[integated_meta$Time==i,])
  
  # Correlation of genes and E
  time_genes_results[index,1] <- cor(Integrated_E_mat_nsprcomp$x[time_subset,1],UMAP_annot[time_subset,]$VIM)
  time_genes_results[index,2] <- cor(Integrated_E_mat_nsprcomp$x[time_subset,1],UMAP_annot[time_subset,]$FN1)
  time_genes_results[index,3] <- cor(Integrated_E_mat_nsprcomp$x[time_subset,1],UMAP_annot[time_subset,]$EPCAM)
  
  # Correlation of genes and M1
  time_genes_results[index,4] <- cor(Integrated_M_mat_nsprcomp$x[time_subset,1],UMAP_annot[time_subset,]$VIM)
  time_genes_results[index,5] <- cor(Integrated_M_mat_nsprcomp$x[time_subset,1],UMAP_annot[time_subset,]$FN1)
  time_genes_results[index,6] <- cor(Integrated_M_mat_nsprcomp$x[time_subset,1],UMAP_annot[time_subset,]$EPCAM)
  
  # Correlation of genes and M2
  time_genes_results[index,7] <- cor(Integrated_M_mat_nsprcomp$x[time_subset,2],UMAP_annot[time_subset,]$VIM)
  time_genes_results[index,8] <- cor(Integrated_M_mat_nsprcomp$x[time_subset,2],UMAP_annot[time_subset,]$FN1)
  time_genes_results[index,9] <- cor(Integrated_M_mat_nsprcomp$x[time_subset,2],UMAP_annot[time_subset,]$EPCAM)
  index = index+1
}

### Line graps
# VIM
plot(c(0,1,2,3,4,8),time_genes_results[,1],col="red",ylim=c(-0.65,0.65),type="l",main="VIM Time",xlab="days",ylab="correlation")
lines(c(0,1,2,3,4,8),time_genes_results[,4],col="green")
lines(c(0,1,2,3,4,8),time_genes_results[,7],col="blue")

# FN1
plot(c(0,1,2,3,4,8),time_genes_results[,2],col="red",ylim=c(-0.65,0.65),type="l",main="FN1 Time",xlab="days",ylab="correlation")
lines(c(0,1,2,3,4,8),time_genes_results[,5],col="green")
lines(c(0,1,2,3,4,8),time_genes_results[,8],col="blue")

# EPCAM
plot(c(0,1,2,3,4,8),time_genes_results[,3],col="red",ylim=c(-0.65,0.65),type="l",main="EPCAM Time",xlab="days",ylab="correlation")
lines(c(0,1,2,3,4,8),time_genes_results[,6],col="green")
lines(c(0,1,2,3,4,8),time_genes_results[,9],col="blue")

# UMAPs
ggplot(UMAP_annot, aes(x = UMAP_1, y = UMAP_2, colour=FN1)) +
  geom_point(size=0.5) + theme_bw() + scale_color_gradient2(midpoint=0, low="darkblue", mid="white",
                                                            high="red", space ="Lab",limits=c(-3,3)) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 18),aspect.ratio=1)

ggplot(UMAP_annot[time_samples,], aes(x = UMAP_1, y = UMAP_2, colour=FN1)) +
  geom_point(size=0.5) + theme_bw() + scale_color_gradient2(midpoint=0, low="darkblue", mid="white",
                                                            high="red", space ="Lab",limits=c(-3,3)) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 18),aspect.ratio=1)

ggplot(UMAP_annot[dose_samples,], aes(x = UMAP_1, y = UMAP_2, colour=FN1)) +
  geom_point(size=0.5) + theme_bw() + scale_color_gradient2(midpoint=0, low="darkblue", mid="white",
                                                            high="red", space ="Lab",limits=c(-3,3)) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 18),aspect.ratio=1)

### Dose Samples vs Genes ###
dose_genes_results <- matrix(,nrow=8,ncol=9)
index = 1
for (i in c(1,2,3,4,8,9,12,13)){
  # Get each dose
  dose_subset <- row.names(integated_meta[integated_meta$GBC_values==i,])
  
  # Correlation of genes and E
  dose_genes_results[index,1] <- cor(Integrated_E_mat_nsprcomp$x[dose_subset,1],UMAP_annot[dose_subset,]$VIM)
  dose_genes_results[index,2] <- cor(Integrated_E_mat_nsprcomp$x[dose_subset,1],UMAP_annot[dose_subset,]$FN1)
  dose_genes_results[index,3] <- cor(Integrated_E_mat_nsprcomp$x[dose_subset,1],UMAP_annot[dose_subset,]$EPCAM)
  
  # Correlation of genes and M1
  dose_genes_results[index,4] <- cor(Integrated_M_mat_nsprcomp$x[dose_subset,1],UMAP_annot[dose_subset,]$VIM)
  dose_genes_results[index,5] <- cor(Integrated_M_mat_nsprcomp$x[dose_subset,1],UMAP_annot[dose_subset,]$FN1)
  dose_genes_results[index,6] <- cor(Integrated_M_mat_nsprcomp$x[dose_subset,1],UMAP_annot[dose_subset,]$EPCAM)
  
  # Correlation of genes and M2
  dose_genes_results[index,7] <- cor(Integrated_M_mat_nsprcomp$x[dose_subset,2],UMAP_annot[dose_subset,]$VIM)
  dose_genes_results[index,8] <- cor(Integrated_M_mat_nsprcomp$x[dose_subset,2],UMAP_annot[dose_subset,]$FN1)
  dose_genes_results[index,9] <- cor(Integrated_M_mat_nsprcomp$x[dose_subset,2],UMAP_annot[dose_subset,]$EPCAM)
  index = index+1
}

### Line graphs

# VIM
plot(c(0,12.5,25,50,100,200,400,800),dose_genes_results[,1],col="red",ylim=c(-0.65,0.65),type="l",main="VIM dose",xlab="days",ylab="correlation")
lines(c(0,12.5,25,50,100,200,400,800),dose_genes_results[,4],col="green")
lines(c(0,12.5,25,50,100,200,400,800),dose_genes_results[,7],col="blue")

# FN1
plot(c(0,12.5,25,50,100,200,400,800),dose_genes_results[,2],col="red",ylim=c(-0.65,0.65),type="l",main="FN1 dose",xlab="days",ylab="correlation")
lines(c(0,12.5,25,50,100,200,400,800),dose_genes_results[,5],col="green")
lines(c(0,12.5,25,50,100,200,400,800),dose_genes_results[,8],col="blue")

# EPCAM
plot(c(0,12.5,25,50,100,200,400,800),dose_genes_results[,3],col="red",ylim=c(-0.65,0.65),type="l",main="EPCAM dose",xlab="days",ylab="correlation")
lines(c(0,12.5,25,50,100,200,400,800),dose_genes_results[,6],col="green")
lines(c(0,12.5,25,50,100,200,400,800),dose_genes_results[,9],col="blue")

#############################
### Without Doublets Time ###
#############################

time_doublets <- readRDS("time_doublets.rds")
time_wo_doublets <- setdiff(time_samples,time_doublets)

integated_meta_drop_doublet <- integated_meta[time_wo_doublets,]

### Time Samples ###
time_E1_vs_M1_drop <- matrix(,nrow=6,ncol=1)
time_E1_vs_M2_drop <- matrix(,nrow=6,ncol=1)
index = 1
for (i in c(0,1,2,3,4,8)){
  # Get time samples
  time_subset <- row.names(integated_meta_drop_doublet[integated_meta_drop_doublet$Time==i,])
  
  # E1 vs M1 and M2 correlation
  e1m1 <- cor(Integrated_E_mat_nsprcomp$x[time_subset,1],Integrated_M_mat_nsprcomp$x[time_subset,1])
  e1m2 <- cor(Integrated_E_mat_nsprcomp$x[time_subset,1],Integrated_M_mat_nsprcomp$x[time_subset,2])
  time_E1_vs_M1_drop[index,1] <- e1m1
  time_E1_vs_M2_drop[index,1] <- e1m2
  index = index + 1
}
cor(time_E1_vs_M1_drop,time_E1_vs_M1)
cor(time_E1_vs_M2_drop,time_E1_vs_M2)

time_genes_results_drop <- matrix(,nrow=6,ncol=9)
index = 1
for (i in c(0,1,2,3,4,8)){
  
  # Get time samples
  time_subset <- row.names(integated_meta_drop_doublet[integated_meta_drop_doublet$Time==i,])
  
  # # Correlation of genes and E
  time_genes_results_drop[index,1] <- cor(Integrated_E_mat_nsprcomp$x[time_subset,1],UMAP_annot[time_subset,]$VIM)
  time_genes_results_drop[index,2] <- cor(Integrated_E_mat_nsprcomp$x[time_subset,1],UMAP_annot[time_subset,]$FN1)
  time_genes_results_drop[index,3] <- cor(Integrated_E_mat_nsprcomp$x[time_subset,1],UMAP_annot[time_subset,]$EPCAM)
  
  # # Correlation of genes and M1
  time_genes_results_drop[index,4] <- cor(Integrated_M_mat_nsprcomp$x[time_subset,1],UMAP_annot[time_subset,]$VIM)
  time_genes_results_drop[index,5] <- cor(Integrated_M_mat_nsprcomp$x[time_subset,1],UMAP_annot[time_subset,]$FN1)
  time_genes_results_drop[index,6] <- cor(Integrated_M_mat_nsprcomp$x[time_subset,1],UMAP_annot[time_subset,]$EPCAM)
  
  # # Correlation of genes and M2
  time_genes_results_drop[index,7] <- cor(Integrated_M_mat_nsprcomp$x[time_subset,2],UMAP_annot[time_subset,]$VIM)
  time_genes_results_drop[index,8] <- cor(Integrated_M_mat_nsprcomp$x[time_subset,2],UMAP_annot[time_subset,]$FN1)
  time_genes_results_drop[index,9] <- cor(Integrated_M_mat_nsprcomp$x[time_subset,2],UMAP_annot[time_subset,]$EPCAM)
  index = index+1
}
cor(time_genes_results_drop[,1],time_genes_results[,1])
cor(time_genes_results_drop[,2],time_genes_results[,2])
cor(time_genes_results_drop[,3],time_genes_results[,3])
cor(time_genes_results_drop[,4],time_genes_results[,4])
cor(time_genes_results_drop[,5],time_genes_results[,5])
cor(time_genes_results_drop[,6],time_genes_results[,6])
cor(time_genes_results_drop[,7],time_genes_results[,7])
cor(time_genes_results_drop[,8],time_genes_results[,8])
cor(time_genes_results_drop[,9],time_genes_results[,9])
