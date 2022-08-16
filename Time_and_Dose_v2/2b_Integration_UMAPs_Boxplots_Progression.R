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

# Progression Plot
gg_means_plot <- function(plot_data,color_n,x_down,x_up,y_down,y_up){
  manual_colors <- hue_pal()(color_n)
  
  # Base plot of E/M scores
  base_plot <- ggplot(plot_data, aes(x=Egenes, y=Mgenes)) + geom_point(size=0.1, alpha=0.75, color='lightgrey') + stat_density_2d(aes(alpha=(..level..)^1, fill=key_label),geom="polygon",show.legend=F) + scale_fill_manual(values=manual_colors) + scale_alpha_continuous(range=c(0,0.4)) + theme_bw() + theme(axis.text = element_text(size = 10), axis.title=element_text(size=11), aspect.ratio=1)
  #base_plot
  
  # Define means and sds
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
  
  # Use these values to size the head of the error bar
  E_20th = (x_up - x_down)/20
  M_20th  = (y_up - y_down)/20
  
  # Add Error Bars
  base_plot_centers <- base_plot + geom_point(data=data_ceneters, aes(x=Emean, y=Mmean),fill=manual_colors,pch=21,cex=4.0) + geom_errorbar(data=data_ceneters, aes(x=Emean, y=Mmean, ymin=Mmean-Msd, ymax=Mmean+Msd), colour="black", width=E_20th) + geom_errorbarh(data=data_ceneters, aes(x=Emean, y=Mmean, xmin=Emean-Esd, xmax=Emean+Esd), colour="black", height=M_20th) + xlim(x_down,x_up) + ylim(y_down,y_up) + xlab("E PC") + ylab("M PC")
  #base_plot_centers
  return(base_plot_centers)
}

# Progression Plot Means 
gg_means <- function(plot_data,color_n,x_down,x_up,y_down,y_up){
  
  # Define means and sds
  data_Emeans <- r1<-with(plot_data, tapply(Egenes, key, mean))
  data_Mmeans <- r1<-with(plot_data, tapply(Mgenes, key, mean))
  data_Esds <- r1<-with(plot_data, tapply(Egenes, key, sd))
  data_Msds <- r1<-with(plot_data, tapply(Mgenes, key, sd))
  
  # Return Data
  data_ceneters <- as.data.frame(t(rbind(data_Emeans,data_Mmeans,data_Esds,data_Msds)))
  names(data_ceneters) <- c("Emean","Mmean","Esd","Msd")
  return(data_ceneters)
}

########################
### Generate Figures ###
########################

### Load Processed Data ###

integated_meta <- readRDS("Integrated_MetaData.rds")
UMAP_values <- readRDS("Integrated_UMAP_values.rds")

dose_samples <- readRDS("DoseSamples.rds")
time_samples <- readRDS("TimeSamples.rds")

Integrated_E_mat_nsprcomp <- readRDS("_oldRDS/Integrated_E_mat_nsprcomp.rds")
Integrated_M_mat_nsprcomp <- readRDS("_oldRDS/Integrated_M_mat_nsprcomp.rds")

### PC UMAP Plots ###

# Add annotation data to UMAP_values
UMAP_annot <- merge(UMAP_values,integated_meta,by=0)
row.names(UMAP_annot) <- UMAP_annot$Row.names
UMAP_annot$GBC_pM_factor <- factor(UMAP_annot$GBC_pM)
UMAP_annot$Time_factor <- factor(UMAP_annot$Time)

# Add PC values
UMAP_samples <- row.names(UMAP_annot)
UMAP_annot$E_PC1 <- Integrated_E_mat_nsprcomp$x[UMAP_samples,1]

UMAP_annot$M_PC1 <- Integrated_M_mat_nsprcomp$x[UMAP_samples,1]
UMAP_annot$M_PC2 <- Integrated_M_mat_nsprcomp$x[UMAP_samples,2]

### Correlation within dose and time
cor(UMAP_annot[time_samples,]$Time,UMAP_annot[time_samples,]$M_PC1) # 0.477845
cor(UMAP_annot[time_samples,]$Time,UMAP_annot[time_samples,]$M_PC2) # 0.1552711

cor(UMAP_annot[dose_samples,]$GBC_pM,UMAP_annot[dose_samples,]$M_PC1) # 0.416201
cor(UMAP_annot[dose_samples,]$GBC_pM,UMAP_annot[dose_samples,]$M_PC2) # 0.0558633


cor(UMAP_annot[time_samples,]$E_PC1,UMAP_annot[time_samples,]$M_PC1) # -0.3066327
cor(UMAP_annot[time_samples,]$E_PC1,UMAP_annot[time_samples,]$M_PC2) # -0.6662105

cor(UMAP_annot[dose_samples,]$E_PC1,UMAP_annot[dose_samples,]$M_PC1) # -0.4948821
cor(UMAP_annot[dose_samples,]$E_PC1,UMAP_annot[dose_samples,]$M_PC2) # -0.7645159

### Basic plots

# Dosage and Time Plot
ggplot(UMAP_annot, aes(x = UMAP_1, y = UMAP_2, colour=GBC_pM_factor)) +
  geom_point(size=0.5) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot(UMAP_annot, aes(x = UMAP_1, y = UMAP_2, colour=Time_factor)) +
  geom_point(size=0.5) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot(UMAP_annot, aes(x = UMAP_1, y = UMAP_2, colour=seurat_clusters)) +
  geom_point(size=0.5) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# E PCs
ggplot(UMAP_annot, aes(x = UMAP_1, y = UMAP_2, colour=E_PC1)) +
  geom_point(size=0.5) + theme_bw() + scale_color_gradient2(midpoint=0, low="darkblue", mid="white",
                                                            high="red", space ="Lab") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# M PCs
ggplot(UMAP_annot, aes(x = UMAP_1, y = UMAP_2, colour=M_PC1)) +
  geom_point(size=0.5) + theme_bw() + scale_color_gradient2(midpoint=0, low="darkblue", mid="white",
                                                            high="red", space ="Lab") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot(UMAP_annot, aes(x = UMAP_1, y = UMAP_2, colour=M_PC2)) +
  geom_point(size=0.5) + theme_bw() + scale_color_gradient2(midpoint=0, low="darkblue", mid="white",
                                                            high="red", space ="Lab") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### Box Plots ###
boxplot(UMAP_annot$E_PC1 ~ UMAP_annot$Time)
boxplot(UMAP_annot$E_PC1 ~ UMAP_annot$GBC_pM)

boxplot(UMAP_annot$M_PC1 ~ UMAP_annot$Time)
boxplot(UMAP_annot$M_PC1 ~ UMAP_annot$GBC_pM)

boxplot(UMAP_annot$M_PC2 ~ UMAP_annot$Time)
boxplot(UMAP_annot$M_PC2 ~ UMAP_annot$GBC_pM)

################################################
### Seperate Plots for Time and Dose samples ###
################################################

UMAP_annot_dose <- UMAP_annot[dose_samples,]
UMAP_annot_time <- UMAP_annot[time_samples,]

# UMAPs
UMAP_batch <- ggplot(UMAP_annot, aes(x = UMAP_1, y = UMAP_2, colour=Batch)) +
  geom_point(size=0.25) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1,legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 3)),guide_colorbar(direction="horizontal"))

UMAP_annot_time$Time <- as.factor(UMAP_annot_time$Time)
UMAP_Time <- ggplot(UMAP_annot_time, aes(x = UMAP_1, y = UMAP_2, colour=Time)) +
  geom_point(size=0.25) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1,legend.position = "bottom") + 
  guides(color = guide_legend(override.aes = list(size = 3)),guide_colorbar(direction="horizontal"))

# Supplemental
ggplot(UMAP_annot_time, aes(x = UMAP_1, y = UMAP_2, colour=seurat_clusters)) +
  geom_point(size=0.25) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1,legend.position = "bottom") + 
  guides(color = guide_legend(override.aes = list(size = 3)),guide_colorbar(direction="horizontal"))

UMAP_annot_dose$GBC_pM <- as.factor(UMAP_annot_dose$GBC_pM)
UMAP_Dose <- ggplot(UMAP_annot_dose, aes(x = UMAP_1, y = UMAP_2, colour=GBC_pM)) +
  geom_point(size=0.25) + theme_bw() + labs(color="pM") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1,legend.position = "bottom") + 
  guides(color = guide_legend(override.aes = list(size = 3)),guide_colorbar(direction="horizontal"))

# Supplemental
ggplot(UMAP_annot_dose, aes(x = UMAP_1, y = UMAP_2, colour=seurat_clusters)) +
  geom_point(size=0.25) + theme_bw() + labs(color="pM") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1,legend.position = "bottom") + 
  guides(color = guide_legend(override.aes = list(size = 3)),guide_colorbar(direction="horizontal"))


par(mfrow = c(1, 3))
grid.arrange(UMAP_batch,UMAP_Time,UMAP_Dose,nrow=1)
# Figure 2, Row 1

# Progression
plot_data <- as.data.frame(cbind(UMAP_annot_dose$E_PC1,UMAP_annot_dose$M_PC1,UMAP_annot_dose$GBC_pM))
colnames(plot_data) <- c("Egenes","Mgenes","key")
plot_data$key_label <- factor(plot_data$key)

manual_colors_dose <- hue_pal()(8)
base_plot_centers_dose <- gg_means_plot(plot_data,8,-7.5,10,-7.5,10)
base_plot_centers_dose

plot_data <- as.data.frame(cbind(UMAP_annot_time$E_PC1,UMAP_annot_time$M_PC1,UMAP_annot_time$Time))
colnames(plot_data) <- c("Egenes","Mgenes","key")
plot_data$key_label <- factor(plot_data$key)

manual_colors_time <- hue_pal()(6)
base_plot_centers_time <- gg_means_plot(plot_data,6,-7.5,10,-7.5,10)
base_plot_centers_time

par(mfrow = c(1, 2))
grid.arrange(base_plot_centers_dose,base_plot_centers_time,nrow=1)
# Figure 2, Row 2

# Box plots
boxplot(UMAP_annot_time$E_PC1 ~ UMAP_annot_time$Time,col=hue_pal()(6),xlab="Time",ylab="E-score")
boxplot(UMAP_annot_time$M_PC1 ~ UMAP_annot_time$Time,col=hue_pal()(6),xlab="TIme",ylab="E-score")

boxplot(UMAP_annot_dose$E_PC1 ~ UMAP_annot_dose$GBC_pM,col=hue_pal()(8),xlab="Dosage",ylab="E-score")
boxplot(UMAP_annot_dose$M_PC1 ~ UMAP_annot_dose$GBC_pM,col=hue_pal()(8),xlab="Dosage",ylab="M-score")

### Seperation Tables ### 
mwu_greater <- function(A, B){
  mwu_test <- wilcox.test(A, B,alternative='greater')
  auc_roc <- mwu_test$statistic/(length(A) * length(B))
  #print(mwu_test$p.value)
  #print(auc_roc)
  return (c(mwu_test$p.value,auc_roc))
}

# Go through pairs of dose samples
dose_results <- matrix(,nrow=7,ncol=8)
dose_values <- c(1,2,3,4,8,9,12,13)
for (index in 1:7){
  
  # Get each value and the next
  v1 = dose_values[index]
  v2 = dose_values[index+1]
  
  # Get samples
  lesser_treatmeant = row.names(integated_meta[integated_meta$GBC_values == v1,])
  greater_treatment = row.names(integated_meta[integated_meta$GBC_values == v2,])
  
  # Check
  #print(paste0(v1," ",v2))
  #print(table(dose_meta[lesser_treatmeant,]$GBC_pM))
  #print(table(dose_meta[greater_treatment,]$GBC_pM))
  
  dose_results[index,1] <- names(table(integated_meta[lesser_treatmeant,]$GBC_pM))
  dose_results[index,2] <- names(table(integated_meta[greater_treatment,]$GBC_pM))
  
  # Run Tests for E1, M1, and M2
  dose_results[index,3:4] <- mwu_greater(Integrated_E_mat_nsprcomp$x[lesser_treatmeant,1],Integrated_E_mat_nsprcomp$x[greater_treatment,1]) # use lesser for E
  dose_results[index,5:6] <- mwu_greater(Integrated_M_mat_nsprcomp$x[greater_treatment,1],Integrated_M_mat_nsprcomp$x[lesser_treatmeant,1])
  dose_results[index,7:8] <- mwu_greater(Integrated_M_mat_nsprcomp$x[greater_treatment,2],Integrated_M_mat_nsprcomp$x[lesser_treatmeant,2])
}
dose_results

# Go through pairs of time samples
time_results <- matrix(,nrow=5,ncol=8)
time_values <- c(0,1,2,3,4,8)
for (index in 1:5){
  
  # Get each value and the next
  v1 = time_values[index]
  v2 = time_values[index+1]
  
  # Get samples
  lesser_treatmeant = row.names(integated_meta[integated_meta$Time == v1,])
  greater_treatment = row.names(integated_meta[integated_meta$Time == v2,])
  
  # Check
  #print(paste0(v1," ",v2))
  #print(table(integated_meta[lesser_treatmeant,]$GBC_pM))
  #print(table(integated_meta[greater_treatment,]$GBC_pM))
  
  time_results[index,1] <- names(table(integated_meta[lesser_treatmeant,]$Time))
  time_results[index,2] <- names(table(integated_meta[greater_treatment,]$Time))
  
  # Run Tests for E1, M1, and M2
  time_results[index,3:4] <- mwu_greater(Integrated_E_mat_nsprcomp$x[lesser_treatmeant,1],Integrated_E_mat_nsprcomp$x[greater_treatment,1]) # use lesser for E
  time_results[index,5:6] <- mwu_greater(Integrated_M_mat_nsprcomp$x[greater_treatment,1],Integrated_M_mat_nsprcomp$x[lesser_treatmeant,1])
  time_results[index,7:8] <- mwu_greater(Integrated_M_mat_nsprcomp$x[greater_treatment,2],Integrated_M_mat_nsprcomp$x[lesser_treatmeant,2])
}
time_results

### Run Time without doublets ###
time_doublets <- readRDS("time_doublets.rds")
integated_meta_drop_doublet <- integated_meta[setdiff(row.names(integated_meta),time_doublets),]

time_drop_results <- matrix(,nrow=5,ncol=8)
time_values <- c(0,1,2,3,4,8)
for (index in 1:5){
  v1 = time_values[index]
  v2 = time_values[index+1]
  lesser_treatmeant = row.names(integated_meta_drop_doublet[integated_meta_drop_doublet$Time == v1,])
  greater_treatment = row.names(integated_meta_drop_doublet[integated_meta_drop_doublet$Time == v2,])
  
  # Check
  #print(paste0(v1," ",v2))
  #print(table(dose_meta[lesser_treatmeant,]$GBC_pM))
  #print(table(dose_meta[greater_treatment,]$GBC_pM))
  
  time_drop_results[index,1] <- names(table(integated_meta_drop_doublet[lesser_treatmeant,]$Time))
  time_drop_results[index,2] <- names(table(integated_meta_drop_doublet[greater_treatment,]$Time))
  
  time_drop_results[index,3:4] <- mwu_greater(Integrated_E_mat_nsprcomp$x[lesser_treatmeant,1],Integrated_E_mat_nsprcomp$x[greater_treatment,1]) # use lesser for E
  time_drop_results[index,5:6] <- mwu_greater(Integrated_M_mat_nsprcomp$x[greater_treatment,1],Integrated_M_mat_nsprcomp$x[lesser_treatmeant,1])
  time_drop_results[index,7:8] <- mwu_greater(Integrated_M_mat_nsprcomp$x[greater_treatment,2],Integrated_M_mat_nsprcomp$x[lesser_treatmeant,2])
}
time_drop_results

### Correlation ###
time_wo_doublets <- setdiff(time_samples,time_doublets)
cor(UMAP_annot[dose_samples,]$E_PC1,UMAP_annot[dose_samples,]$M_PC1)
cor(UMAP_annot[time_samples,]$E_PC1,UMAP_annot[time_samples,]$M_PC1)
cor(UMAP_annot[time_wo_doublets,]$E_PC1,UMAP_annot[time_wo_doublets,]$M_PC1)

cor(UMAP_annot[dose_samples,]$E_PC1,UMAP_annot[dose_samples,]$M_PC2)
cor(UMAP_annot[time_samples,]$E_PC1,UMAP_annot[time_samples,]$M_PC2)
cor(UMAP_annot[time_wo_doublets,]$E_PC1,UMAP_annot[time_wo_doublets,]$M_PC2)