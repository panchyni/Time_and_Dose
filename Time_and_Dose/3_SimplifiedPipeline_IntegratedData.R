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

#################
### Functions ###
#################

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
  
  base_plot_centers <- base_plot + geom_point(data=data_ceneters, aes(x=Emean, y=Mmean),fill=manual_colors,pch=21,cex=4.0) + geom_errorbar(data=data_ceneters, aes(x=Emean, y=Mmean, ymin=Mmean-Msd, ymax=Mmean+Msd), colour="black", width=E_20th) + geom_errorbarh(data=data_ceneters, aes(x=Emean, y=Mmean, xmin=Emean-Esd, xmax=Emean+Esd), colour="black", height=M_20th) + xlim(x_down,x_up) + ylim(y_down,y_up)
  #base_plot_centers
  return(base_plot_centers)
}

gg_means <- function(plot_data){
  data_Emeans <- r1<-with(plot_data, tapply(Egenes, key, mean))
  data_Mmeans <- r1<-with(plot_data, tapply(Mgenes, key, mean))
  data_Esds <- r1<-with(plot_data, tapply(Egenes, key, sd))
  data_Msds <- r1<-with(plot_data, tapply(Mgenes, key, sd))
  
  data_ceneters <- as.data.frame(t(rbind(data_Emeans,data_Mmeans,data_Esds,data_Msds)))
  names(data_ceneters) <- c("Emean","Mmean","Esd","Msd")
  
  return(data_ceneters)
}

######################
### Integrated All ###
######################

integration.combined.All <- readRDS('SeuratData/Seurat_integrated_TimeAndCon_CellCycle_hg38.rds')

### Fix annotation ###

# Drop NAs 
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
# 10 samples are present in barcodes, these are the NAs dropped during processning

p1 <- DimPlot(integration.combined.All, reduction = "umap",group.by = "Batch")
p2 <- DimPlot(integration.combined.All, reduction = "umap",group.by = "Time")
p3 <- DimPlot(integration.combined.All, reduction = "umap", label = TRUE, repel = TRUE)
p4 <- DimPlot(integration.combined.All, reduction = "umap",group.by = "Phase")

p1
p2
p3
p4

### Figure 2 ###
pdf(file="Figures/Panels/Fig2_Integration_UMAP_Batch.pdf")
p1
dev.off()

# Custom UMAP Plots #
UMAP_values <- integration.combined.All[["umap"]]@cell.embeddings
UMAP_annot <- merge(UMAP_values,integration.combined.All@meta.data,by=0)

ggplot(UMAP_annot, aes(x = UMAP_1, y = UMAP_2, colour=Batch)) +
  geom_density_2d(bins=6,alpha=0.5,size=1) + theme_bw()

UMAP_annot_Time <- UMAP_annot[UMAP_annot$Time < 100,]
UMAP_annot_Con <- UMAP_annot[UMAP_annot$GBC_pM > -1,]
UMAP_annot_Time$Time <- as.factor(UMAP_annot_Time$Time)
UMAP_annot_Con$GBC_pM <- as.factor(UMAP_annot_Con$GBC_pM)
dim(UMAP_annot_Time)
dim(UMAP_annot_Con)

### Figure 2 ###
pdf(file="Figures/Panels/Fig2_Integration_UMAP_Time.pdf")
ggplot(UMAP_annot_Time, aes(x = UMAP_1, y = UMAP_2, colour=Time)) +
  geom_point(size=0.5) + theme_bw()
dev.off()

### Figure 2 ###
pdf(file="Figures/Panels/Fig2_Integration_UMAP_Con.pdf")
ggplot(UMAP_annot_Con, aes(x = UMAP_1, y = UMAP_2, colour=GBC_pM)) +
  geom_point(size=0.5) + theme_bw()
dev.off()

###############
### Scoring ###
###############

seurat_combat_integratd_norm <- integration.combined.All@assays$integrated@data
seurat_combat_integratd_norm_SampleTotal <- colSums(seurat_combat_integratd_norm)
integration.combined.All$IntegratedTotal <- seurat_combat_integratd_norm_SampleTotal

seurat_combo_All_mat <- GetAssayData(object = integration.combined.All, slot = "scale.data")
dim(seurat_combo_All_mat)

# Genes Sets #
EMT_genes <- read.table("EMTGenesUpdateV2.txt",header=TRUE)
E_genes <- subset(EMT_genes, EMT_genes$Annotation == 'E')
M_genes <- subset(EMT_genes, EMT_genes$Annotation == 'M')
Gene_Sets <- list('Egenes'=E_genes$Gene,'Mgenes'=M_genes$Gene)

seurat_combo_E <- ScaleData(integration.combined.All, vars.to.regress=c("S.Score","G2M.Score"),verbose = FALSE, features=E_genes$Gene)
seurat_combo_M <- ScaleData(integration.combined.All, vars.to.regress=c("S.Score","G2M.Score"), verbose = FALSE, features=M_genes$Gene)

#seurat_combo_E <- ScaleData(integration.combined.All, vars.to.regress=c("S.Score","G2M.Score","IntegratedTotal"),verbose = FALSE, features=E_genes$Gene)
#seurat_combo_M <- ScaleData(integration.combined.All, vars.to.regress=c("S.Score","G2M.Score","IntegratedTotal"), verbose = FALSE, features=M_genes$Gene)

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

# Time #
plot_data <- as.data.frame(cbind(seurat_combo_E_mat_nsprcomp_annot[Time_Samples,]$PC1,seurat_combo_M_mat_nsprcomp_annot[Time_Samples,]$PC1,seurat_combo_M_mat_nsprcomp_annot[Time_Samples,]$Time))
colnames(plot_data) <- c("Egenes","Mgenes","key")
plot_data$key_label <- factor(plot_data$key)

manual_colors <- hue_pal()(6)
base_plot_centers <- gg_means_plot(plot_data,6,-7.5,10,-7.5,10)
base_plot_centers

### Figure 2 ###
pdf(file="Figures/Panels/Fig2_Integration_EMSpace_Progression_Time.pdf")
base_plot_centers
dev.off()

### Figure 2 ###
pdf(file="Figures/Panels/Fig2_Integration_EMSpace_Boxplots_Time.pdf")
par(mfrow=c(1,2))
boxplot(seurat_combo_E_mat_nsprcomp_annot[Time_Samples,]$PC1 ~ seurat_combo_E_mat_nsprcomp_annot[Time_Samples,]$Time,col=manual_colors,xlab="Dosage",ylab="E-score")
boxplot(seurat_combo_M_mat_nsprcomp_annot[Time_Samples,]$PC1 ~ seurat_combo_M_mat_nsprcomp_annot[Time_Samples,]$Time,col=manual_colors,xlab="Dosage",ylab="M-score")
dev.off()

# Concentration #
plot_data <- as.data.frame(cbind(seurat_combo_E_mat_nsprcomp_annot[ConGood_Samples,]$PC1,seurat_combo_M_mat_nsprcomp_annot[ConGood_Samples,]$PC1,seurat_combo_M_mat_nsprcomp_annot[ConGood_Samples,]$GBC_pM))
colnames(plot_data) <- c("Egenes","Mgenes","key")
plot_data$key_label <- factor(plot_data$key)

manual_colors <- hue_pal()(8)
base_plot_centers <- gg_means_plot(plot_data,8,-7.5,10,-7.5,10)
base_plot_centers

### Figure 2 ###
pdf(file="Figures/Panels/Fig2_Integration_EMSpace_Progression_Con.pdf")
base_plot_centers
dev.off()

### Figure 2 ###
pdf(file="Figures/Panels/Fig2_Integration_EMSpace_Boxplots_Con.pdf")
par(mfrow=c(1,2))
boxplot(seurat_combo_E_mat_nsprcomp_annot[ConGood_Samples,]$PC1 ~ seurat_combo_E_mat_nsprcomp_annot[ConGood_Samples,]$GBC_pM,col=manual_colors,xlab="Dosage",ylab="E-score")
boxplot(seurat_combo_M_mat_nsprcomp_annot[ConGood_Samples,]$PC1 ~ seurat_combo_M_mat_nsprcomp_annot[ConGood_Samples,]$GBC_pM,col=manual_colors,xlab="Dosage",ylab="M-score")
dev.off()

### Correlation of E/M Scores ###
cor(seurat_combo_E_mat_nsprcomp_annot[Time_Samples,]$PC1,seurat_combo_M_mat_nsprcomp_annot[Time_Samples,]$PC1)
cor(seurat_combo_E_mat_nsprcomp_annot[ConGood_Samples,]$PC1,seurat_combo_M_mat_nsprcomp_annot[ConGood_Samples,]$PC1)
# -0.2881375
# -0.5316009

plot_data <- as.data.frame(cbind(seurat_combo_E_mat_nsprcomp_annot[Time_Samples,]$PC1,seurat_combo_M_mat_nsprcomp_annot[Time_Samples,]$PC1,seurat_combo_M_mat_nsprcomp_annot[Time_Samples,]$Time))
colnames(plot_data) <- c("Egenes","Mgenes","key")
plot_data$key_label <- factor(plot_data$key)
gg_means(plot_data)

plot_data <- as.data.frame(cbind(seurat_combo_E_mat_nsprcomp_annot[ConGood_Samples,]$PC1,seurat_combo_M_mat_nsprcomp_annot[ConGood_Samples,]$PC1,seurat_combo_M_mat_nsprcomp_annot[ConGood_Samples,]$GBC_pM))
colnames(plot_data) <- c("Egenes","Mgenes","key")
plot_data$key_label <- factor(plot_data$key)
gg_means(plot_data)

### Mann Whitney U-test ###

# nnPCA E Time
nnPCA_Egenes_0 <- seurat_combo_E_mat_nsprcomp_annot[seurat_combo_E_mat_nsprcomp_annot$Time == 0, ]
nnPCA_Egenes_1 <- seurat_combo_E_mat_nsprcomp_annot[seurat_combo_E_mat_nsprcomp_annot$Time == 1, ]
nnPCA_Egenes_2 <- seurat_combo_E_mat_nsprcomp_annot[seurat_combo_E_mat_nsprcomp_annot$Time == 2, ]
nnPCA_Egenes_3 <- seurat_combo_E_mat_nsprcomp_annot[seurat_combo_E_mat_nsprcomp_annot$Time == 3, ]
nnPCA_Egenes_4 <- seurat_combo_E_mat_nsprcomp_annot[seurat_combo_E_mat_nsprcomp_annot$Time == 4, ]
nnPCA_Egenes_8 <- seurat_combo_E_mat_nsprcomp_annot[seurat_combo_E_mat_nsprcomp_annot$Time == 8, ]

dim(nnPCA_Egenes_0)
dim(nnPCA_Egenes_1)
dim(nnPCA_Egenes_2)
dim(nnPCA_Egenes_3)
dim(nnPCA_Egenes_4)
dim(nnPCA_Egenes_8)
table(integration.combined.All@meta.data$Time)

nnPCA_E_0vs1 <- wilcox.test(nnPCA_Egenes_0$PC1,nnPCA_Egenes_1$PC1,alternative='greater')
nnPCA_E_0vs1$statistic/(length(nnPCA_Egenes_0$PC1)*length(nnPCA_Egenes_1$PC1))
nnPCA_E_0vs1$p.value

nnPCA_E_1vs2 <- wilcox.test(nnPCA_Egenes_1$PC1,nnPCA_Egenes_2$PC1,alternative='greater')
nnPCA_E_1vs2$statistic/(length(nnPCA_Egenes_1$PC1)*length(nnPCA_Egenes_2$PC1))
nnPCA_E_1vs2$p.value

nnPCA_E_2vs3 <- wilcox.test(nnPCA_Egenes_2$PC1,nnPCA_Egenes_3$PC1,alternative='greater')
nnPCA_E_2vs3$statistic/(length(nnPCA_Egenes_2$PC1)*length(nnPCA_Egenes_3$PC1))
nnPCA_E_2vs3$p.value

nnPCA_E_3vs4 <- wilcox.test(nnPCA_Egenes_3$PC1,nnPCA_Egenes_4$PC1,alternative='greater')
nnPCA_E_3vs4$statistic/(length(nnPCA_Egenes_3$PC1)*length(nnPCA_Egenes_4$PC1))
nnPCA_E_3vs4$p.value

nnPCA_E_4vs8 <- wilcox.test(nnPCA_Egenes_4$PC1,nnPCA_Egenes_8$PC1,alternative='greater')
nnPCA_E_4vs8$statistic/(length(nnPCA_Egenes_4$PC1)*length(nnPCA_Egenes_8$PC1))
nnPCA_E_4vs8$p.value


# nnPCA M Time
nnPCA_Mgenes_0 <- seurat_combo_M_mat_nsprcomp_annot[seurat_combo_M_mat_nsprcomp_annot$Time == 0, ]
nnPCA_Mgenes_1 <- seurat_combo_M_mat_nsprcomp_annot[seurat_combo_M_mat_nsprcomp_annot$Time == 1, ]
nnPCA_Mgenes_2 <- seurat_combo_M_mat_nsprcomp_annot[seurat_combo_M_mat_nsprcomp_annot$Time == 2, ]
nnPCA_Mgenes_3 <- seurat_combo_M_mat_nsprcomp_annot[seurat_combo_M_mat_nsprcomp_annot$Time == 3, ]
nnPCA_Mgenes_4 <- seurat_combo_M_mat_nsprcomp_annot[seurat_combo_M_mat_nsprcomp_annot$Time == 4, ]
nnPCA_Mgenes_8 <- seurat_combo_M_mat_nsprcomp_annot[seurat_combo_M_mat_nsprcomp_annot$Time == 8, ]

dim(nnPCA_Mgenes_0)
dim(nnPCA_Mgenes_1)
dim(nnPCA_Mgenes_2)
dim(nnPCA_Mgenes_3)
dim(nnPCA_Mgenes_4)
dim(nnPCA_Mgenes_8)
table(integration.combined.All@meta.data$Time)

nnPCA_M_0vs1 <- wilcox.test(nnPCA_Mgenes_0$PC1,nnPCA_Mgenes_1$PC1,alternative='less')
1-nnPCA_M_0vs1$statistic/(length(nnPCA_Mgenes_0$PC1)*length(nnPCA_Mgenes_1$PC1))
nnPCA_M_0vs1$p.value

nnPCA_M_1vs2 <- wilcox.test(nnPCA_Mgenes_1$PC1,nnPCA_Mgenes_2$PC1,alternative='less')
1-nnPCA_M_1vs2$statistic/(length(nnPCA_Mgenes_1$PC1)*length(nnPCA_Mgenes_2$PC1))
nnPCA_M_1vs2$p.value

nnPCA_M_2vs3 <- wilcox.test(nnPCA_Mgenes_2$PC1,nnPCA_Mgenes_3$PC1,alternative='less')
1-nnPCA_M_2vs3$statistic/(length(nnPCA_Mgenes_2$PC1)*length(nnPCA_Mgenes_3$PC1))
nnPCA_M_2vs3$p.value

nnPCA_M_3vs4 <- wilcox.test(nnPCA_Mgenes_3$PC1,nnPCA_Mgenes_4$PC1,alternative='less')
1-nnPCA_M_3vs4$statistic/(length(nnPCA_Mgenes_3$PC1)*length(nnPCA_Mgenes_4$PC1))
nnPCA_M_3vs4$p.value

nnPCA_M_4vs8 <- wilcox.test(nnPCA_Mgenes_4$PC1,nnPCA_Mgenes_8$PC1,alternative='less')
1-nnPCA_M_4vs8$statistic/(length(nnPCA_Mgenes_4$PC1)*length(nnPCA_Mgenes_8$PC1))
nnPCA_M_4vs8$p.value

# nnPCA E GBC 
nnPCA_Egenes_GBC_1 <- seurat_combo_E_mat_nsprcomp_annot[seurat_combo_E_mat_nsprcomp_annot$GBC_values == 1, ]
nnPCA_Egenes_GBC_2 <- seurat_combo_E_mat_nsprcomp_annot[seurat_combo_E_mat_nsprcomp_annot$GBC_values == 2, ]
nnPCA_Egenes_GBC_3 <- seurat_combo_E_mat_nsprcomp_annot[seurat_combo_E_mat_nsprcomp_annot$GBC_values == 3, ]
nnPCA_Egenes_GBC_4 <- seurat_combo_E_mat_nsprcomp_annot[seurat_combo_E_mat_nsprcomp_annot$GBC_values == 4, ]
nnPCA_Egenes_GBC_8 <- seurat_combo_E_mat_nsprcomp_annot[seurat_combo_E_mat_nsprcomp_annot$GBC_values == 8, ]
nnPCA_Egenes_GBC_9 <- seurat_combo_E_mat_nsprcomp_annot[seurat_combo_E_mat_nsprcomp_annot$GBC_values == 9, ]
nnPCA_Egenes_GBC_12 <- seurat_combo_E_mat_nsprcomp_annot[seurat_combo_E_mat_nsprcomp_annot$GBC_values == 12, ]
nnPCA_Egenes_GBC_13 <- seurat_combo_E_mat_nsprcomp_annot[seurat_combo_E_mat_nsprcomp_annot$GBC_values == 13, ]

dim(nnPCA_Egenes_GBC_1)
dim(nnPCA_Egenes_GBC_2)
dim(nnPCA_Egenes_GBC_3)
dim(nnPCA_Egenes_GBC_4)
dim(nnPCA_Egenes_GBC_8)
dim(nnPCA_Egenes_GBC_9)
dim(nnPCA_Egenes_GBC_12)
dim(nnPCA_Egenes_GBC_13)
table(integration.combined.All@meta.data$GBC)

nnPCA_E_GBC_1vs2 <- wilcox.test(nnPCA_Egenes_GBC_1$PC1,nnPCA_Egenes_GBC_2$PC1,alternative="greater")
nnPCA_E_GBC_1vs2$statistic/(length(nnPCA_Egenes_GBC_1$PC1)*length(nnPCA_Egenes_GBC_2$PC1))
nnPCA_E_GBC_1vs2$p.value

nnPCA_E_GBC_2vs3 <- wilcox.test(nnPCA_Egenes_GBC_2$PC1,nnPCA_Egenes_GBC_3$PC1,alternative="greater")
nnPCA_E_GBC_2vs3$statistic/(length(nnPCA_Egenes_GBC_2$PC1)*length(nnPCA_Egenes_GBC_3$PC1))
nnPCA_E_GBC_2vs3$p.value

nnPCA_E_GBC_3vs4 <- wilcox.test(nnPCA_Egenes_GBC_3$PC1,nnPCA_Egenes_GBC_4$PC1,alternative="greater")
nnPCA_E_GBC_3vs4$statistic/(length(nnPCA_Egenes_GBC_3$PC1)*length(nnPCA_Egenes_GBC_4$PC1))
nnPCA_E_GBC_3vs4$p.value

nnPCA_E_GBC_4vs8 <- wilcox.test(nnPCA_Egenes_GBC_4$PC1,nnPCA_Egenes_GBC_8$PC1,alternative="greater")
nnPCA_E_GBC_4vs8$statistic/(length(nnPCA_Egenes_GBC_4$PC1)*length(nnPCA_Egenes_GBC_8$PC1))
nnPCA_E_GBC_4vs8$p.value

nnPCA_E_GBC_8vs9 <- wilcox.test(nnPCA_Egenes_GBC_8$PC1,nnPCA_Egenes_GBC_9$PC1,alternative="greater")
nnPCA_E_GBC_8vs9$statistic/(length(nnPCA_Egenes_GBC_8$PC1)*length(nnPCA_Egenes_GBC_9$PC1))
nnPCA_E_GBC_8vs9$p.value

nnPCA_E_GBC_9vs12 <- wilcox.test(nnPCA_Egenes_GBC_9$PC1,nnPCA_Egenes_GBC_12$PC1,alternative="greater")
nnPCA_E_GBC_9vs12$statistic/(length(nnPCA_Egenes_GBC_9$PC1)*length(nnPCA_Egenes_GBC_12$PC1))
nnPCA_E_GBC_9vs12$p.value

nnPCA_E_GBC_12vs13 <- wilcox.test(nnPCA_Egenes_GBC_12$PC1,nnPCA_Egenes_GBC_13$PC1,alternative="greater")
nnPCA_E_GBC_12vs13$statistic/(length(nnPCA_Egenes_GBC_12$PC1)*length(nnPCA_Egenes_GBC_13$PC1))
nnPCA_E_GBC_12vs13$p.value

# nnPCA M GBC 
nnPCA_Mgenes_GBC_1 <- seurat_combo_M_mat_nsprcomp_annot[seurat_combo_M_mat_nsprcomp_annot$GBC_values == 1, ]
nnPCA_Mgenes_GBC_2 <- seurat_combo_M_mat_nsprcomp_annot[seurat_combo_M_mat_nsprcomp_annot$GBC_values == 2, ]
nnPCA_Mgenes_GBC_3 <- seurat_combo_M_mat_nsprcomp_annot[seurat_combo_M_mat_nsprcomp_annot$GBC_values == 3, ]
nnPCA_Mgenes_GBC_4 <- seurat_combo_M_mat_nsprcomp_annot[seurat_combo_M_mat_nsprcomp_annot$GBC_values == 4, ]
nnPCA_Mgenes_GBC_8 <- seurat_combo_M_mat_nsprcomp_annot[seurat_combo_M_mat_nsprcomp_annot$GBC_values == 8, ]
nnPCA_Mgenes_GBC_9 <- seurat_combo_M_mat_nsprcomp_annot[seurat_combo_M_mat_nsprcomp_annot$GBC_values == 9, ]
nnPCA_Mgenes_GBC_12 <- seurat_combo_M_mat_nsprcomp_annot[seurat_combo_M_mat_nsprcomp_annot$GBC_values == 12, ]
nnPCA_Mgenes_GBC_13 <- seurat_combo_M_mat_nsprcomp_annot[seurat_combo_M_mat_nsprcomp_annot$GBC_values == 13, ]

dim(nnPCA_Mgenes_GBC_1)
dim(nnPCA_Mgenes_GBC_2)
dim(nnPCA_Mgenes_GBC_3)
dim(nnPCA_Mgenes_GBC_4)
dim(nnPCA_Mgenes_GBC_8)
dim(nnPCA_Mgenes_GBC_9)
dim(nnPCA_Mgenes_GBC_12)
dim(nnPCA_Mgenes_GBC_13)
table(integration.combined.All@meta.data$GBC)

nnPCA_M_GBC_1vs2 <- wilcox.test(nnPCA_Mgenes_GBC_1$PC1,nnPCA_Mgenes_GBC_2$PC1,alternative='less')
1-nnPCA_M_GBC_1vs2$statistic/(length(nnPCA_Mgenes_GBC_1$PC1)*length(nnPCA_Mgenes_GBC_2$PC1))
nnPCA_M_GBC_1vs2$p.value

nnPCA_M_GBC_2vs3 <- wilcox.test(nnPCA_Mgenes_GBC_2$PC1,nnPCA_Mgenes_GBC_3$PC1,alternative='less')
1-nnPCA_M_GBC_2vs3$statistic/(length(nnPCA_Mgenes_GBC_2$PC1)*length(nnPCA_Mgenes_GBC_3$PC1))
nnPCA_M_GBC_2vs3$p.value

nnPCA_M_GBC_3vs4 <- wilcox.test(nnPCA_Mgenes_GBC_3$PC1,nnPCA_Mgenes_GBC_4$PC1,alternative='less')
1-nnPCA_M_GBC_3vs4$statistic/(length(nnPCA_Mgenes_GBC_3$PC1)*length(nnPCA_Mgenes_GBC_4$PC1))
nnPCA_M_GBC_3vs4$p.value

nnPCA_M_GBC_4vs8 <- wilcox.test(nnPCA_Mgenes_GBC_4$PC1,nnPCA_Mgenes_GBC_8$PC1,alternative='less')
1-nnPCA_M_GBC_4vs8$statistic/(length(nnPCA_Mgenes_GBC_4$PC1)*length(nnPCA_Mgenes_GBC_8$PC1))
nnPCA_M_GBC_4vs8$p.value

nnPCA_M_GBC_8vs9 <- wilcox.test(nnPCA_Mgenes_GBC_8$PC1,nnPCA_Mgenes_GBC_9$PC1,alternative='less')
1-nnPCA_M_GBC_8vs9$statistic/(length(nnPCA_Mgenes_GBC_8$PC1)*length(nnPCA_Mgenes_GBC_9$PC1))
nnPCA_M_GBC_8vs9$p.value

nnPCA_M_GBC_9vs12 <- wilcox.test(nnPCA_Mgenes_GBC_9$PC1,nnPCA_Mgenes_GBC_12$PC1,alternative='less')
1-nnPCA_M_GBC_9vs12$statistic/(length(nnPCA_Mgenes_GBC_9$PC1)*length(nnPCA_Mgenes_GBC_12$PC1))
nnPCA_M_GBC_9vs12$p.value

nnPCA_M_GBC_12vs13 <- wilcox.test(nnPCA_Mgenes_GBC_12$PC1,nnPCA_Mgenes_GBC_13$PC1,alternative='less')
1-nnPCA_M_GBC_12vs13$statistic/(length(nnPCA_Mgenes_GBC_12$PC1)*length(nnPCA_Mgenes_GBC_13$PC1))
nnPCA_M_GBC_12vs13$p.value

#######################################################
### Correlation of NONEMT Genes with E and M-scores ###
#######################################################

dim(seurat_combo_All_mat)
seurat_combo_All_mat_NOEMT <- seurat_combo_All_mat[-which(row.names(seurat_combo_All_mat) %in% EMT_genes$Gene),]
dim(seurat_combo_All_mat_NOEMT)

seurat_combo_All_mat_NOEMT_TimeT <- t(seurat_combo_All_mat_NOEMT)[Time_Samples,]
seurat_combo_All_mat_NOEMT_ConT <- t(seurat_combo_All_mat_NOEMT)[ConGood_Samples,]
dim(seurat_combo_All_mat_NOEMT_TimeT)
dim(seurat_combo_All_mat_NOEMT_ConT)

correl_TimeE <- cor(seurat_combo_All_mat_NOEMT_TimeT,seurat_combo_E_mat_nsprcomp_annot[Time_Samples,]$PC1)
hist(correl_TimeE)
correl_TimeM <- cor(seurat_combo_All_mat_NOEMT_TimeT,seurat_combo_M_mat_nsprcomp_annot[Time_Samples,]$PC1)
hist(correl_TimeM)

correl_ConE <- cor(seurat_combo_All_mat_NOEMT_ConT,seurat_combo_E_mat_nsprcomp_annot[ConGood_Samples,]$PC1)
hist(correl_ConE)
correl_ConM <- cor(seurat_combo_All_mat_NOEMT_ConT,seurat_combo_M_mat_nsprcomp_annot[ConGood_Samples,]$PC1)
hist(correl_ConM)

########################## Figure 3 ########################## 

# Correlated with M-scores
correl_TimeE_Bot25 <- correl_TimeE[correl_TimeE < -0.25,]
length(correl_TimeE_Bot25)

correl_TimeM_Top25 <- correl_TimeM[correl_TimeM > 0.25,]
length(correl_TimeM_Top25)
length(intersect(names(correl_TimeE_Bot25),names(correl_TimeM_Top25)))

correl_ConE_Bot25 <- correl_ConE[correl_ConE < -0.25,]
length(correl_ConE_Bot25)

correl_ConM_Top25 <- correl_ConM[correl_ConM > 0.25,]
length(correl_ConM_Top25)
length(intersect(names(correl_ConE_Bot25),names(correl_ConM_Top25)))

# Correlated with E-scores
correl_TimeE_Top25 <- correl_TimeE[correl_TimeE > 0.25,]
length(correl_TimeE_Top25)

correl_TimeM_Bot25 <- correl_TimeM[correl_TimeM < -0.25,]
length(correl_TimeM_Bot25)
length(intersect(names(correl_TimeE_Top25),names(correl_TimeM_Bot25)))

correl_ConE_Top25 <- correl_ConE[correl_ConE > 0.25,]
length(correl_ConE_Top25)

correl_ConM_Bot25 <- correl_ConM[correl_ConM < -0.25,]
length(correl_ConM_Bot25)
length(intersect(names(correl_ConE_Top25),names(correl_ConM_Bot25)))

#
Time_25_genes <- intersect(names(correl_TimeE_Bot25),names(correl_TimeM_Top25))
Time_25_dn_genes <- intersect(names(correl_TimeE_Top25),names(correl_TimeM_Bot25))

pdf(file="Figures/Panels/Fig3_TimeCorr.pdf")
plot(correl_TimeE[,1],correl_TimeM[,1],pch=16)
points(correl_TimeE[Time_25_genes,1],correl_TimeM[Time_25_genes,1],pch=16,col='red')
points(correl_TimeE[Time_25_dn_genes,1],correl_TimeM[Time_25_dn_genes,1],pch=16,col='blue')
cor(correl_TimeE[,1],correl_TimeM[,1])
dev.off()

#
Con_25_genes <- intersect(names(correl_ConE_Bot25),names(correl_ConM_Top25))
Con_25_dn_genes <- intersect(names(correl_ConE_Top25),names(correl_ConM_Bot25))

pdf(file="Figures/Panels/Fig3_ConCorr.pdf")
plot(correl_ConE[,1],correl_ConM[,1],pch=16)
points(correl_ConE[Con_25_genes,1],correl_ConM[Con_25_genes,1],pch=16,col='red')
points(correl_ConE[Con_25_dn_genes,1],correl_ConM[Con_25_dn_genes,1],pch=16,col='blue')
cor(correl_ConE[,1],correl_ConM[,1])
dev.off()


# UP 25
Time_genes_up <- intersect(names(correl_TimeE_Bot25),names(correl_TimeM_Top25))
Con_genes_up <- intersect(names(correl_ConE_Bot25),names(correl_ConM_Top25))

length(Time_genes_up) # 18
length(Con_genes_up) # 83

Con_only_up <- setdiff(Con_genes_up,Time_genes_up)
Con_only_up # 67

Time_only_up <- setdiff(Time_genes_up,Con_genes_up)
Time_only_up # 2

# DOWN 25
Time_genes_dn <- intersect(names(correl_TimeE_Top25),names(correl_TimeM_Bot25))
Con_genes_dn <- intersect(names(correl_ConE_Top25),names(correl_ConM_Bot25))

length(Time_genes_dn) # 26
length(Con_genes_dn) # 108

Con_only_dn <- setdiff(Con_genes_dn,Time_genes_dn)
Con_only_dn # 82

Time_only_dn <- setdiff(Time_genes_dn,Con_genes_dn)
Time_only_dn # 0

library(clusterProfiler)
library(GOSemSim)
library(GO.db)
library(HGNChelper)

hsGO <- godata('org.Hs.eg.db', ont="MF")
hsGO_mf <- godata('org.Hs.eg.db', ont="MF")
hsGO_bp <- godata('org.Hs.eg.db', ont="BP")

# Universe is all non EMT genes retained post integration
all_genes <- setdiff(row.names(integration.combined.All@assays$integrated),EMT_genes$Gene)

Con_only_up_bp <- enrichGO(gene = Con_only_up,
                           universe      = all_genes,
                           OrgDb         = org.Hs.eg.db,
                           keyType = "SYMBOL",
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 0.05,
                           readable      = FALSE)
Con_only_up_bp
# none

Con_only_dn_bp <- enrichGO(gene = Con_only_dn,
                           universe      = all_genes,
                           OrgDb         = org.Hs.eg.db,
                           keyType = "SYMBOL",
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 0.05,
                           readable      = FALSE)
Con_only_dn_bp
Con_only_dn_bp@result[1:5,] # of 34


Con_only_up_mf <- enrichGO(gene = Con_only_up,
                           universe      = all_genes,
                           OrgDb         = org.Hs.eg.db,
                           keyType = "SYMBOL",
                           ont           = "MF",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 0.05,
                           readable      = FALSE)
Con_only_up_mf
Con_only_up_mf@result[1:2,]

Con_only_dn_mf <- enrichGO(gene = Con_only_dn,
                           universe      = all_genes,
                           OrgDb         = org.Hs.eg.db,
                           keyType = "SYMBOL",
                           ont           = "MF",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 0.05,
                           readable      = FALSE)
Con_only_dn_mf
Con_only_dn_mf@result[1:5,] # of 9

Con_only_up_cc <- enrichGO(gene = Con_only_up,
                           universe      = all_genes,
                           OrgDb         = org.Hs.eg.db,
                           keyType = "SYMBOL",
                           ont           = "CC",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 0.05,
                           readable      = FALSE)
Con_only_up_cc
Con_only_up_cc@result[1:1,]

Con_only_dn_cc <- enrichGO(gene = Con_only_dn,
                           universe      = all_genes,
                           OrgDb         = org.Hs.eg.db,
                           keyType = "SYMBOL",
                           ont           = "CC",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 0.05,
                           readable      = FALSE)
Con_only_dn_cc
Con_only_dn_cc@result[1:5,] # of 10

Combined_GO_Results <- rbind( Con_only_dn_bp@result[1:5,1:7],
                              Con_only_up_mf@result[1:2,1:7],
                              Con_only_dn_mf@result[1:5,1:7],
                              Con_only_up_cc@result[1:1,1:7],
                              Con_only_dn_cc@result[1:5,1:7])

Combined_GO_Results$Source <- "Con_dn"
Combined_GO_Results$Source[c(6:7,13)] <- "Con_up"

OddRatios <- matrix(, nrow = dim(Combined_GO_Results)[1], ncol = 4)

library(stringr)
for(i in 1:dim(Combined_GO_Results)[1]){
  current_geneRatio <- Combined_GO_Results$GeneRatio[i]
  current_bgRatio <- Combined_GO_Results$BgRatio[i]
  
  geneRatio_list <- str_split(current_geneRatio, "/", n = Inf, simplify = FALSE)
  bgRation_list <- str_split(current_bgRatio, "/", n = Inf, simplify = FALSE)
  
  numeric_geneRatio <- as.numeric(geneRatio_list[[1]][1]) /as.numeric(geneRatio_list[[1]][2] )
  numeric_bgRatio <- as.numeric(bgRation_list[[1]][1]) /as.numeric(bgRation_list[[1]][2] )
  
  OddRatio <- numeric_geneRatio/numeric_bgRatio
  
  OddRatios[i,1] <- numeric_geneRatio
  OddRatios[i,2] <- numeric_bgRatio
  OddRatios[i,3] <- OddRatio
  OddRatios[i,4] <- log(OddRatio)
  
}
Combined_GO_Results$LogOdds <- OddRatios[,4]

attach(Combined_GO_Results)
Combined_GO_Results_Sort <- Combined_GO_Results[order(Source,LogOdds),]
detach(Combined_GO_Results)

Combined_GO_Results_Sort$Description <- factor(Combined_GO_Results_Sort$Description,levels=Combined_GO_Results_Sort$Description[order(Combined_GO_Results_Sort$Source,Combined_GO_Results_Sort$LogOdds)])

p <- ggplot(Combined_GO_Results_Sort, aes(x = Description, y = LogOdds))+
  geom_col(aes(fill = Source), width = 0.7)

pdf(file="Figures/Panels/Fig3_GO_Enrichment_Bar.pdf")
p + coord_flip()
dev.off()

# High M in only Con #
Con_only_up_AllMCorr <- as.data.frame(cbind(correl_ConM[Con_only_up,],correl_TimeM[Con_only_up,],correl_ConM[Con_only_up,]-correl_TimeM[Con_only_up,]))

##########################
### Overlapping Scores ###
##########################

overlapping_Escores <- seurat_combo_E_mat_nsprcomp_annot
overlapping_Mscores <- seurat_combo_M_mat_nsprcomp_annot

BatchA <- row.names(overlapping_Escores[overlapping_Escores$Batch == "A",])
BatchB <- row.names(overlapping_Escores[overlapping_Escores$Batch == "B",])
BatchC <- row.names(overlapping_Escores[overlapping_Escores$Batch == "C",])

### EM Contour ###
EM_contour_df <- as.data.frame(cbind(overlapping_Escores$PC1,overlapping_Mscores$PC1,overlapping_Escores$Batch,overlapping_Mscores$Batch))
names(EM_contour_df) <- c("Escores","Mscores","Batch","Expr")
EM_contour_df[EM_contour_df$Expr=='A',]$Expr <- "Time"
EM_contour_df[EM_contour_df$Expr=='B',]$Expr <- "Time"
EM_contour_df[EM_contour_df$Expr=='C',]$Expr <- "Dose"

EM_contour_df$Escores <- as.numeric(as.character(EM_contour_df$Escores))
EM_contour_df$Mscores <- as.numeric(as.character(EM_contour_df$Mscores))

ggplot(EM_contour_df, aes(x = Escores, y = Mscores, colour=Expr)) +
  geom_density_2d(bins=10,alpha=0.5,size=1) + theme_bw()

### Bubble Plot Sample Sets ###
Con_Samples <- ConGood_Samples
ConAll_samples <- rownames(integration.combined.All@meta.data[integration.combined.All@meta.data$Batch == "C",])
length(Time_Samples)
length(Con_Samples)
length(ConAll_samples)

Labeled_Samples <- c(Time_Samples,Con_Samples)
length(Labeled_Samples)

### Bubble Plot EMscore ###

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

EMCount_Heatmap <- function(bonds_df){
  plot(overlapping_Escores$PC1, overlapping_Mscores$PC1,pch=16,cex=0.5)
  
  odds_mat <- matrix(, nrow = dim(bonds_df)[1]+1, ncol = dim(bonds_df)[1]+1)
  pv_mat <- matrix(, nrow = dim(bonds_df)[1]+1, ncol = dim(bonds_df)[1]+1)
  counts_mat <- matrix(, nrow = dim(bonds_df)[1]+1, ncol = dim(bonds_df)[1]+1)
  
  for(i in 1:(dim(bonds_df)[1]+1)){
    for(j in 1:(dim(bonds_df)[1]+1)){
      E_boundary_samples <- GetBoundarySamples(i,1,bonds_df)
      M_boundary_sampels <- GetBoundarySamples(j,2,bonds_df)
      intersect_samples = intersect(E_boundary_samples, M_boundary_sampels)
      
      # Basic Counts
      print(i)
      print(j)
      print(length(intersect_samples))
      
      # Con and Time Counts
      cluster_time_samples <- length(intersect(intersect_samples,Time_Samples))
      cluster_con_samples <- length(intersect(intersect_samples,Con_Samples))
      cluster_allcon_samples <- length(intersect(intersect_samples,ConAll_samples))
      
      outcluster_time_samples <- length(Time_Samples) - cluster_time_samples
      outcluster_con_samples <- length(Con_Samples) - cluster_con_samples
      outcluster_allcon_smaples <- length(ConAll_samples) - cluster_allcon_samples
      
      #print(cluster_time_samples)
      #print(cluster_con_samples)
      #print(cluster_allcon_samples)
      #print(outcluster_time_samples)
      #print(outcluster_con_samples)
      #print(outcluster_allcon_smaples)
      
      # Fisher Exact
      Time_vs_Con <- matrix(c(cluster_time_samples, outcluster_time_samples, cluster_con_samples, outcluster_con_samples), nrow = 2,
                            dimnames =
                              list(c("In Cluster", "Out Cluster"),
                                   c("Time ", "Con")))
      test <- fisher.test(Time_vs_Con)
      print(Time_vs_Con)
      print(test$p.value)
      print(test$estimate)
      
      odds_mat[j,i] = test$estimate
      pv_mat[j,i] = test$p.value
      counts_mat[j,i] = cluster_time_samples + cluster_con_samples
      
    }
  }
  
  return(list(odds_mat,pv_mat,counts_mat))
}

#results <- EMCount_Heatmap(Bounds_IQR)
results <- EMCount_Heatmap(Bounds_IQR_med)
odd_mat = results[[1]]
pv_mat = results[[2]]
counts_mat = results[[3]]

#log2(odd_mat[nrow(odd_mat):1, ])
#pv_mat[nrow(pv_mat):1, ]
#counts_mat[nrow(pv_mat):1, ]

# 4x4 --- Bubble Plot
odd_mat_long <- matrix(, nrow = 16, ncol = 4)
for (i in 1:4){
  for (j in 1:4){
    row = 4*(i-1) + j
    odd_mat_long[row,1] = j
    odd_mat_long[row,2] = i
    odd_mat_long[row,3] = log2(odd_mat[i,j])
    odd_mat_long[row,4] = counts_mat[i,j]
  }
}
odd_mat_long <- as.data.frame(odd_mat_long)
colnames(odd_mat_long) <- c("X","Y","Odds","Count")

### Figure 4 ###
#pdf(file="Figures/Panels/Fig4_OddsBubble_v2.pdf")
bubble_p <- ggplot(odd_mat_long, aes(x = X, y = Y,size=Count,fill=Odds)) + 
  geom_point(alpha=1.0, shape=21, color="black") +
  scale_size(range = c(5, 15)) +
  scale_fill_gradient2(low='blue',high='red',mid='white',midpoint=0,limits=c(-3,3),oob=squish) +
  theme_bw() + coord_fixed(ratio = 1) + ylim(c(0.75,4.25)) + xlim(c(0.75,4.25))
bubble_p
#dev.off()

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

# Use unintegrated data  for Marker ID
integration.combined.All.SpoofQuart <- integration.combined.All
DefaultAssay(integration.combined.All.SpoofQuart) <- "RNA"

# Drop Unannoated Con Samples
integration.combined.All.SpoofQuart <- subset(integration.combined.All.SpoofQuart, subset = GBC != 'GBC0' | Batch == 'A' | Batch == 'B')
integration.combined.All.SpoofQuart

integration.combined.All.SpoofQuart$quart_clusters <- 0

integration.combined.All.SpoofQuart@meta.data[quart_clusters[[1]],]$quart_clusters <- 1
integration.combined.All.SpoofQuart@meta.data[quart_clusters[[2]],]$quart_clusters <- 2
integration.combined.All.SpoofQuart@meta.data[quart_clusters[[3]],]$quart_clusters <- 3
integration.combined.All.SpoofQuart@meta.data[quart_clusters[[4]],]$quart_clusters <- 4
integration.combined.All.SpoofQuart@meta.data[quart_clusters[[5]],]$quart_clusters <- 5
integration.combined.All.SpoofQuart@meta.data[quart_clusters[[6]],]$quart_clusters <- 6
integration.combined.All.SpoofQuart@meta.data[quart_clusters[[7]],]$quart_clusters <- 7
integration.combined.All.SpoofQuart@meta.data[quart_clusters[[8]],]$quart_clusters <- 8
integration.combined.All.SpoofQuart@meta.data[quart_clusters[[9]],]$quart_clusters <- 9
integration.combined.All.SpoofQuart@meta.data[quart_clusters[[10]],]$quart_clusters <- 10
integration.combined.All.SpoofQuart@meta.data[quart_clusters[[11]],]$quart_clusters <- 11
integration.combined.All.SpoofQuart@meta.data[quart_clusters[[12]],]$quart_clusters <- 12
integration.combined.All.SpoofQuart@meta.data[quart_clusters[[13]],]$quart_clusters <- 13
integration.combined.All.SpoofQuart@meta.data[quart_clusters[[14]],]$quart_clusters <- 14
integration.combined.All.SpoofQuart@meta.data[quart_clusters[[15]],]$quart_clusters <- 15
integration.combined.All.SpoofQuart@meta.data[quart_clusters[[16]],]$quart_clusters <- 16

integration.combined.All.SpoofQuart@meta.data$quart_clusters <- as.factor(integration.combined.All.SpoofQuart@meta.data$quart_clusters)

table(integration.combined.All.SpoofQuart@meta.data$quart_clusters)

# Remove borderline samples (6 of 22765)
integration.combined.All.SpoofQuart <- subset(integration.combined.All.SpoofQuart, subset = quart_clusters != 0)
integration.combined.All.SpoofQuart

table(integration.combined.All.SpoofQuart@meta.data$quart_clusters)

# Check Assignments
Escores <- as.data.frame(seurat_combo_E_mat_nsprcomp$x[union(Time_Samples,Con_Samples),1])
Mscores <- as.data.frame(seurat_combo_M_mat_nsprcomp$x[union(Time_Samples,Con_Samples),1])
MetaData <- integration.combined.All.SpoofQuart@meta.data

MergedScores <- base::merge(Escores, Mscores, by="row.names")
row.names(MergedScores) <- MergedScores$Row.names
MergedScores <- MergedScores[,2:3]
colnames(MergedScores) <- c("Escores","Mscores")
MergedData<- base::merge(MergedScores, MetaData, by="row.names")
  
ggplot(data = MergedData, aes(x=Escores, y=Mscores, color=quart_clusters)) + geom_point()

# Change active identy for marker finding
integration.combined.All.SpoofQuart <- SetIdent(integration.combined.All.SpoofQuart, value = integration.combined.All.SpoofQuart@meta.data$quart_clusters)
table(integration.combined.All.SpoofQuart@active.ident)

# Find Markers
AllClusterMarkers_QuartClusters <- FindAllMarkers(integration.combined.All.SpoofQuart)
saveRDS(AllClusterMarkers_QuartClusters,"Integration_AllClusterMarkers_QuartAlt.rds")

# Subset for GO analysis
AllClusterMarkers_QuartClusters$rawFC <- 2^AllClusterMarkers_QuartClusters$avg_log2FC

RawFC_1stquart <- quantile(AllClusterMarkers_QuartClusters$rawFC,0.25)
RawFC_3rdquart <- quantile(AllClusterMarkers_QuartClusters$rawFC,0.75)

AllClusterMarkers_QuartClusters_FC3rd <- AllClusterMarkers_QuartClusters[AllClusterMarkers_QuartClusters$rawFC > RawFC_3rdquart,]
AllClusterMarkers_QuartClusters_FC3rd_adjPV05 <- AllClusterMarkers_QuartClusters_FC3rd[AllClusterMarkers_QuartClusters_FC3rd$p_val_adj < 0.05,]

AllClusterMarkers_QuartClusters_FC1st <- AllClusterMarkers_QuartClusters[AllClusterMarkers_QuartClusters$rawFC < RawFC_1stquart,]
AllClusterMarkers_QuartClusters_FC1st_adjPV05 <- AllClusterMarkers_QuartClusters_FC1st[AllClusterMarkers_QuartClusters_FC1st$p_val_adj < 0.05,]

UpGenes <- table(AllClusterMarkers_QuartClusters_FC3rd_adjPV05$cluster)
DnGenes <- table(AllClusterMarkers_QuartClusters_FC1st_adjPV05$cluster)

UpGenes_df <- as.data.frame(UpGenes)
UpGenes_df$Var1 <- as.numeric(as.character(UpGenes_df$Var1))

DnGenes_df <- as.data.frame(DnGenes)
DnGenes_df$Var1 <- as.numeric(as.character(DnGenes_df$Var1))

# 4x4 --- Bubble Plot
odd_mat_long <- matrix(, nrow = 16, ncol =9)
index = 1 
for (i in 1:4){
  for (j in 1:4){
    odd_mat_long[index,1] = i
    odd_mat_long[index,2] = j
    odd_mat_long[index,3] <- length(quart_clusters[[index]])
    odd_mat_long[index,4] = UpGenes_df[UpGenes_df$Var1==index,]$Freq
    odd_mat_long[index,5] = DnGenes_df[DnGenes_df$Var1==index,]$Freq
    
    VIM_FC <- AllClusterMarkers_QuartClusters[AllClusterMarkers_QuartClusters$gene=="VIM",]
    if (index %in% VIM_FC$cluster){
      odd_mat_long[index,6] <- VIM_FC[VIM_FC$cluster==index,]$avg_log2FC
    }
    else{
      odd_mat_long[index,6] <- NA
    }
    
    FN1_FC <- AllClusterMarkers_QuartClusters[AllClusterMarkers_QuartClusters$gene=="FN1",]
    if (index %in% FN1_FC$cluster){
      odd_mat_long[index,7] <- FN1_FC[FN1_FC$cluster==index,]$avg_log2FC
    }
    else{
      odd_mat_long[index,7] <- NA
    }
    
    CDH1_FC <- AllClusterMarkers_QuartClusters[AllClusterMarkers_QuartClusters$gene=="CDH1",]
    if (index %in% CDH1_FC$cluster){
      odd_mat_long[index,8] <- CDH1_FC[CDH1_FC$cluster==index,]$avg_log2FC
    }
    else{
      odd_mat_long[index,8] <- NA
    }
    
    EPCAM_FC <- AllClusterMarkers_QuartClusters[AllClusterMarkers_QuartClusters$gene=="EPCAM",]
    if (index %in% EPCAM_FC$cluster){
      odd_mat_long[index,9] <- EPCAM_FC[EPCAM_FC$cluster==index,]$avg_log2FC
    }
    else{
      odd_mat_long[index,9] <- NA
    }
    index = index + 1
  }
}
odd_mat_long <- as.data.frame(odd_mat_long)
colnames(odd_mat_long) <- c("X","Y","Count","Up","Dn","VIM","FN1","CDH1",'EPCAM')

### Figure 4 ###
# Don't aut-save figures here because of plot dimension affects 

#pdf(file="Figures/Panels/Fig4_VIMBubble_v2.pdf")
ggplot(odd_mat_long, aes(x = X, y = Y,size=Count,fill=VIM)) + 
  geom_point(alpha=1.0, shape=21, color="black") +
  scale_size(range = c(5, 15)) +
  scale_fill_gradient2(low='blue',high='red',mid='white',midpoint=0,limits=c(-2,1),oob=squish) +
  theme_bw() + coord_fixed(ratio = 1) + ylim(c(0.75,4.25)) + xlim(c(0.75,4.25))
#dev.off()

#pdf(file="Figures/Panels/Fig4_FN1Bubble_v2.pdf")
ggplot(odd_mat_long, aes(x = X, y = Y,size=Count,fill=FN1)) + 
  geom_point(alpha=1.0, shape=21, color="black") +
  scale_size(range = c(5, 15)) +
  scale_fill_gradient2(low='blue',high='red',mid='white',midpoint=0,limits=c(-3,3),oob=squish) +
  theme_bw() + coord_fixed(ratio = 1) + ylim(c(0.75,4.25)) + xlim(c(0.75,4.25))
#dev.off()

#pdf(file="Figures/Panels/Fig4_CDH1Bubble_v2.pdf")
ggplot(odd_mat_long, aes(x = X, y = Y,size=Count,fill=CDH1)) + 
  geom_point(alpha=1.0, shape=21, color="black") +
  scale_size(range = c(5, 15)) +
  scale_fill_gradient2(low='blue',high='red',mid='white',midpoint=0,limits=c(0,0.5),oob=squish) +
  theme_bw() + coord_fixed(ratio = 1) + ylim(c(0.75,4.25)) + xlim(c(0.75,4.25))
#dev.off()

ggplot(odd_mat_long, aes(x = X, y = Y,size=Count,fill=EPCAM)) + 
  geom_point(alpha=1.0, shape=21, color="black") +
  scale_size(range = c(5, 15)) +
  scale_fill_gradient2(low='blue',high='red',mid='white',midpoint=0,limits=c(-0.5,0.5),oob=squish) +
  theme_bw() + coord_fixed(ratio = 1) + ylim(c(0.75,4.25)) + xlim(c(0.75,4.25))

#pdf(file="Figures/Panels/Fig4_UpRegulatedGenes_v2.pdf")
ggplot(odd_mat_long, aes(x = X, y = Y,size=Count,fill=Up)) + 
  geom_point(alpha=1.0, shape=21, color="black") +
  scale_size(range = c(5, 15)) +
  scale_fill_gradient2(low='black',high='red',mid='yellow',midpoint=80,limits=c(0,250),oob=squish) +
  theme_bw() + coord_fixed(ratio = 1) + ylim(c(0.75,4.25)) + xlim(c(0.75,4.25))
#dev.off()

### Figure 4 ###

##################################
### Time/Concentration Markers ###
##################################

### Quarts ###
integration.combined.All.SpoofQuart.Time <- subset(integration.combined.All.SpoofQuart, subset = Batch %in% c("A","B"))
integration.combined.All.SpoofQuart.Con <- subset(integration.combined.All.SpoofQuart, subset = Batch %in% c("C"))

table(integration.combined.All.SpoofQuart.Time@active.ident)
AllClusterMarkers_TimeOnly_Quart <- FindAllMarkers(integration.combined.All.SpoofQuart.Time)
saveRDS(AllClusterMarkers_TimeOnly_Quart,"AllClusterMarkers_Integration_Quart_TimeOnly_alt.rds")

table(integration.combined.All.SpoofQuart.Con@active.ident)
AllClusterMarkers_ConOnly_Quart <- FindAllMarkers(integration.combined.All.SpoofQuart.Con)
saveRDS(AllClusterMarkers_ConOnly_Quart,"AllClusterMarkers_Integration_Quart_ConOnly_alt.rds")

##########################
### GO Annotation in R ###
##########################

#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#
#BiocManager::install("clusterProfiler")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("GO.db")
#BiocManager::install("GOSemSim")
#BiocManager::install("org.Hs.eg.db")

library(clusterProfiler)
library(GOSemSim)
library(GO.db)
library(HGNChelper)

# Check input format
data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]

hsGO <- godata('org.Hs.eg.db', ont="MF")
hsGO_mf <- godata('org.Hs.eg.db', ont="MF")
hsGO_bp <- godata('org.Hs.eg.db', ont="BP")

integration.combined.All.SpoofQuart.mat <- GetAssayData(object = integration.combined.All.SpoofQuart, slot = "data")
all_genes <- row.names(integration.combined.All.SpoofQuart.mat)

# GO Enrichment
cluster1_up_genes <- AllClusterMarkers_QuartClusters_FC3rd_adjPV05[AllClusterMarkers_QuartClusters_FC3rd_adjPV05$cluster==1,]$gene
ego_cluster1 <- enrichGO(gene = cluster1_up_genes,
                         universe      = all_genes,
                         OrgDb         = org.Hs.eg.db,
                         keyType = "SYMBOL",
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05,
                         readable      = FALSE)

# Loop through all clusters
GO_UP_MF_all <- list()
GO_DN_MF_all <- list()
GO_UP_BP_all <- list()
GO_DN_BP_all <- list()

for(i in 1:16) {                                    
  cluster_up_genes <- AllClusterMarkers_QuartClusters_FC3rd_adjPV05[AllClusterMarkers_QuartClusters_FC3rd_adjPV05$cluster==i,]$gene
  cluster_dn_genes <- AllClusterMarkers_QuartClusters_FC1st_adjPV05[AllClusterMarkers_QuartClusters_FC1st_adjPV05$cluster==i,]$gene
  print(i)
  
  ego_cluster_up_mf <- enrichGO(gene = cluster_up_genes,
                                universe      = all_genes,
                                OrgDb         = org.Hs.eg.db,
                                keyType = "SYMBOL",
                                ont           = "MF",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.05,
                                qvalueCutoff  = 0.05,
                                readable      = FALSE)
  
  ego_cluster_dn_mf <- enrichGO(gene = cluster_dn_genes,
                                universe      = all_genes,
                                OrgDb         = org.Hs.eg.db,
                                keyType = "SYMBOL",
                                ont           = "MF",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.05,
                                qvalueCutoff  = 0.05,
                                readable      = FALSE)
  
  ego_cluster_up_bp <- enrichGO(gene = cluster_up_genes,
                                universe      = all_genes,
                                OrgDb         = org.Hs.eg.db,
                                keyType = "SYMBOL",
                                ont           = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.05,
                                qvalueCutoff  = 0.05,
                                readable      = FALSE)
  
  ego_cluster_dn_bp <- enrichGO(gene = cluster_dn_genes,
                                universe      = all_genes,
                                OrgDb         = org.Hs.eg.db,
                                keyType = "SYMBOL",
                                ont           = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.05,
                                qvalueCutoff  = 0.05,
                                readable      = FALSE)
  
  if(length(cluster_up_genes) < 1){
    ego_cluster_up_mf <- "NA"
    ego_cluster_up_bp <- "NA"
  }
  if(length(cluster_dn_genes) < 1){
    ego_cluster_dn_mf <- "NA"
    ego_cluster_dn_bp <- "NA"
  }
  
  GO_UP_MF_all[[i]] <- ego_cluster_up_mf
  GO_DN_MF_all[[i]] <- ego_cluster_dn_mf
  GO_UP_BP_all[[i]] <- ego_cluster_up_bp
  GO_DN_BP_all[[i]] <- ego_cluster_dn_bp
}

saveRDS(GO_UP_MF_all,"GO_UP_MF_all_Integrated_setseed5alt.rds")
saveRDS(GO_DN_MF_all,"GO_DN_MF_all_Integrated_setseed5alt.rds")
saveRDS(GO_UP_BP_all,"GO_UP_BP_all_Integrated_setseed5alt.rds")
saveRDS(GO_DN_BP_all,"GO_DN_BP_all_Integrated_setseed5alt.rds")

### Molecular Function Up ###
molecularfunction_Up <- matrix(nrow=10000,ncol=2)
n = 1
for (i in c(1:8,10:16)){
  # No sig MF terms for 9
  GO_Result <- GO_UP_MF_all[[i]]@result
  GO_Result_padj95 <- GO_Result[GO_Result$p.adjust < 0.05,]
  GO_Description <- GO_Result_padj95$Description
  
  for (j in 1:length(GO_Description)){
    molecularfunction_Up[n,1] <- GO_Description[j]
    molecularfunction_Up[n,2] <- i
    n = n + 1
  }
}
molecularfunction_Up <- molecularfunction_Up[complete.cases(molecularfunction_Up), ]
molecularfunction_Up_table_df <- as.data.frame(table(molecularfunction_Up[,1]))
dim(molecularfunction_Up_table_df)

length(unique(molecularfunction_Up_table_df$Var1))

#MF Upregualted Time
molecularfunction_Up_Time <- matrix(nrow=4000,ncol=2)
n=1
for (i in c(1,2,3,5,14)){
  GO_Result <- GO_UP_MF_all[[i]]@result
  GO_Result_padj95 <- GO_Result[GO_Result$p.adjust < 0.05,]
  GO_Description <- GO_Result_padj95$Description
  
  for (j in 1:length(GO_Description)){
    molecularfunction_Up_Time[n,1] <- GO_Description[j]
    molecularfunction_Up_Time[n,2] <- i
    n = n + 1
  }
}
molecularfunction_Up_Time <- molecularfunction_Up_Time[complete.cases(molecularfunction_Up_Time), ]
molecularfunction_Up_Time_table_df <- as.data.frame(table(molecularfunction_Up_Time[,1]))
dim(molecularfunction_Up_Time_table_df)

#MF Upregualted Con
molecularfunction_Up_Con <- matrix(nrow=4000,ncol=2)
n=1
for (i in c(4,6,7,8,11)){
  GO_Result <- GO_UP_MF_all[[i]]@result
  GO_Result_padj95 <- GO_Result[GO_Result$p.adjust < 0.05,]
  GO_Description <- GO_Result_padj95$Description
  
  for (j in 1:length(GO_Description)){
    molecularfunction_Up_Con[n,1] <- GO_Description[j]
    molecularfunction_Up_Con[n,2] <- i
    n = n + 1
  }
}
molecularfunction_Up_Con <- molecularfunction_Up_Con[complete.cases(molecularfunction_Up_Con), ]
molecularfunction_Up_Con_table_df <- as.data.frame(table(molecularfunction_Up_Con[,1]))
dim(molecularfunction_Up_Con_table_df)

# Combine count between Time and Dose
Combined_GOTerms_MF_UP <- c(as.character(molecularfunction_Up_Time_table_df$Var1),as.character(molecularfunction_Up_Con_table_df$Var1))

Combined_GOTerms_MF_UP_mat <- matrix(nrow=length(Combined_GOTerms_MF_UP),ncol=3)
n = 1
for (term in Combined_GOTerms_MF_UP){
  
  Combined_GOTerms_MF_UP_mat[n,1] = as.character(term)
  
  if (term %in% molecularfunction_Up_Time_table_df$Var1){
    Combined_GOTerms_MF_UP_mat[n,2] <- as.numeric(molecularfunction_Up_Time_table_df[molecularfunction_Up_Time_table_df$Var1 == term,2])
  }
  else{
    Combined_GOTerms_MF_UP_mat[n,2] <- 0 
  }
  
  if (term %in% molecularfunction_Up_Con_table_df$Var1){
    Combined_GOTerms_MF_UP_mat[n,3] <- as.numeric(molecularfunction_Up_Con_table_df[molecularfunction_Up_Con_table_df$Var1 == term,2])
  }
  else{
    Combined_GOTerms_MF_UP_mat[n,3] <- 0 
  }
  
  n = n + 1
}
Combined_GOTerms_MF_UP_mat_df <- as.data.frame(Combined_GOTerms_MF_UP_mat)
Combined_GOTerms_MF_UP_mat_df$V1 <- as.character(Combined_GOTerms_MF_UP_mat_df$V1)
Combined_GOTerms_MF_UP_mat_df$V2 <- as.numeric(Combined_GOTerms_MF_UP_mat_df$V2)
Combined_GOTerms_MF_UP_mat_df$V3 <- as.numeric(Combined_GOTerms_MF_UP_mat_df$V3)

colnames(Combined_GOTerms_MF_UP_mat_df) <- c("Name","Time","Con")
Combined_GOTerms_MF_UP_mat_df$Ratio <- Combined_GOTerms_MF_UP_mat_df$Time/Combined_GOTerms_MF_UP_mat_df$Con

Combined_GOTerms_MF_UP_mat_df_TimeTerms <- Combined_GOTerms_MF_UP_mat_df[Combined_GOTerms_MF_UP_mat_df$Ratio >= 2.0,]
Combined_GOTerms_MF_UP_mat_df_TimeTerms <- Combined_GOTerms_MF_UP_mat_df_TimeTerms[Combined_GOTerms_MF_UP_mat_df_TimeTerms$Time >= 3,]
dim(Combined_GOTerms_MF_UP_mat_df_TimeTerms)
Combined_GOTerms_MF_UP_mat_df_ConTerms <- Combined_GOTerms_MF_UP_mat_df[Combined_GOTerms_MF_UP_mat_df$Ratio <= 0.5,]
Combined_GOTerms_MF_UP_mat_df_ConTerms <- Combined_GOTerms_MF_UP_mat_df_ConTerms[Combined_GOTerms_MF_UP_mat_df_ConTerms$Con >= 3,]
dim(Combined_GOTerms_MF_UP_mat_df_ConTerms)

### Biological Process Up ###
biologicalprocess_Up <- matrix(nrow=20000,ncol=2)
n = 1
for (i in c(1:4,6:16)){
  # No BP terms for 5
  GO_Result <- GO_UP_BP_all[[i]]@result
  GO_Result_padj95 <- GO_Result[GO_Result$p.adjust < 0.05,]
  GO_Description <- GO_Result_padj95$Description
  
  for (j in 1:length(GO_Description)){
    biologicalprocess_Up[n,1] <- GO_Description[j]
    biologicalprocess_Up[n,2] <- i
    n = n + 1
  }
}
biologicalprocess_Up <- biologicalprocess_Up[complete.cases(biologicalprocess_Up), ]
biologicalprocess_Up_table_df <- as.data.frame(table(biologicalprocess_Up[,1]))

dim(biologicalprocess_Up_table_df)
length(unique(biologicalprocess_Up_table_df$Var1))

#MF Upregualted Time
biologicalprocess_Up_Time <- matrix(nrow=4000,ncol=2)
n=1
for (i in c(1,2,3,9,14)){
  GO_Result <- GO_UP_BP_all[[i]]@result
  GO_Result_padj95 <- GO_Result[GO_Result$p.adjust < 0.05,]
  GO_Description <- GO_Result_padj95$Description
  
  for (j in 1:length(GO_Description)){
    biologicalprocess_Up_Time[n,1] <- GO_Description[j]
    biologicalprocess_Up_Time[n,2] <- i
    n = n + 1
  }
}
biologicalprocess_Up_Time <- biologicalprocess_Up_Time[complete.cases(biologicalprocess_Up_Time), ]
biologicalprocess_Up_Time_table_df <- as.data.frame(table(biologicalprocess_Up_Time[,1]))
dim(biologicalprocess_Up_Time_table_df)

#MF Upregualted Con
biologicalprocess_Up_Con <- matrix(nrow=4000,ncol=2)
n=1
for (i in c(4,6,7,8,11)){
  GO_Result <- GO_UP_BP_all[[i]]@result
  GO_Result_padj95 <- GO_Result[GO_Result$p.adjust < 0.05,]
  GO_Description <- GO_Result_padj95$Description
  
  for (j in 1:length(GO_Description)){
    biologicalprocess_Up_Con[n,1] <- GO_Description[j]
    biologicalprocess_Up_Con[n,2] <- i
    n = n + 1
  }
}
biologicalprocess_Up_Con <- biologicalprocess_Up_Con[complete.cases(biologicalprocess_Up_Con), ]
biologicalprocess_Up_Con_table_df <- as.data.frame(table(biologicalprocess_Up_Con[,1]))
dim(biologicalprocess_Up_Con)


# Combine count between Time and Dose
Combined_GOTerms_BP_UP <- union(as.character(biologicalprocess_Up_Time_table_df$Var1),as.character(biologicalprocess_Up_Con_table_df$Var1))

Combined_GOTerms_BP_UP_mat <- matrix(nrow=length(Combined_GOTerms_BP_UP),ncol=3)
n = 1
for (term in Combined_GOTerms_BP_UP){
  
  Combined_GOTerms_BP_UP_mat[n,1] = as.character(term)
  
  if (term %in% biologicalprocess_Up_Time_table_df$Var1){
    Combined_GOTerms_BP_UP_mat[n,2] <- as.numeric(biologicalprocess_Up_Time_table_df[biologicalprocess_Up_Time_table_df$Var1 == term,2])
  }
  else{
    Combined_GOTerms_BP_UP_mat[n,2] <- 0 
  }
  
  if (term %in% biologicalprocess_Up_Con_table_df$Var1){
    Combined_GOTerms_BP_UP_mat[n,3] <- as.numeric(biologicalprocess_Up_Con_table_df[biologicalprocess_Up_Con_table_df$Var1 == term,2])
  }
  else{
    Combined_GOTerms_BP_UP_mat[n,3] <- 0 
  }
  
  n = n + 1
}
Combined_GOTerms_BP_UP_mat_df <- as.data.frame(Combined_GOTerms_BP_UP_mat)
Combined_GOTerms_BP_UP_mat_df$V1 <- as.character(Combined_GOTerms_BP_UP_mat_df$V1)
Combined_GOTerms_BP_UP_mat_df$V2 <- as.numeric(Combined_GOTerms_BP_UP_mat_df$V2)
Combined_GOTerms_BP_UP_mat_df$V3 <- as.numeric(Combined_GOTerms_BP_UP_mat_df$V3)

colnames(Combined_GOTerms_BP_UP_mat_df) <- c("Name","Time","Con")
Combined_GOTerms_BP_UP_mat_df$Ratio <- Combined_GOTerms_BP_UP_mat_df$Time/Combined_GOTerms_BP_UP_mat_df$Con

Combined_GOTerms_BP_UP_mat_df_TimeTerms <- Combined_GOTerms_BP_UP_mat_df[Combined_GOTerms_BP_UP_mat_df$Ratio >= 2.0,]
Combined_GOTerms_BP_UP_mat_df_TimeTerms <- Combined_GOTerms_BP_UP_mat_df_TimeTerms[Combined_GOTerms_BP_UP_mat_df_TimeTerms$Time >= 3,]
dim(Combined_GOTerms_BP_UP_mat_df_TimeTerms)
Combined_GOTerms_BP_UP_mat_df_ConTerms <- Combined_GOTerms_BP_UP_mat_df[Combined_GOTerms_BP_UP_mat_df$Ratio <= 0.5,]
Combined_GOTerms_BP_UP_mat_df_ConTerms <- Combined_GOTerms_BP_UP_mat_df_ConTerms[Combined_GOTerms_BP_UP_mat_df_ConTerms$Con >= 3,]
dim(Combined_GOTerms_BP_UP_mat_df_ConTerms)

### Figure 4 V2 ###
AllFiltered_GO_Terms <- rbind(Combined_GOTerms_BP_UP_mat_df_TimeTerms,Combined_GOTerms_BP_UP_mat_df_ConTerms,Combined_GOTerms_MF_UP_mat_df_TimeTerms,Combined_GOTerms_MF_UP_mat_df_ConTerms)
AllFiltered_GO_Terms$Diff <- AllFiltered_GO_Terms$Con - AllFiltered_GO_Terms$Time
dim(AllFiltered_GO_Terms)

#View(as.data.frame(unique(AllFiltered_GO_Terms$Name)))
length(unique(AllFiltered_GO_Terms$Name))
# 102

# Select 5 for each
SelectCluster_GO <- rbind(
  
  # Con --- Best 5 Differential ---
  c('stress response to metal ion',0,4,4),
  c('extracellular matrix structural constituent',0,3,3),
  c('positive regulation of cell motility',0,3,3),
  c('integrin binding',0,3,3),
  c('endoderm development',0,3,3),
  
  # Time --- Best 5 Differential ---
  c('oxidative phosphorylation',3,0,-3),
  c('neuron death',3,0,-3),
  c('neutrophil activation',4,1,-3),
  c('cellular detoxification',4,1,-3),
  c('generation of precursor metabolites and energy',4,0,-4)
)

SelectCluster_GO_df <- as.data.frame(SelectCluster_GO)
SelectCluster_GO_df$V1 <- as.character(SelectCluster_GO_df$V1)
SelectCluster_GO_df$V2 <- as.numeric(SelectCluster_GO_df$V2)
SelectCluster_GO_df$V3 <- as.numeric(SelectCluster_GO_df$V3)

SelectCluster_GO_long <- rbind(cbind(SelectCluster_GO_df$V1,SelectCluster_GO_df$V2),cbind(SelectCluster_GO_df$V1,SelectCluster_GO_df$V3))
SelectCluster_GO_long_df <- as.data.frame(SelectCluster_GO_long)
SelectCluster_GO_long_df$V1 <- as.character(SelectCluster_GO_long_df$V1)
SelectCluster_GO_long_df$V2 <- as.numeric(SelectCluster_GO_long_df$V2)
SelectCluster_GO_long_df$V3 <- "Time"
SelectCluster_GO_long_df$V3[11:20] <- "Con"

# Order
SelectCluster_GO_long_df$V1 <- factor(SelectCluster_GO_long_df$V1,levels=c('stress response to metal ion',
                                                                           'extracellular matrix structural constituent',
                                                                           'positive regulation of cell motility',
                                                                           'integrin binding',
                                                                           'endoderm development',
                                                                           'oxidative phosphorylation',
                                                                           'neuron death',
                                                                           'neutrophil activation',
                                                                           'cellular detoxification',
                                                                           'generation of precursor metabolites and energy'))

p <- ggplot(SelectCluster_GO_long_df, aes(x=V1, y=V2, fill=V3)) + 
  geom_bar(stat="identity", position=position_dodge())

p + theme_minimal() + coord_flip()
