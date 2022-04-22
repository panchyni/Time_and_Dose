### Folder
setwd("<path_to_folder>/Time_and_Dose")

### IMPORTANT NOTE ###
# Some wonkiness is expected with UMAP across systems, even with set seed
# https://stackoverflow.com/questions/67101829/seurat-umap-visualization-result-is-mirrored-after-running-in-two-identical-envi

# Dependencies
# install.packages('https://cran.r-project.org/src/contrib/Archive/spatstat/spatstat_1.64-1.tar.gz', repos=NULL,type="source", INSTALL_opts = "--no-lock") 
library(Seurat) #v3.x
library(viridis)
library(cowplot)
library(ggplot2)
library(scales)
library(HGNChelper)

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
  
  base_plot_centers <- base_plot + geom_point(data=data_ceneters, aes(x=Emean, y=Mmean),fill=manual_colors,pch=21,cex=4.0) + geom_errorbar(data=data_ceneters, aes(x=Emean, y=Mmean, ymin=Mmean-Msd, ymax=Mmean+Msd), colour="black", width=E_20th) + geom_errorbarh(data=data_ceneters, aes(x=Emean, y=Mmean, xmin=Emean-Esd, xmax=Emean+Esd), colour="black", height=M_20th) + xlim(x_down,x_up) + ylim(y_down,y_up) + xlab("Egenes") + ylab("Mgenes")
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

#####################
### Concentration ###
#####################

ConCC_integration.combined <- readRDS("SeuratData/Seurat_concentration_DoseOnly_cellcycle_correct_hg38.rds")

# Fix annotation
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
# The missing ten are filtered (these are the 10 NA sample post stringent filtering)

p1 <- DimPlot(ConCC_integration.combined, reduction = "umap",group.by = "Batch")
p2 <- DimPlot(ConCC_integration.combined, reduction = "umap",group.by = "Time")
p3 <- DimPlot(ConCC_integration.combined, reduction = "umap", label = TRUE, repel = TRUE)
p4 <- DimPlot(ConCC_integration.combined, reduction = "umap",group.by = "Phase")

p1
p2
p3
p4

# GBC Labels #

ConCC_integration.combined.labeled <- subset(x = ConCC_integration.combined, subset = GBC != 'GBC0')
ConCC_integration.combined.unlabeled <- subset(x = ConCC_integration.combined, subset = GBC == 'GBC0')

ConCC_integration.combined.labeled@meta.data$GBC_pM_factor <- factor(ConCC_integration.combined.labeled@meta.data$GBC_pM)
ConCC_integration.combined.unlabeled@meta.data$GBC_pM_factor <- factor(ConCC_integration.combined.unlabeled@meta.data$GBC_pM)

UMAP_GBC_Labeled <- DimPlot(ConCC_integration.combined.labeled, reduction = "umap",group.by = "GBC_pM_factor")
UMAP_GBC_UnLabeled <- DimPlot(ConCC_integration.combined.unlabeled, reduction = "umap",group.by = "GBC_pM_factor")

### Figure 1 ###
pdf(file="Figures/Panels/Fig1_UMAP_GBC_Labeled.pdf")
UMAP_GBC_Labeled
dev.off()

pdf(file="Figures/Panels/FigS1_UMAP_GBC_UnLabeled.pdf")
UMAP_GBC_UnLabeled
dev.off()

# UMAP contour #
UMAP_values <- ConCC_integration.combined[["umap"]]@cell.embeddings
UMAP_annot <- merge(UMAP_values,ConCC_integration.combined@meta.data,by=0)

ggplot(UMAP_annot, aes(x = UMAP_1, y = UMAP_2, colour=Batch)) +
  geom_density_2d(bins=6,alpha=0.5,size=1) + theme_bw()

UMAP_annot_Con <- UMAP_annot[UMAP_annot$GBC_pM > -1,]
UMAP_annot_Else <- UMAP_annot[UMAP_annot$GBC_pM == -100,]
UMAP_annot_Con$GBC_pM <- as.factor(UMAP_annot_Con$GBC_pM)
dim(UMAP_annot_Con)

ggplot(UMAP_annot_Con, aes(x = UMAP_1, y = UMAP_2, colour=GBC_pM)) +
  geom_point(size=0.5) + theme_bw()

ggplot(UMAP_annot_Con, aes(x = UMAP_1, y = UMAP_2, colour=GBC_pM)) +
  geom_density_2d(bins=8,alpha=0.5,size=1) + theme_bw()

# Else
plot(UMAP_annot_Con$UMAP_1,UMAP_annot_Con$UMAP_2, col='black',pch=16,cex=0.25) 
points(UMAP_annot_Else$UMAP_1,UMAP_annot_Else$UMAP_2, col='red',pch=16,cex=0.25) 

seurat_combo_All_mat <- GetAssayData(object = ConCC_integration.combined, slot = "scale.data")

UMAP_annot$CDH1 <- seurat_combo_All_mat['CDH1',UMAP_annot$Row.names]
UMAP_annot$VIM <- seurat_combo_All_mat['VIM',UMAP_annot$Row.names]
UMAP_annot$FN1 <- seurat_combo_All_mat['FN1',UMAP_annot$Row.names]

### Figure 1 ###
pdf(file="Figures/Panels/Fig1_UMAP_CDH1.pdf")
ggplot(UMAP_annot, aes(x = UMAP_1, y = UMAP_2, colour=CDH1)) +
  geom_point(size=0.5) + theme_bw() + scale_color_gradient2(midpoint=0, low="darkblue", mid="white",
                                                            high="red", space ="Lab",limits=c(-3,3))
dev.off()

pdf(file="Figures/Panels/Fig1_UMAP_VIM.pdf")
ggplot(UMAP_annot, aes(x = UMAP_1, y = UMAP_2, colour=VIM)) +
  geom_point(size=0.5) + theme_bw() + scale_color_gradient2(midpoint=0, low="darkblue", mid="white",
                                                            high="red", space ="Lab",limits=c(-3,3))
dev.off()

pdf(file="Figures/Panels/Fig1_UMAP_FN1.pdf")
ggplot(UMAP_annot, aes(x = UMAP_1, y = UMAP_2, colour=FN1)) +
  geom_point(size=0.5) + theme_bw() + scale_color_gradient2(midpoint=0, low="darkblue", mid="white",
                                                            high="red", space ="Lab",limits=c(-3,3))
dev.off()

# Genes Sets #
EMT_genes <- read.table("EMTGenesUpdateV2.txt",header=TRUE)
E_genes <- subset(EMT_genes, EMT_genes$Annotation == 'E')
M_genes <- subset(EMT_genes, EMT_genes$Annotation == 'M')
Gene_Sets <- list('Egenes'=E_genes$Gene,'Mgenes'=M_genes$Gene)

### Z-score ###
#seurat_combo_All <- ScaleData(TimeConCC_integration.combined, verbose = FALSE, vars.to.regress=c("S.Score","G2M.Score"))

seurat_combo_Con_mat <- GetAssayData(object = ConCC_integration.combined, slot = "scale.data")
dim(seurat_combo_Con_mat)

library(GSVA)
seurat_combo_Con_zscore <- gsva(as.matrix(seurat_combo_Con_mat),Gene_Sets,method='zscore',kcdf='Gaussian')
seurat_combo_Con_zscore_t <- t(seurat_combo_Con_zscore)
seurat_combo_Con_zscore_t_annot <- cbind(seurat_combo_Con_zscore_t,ConCC_integration.combined@meta.data)
seurat_combo_Con_zscore_t_annot_GBCGood <- seurat_combo_Con_zscore_t_annot[seurat_combo_Con_zscore_t_annot$GBC_pM >= 0,]

cor(seurat_combo_Con_zscore_t_annot_GBCGood$Egenes,seurat_combo_Con_zscore_t_annot_GBCGood$GBC_pM)
cor(seurat_combo_Con_zscore_t_annot_GBCGood$Mgenes,seurat_combo_Con_zscore_t_annot_GBCGood$GBC_pM)

model <- lm(seurat_combo_Con_zscore_t_annot_GBCGood$GBC_pM ~ seurat_combo_Con_zscore_t_annot_GBCGood$Egenes + seurat_combo_Con_zscore_t_annot_GBCGood$Mgenes)
summary(model)

plot(seurat_combo_Con_zscore_t_annot$Egenes,seurat_combo_Con_zscore_t_annot$Mgenes)

plot_data <- as.data.frame(cbind(seurat_combo_Con_zscore_t_annot_GBCGood$Egenes,seurat_combo_Con_zscore_t_annot_GBCGood$Mgenes,seurat_combo_Con_zscore_t_annot_GBCGood$GBC_pM))
colnames(plot_data) <- c("Egenes","Mgenes","key")
plot_data$key_label <- factor(plot_data$key)

manual_colors <- hue_pal()(8)
base_plot_centers <- gg_means_plot(plot_data,8,-7.5,10,-7.5,10)
base_plot_centers

### Figure 1 ###
pdf(file="Figures/Panels/Fig1_Zscore_Progression.pdf")
base_plot_centers
dev.off()

boxplot(seurat_combo_Con_zscore_t_annot_GBCGood$Egenes ~ seurat_combo_Con_zscore_t_annot_GBCGood$GBC_pM,col=manual_colors,xlab="Egenes",y="Freq")
boxplot(seurat_combo_Con_zscore_t_annot_GBCGood$Mgenes ~ seurat_combo_Con_zscore_t_annot_GBCGood$GBC_pM,col=manual_colors,xlab="Mgenes",y="Freq")

### NN-PCA ###
library(nsprcomp)
set.seed(5) # Controls randomness for nnPC, minor variations , mainly to keep enriched gene lists consistent (+/-1 genes)

# For annotation data
seuratCon_combo_E <- ScaleData(ConCC_integration.combined, verbose = FALSE, vars.to.regress=c("S.Score","G2M.Score"), features=E_genes$Gene)
seuratCon_combo_M <- ScaleData(ConCC_integration.combined, verbose = FALSE, vars.to.regress=c("S.Score","G2M.Score"), features=M_genes$Gene)

#seuratCon_combo_E_mat <- GetAssayData(object = seuratCon_combo_E, slot = "scale.data")
#seuratCon_combo_M_mat <- GetAssayData(object = seuratCon_combo_M, slot = "scale.data")

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

# Check #
dim(seuratCon_combo_E_mat_mat_nsprcomp_annot)
dim(seuratCon_combo_M_mat_mat_nsprcomp_annot)

seuratCon_combo_E_mat_mat_nsprcomp_annot_GBCGood <- seuratCon_combo_E_mat_mat_nsprcomp_annot[seuratCon_combo_E_mat_mat_nsprcomp_annot$GBC_pM >= 0,]
seuratCon_combo_M_mat_mat_nsprcomp_annot_GBCGood <- seuratCon_combo_M_mat_mat_nsprcomp_annot[seuratCon_combo_M_mat_mat_nsprcomp_annot$GBC_pM >= 0,]

cor(seuratCon_combo_E_mat_mat_nsprcomp_annot_GBCGood$PC1,as.numeric(seuratCon_combo_E_mat_mat_nsprcomp_annot_GBCGood$GBC_pM))
cor(seuratCon_combo_M_mat_mat_nsprcomp_annot_GBCGood$PC1,as.numeric(seuratCon_combo_M_mat_mat_nsprcomp_annot_GBCGood$GBC_pM))

# Model
lm_fit <- lm(as.numeric(seuratCon_combo_E_mat_mat_nsprcomp_annot_GBCGood$GBC_pM) ~ seuratCon_combo_E_mat_mat_nsprcomp_annot_GBCGood$PC1 + seuratCon_combo_M_mat_mat_nsprcomp_annot_GBCGood$PC1)
summary(lm_fit)

### Progression ###
plot_data <- as.data.frame(cbind(seuratCon_combo_E_mat_mat_nsprcomp_annot_GBCGood$PC1,seuratCon_combo_M_mat_mat_nsprcomp_annot_GBCGood$PC1,seuratCon_combo_E_mat_mat_nsprcomp_annot_GBCGood$GBC_pM))
colnames(plot_data) <- c("Egenes","Mgenes","key")
plot_data$key_label <- factor(plot_data$key)

manual_colors <- hue_pal()(8)
base_plot_centers <- gg_means_plot(plot_data,8,-7.5,10,-7.5,10)
base_plot_centers

### Figure 1 ###
pdf(file="Figures/Panels/Fig1_EM_Progression.pdf")
base_plot_centers
dev.off()

boxplot(seuratCon_combo_E_mat_mat_nsprcomp_annot_GBCGood$PC1 ~ seuratCon_combo_E_mat_mat_nsprcomp_annot_GBCGood$GBC_pM,col=manual_colors,xlab="Egenes",y="Freq")
boxplot(seuratCon_combo_M_mat_mat_nsprcomp_annot_GBCGood$PC1 ~ seuratCon_combo_M_mat_mat_nsprcomp_annot_GBCGood$GBC_pM,col=manual_colors,xlab="Mgenes",y="Freq")

### Supplemental Figure 1 ###
plot_data <- as.data.frame(cbind(seuratCon_combo_E_mat_mat_nsprcomp_annot_GBCGood$PC1,seuratCon_combo_M_mat_mat_nsprcomp_annot_GBCGood$PC1,seuratCon_combo_E_mat_mat_nsprcomp_annot_GBCGood$seurat_clusters))
colnames(plot_data) <- c("Egenes","Mgenes","key")
plot_data$key_label <- factor(plot_data$key)

manual_colors <- hue_pal()(7)
base_plot_centers <- gg_means_plot(plot_data,7,-7.5,10,-7.5,10)
base_plot_centers

### Log Models ###

# Z-score
model <- lm(log(seurat_combo_Con_zscore_t_annot_GBCGood$GBC_pM+1) ~ seurat_combo_Con_zscore_t_annot_GBCGood$Egenes + seurat_combo_Con_zscore_t_annot_GBCGood$Mgenes)
summary(model)

# nnPCA
timecourse_all_lm_fit <- lm(log(as.numeric(seuratCon_combo_E_mat_mat_nsprcomp_annot_GBCGood$GBC_pM+1)) ~ seuratCon_combo_E_mat_mat_nsprcomp_annot_GBCGood$PC1 + seuratCon_combo_M_mat_mat_nsprcomp_annot_GBCGood$PC1)
summary(timecourse_all_lm_fit)

### Mann Whitney U-test ###

#  Zscore E GBC #
ZScore_GBC_1 <- seurat_combo_Con_zscore_t_annot[seurat_combo_Con_zscore_t_annot$GBC_values == 1, ]
ZScore_GBC_2 <- seurat_combo_Con_zscore_t_annot[seurat_combo_Con_zscore_t_annot$GBC_values == 2, ]
ZScore_GBC_3 <- seurat_combo_Con_zscore_t_annot[seurat_combo_Con_zscore_t_annot$GBC_values == 3, ]
ZScore_GBC_4 <- seurat_combo_Con_zscore_t_annot[seurat_combo_Con_zscore_t_annot$GBC_values == 4, ]
ZScore_GBC_8 <- seurat_combo_Con_zscore_t_annot[seurat_combo_Con_zscore_t_annot$GBC_values == 8, ]
ZScore_GBC_9 <- seurat_combo_Con_zscore_t_annot[seurat_combo_Con_zscore_t_annot$GBC_values == 9, ]
ZScore_GBC_12 <- seurat_combo_Con_zscore_t_annot[seurat_combo_Con_zscore_t_annot$GBC_values == 12, ]
ZScore_GBC_13 <- seurat_combo_Con_zscore_t_annot[seurat_combo_Con_zscore_t_annot$GBC_values == 13, ]

dim(ZScore_GBC_1)
dim(ZScore_GBC_2)
dim(ZScore_GBC_3)
dim(ZScore_GBC_4)
dim(ZScore_GBC_8)
dim(ZScore_GBC_9)
dim(ZScore_GBC_12)
dim(ZScore_GBC_13)

Zscore_E_GBC_1vs2 <- wilcox.test(ZScore_GBC_1$Egenes,ZScore_GBC_2$Egenes,alternative='greater')
Zscore_E_GBC_1vs2$statistic/(length(ZScore_GBC_1$Egenes)*length(ZScore_GBC_2$Egenes))
Zscore_E_GBC_1vs2$p.value

Zscore_E_GBC_2vs3 <- wilcox.test(ZScore_GBC_2$Egenes,ZScore_GBC_3$Egenes,alternative='greater')
Zscore_E_GBC_2vs3$statistic/(length(ZScore_GBC_2$Egenes)*length(ZScore_GBC_3$Egenes))
Zscore_E_GBC_2vs3$p.value

Zscore_E_GBC_3vs4 <- wilcox.test(ZScore_GBC_3$Egenes,ZScore_GBC_4$Egenes,alternative='greater')
Zscore_E_GBC_3vs4$statistic/(length(ZScore_GBC_3$Egenes)*length(ZScore_GBC_4$Egenes))
Zscore_E_GBC_3vs4$p.value

Zscore_E_GBC_4vs8 <- wilcox.test(ZScore_GBC_4$Egenes,ZScore_GBC_8$Egenes,alternative='greater')
Zscore_E_GBC_4vs8$statistic/(length(ZScore_GBC_4$Egenes)*length(ZScore_GBC_8$Egenes))
Zscore_E_GBC_4vs8$p.value

Zscore_E_GBC_8vs9 <- wilcox.test(ZScore_GBC_8$Egenes,ZScore_GBC_9$Egenes,alternative='greater')
Zscore_E_GBC_8vs9$statistic/(length(ZScore_GBC_8$Egenes)*length(ZScore_GBC_9$Egenes))
Zscore_E_GBC_8vs9$p.value

Zscore_E_GBC_9vs12 <- wilcox.test(ZScore_GBC_9$Egenes,ZScore_GBC_12$Egenes,alternative='greater')
Zscore_E_GBC_9vs12$statistic/(length(ZScore_GBC_9$Egenes)*length(ZScore_GBC_12$Egenes))
Zscore_E_GBC_9vs12$p.value

Zscore_E_GBC_12vs13 <- wilcox.test(ZScore_GBC_12$Egenes,ZScore_GBC_13$Egenes,alternative='greater')
Zscore_E_GBC_12vs13$statistic/(length(ZScore_GBC_12$Egenes)*length(ZScore_GBC_13$Egenes))
Zscore_E_GBC_12vs13$p.value

#  Zscore M GBC #
Zscore_M_GBC_1vs2 <- wilcox.test(ZScore_GBC_1$Mgenes,ZScore_GBC_2$Mgenes,alternative='less')
1-Zscore_M_GBC_1vs2$statistic/(length(ZScore_GBC_1$Mgenes)*length(ZScore_GBC_2$Mgenes))
Zscore_M_GBC_1vs2$p.value

Zscore_M_GBC_2vs3 <- wilcox.test(ZScore_GBC_2$Mgenes,ZScore_GBC_3$Mgenes,alternative='less')
1-Zscore_M_GBC_2vs3$statistic/(length(ZScore_GBC_2$Mgenes)*length(ZScore_GBC_3$Mgenes))
Zscore_M_GBC_2vs3$p.value

Zscore_M_GBC_3vs4 <- wilcox.test(ZScore_GBC_3$Mgenes,ZScore_GBC_4$Mgenes,alternative='less')
1-Zscore_M_GBC_3vs4$statistic/(length(ZScore_GBC_3$Mgenes)*length(ZScore_GBC_4$Mgenes))
Zscore_M_GBC_3vs4$p.value

Zscore_M_GBC_4vs8 <- wilcox.test(ZScore_GBC_4$Mgenes,ZScore_GBC_8$Mgenes,alternative='less')
1-Zscore_M_GBC_4vs8$statistic/(length(ZScore_GBC_4$Mgenes)*length(ZScore_GBC_8$Mgenes))
Zscore_M_GBC_4vs8$p.value

Zscore_M_GBC_8vs9 <- wilcox.test(ZScore_GBC_8$Mgenes,ZScore_GBC_9$Mgenes,alternative='less')
1-Zscore_M_GBC_8vs9$statistic/(length(ZScore_GBC_8$Mgenes)*length(ZScore_GBC_9$Mgenes))
Zscore_M_GBC_8vs9$p.value

Zscore_M_GBC_9vs12 <- wilcox.test(ZScore_GBC_9$Mgenes,ZScore_GBC_12$Mgenes,alternative='less')
1-Zscore_M_GBC_9vs12$statistic/(length(ZScore_GBC_9$Mgenes)*length(ZScore_GBC_12$Mgenes))
Zscore_M_GBC_9vs12$p.value

Zscore_M_GBC_12vs13 <- wilcox.test(ZScore_GBC_12$Mgenes,ZScore_GBC_13$Mgenes,alternative='less')
1-Zscore_M_GBC_12vs13$statistic/(length(ZScore_GBC_12$Mgenes)*length(ZScore_GBC_13$Mgenes))
Zscore_M_GBC_12vs13$p.value

# nnPCA E GBC 
nnPCA_Egenes_GBC_1 <- seuratCon_combo_E_mat_mat_nsprcomp_annot[seuratCon_combo_E_mat_mat_nsprcomp_annot$GBC_values == 1, ]
nnPCA_Egenes_GBC_2 <- seuratCon_combo_E_mat_mat_nsprcomp_annot[seuratCon_combo_E_mat_mat_nsprcomp_annot$GBC_values == 2, ]
nnPCA_Egenes_GBC_3 <- seuratCon_combo_E_mat_mat_nsprcomp_annot[seuratCon_combo_E_mat_mat_nsprcomp_annot$GBC_values == 3, ]
nnPCA_Egenes_GBC_4 <- seuratCon_combo_E_mat_mat_nsprcomp_annot[seuratCon_combo_E_mat_mat_nsprcomp_annot$GBC_values == 4, ]
nnPCA_Egenes_GBC_8 <- seuratCon_combo_E_mat_mat_nsprcomp_annot[seuratCon_combo_E_mat_mat_nsprcomp_annot$GBC_values == 8, ]
nnPCA_Egenes_GBC_9 <- seuratCon_combo_E_mat_mat_nsprcomp_annot[seuratCon_combo_E_mat_mat_nsprcomp_annot$GBC_values == 9, ]
nnPCA_Egenes_GBC_12 <- seuratCon_combo_E_mat_mat_nsprcomp_annot[seuratCon_combo_E_mat_mat_nsprcomp_annot$GBC_values == 12, ]
nnPCA_Egenes_GBC_13 <- seuratCon_combo_E_mat_mat_nsprcomp_annot[seuratCon_combo_E_mat_mat_nsprcomp_annot$GBC_values == 13, ]

dim(nnPCA_Egenes_GBC_1)
dim(nnPCA_Egenes_GBC_2)
dim(nnPCA_Egenes_GBC_3)
dim(nnPCA_Egenes_GBC_4)
dim(nnPCA_Egenes_GBC_8)
dim(nnPCA_Egenes_GBC_9)
dim(nnPCA_Egenes_GBC_12)
dim(nnPCA_Egenes_GBC_13)

nnPCA_E_GBC_1vs2 <- wilcox.test(nnPCA_Egenes_GBC_1$PC1,nnPCA_Egenes_GBC_2$PC1,alternative='greater')
nnPCA_E_GBC_1vs2$statistic/(length(nnPCA_Egenes_GBC_1$PC1)*length(nnPCA_Egenes_GBC_2$PC1))
nnPCA_E_GBC_1vs2$p.value

nnPCA_E_GBC_2vs3 <- wilcox.test(nnPCA_Egenes_GBC_2$PC1,nnPCA_Egenes_GBC_3$PC1,alternative='greater')
nnPCA_E_GBC_2vs3$statistic/(length(nnPCA_Egenes_GBC_2$PC1)*length(nnPCA_Egenes_GBC_3$PC1))
nnPCA_E_GBC_2vs3$p.value
 
nnPCA_E_GBC_3vs4 <- wilcox.test(nnPCA_Egenes_GBC_3$PC1,nnPCA_Egenes_GBC_4$PC1,alternative='greater')
nnPCA_E_GBC_3vs4$statistic/(length(nnPCA_Egenes_GBC_3$PC1)*length(nnPCA_Egenes_GBC_4$PC1))
nnPCA_E_GBC_3vs4$p.value

nnPCA_E_GBC_4vs8 <- wilcox.test(nnPCA_Egenes_GBC_4$PC1,nnPCA_Egenes_GBC_8$PC1,alternative='greater')
nnPCA_E_GBC_4vs8$statistic/(length(nnPCA_Egenes_GBC_4$PC1)*length(nnPCA_Egenes_GBC_8$PC1))
nnPCA_E_GBC_4vs8$p.value

nnPCA_E_GBC_8vs9 <- wilcox.test(nnPCA_Egenes_GBC_8$PC1,nnPCA_Egenes_GBC_9$PC1,alternative='greater')
nnPCA_E_GBC_8vs9$statistic/(length(nnPCA_Egenes_GBC_8$PC1)*length(nnPCA_Egenes_GBC_9$PC1))
nnPCA_E_GBC_8vs9$p.value

nnPCA_E_GBC_9vs12 <- wilcox.test(nnPCA_Egenes_GBC_9$PC1,nnPCA_Egenes_GBC_12$PC1,alternative='greater')
nnPCA_E_GBC_9vs12$statistic/(length(nnPCA_Egenes_GBC_9$PC1)*length(nnPCA_Egenes_GBC_12$PC1))
nnPCA_E_GBC_9vs12$p.value

nnPCA_E_GBC_12vs13 <- wilcox.test(nnPCA_Egenes_GBC_12$PC1,nnPCA_Egenes_GBC_13$PC1,alternative='greater')
nnPCA_E_GBC_12vs13$statistic/(length(nnPCA_Egenes_GBC_12$PC1)*length(nnPCA_Egenes_GBC_13$PC1))
nnPCA_E_GBC_12vs13$p.value

# nnPCA M GBC 
nnPCA_Mgenes_GBC_1 <- seuratCon_combo_M_mat_mat_nsprcomp_annot[seuratCon_combo_M_mat_mat_nsprcomp_annot$GBC_values == 1, ]
nnPCA_Mgenes_GBC_2 <- seuratCon_combo_M_mat_mat_nsprcomp_annot[seuratCon_combo_M_mat_mat_nsprcomp_annot$GBC_values == 2, ]
nnPCA_Mgenes_GBC_3 <- seuratCon_combo_M_mat_mat_nsprcomp_annot[seuratCon_combo_M_mat_mat_nsprcomp_annot$GBC_values == 3, ]
nnPCA_Mgenes_GBC_4 <- seuratCon_combo_M_mat_mat_nsprcomp_annot[seuratCon_combo_M_mat_mat_nsprcomp_annot$GBC_values == 4, ]
nnPCA_Mgenes_GBC_8 <- seuratCon_combo_M_mat_mat_nsprcomp_annot[seuratCon_combo_M_mat_mat_nsprcomp_annot$GBC_values == 8, ]
nnPCA_Mgenes_GBC_9 <- seuratCon_combo_M_mat_mat_nsprcomp_annot[seuratCon_combo_M_mat_mat_nsprcomp_annot$GBC_values == 9, ]
nnPCA_Mgenes_GBC_12 <- seuratCon_combo_M_mat_mat_nsprcomp_annot[seuratCon_combo_M_mat_mat_nsprcomp_annot$GBC_values == 12, ]
nnPCA_Mgenes_GBC_13 <- seuratCon_combo_M_mat_mat_nsprcomp_annot[seuratCon_combo_M_mat_mat_nsprcomp_annot$GBC_values == 13, ]

dim(nnPCA_Mgenes_GBC_1)
dim(nnPCA_Mgenes_GBC_2)
dim(nnPCA_Mgenes_GBC_3)
dim(nnPCA_Mgenes_GBC_4)
dim(nnPCA_Mgenes_GBC_8)
dim(nnPCA_Mgenes_GBC_9)
dim(nnPCA_Mgenes_GBC_12)
dim(nnPCA_Mgenes_GBC_13)

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

##################
### GMM Models ###
##################
library(mclust)

conE <- Mclust(seuratCon_combo_E_mat_mat_nsprcomp_annot$PC1,modelNames = "V")
conM <- Mclust(seuratCon_combo_M_mat_mat_nsprcomp_annot$PC1,modelNames = "V")

conE$BIC
# V
#1 -58129.89
#2 -56945.93
#3 -56814.33
#4 -56833.76
#5 -56854.87
#6 -56873.50
#7 -56891.79
#8 -56918.20
#9 -56942.48

conM$BIC
#1 -54484.44
#2 -53295.82
#3 -53216.51
#4 -53234.18
#5 -53242.48
#6 -53258.35
#7 -53279.84
#8 -53307.02
#9 -53337.72
library("MineICA")
plotMix(conE,seuratCon_combo_E_mat_mat_nsprcomp_annot$PC1,36,title="Dose E")
plotMix(conM,seuratCon_combo_M_mat_mat_nsprcomp_annot$PC1,36,title="Dose M")

###############
### Overlap ###
###############

# Integrated Time #
integration.combined.TimeOnly <- readRDS('SeuratData/Seurat_integrated_TimeOnly_CellCycle_hg38.rds')

p1 <- DimPlot(integration.combined.TimeOnly, reduction = "umap",group.by = "Batch")
p2 <- DimPlot(integration.combined.TimeOnly, reduction = "umap",group.by = "Time")
p3 <- DimPlot(integration.combined.TimeOnly, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2 + p3

seurat_combo_E <- ScaleData(integration.combined.TimeOnly, ,vars.to.regress=c("S.Score","G2M.Score"),verbose = FALSE, features=E_genes$Gene)
seurat_combo_M <- ScaleData(integration.combined.TimeOnly, ,vars.to.regress=c("S.Score","G2M.Score"), verbose = FALSE, features=M_genes$Gene)

### NN-PCA ###
set.seed(5) # Controls randomness for nnPC, minor variations , mainly to keep enriched gene lists consistent (+/-1 genes)

# E
seurat_combo_E_mat <- GetAssayData(object = seurat_combo_E, slot = "scale.data")
dim(seurat_combo_E_mat)

seurat_combo_E_mat_nsprcomp <- nsprcomp(t(seurat_combo_E_mat),nneg=TRUE,ncomp=5)
seurat_combo_E_mat_nsprcomp_annot <- cbind(seurat_combo_E_mat_nsprcomp$x,seurat_combo_E@meta.data)
cor(seurat_combo_E_mat_nsprcomp_annot$PC1,seurat_combo_E_mat_nsprcomp_annot$Time)

# M 
seurat_combo_M_mat <- GetAssayData(object = seurat_combo_M, slot = "scale.data")
dim(seurat_combo_M_mat)

seurat_combo_M_mat_nsprcomp <- nsprcomp(t(seurat_combo_M_mat),nneg=TRUE,ncomp=5)
seurat_combo_M_mat_nsprcomp_annot <- cbind(seurat_combo_M_mat_nsprcomp$x,seurat_combo_M@meta.data)
cor(seurat_combo_M_mat_nsprcomp_annot$PC1,seurat_combo_M_mat_nsprcomp_annot$Time)

model <- lm(seurat_combo_E_mat_nsprcomp_annot$Time ~ seurat_combo_E_mat_nsprcomp_annot$PC1 + seurat_combo_M_mat_nsprcomp_annot$PC1)
summary(model)

boxplot(seurat_combo_E_mat_nsprcomp_annot$PC1 ~ seurat_combo_E_mat_nsprcomp_annot$Time,col=manual_colors,xlab="Egenes",y="Freq")
boxplot(seurat_combo_M_mat_nsprcomp_annot$PC1 ~ seurat_combo_M_mat_nsprcomp_annot$Time,col=manual_colors,xlab="Mgenes",y="Freq")

##########################
### Overlapping Scores ###
##########################

overlapping_Escores <- rbind(seurat_combo_E_mat_nsprcomp_annot[1:14],seuratCon_combo_E_mat_mat_nsprcomp_annot_GBCGood[,1:14])
overlapping_Mscores <- rbind(seurat_combo_M_mat_nsprcomp_annot[1:14],seuratCon_combo_M_mat_mat_nsprcomp_annot_GBCGood[,1:14])

#overlapping_Escores <- rbind(seurat_combo_E_mat_nsprcomp_annot[1:14],seuratCon_combo_E_mat_mat_nsprcomp_annot[,1:14])
#overlapping_Mscores <- rbind(seurat_combo_M_mat_nsprcomp_annot[1:14],seuratCon_combo_M_mat_mat_nsprcomp_annot[,1:14])


BatchA <- row.names(overlapping_Escores[overlapping_Escores$Batch == "A",])
BatchB <- row.names(overlapping_Escores[overlapping_Escores$Batch == "B",])
BatchC <- row.names(overlapping_Escores[overlapping_Escores$Batch == "C",])

### EM Scores ###
plot(overlapping_Escores[BatchB,]$PC1,overlapping_Mscores[BatchB,]$PC1,pch=16,cex=0.5,col=rgb(red=0.0, green=1.0, blue=0.0, alpha=0.33))
points(overlapping_Escores[BatchA,]$PC1,overlapping_Mscores[BatchA,]$PC1,pch=16,cex=0.5,col=rgb(red=1.0, green=0.0, blue=0.0, alpha=0.33))
points(overlapping_Escores[BatchC,]$PC1,overlapping_Mscores[BatchC,]$PC1,pch=16,cex=0.5,col=rgb(red=0.0, green=0.0, blue=1.0, alpha=0.33))
# Supplemental Overlap Figure

##################################
### Supplemental --- Unlabeled ###
##################################

### Labeled Points ###
dim(seuratCon_combo_E_mat_mat_nsprcomp_annot_GBCGood)
dim(seuratCon_combo_M_mat_mat_nsprcomp_annot_GBCGood)
table(seuratCon_combo_M_mat_mat_nsprcomp_annot_GBCGood$seurat_clusters)

plot_data <- as.data.frame(cbind(seuratCon_combo_E_mat_mat_nsprcomp_annot_GBCGood$PC1,seuratCon_combo_M_mat_mat_nsprcomp_annot_GBCGood$PC1,seuratCon_combo_E_mat_mat_nsprcomp_annot_GBCGood$seurat_clusters))
colnames(plot_data) <- c("Egenes","Mgenes","key")
plot_data$key_label <- factor(plot_data$key)

manual_colors <- hue_pal()(7)
dose_means_plot <- gg_means_plot(plot_data,7,-7.5,10.5,-7.5,10)
dose_means_plot

dose_means <- gg_means(plot_data)

seuratCon_combo_E_mat_mat_nsprcomp_annot_Unlabeled <- seuratCon_combo_E_mat_mat_nsprcomp_annot[seuratCon_combo_E_mat_mat_nsprcomp_annot$GBC_pM == -100,]
seuratCon_combo_M_mat_mat_nsprcomp_annot_Unlabeled <- seuratCon_combo_M_mat_mat_nsprcomp_annot[seuratCon_combo_M_mat_mat_nsprcomp_annot$GBC_pM == -100,]
table(seuratCon_combo_E_mat_mat_nsprcomp_annot_Unlabeled$seurat_clusters)
#   0   1   2   3   4   5   6 
#   785 485 352 198 230 271  45 

dim(seuratCon_combo_E_mat_mat_nsprcomp_annot_Unlabeled)
dim(seuratCon_combo_M_mat_mat_nsprcomp_annot_Unlabeled)

plot_data_unlabeled <- as.data.frame(cbind(seuratCon_combo_E_mat_mat_nsprcomp_annot_Unlabeled$PC1,seuratCon_combo_M_mat_mat_nsprcomp_annot_Unlabeled$PC1,seuratCon_combo_E_mat_mat_nsprcomp_annot_Unlabeled$seurat_clusters))
colnames(plot_data_unlabeled) <- c("Egenes","Mgenes","key")
plot_data_unlabeled$key_label <- factor(plot_data_unlabeled$key)

manual_colors <- hue_pal()(7)
dose_means_unlabeled_plot <- gg_means_plot(plot_data_unlabeled,7,-7.5,10.5,-7.5,10)
dose_means_unlabeled_plot

dose_means_unlabeled <- gg_means(plot_data_unlabeled)

plot_data <- as.data.frame(cbind(seuratCon_combo_E_mat_mat_nsprcomp_annot$PC1,seuratCon_combo_M_mat_mat_nsprcomp_annot$PC1))
colnames(plot_data) <- c("Egenes","Mgenes") 
base_plot <- ggplot(plot_data, aes(x=Egenes, y=Mgenes)) + geom_point(size=0.1, alpha=0.75, color='lightgrey') + theme_bw()
base_plot_centers <- base_plot + geom_point(data=dose_means, aes(x=Emean, y=Mmean),fill='black',pch=21,cex=4.0) + geom_errorbar(data=dose_means, aes(x=Emean, y=Mmean, ymin=Mmean-Msd, ymax=Mmean+Msd), colour="black") + geom_errorbarh(data=dose_means, aes(x=Emean, y=Mmean, xmin=Emean-Esd, xmax=Emean+Esd), colour="black")
base_plot_centers <- base_plot_centers + geom_point(data=dose_means_unlabeled, aes(x=Emean, y=Mmean),fill='red',pch=21,cex=4.0) + geom_errorbar(data=dose_means_unlabeled, aes(x=Emean, y=Mmean, ymin=Mmean-Msd, ymax=Mmean+Msd), colour="black") + geom_errorbarh(data=dose_means_unlabeled, aes(x=Emean, y=Mmean, xmin=Emean-Esd, xmax=Emean+Esd), colour="black")
base_plot_centers

