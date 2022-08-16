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
library(openxlsx)

set.seed(5) # Controls randomness for nnPC, minor variations , mainly to keep enriched gene lists consistent (+/-1 genes)

#################
### Load Data ###
#################

### Load Processed Data ###

integated_meta <- readRDS("Integrated_MetaData.rds")
UMAP_values <- readRDS("Integrated_UMAP_values.rds")

dose_samples <- readRDS("DoseSamples.rds")
time_samples <- readRDS("TimeSamples.rds")

Integrated_E_mat_nsprcomp <- readRDS("Integrated_E_mat_nsprcomp.rds")
Integrated_M_mat_nsprcomp <- readRDS("Integrated_M_mat_nsprcomp.rds")

### Load Genes Sets ###
EMT_genes <- read.table("EMTGeneUpdate3.txt",header=TRUE)
E_genes <- subset(EMT_genes, EMT_genes$Annotation == 'E')
M_genes <- subset(EMT_genes, EMT_genes$Annotation == 'M')
Gene_Sets <- list('Egenes'=E_genes$Gene,'Mgenes'=M_genes$Gene)

### Load Correlation Data ###

# Calcluate Correlation with Time Samples
correl_TimeE <- readRDS("correl_TimeE.rds")
#correl_TimeM <- readRDS("correl_TimeM.rds")
correl_TimeM <- readRDS("correl_TimeM2.rds")

# correl_TimeE <- readRDS("correl_TimeE_wodoublets.rds")
# correl_TimeM <- readRDS("correl_TimeM2_wodoublts.rds")

# Calcualte Correlation with Dose Samples
correl_ConE <- readRDS("correl_ConE.rds")
#correl_ConM <- readRDS("correl_ConM.rds")
correl_ConM <- readRDS("correl_ConM2.rds")

all_genes <- readRDS("AnitCorr_all_genes.rds")

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

# Time Scatter Plot
Time_25_genes <- intersect(names(correl_TimeE_Bot25),names(correl_TimeM_Top25))
Time_25_dn_genes <- intersect(names(correl_TimeE_Top25),names(correl_TimeM_Bot25))

#pdf(file="Figures/Panels/Fig3_TimeCorr.pdf")
par(pty="s")
plot(correl_TimeE[,1],correl_TimeM[,1],pch=16,xlab ="Corr. with EPC1", ylab='Corr. with M PC2')
points(correl_TimeE[Time_25_genes,1],correl_TimeM[Time_25_genes,1],pch=16,col='red')
points(correl_TimeE[Time_25_dn_genes,1],correl_TimeM[Time_25_dn_genes,1],pch=16,col='blue')
cor(correl_TimeE[,1],correl_TimeM[,1])
#dev.off()

# Dose Scatter Plot
Con_25_genes <- intersect(names(correl_ConE_Bot25),names(correl_ConM_Top25))
Con_25_dn_genes <- intersect(names(correl_ConE_Top25),names(correl_ConM_Bot25))

#pdf(file="Figures/Panels/Fig3_ConCorr.pdf")
par(pty="s")
plot(correl_ConE[,1],correl_ConM[,1],pch=16,xlab ="Corr. with EPC1", ylab='Corr. with M PC2')
points(correl_ConE[Con_25_genes,1],correl_ConM[Con_25_genes,1],pch=16,col='red')
points(correl_ConE[Con_25_dn_genes,1],correl_ConM[Con_25_dn_genes,1],pch=16,col='blue')
cor(correl_ConE[,1],correl_ConM[,1])
#dev.off()


# UP 25
Time_genes_up <- intersect(names(correl_TimeE_Bot25),names(correl_TimeM_Top25))
Con_genes_up <- intersect(names(correl_ConE_Bot25),names(correl_ConM_Top25))

# Overlapping Genes
length(Time_genes_up) # 57
length(Con_genes_up) # 114
length(intersect(Time_genes_up,Con_genes_up)) # 39

# Unique Genes
Con_only_up <- setdiff(Con_genes_up,Time_genes_up)
Con_only_up # 75

Time_only_up <- setdiff(Time_genes_up,Con_genes_up)
Time_only_up # 18

# DOWN 25
Time_genes_dn <- intersect(names(correl_TimeE_Top25),names(correl_TimeM_Bot25))
Con_genes_dn <- intersect(names(correl_ConE_Top25),names(correl_ConM_Bot25))

# Overlapping genes
length(Time_genes_dn) # 45
length(Con_genes_dn) # 97
length(intersect(Time_genes_dn,Con_genes_dn)) # 38

# Unique ge
Con_only_dn <- setdiff(Con_genes_dn,Time_genes_dn)
Con_only_dn # 59

Time_only_dn <- setdiff(Time_genes_dn,Con_genes_dn)
Time_only_dn # 7

# BiocManager::install("clusterProfiler")
# BiocManager::install("GO.db")
# BiocManager::install("GO.db")
# BiocManager::install("org.Hs.eg.db")
library(clusterProfiler)
library(GOSemSim)
library(GO.db)
library(org.Hs.eg.db)
library(HGNChelper)

#######################
### Dose Only Genes ###
#######################

# Con Only UP Biological Process
Con_only_up_GO <- enrichGO(gene = Con_only_up,
                           universe      = all_genes,
                           OrgDb         = org.Hs.eg.db,
                           keyType = "SYMBOL",
                           ont           = "ALL",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 0.05,
                           readable      = FALSE)
Con_only_up_GO
saveRDS(Con_only_up_GO,"Con_only_up_GO_M2.rds")
write.xlsx(Con_only_up_GO@result,"Dose_M2_Mprogram_GOTerms.xlsx")

# Con Only DN Biological Process
Con_only_dn_GO <- enrichGO(gene = Con_only_dn,
                           universe      = all_genes,
                           OrgDb         = org.Hs.eg.db,
                           keyType = "SYMBOL",
                           ont           = "ALL",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 0.05,
                           readable      = FALSE)
Con_only_dn_GO
saveRDS(Con_only_dn_GO,"Con_only_dn_GO_M2.rds")
write.xlsx(Con_only_dn_GO@result,"Dose_M2_Eprogram_GOTerms.xlsx")


#######################
### Time Only Genes ###
#######################

# Time ONly UP Biological Process
Time_only_up_GO <- enrichGO(gene = Time_only_up,
                           universe      = all_genes,
                           OrgDb         = org.Hs.eg.db,
                           keyType = "SYMBOL",
                           ont           = "ALL",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 0.05,
                           readable      = FALSE)

saveRDS(Time_only_up_GO,"Time_only_up_GO_M2.rds")
write.xlsx(Time_only_up_GO@result,"Time_M2_Mprogram_GOTerms.xlsx")


# Time ONly DN Biological Process
Time_only_dn_GO <- enrichGO(gene = Time_only_dn,
                           universe      = all_genes,
                           OrgDb         = org.Hs.eg.db,
                           keyType = "SYMBOL",
                           ont           = "ALL",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 0.05,
                           readable      = FALSE)
saveRDS(Time_only_dn_GO,"Time_only_dn_GO_M2.rds")
write.xlsx(Time_only_dn_GO@result,"Time_M2_Eprogram_GOTerms.xlsx")

### Look at Genes with large difference in Mcorr between Data sets ###

Dose_AllMCorr <- as.data.frame(cbind(correl_ConM[Con_only_up,],correl_TimeM[Con_only_up,],correl_ConM[Con_only_up,]-correl_TimeM[Con_only_up,]))
## Dose
# CADM3 -- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3622948/, correlated change with VIM and FN1, but weird paper, not sure which line is supposed to be M
# SYN1 --- obliqued reference to EMT https://www.assaygenie.com/content/Antibody%20Genie/PACO/PACO19833.pdf
# NKAIN1 --- cancer marker gene (https://academic.oup.com/noa/article/3/1/vdab144/6386379,https://academic.oup.com/noa/article/3/1/vdab144/6386379)
#        --- but part of a gene set anti-correlated with MEK activated EMT (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8131078/)
Time_AllMCorr <- as.data.frame(cbind(correl_ConM[Time_only_up,],correl_TimeM[Time_only_up,],correl_TimeM[Time_only_up,]-correl_ConM[Time_only_up,]))
## Time
#SLC20A1 -- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6755178/ (Down Regulation?)
#AKR1C1 -- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5048132/
#HPCAL1 -- Wnt signaling (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6484330/), but also downregulated in EMT in KK47 cancer cells (https://www.mdpi.com/1420-3049/21/1/84/htm)
