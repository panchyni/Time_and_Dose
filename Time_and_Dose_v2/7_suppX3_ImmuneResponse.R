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
UMAP_values <- readRDS("Integrated_UMAP_values.rds")

dose_samples <- readRDS("DoseSamples.rds")
time_samples <- readRDS("TimeSamples.rds")

Integrated_E_mat_nsprcomp <- readRDS("Integrated_E_mat_nsprcomp.rds")
Integrated_M_mat_nsprcomp <- readRDS("Integrated_M_mat_nsprcomp.rds")

### PC UMAP Plots ###

# Add annotation data to UMAP_values
UMAP_annot <- merge(UMAP_values,integated_meta,by=0)
row.names(UMAP_annot) <- UMAP_annot$Row.names
UMAP_annot$GBC_pM_factor <- factor(UMAP_annot$GBC_pM)
UMAP_annot$Time_factor <- factor(UMAP_annot$Time)

UMAP_samples <- row.names(UMAP_annot)

### Do Immune Analysis

# From MSigDB
PD_UPvsControl <- read.delim("GSE24026_PD1_LIGATION_VS_CTRL_IN_ACT_TCELL_LINE_UP.v7.5.1.gmt",header=FALSE)
# Double checked file by redownloading, ignoring warning
PD_UPvsControl_genes <- PD_UPvsControl[3:199]

# Check against last line from tsv meta data file
#PD_UPvsControl_check <- read.csv("GSE24026_PD1_LIGATION_VS_CTRL_IN_ACT_TCELL_LINE_UP.v7.5.1.tsv.MAPPEDGENES",header=FALSE)
#length(intersect(PD_UPvsControl_genes,PD_UPvsControl_check))
#setdiff(PD_UPvsControl_check,PD_UPvsControl_genes)
#setdiff(PD_UPvsControl_genes,PD_UPvsControl_check)
# Only diff is some unmapped genes

# Check genes
library(HGNChelper)
source("NonOverlappingReplacement.R")
currentmamp <- getCurrentHumanMap()
check_PD1_up <- checkGeneSymbols(PD_UPvsControl_genes, map = currentmamp, unmapped.as.na = FALSE,species='human') # check
# One false with a suggestion
check_PD1_NoDup <- NonOverlappingReplacement(check_PD1_up)
PD_UPvsControl_genes <- check_PD1_NoDup$Suggested.Symbol

# Get and correct scaled data
integated_scale <- readRDS("Integrated_ScaledData.rds")
dim(integated_scale)
# Note: I checked EMT genes at the time of processing the raw data and all but one EMT Gene is present
# in the gene annotation used at the time (the 2019 one stored with HGNChelper). The missing is OR7E14P
# which is not present in the filtered CellRanger genome index as of the time we processed the data
# (Its present in unfiltered, full sett of genome features including what I am guess are low conf. pseduogenes
# and other features )
check_data <- checkGeneSymbols(row.names(integated_scale), map = currentmamp, unmapped.as.na = FALSE,species='human') # check
check_data_NoDup <- NonOverlappingReplacement(check_data)
row.names(integated_scale) <- check_data_NoDup$Suggested.Symbol

# Make data
Scaled_Integrated_data_PD1up <- integated_scale[row.names(integated_scale) %in% PD_UPvsControl_genes,]
rm(integated_scale)
dim(Scaled_Integrated_data_PD1up)
# Covers 139 of 197 genes

### Run NN-PCA ###
library(nsprcomp)
set.seed(5) # Controls randomness for nnPC, minor variations, mainly to keep enriched gene lists consistent (+/-1 genes)

Integrated_PD1up_mat_nsprcomp <- nsprcomp(t(Scaled_Integrated_data_PD1up),nneg=TRUE,ncomp=25,em_maxiter = 10000,em_tol = 0.00001)

saveRDS(Integrated_PD1up_mat_nsprcomp,"Integrated_PD1up_mat_nsprcomp.rds")
Integrated_PD1up_mat_nsprcomp <- readRDS("Integrated_PD1up_mat_nsprcomp.rds")

UMAP_annot$PD1_up <- Integrated_PD1up_mat_nsprcomp$x[UMAP_samples,1]

# DUMAP Plot
ggplot(UMAP_annot, aes(x = UMAP_1, y = UMAP_2, colour=PD1_up)) +
  geom_point(size=0.5) + theme_bw() + scale_color_gradient2(midpoint=0, low="darkblue", mid="white",
                                                            high="red", space ="Lab",guide=guide_colorbar(direction="horizontal")) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1,legend.position = "bottom")

cor(UMAP_annot[dose_samples,]$GBC_pM,UMAP_annot[dose_samples,]$PD1_up) # 0.347414
cor(UMAP_annot[time_samples,]$Time,UMAP_annot[time_samples,]$PD1_up) # 0.3441808
