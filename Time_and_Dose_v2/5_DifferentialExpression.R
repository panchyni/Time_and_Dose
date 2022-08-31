########################
### Setup Enviroment ###
########################

### Folder
setwd("/Users/nicholaspanchy/Documents/Work_UTK/DoseTime_FromScratch2/")
#setwd("/home/panchy/DoseTIme/_fromscratch_server")

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

###########################
### Load Processed Data ###
###########################

integated_meta <- readRDS("Integrated_MetaData.rds")
UMAP_values <- readRDS("Integrated_UMAP_values.rds")

dose_samples <- readRDS("DoseSamples.rds")
time_samples <- readRDS("TimeSamples.rds")

Integrated_E_mat_nsprcomp <- readRDS("Integrated_E_mat_nsprcomp.rds")
Integrated_M_mat_nsprcomp <- readRDS("Integrated_M_mat_nsprcomp.rds")

# Bubble Data
Con_Samples <- dose_samples
Time_Samples <- time_samples
Labeled_Samples <- c(Time_Samples,Con_Samples)

Bounds_IQR_med <- readRDS("Bounds_IQR_med.rds")
odd_mat_long <- readRDS("BubbleOddsMat.rds")
quart_clusters <- readRDS("EMQuartileClusters.rds")

# Use unintegrated data  for Marker ID
integration.combined.All <- readRDS("SeuratData/Seurat_integrated_TimeAndCon_FullCorrection_hg38_Processed.rds")
integration.combined.All <- subset(integration.combined.All, subset = GBC != 'GBC0' | Batch == 'A' | Batch == 'B')
integration.combined.All@meta.data$label <- 0

# Replace label with quart_clusters
integration.combined.All@meta.data[quart_clusters[[1]],]$label <- 1
integration.combined.All@meta.data[quart_clusters[[2]],]$label <- 2
integration.combined.All@meta.data[quart_clusters[[3]],]$label <- 3
integration.combined.All@meta.data[quart_clusters[[4]],]$label <- 4
integration.combined.All@meta.data[quart_clusters[[5]],]$label <- 5
integration.combined.All@meta.data[quart_clusters[[6]],]$label <- 6
integration.combined.All@meta.data[quart_clusters[[7]],]$label <- 7
integration.combined.All@meta.data[quart_clusters[[8]],]$label <- 8
integration.combined.All@meta.data[quart_clusters[[9]],]$label <- 9
integration.combined.All@meta.data[quart_clusters[[10]],]$label <- 10
integration.combined.All@meta.data[quart_clusters[[11]],]$label <- 11
integration.combined.All@meta.data[quart_clusters[[12]],]$label <- 12
integration.combined.All@meta.data[quart_clusters[[13]],]$label <- 13
integration.combined.All@meta.data[quart_clusters[[14]],]$label <- 14
integration.combined.All@meta.data[quart_clusters[[15]],]$label <- 15
integration.combined.All@meta.data[quart_clusters[[16]],]$label <- 16

table(integration.combined.All@meta.data$label)
saveRDS(integration.combined.All,"SeuratData/Seurat_integrated_TimeAndCon_FullCorrection_hg38_ProcessedforMarkers.rds")
integration.combined.All <- readRDS("SeuratData/Seurat_integrated_TimeAndCon_FullCorrection_hg38_ProcessedforMarkers.rds")
DefaultAssay(integration.combined.All) <- "RNA"

# Check Assignments
Escores <- as.data.frame(Integrated_E_mat_nsprcomp$x[union(Time_Samples,Con_Samples),1])
Mscores <- as.data.frame(Integrated_M_mat_nsprcomp$x[union(Time_Samples,Con_Samples),1])
MetaData <- integration.combined.All@meta.data

MergedScores <- base::merge(Escores, Mscores, by="row.names")
row.names(MergedScores) <- MergedScores$Row.names
MergedScores <- MergedScores[,2:3]
colnames(MergedScores) <- c("Escores","Mscores")
MergedData<- base::merge(MergedScores, MetaData, by="row.names")

check_figure_gg <- ggplot(data = MergedData, aes(x=Escores, y=Mscores, color=label)) + geom_point()
check_figure_gg

# Find All Markers is one vs.rest, see https://www.biostars.org/p/409790/
integration.combined.All@active.ident <- as.factor(integration.combined.All$label)

AllClusterMarkers_Wilcox_QuartClusters <- FindAllMarkers(integration.combined.All,test.use="LR",latent.vars=c("nCount_RNA","percent.mito","S.Score","G2M.Score","Batch"))
saveRDS(AllClusterMarkers_Wilcox_QuartClusters,"AllClusterMarkers_SeuratLRBatchCorr_AllData.rds")

### Alternate Quarts ###
integration.combined.All.Time <- subset(integration.combined.All, subset = Batch %in% c("A","B"))
integration.combined.All.Con <- subset(integration.combined.All, subset = Batch %in% c("C"))

AllClusterMarkers_TimeOnly_Quart <- FindAllMarkers(integration.combined.All.Time,test.use="LR",latent.vars=c("nCount_RNA","percent.mito","S.Score","G2M.Score","Batch"))
saveRDS(AllClusterMarkers_TimeOnly_Quart,"AllClusterMarkers_SeuratLRBatchCorr_TimeData.rds")
AllClusterMarkers_ConOnly_Quart <- FindAllMarkers(integration.combined.All.Con,test.use="LR",latent.vars=c("nCount_RNA","percent.mito","S.Score","G2M.Score"))
saveRDS(AllClusterMarkers_ConOnly_Quart,"AllClusterMarkers_SeuratLRBatchCorr_ConData.rds")

# All all results for CDH1, EPCAM, VIM, and FN1 for plots
AllClusterMarkers_EMTGenes_AllResults <- FindAllMarkers(integration.combined.All,test.use="LR",latent.vars=c("nCount_RNA","percent.mito","S.Score","G2M.Score","Batch"),logfc.threshold=0,min.pct=0,return.thresh=1,features=c("VIM","FN1","CDH1","EPCAM"))
saveRDS(AllClusterMarkers_EMTGenes_AllResults,"AllClusterMarkers_EMTGenes_AllResults.rds")

# Middle samples
integration.combined.Middle <- subset(integration.combined.All, subset = label %in% c(6,7,10,11))
print(dim(integration.combined.Middle))
print(table(integration.combined.Middle$label))

integration.combined.Middle@meta.data$Treatment <- "Time"
integration.combined.Middle@meta.data[integration.combined.Middle@meta.data$Batch == "C",]$Treatment <- "Dose"
integration.combined.Middle@active.ident <- as.factor(integration.combined.Middle$Treatment)

MiddleMarkers_EMTGenes_AllResults <- FindAllMarkers(integration.combined.Middle,test.use="LR",latent.vars=c("nCount_RNA","percent.mito","S.Score","G2M.Score"),logfc.threshold=0,min.pct=0,return.thresh=1,features=c("VIM","FN1","CDH1","EPCAM"))
saveRDS(MiddleMarkers_EMTGenes_AllResults,"MiddleMarkers_EMTGenes_AllResults.rds")


#######################
### Compare Results ###
#######################

AllMarkers_Norm <- readRDS("AllClusterMarkers_SeuratLRBatchCorr_AllData.rds")
AllMarkers_Time <- readRDS("AllClusterMarkers_SeuratLRBatchCorr_TimeData.rds")
AllMarkers_Dose <- readRDS("AllClusterMarkers_SeuratLRBatchCorr_ConData.rds")

### Correct Marker for Gene Dropout Scores ### 

# Correct p-values for total number of non-zero following Seurat recommendations (Bonferoni)
# Totals tests = # genes in each data set + 351 rested EMT genes * 16 clusters
AllMarkers_Norm$p_val_adj2 <- pmin(AllMarkers_Norm$p_val*(20085*16 + 351*16),1)
AllMarkers_Time$p_val_adj2 <- pmin(AllMarkers_Time$p_val*(18570*16 + 351*16),1)
AllMarkers_Dose$p_val_adj2 <- pmin(AllMarkers_Dose$p_val*(19186*16 + 351*16),1)

AllMarkers_Norm_adpv05_filter <- AllMarkers_Norm[AllMarkers_Norm$p_val_adj2 < 0.05,]
AllMarkers_Time_adpv05_filter <- AllMarkers_Time[AllMarkers_Time$p_val_adj2 < 0.05,]
AllMarkers_Dose_adpv05_filter <- AllMarkers_Dose[AllMarkers_Dose$p_val_adj2 < 0.05,]

dim(AllMarkers_Norm_adpv05_filter)
dim(AllMarkers_Time_adpv05_filter)
dim(AllMarkers_Dose_adpv05_filter)

drop_files_all <- "DropOutGeneSets/DropFilesAll.txt"
drop_files_time <- "DropOutGeneSets/DropFilesTime.txt"
drop_files_dose <- "DropOutGeneSets/DropFilesDose.txt"

CorrectEMTMarkers <- function(marker_df,drop_files,number_of_genes,logFC=FALSE){
  total_EMT_markers = 0
  drop_EMT_markers = 0
  
  gene_drop_files <- read.delim(drop_files,header=FALSE)
  for (f in gene_drop_files$V1){
    #print(f)
    gene_drop_markers <- readRDS(paste0("DropOutGeneSets/",f))
    
    # Account for number of tests
    gene_drop_markers$p_val_adj2 <- pmin(gene_drop_markers$p_val*(number_of_genes*16 + 351*16),1)
    
    # Select genes by p.value (and LogFC is wanted)
    if(logFC==TRUE){
      gene_drop_select <- gene_drop_markers[gene_drop_markers$p_val_adj2 < 0.05 & abs(gene_drop_markers$avg_log2FC) > 0.25,]
    }
    else{
      gene_drop_select <- gene_drop_markers[gene_drop_markers$p_val_adj2 < 0.05,]
    }
    
    # Compare clusters with marker
    drop_clusters <- gene_drop_select$cluster
    all_clusters <- marker_df[marker_df$gene==unique(gene_drop_markers$gene),]$cluster
    diff <- setdiff(all_clusters,drop_clusters)
    
    total_EMT_markers = total_EMT_markers + length(all_clusters)
    drop_EMT_markers = drop_EMT_markers + length(diff)
    
    for(cluster in diff){
      row.remove <- row.names(marker_df[marker_df$gene==unique(gene_drop_markers$gene) & marker_df$cluster==cluster,])
      #print(marker_df[row.remove,])
      marker_df <- marker_df[!(row.names(marker_df) %in% row.remove),]
    }
  }
  print(paste("Total EMT Markers:",total_EMT_markers))
  print(paste("Droped EMT Markers:",drop_EMT_markers))
  print(paste("Remaiing EMT Markers:", total_EMT_markers - drop_EMT_markers))
  print(paste("Percent Remaining:",(total_EMT_markers - drop_EMT_markers)/(total_EMT_markers)))
  
  return(marker_df)
}

AllMarkers_Norm_adpv05_filter_DropCorrection <- CorrectEMTMarkers(AllMarkers_Norm_adpv05_filter,drop_files_all,20085)
# [1] "Total EMT Markers: 709"
# [1] "Droped EMT Markers: 12"
# [1] "Remaiing EMT Markers: 697"
# [1] "Percent Remaining: 0.983074753173484"
AllMarkers_Time_adpv05_filter_DropCorrection <- CorrectEMTMarkers(AllMarkers_Time_adpv05_filter,drop_files_time,18570)
# [1] "Total EMT Markers: 757"
# [1] "Droped EMT Markers: 14"
# [1] "Remaiing EMT Markers: 743"
# [1] "Percent Remaining: 0.981505944517834"
AllMarkers_Dose_adpv05_filter_DropCorrection <- CorrectEMTMarkers(AllMarkers_Dose_adpv05_filter,drop_files_dose,19186)
# [1] "Total EMT Markers: 496"
# [1] "Droped EMT Markers: 33"
# [1] "Remaiing EMT Markers: 463"
# [1] "Percent Remaining: 0.933467741935484"

# Check 
MergeDropGenes <- function(drop_files){
  gene_drop_files <- read.delim(drop_files,header=FALSE)
  merge_drop_files <- readRDS(paste0("DropOutGeneSets/",gene_drop_files$V1[1]))
  for (f in gene_drop_files$V1[2:length(gene_drop_files$V1)]){
    merge_drop_files <- rbind(merge_drop_files,readRDS(paste0("DropOutGeneSets/",f)))
  }
  return(merge_drop_files)
}

drop_files_all_merge <- MergeDropGenes(drop_files_all)
drop_files_all_merge$p_val_adj2 <- pmin(drop_files_all_merge$p_val*(20085*16 + 351*16),1)
drop_files_all_merge_pv05 <- drop_files_all_merge[drop_files_all_merge$p_val_adj2 < 0.05,]
drop_files_all_merge_pv05$Pair <- paste0(drop_files_all_merge_pv05$cluster,drop_files_all_merge_pv05$gene)
AllMarkers_Norm_adpv05_filter_DropCorrection$Pair <- paste0(AllMarkers_Norm_adpv05_filter_DropCorrection$cluster,AllMarkers_Norm_adpv05_filter_DropCorrection$gene)
length(intersect(AllMarkers_Norm_adpv05_filter_DropCorrection$Pair,drop_files_all_merge_pv05$Pair ))

drop_files_time_merge <- MergeDropGenes(drop_files_time)
drop_files_time_merge$p_val_adj2 <- pmin(drop_files_time_merge$p_val*(18570*16 + 351*16),1)
drop_files_time_merge_pv05 <- drop_files_time_merge[drop_files_time_merge$p_val_adj2 < 0.05,]
drop_files_time_merge_pv05$Pair <- paste0(drop_files_time_merge_pv05$cluster,drop_files_time_merge_pv05$gene)
AllMarkers_Time_adpv05_filter_DropCorrection$Pair <- paste0(AllMarkers_Time_adpv05_filter_DropCorrection$cluster,AllMarkers_Time_adpv05_filter_DropCorrection$gene)
length(intersect(AllMarkers_Time_adpv05_filter_DropCorrection$Pair,drop_files_time_merge_pv05$Pair ))

drop_files_dose_merge <- MergeDropGenes(drop_files_dose)
drop_files_dose_merge$p_val_adj2 <- pmin(drop_files_dose_merge$p_val*(19186*16 + 351*16),1)
drop_files_dose_merge_pv05 <- drop_files_dose_merge[drop_files_dose_merge$p_val_adj2 < 0.05,]
drop_files_dose_merge_pv05$Pair <- paste0(drop_files_dose_merge_pv05$cluster,drop_files_dose_merge_pv05$gene)
AllMarkers_Dose_adpv05_filter_DropCorrection$Pair <- paste0(AllMarkers_Dose_adpv05_filter_DropCorrection$cluster,AllMarkers_Dose_adpv05_filter_DropCorrection$gene)
length(intersect(AllMarkers_Dose_adpv05_filter_DropCorrection$Pair,drop_files_dose_merge_pv05$Pair ))

# Check if sign matters
drop_files_all_merge_pv05$Pair <- paste0(drop_files_all_merge_pv05$cluster,drop_files_all_merge_pv05$gene,sign(drop_files_all_merge_pv05$avg_log2FC))
AllMarkers_Norm_adpv05_filter_DropCorrection$Pair <- paste0(AllMarkers_Norm_adpv05_filter_DropCorrection$cluster,AllMarkers_Norm_adpv05_filter_DropCorrection$gene,sign(AllMarkers_Norm_adpv05_filter_DropCorrection$avg_log2FC))
length(intersect(AllMarkers_Norm_adpv05_filter_DropCorrection$Pair,drop_files_all_merge_pv05$Pair ))

drop_files_time_merge_pv05$Pair <- paste0(drop_files_time_merge_pv05$cluster,drop_files_time_merge_pv05$gene,sign(drop_files_time_merge_pv05$avg_log2FC))
AllMarkers_Time_adpv05_filter_DropCorrection$Pair <- paste0(AllMarkers_Time_adpv05_filter_DropCorrection$cluster,AllMarkers_Time_adpv05_filter_DropCorrection$gene,sign(AllMarkers_Time_adpv05_filter_DropCorrection$avg_log2FC))
length(intersect(AllMarkers_Time_adpv05_filter_DropCorrection$Pair,drop_files_time_merge_pv05$Pair ))

drop_files_dose_merge_pv05$Pair <- paste0(drop_files_dose_merge_pv05$cluster,drop_files_dose_merge_pv05$gene,sign(drop_files_dose_merge_pv05$avg_log2FC))
AllMarkers_Dose_adpv05_filter_DropCorrection$Pair <- paste0(AllMarkers_Dose_adpv05_filter_DropCorrection$cluster,AllMarkers_Dose_adpv05_filter_DropCorrection$gene,sign(AllMarkers_Dose_adpv05_filter_DropCorrection$avg_log2FC))
length(intersect(AllMarkers_Dose_adpv05_filter_DropCorrection$Pair,drop_files_dose_merge_pv05$Pair ))


# Check Fold Change 
AllMarkers_Norm_adpv05LogFC025_filter_DropCorrection <- CorrectEMTMarkers(AllMarkers_Norm_adpv05_filter,drop_files_all,20085,logFC=TRUE)
# [1] "Total EMT Markers: 709"
# [1] "Droped EMT Markers: 63"
# [1] "Remaiing EMT Markers: 646"
# [1] "Percent Remaining: 0.91114245416079"
AllMarkers_Time_adpv05LogFC025_filter_DropCorrection <- CorrectEMTMarkers(AllMarkers_Time_adpv05_filter,drop_files_time,18570,logFC=TRUE)
# [1] "Total EMT Markers: 757"
# [1] "Droped EMT Markers: 52"
# [1] "Remaiing EMT Markers: 705"
# [1] "Percent Remaining: 0.931307793923382"
AllMarkers_Dose_adpv05LogFC025_filter_DropCorrection <- CorrectEMTMarkers(AllMarkers_Dose_adpv05_filter,drop_files_dose,19186,logFC=TRUE)
# [1] "Total EMT Markers: 496"
# [1] "Droped EMT Markers: 62"
# [1] "Remaiing EMT Markers: 434"
# [1] "Percent Remaining: 0.875"

### Fix gene annotations before exporting ###

# Check Size
dim(AllMarkers_Norm_adpv05_filter_DropCorrection)
dim(AllMarkers_Time_adpv05_filter_DropCorrection)
dim(AllMarkers_Dose_adpv05_filter_DropCorrection)

# Update gene map
library(HGNChelper)
source("NonOverlappingReplacement.R")
currentmamp <- getCurrentHumanMap()

MarkerGenes_all <- AllMarkers_Norm_adpv05_filter_DropCorrection$gene
check_markers_all <- checkGeneSymbols(MarkerGenes_all,unmapped.as.na=FALSE,map=currentmamp)
check_markers_NoOverlap_all <- NonOverlappingReplacement(check_markers_all)
dim(check_markers_NoOverlap_all)
dim(AllMarkers_Norm_adpv05_filter_DropCorrection)
AllMarkers_Norm_adpv05_filter_DropCorrection$gene <- check_markers_NoOverlap_all$Suggested.Symbol
dim(AllMarkers_Norm_adpv05_filter_DropCorrection)

saveRDS(AllMarkers_Norm_adpv05_filter_DropCorrection,"AllMarkers_pv05_DropFiltered.rds")

MarkerGenes_time <- AllMarkers_Time_adpv05_filter_DropCorrection$gene
check_markers_time <- checkGeneSymbols(MarkerGenes_time,unmapped.as.na=FALSE,map=currentmamp)
check_markers_NoOverlap_time <- NonOverlappingReplacement(check_markers_time)
dim(check_markers_NoOverlap_time)
dim(AllMarkers_Time_adpv05_filter_DropCorrection)
AllMarkers_Time_adpv05_filter_DropCorrection$gene <- check_markers_NoOverlap_time$Suggested.Symbol
dim(AllMarkers_Time_adpv05_filter_DropCorrection)

saveRDS(AllMarkers_Time_adpv05_filter_DropCorrection,"TimeMarkers_pv05_DropFiltered.rds")

MarkerGenes_dose <- AllMarkers_Dose_adpv05_filter_DropCorrection$gene
check_markers_dose <- checkGeneSymbols(MarkerGenes_dose,unmapped.as.na=FALSE,map=currentmamp)
check_markers_NoOverlap_dose <- NonOverlappingReplacement(check_markers_dose)
dim(check_markers_NoOverlap_dose)
dim(AllMarkers_Dose_adpv05_filter_DropCorrection)
AllMarkers_Dose_adpv05_filter_DropCorrection$gene <- check_markers_NoOverlap_dose$Suggested.Symbol
dim(AllMarkers_Dose_adpv05_filter_DropCorrection)

saveRDS(AllMarkers_Dose_adpv05_filter_DropCorrection,"DoseMarkers_pv05_DropFiltered.rds")
