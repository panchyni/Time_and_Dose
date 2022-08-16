########################
### Setup Enviroment ###
########################

### Folder
setwd("/Users/nicholaspanchy/Documents/Work_UTK/DoseTime_FullCorr_FromScratch/")

# Dependencies
library(Seurat)
# remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
library(ggplot2)
library(scales)

# Read Time data
time_data.combined <- readRDS("SeuratData/Seurat_integrated_TimeOnly_FullCorr_hg38.rds")

# Make a batch X Time variables to split data
time_data.combined$TimeBatch <- paste0(time_data.combined$Time,time_data.combined$Batch)
TimeIDs <- names(table(time_data.combined$TimeBatch))

# For each time data set
for (i in TimeIDs){
  
  # Subset set the data
  time_sample <- row.names(time_data.combined@meta.data[time_data.combined@meta.data$TimeBatch == i,])
  data.filt <- subset(time_data.combined,cell=time_sample)
  
  # Change to RNA (non-integrated) data for Doublet detection
  DefaultAssay(data.filt) <- "RNA"
  
  # Check Size
  print(dim(data.filt))
  
  # Reprocess Individual Dose for speed
  data.filt = NormalizeData(data.filt)
  data.filt <- CellCycleScoring(object = data.filt, g2m.features = cc.genes$g2m.genes,
                                s.features = cc.genes$s.genes)
  
  suppressMessages(require(DoubletFinder))
  
  # Take top 2000 variable genes for doublet ID
  data.filt = FindVariableFeatures(data.filt, verbose = F)
  data.filt = ScaleData(data.filt, vars.to.regress = c("nCount_RNA","percent.mito","S.Score","G2M.Score"),
                        verbose = F)
  data.filt = RunPCA(data.filt, verbose = F, npcs = 15)
  data.filt = RunUMAP(data.filt, dims = 1:15, verbose = F)
  
  ## pK Identification (no ground-truth)
  sweep.res.list<- paramSweep_v3(data.filt, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn_ <- find.pK(sweep.stats)
  optimal_pk <- as.numeric(as.character(bcmvn_[bcmvn_$BCmetric==max(bcmvn_$BCmetric),]$pK))
  
  # Find Doublets
  nExp <- round(ncol(data.filt) * 0.024)  # expect 2.4% doublets for 3K cell target
  data.filt <- doubletFinder_v3(data.filt, pN = 0.25, pK = optimal_pk, nExp = nExp, PCs = 1:10)
  
  # name of the DF prediction can change, so extract the correct column name.
  DF.name = colnames(data.filt@meta.data)[grepl("DF.classification", colnames(data.filt@meta.data))]
  cowplot::plot_grid(ncol = 2, DimPlot(data.filt, group.by = "orig.ident") + NoAxes(),
                     DimPlot(data.filt, group.by = DF.name) + NoAxes())
  saveRDS(data.filt@meta.data,paste0("Time_Batch_",i,"_Doublets.rds"))
}

##############################################################################
### Visualize Doublets against fully processed data (Requires 1x Scripts) ####
##############################################################################

# Load prepared data
Time_0A_Doublets <- readRDS("Time_Batch_0A_Doublets.rds")
Time_0B_Doublets <- readRDS("Time_Batch_0B_Doublets.rds")
Time_1B_Doublets <- readRDS("Time_Batch_1B_Doublets.rds")
Time_2B_Doublets <- readRDS("Time_Batch_2B_Doublets.rds")
Time_3B_Doublets <- readRDS("Time_Batch_3B_Doublets.rds")
Time_4A_Doublets <- readRDS("Time_Batch_4A_Doublets.rds")
Time_8A_Doublets <- readRDS("Time_Batch_8A_Doublets.rds")

# Load Meta and UMAP
integated_meta <- readRDS("Integrated_MetaData.rds")
UMAP_values <- readRDS("Integrated_UMAP_values.rds")

dose_samples <- readRDS("DoseSamples.rds")
time_samples <- readRDS("TimeSamples.rds")

integated_meta <- integated_meta[time_samples,]

# Add doublet info to meta data
integated_meta$singlet <- "unlabeled"
integated_meta[row.names(Time_0A_Doublets),]$singlet <- Time_0A_Doublets$DF.classifications_0.25_0.05_32
integated_meta[row.names(Time_0B_Doublets),]$singlet <- Time_0B_Doublets$DF.classifications_0.25_0.14_65
integated_meta[row.names(Time_1B_Doublets),]$singlet <- Time_1B_Doublets$DF.classifications_0.25_0.12_55
integated_meta[row.names(Time_2B_Doublets),]$singlet <- Time_2B_Doublets$DF.classifications_0.25_0.08_57
integated_meta[row.names(Time_3B_Doublets),]$singlet <- Time_3B_Doublets$DF.classifications_0.25_0.22_51
integated_meta[row.names(Time_4A_Doublets),]$singlet <- Time_4A_Doublets$DF.classifications_0.25_0.05_27
integated_meta[row.names(Time_8A_Doublets),]$singlet <- Time_8A_Doublets$DF.classifications_0.25_0.08_45

time_doublets <- row.names(integated_meta[integated_meta$singlet=="Doublet",])
saveRDS(time_doublets,"time_doublets.rds")

# Plot Doublets on the UMAP
UMAP_annot <- merge(UMAP_values,integated_meta,by=0)
row.names(UMAP_annot) <- UMAP_annot$Row.names
UMAP_annot$GBC_pM_factor <- factor(UMAP_annot$GBC_pM)
UMAP_annot$Time_factor <- factor(UMAP_annot$Time)

ggplot(UMAP_annot, aes(x = UMAP_1, y = UMAP_2, colour=singlet)) + geom_point(size=0.5)
UMAP_annot$clusters_and_doublets <- as.character(UMAP_annot$seurat_clusters)
UMAP_annot[UMAP_annot$singlet=='Doublet',]$clusters_and_doublets <- "Doublet"

manual_colors <- hue_pal()(13)
manual_colors <- c(manual_colors,"#000000")
ggplot(UMAP_annot, aes(x = UMAP_1, y = UMAP_2, colour=clusters_and_doublets)) + geom_point(size=0.5) +
  scale_colour_manual(values=manual_colors) + geom_point(size=0.25) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1,legend.position = "bottom")
