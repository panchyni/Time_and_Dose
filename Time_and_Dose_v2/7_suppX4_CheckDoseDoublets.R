########################
### Setup Enviroment ###
########################

### Folder
setwd("/Users/nicholaspanchy/Documents/Work_UTK/DoseTime_FromScratch2/")

# Dependencies
library(Seurat)
# remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
library(ggplot2)
library(scales)
library(HGNChelper)

dose_data.combined <- readRDS("SeuratData/Seurat_concentration_DoseOnly_FullCorr_hg38_v2CC.rds")

### Fix annotation ##

# Check against original list
check_table <- read.table("GBCannot/CBC_GBC_summary.txt.alt_names",head=TRUE)
check_table$CBC_alt <- cbind(lapply(check_table$CBC.1, function(x) {sub(".1", "-1", x)}))
row.names(check_table) <- check_table$CBC_alt
overlap <- intersect(check_table$CBC_alt,row.names(dose_data.combined@meta.data))
overlap <- unlist(overlap)

check_table <- as.data.frame(check_table)
check_table$meta_annot <- 0
check_table[overlap,]$meta_annot <- dose_data.combined@meta.data[overlap ,]$GBC

table(check_table$meta_annot,check_table$GBC)
# Note: 0s are samples present in the annotation, but not post filtering
# 149 samples in check table but not in processed data == 149 0s
# 4 samples that are GBC0 in the processed, but not the annotation (these were dropped in a previous filtering)

# Make replacements
#check_table[check_table$GBC == "GBC1" & check_table$meta_annot == "GBC0",]
dose_data.combined@meta.data["CTGAAACAGTTGTCGT-1",]$GBC = "GBC1"
dose_data.combined@meta.data["CTGAAACAGTTGTCGT-1",]$GBC_values = "1"
dose_data.combined@meta.data["CTGAAACAGTTGTCGT-1",]$GBC_pM = 0
dose_data.combined@meta.data["TGTATTCCATTTGCTT-1",]$GBC = "GBC1"
dose_data.combined@meta.data["TGTATTCCATTTGCTT-1",]$GBC_values = "1"
dose_data.combined@meta.data["TGTATTCCATTTGCTT-1",]$GBC_pM = 0
#check_table[check_table$GBC == "GBC4" & check_table$meta_annot == "GBC0",]
dose_data.combined@meta.data["ACGCCGAGTGGTACAG-1",]$GBC = "GBC4"
dose_data.combined@meta.data["ACGCCGAGTGGTACAG-1",]$GBC_values = "4"
dose_data.combined@meta.data["ACGCCGAGTGGTACAG-1",]$GBC_pM = 50
#check_table[check_table$GBC == "GBC9" & check_table$meta_annot == "GBC0",]
dose_data.combined@meta.data["CCGGGATCAGGCGATA-1",]$GBC = "GBC9"
dose_data.combined@meta.data["CCGGGATCAGGCGATA-1",]$GBC_values = "9"
dose_data.combined@meta.data["CCGGGATCAGGCGATA-1",]$GBC_pM = 200

check_table[overlap,]$meta_annot <- dose_data.combined@meta.data[overlap ,]$GBC
table(check_table$meta_annot,check_table$GBC)
# The are 4 are fixed and annotations line up except for 0s, as expected

# Check correspondence post correction
table(dose_data.combined@meta.data$GBC,dose_data.combined@meta.data$GBC_values)
table(dose_data.combined@meta.data$GBC,dose_data.combined@meta.data$GBC_pM)

# Check missing samples
check_table_0 <- check_table[check_table$meta_annot ==0,]
dim(check_table_0)
kazu_barcodes <- read.table("data_KazuConcentration/barcodes.tsv")
length(intersect(kazu_barcodes$V1,check_table_0$CBC_alt))
# All but ten are missing from the filtered barcode set
# The missing ten are filtered (these are the 10 NA sample post stringent filtering)


# Go through each dose treatment seperately, skip unlabeled
current_map <- getCurrentHumanMap()
for (i in c(1,2,3,4,8,9,12,13)){
  dose_sample <- row.names(dose_data.combined@meta.data[dose_data.combined@meta.data$GBC_values == i,])
  data.filt <- subset(dose_data.combined,cell=dose_sample)
  
  print(dim(data.filt))
  
  # Reprocess Individual Dose for speed
  g2m.genes <- cc.genes$g2m.genes
  s.genes <- cc.genes$s.genes
  
  g2m.check <- checkGeneSymbols(g2m.genes,map=current_map, unmapped.as.na=FALSE)
  s.check <- checkGeneSymbols(s.genes,map=current_map, unmapped.as.na=FALSE)
  
  g2m.genes <- g2m.check$Suggested.Symbol
  s.genes <- s.check$Suggested.Symbol
  
  print(length(intersect(g2m.genes,row.names(data.filt@assays$RNA@data))))
  print(length(intersect(s.genes,row.names(data.filt@assays$RNA@data))))
  
  data.filt = NormalizeData(data.filt)
  data.filt <- CellCycleScoring(object = data.filt, g2m.features = g2m.genes,
                                s.features = s.genes)
  
  suppressMessages(require(DoubletFinder))
  
  data.filt = FindVariableFeatures(data.filt, verbose = F)
  data.filt = ScaleData(data.filt, vars.to.regress = c("nCount_RNA","percent.mito","S.Score","G2M.Score"),
                        verbose = F)
  data.filt = RunPCA(data.filt, verbose = F, npcs = 15)
  data.filt = RunUMAP(data.filt, dims = 1:15, verbose = F)
  
  ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
  sweep.res.list<- paramSweep_v3(data.filt, PCs = 1:15, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn_ <- find.pK(sweep.stats)
  optimal_pk <- as.numeric(as.character(bcmvn_[bcmvn_$BCmetric==max(bcmvn_$BCmetric),]$pK))
  
  # Find Doublets
  annot_freq = length(dose_sample)/dim(dose_data.combined@meta.data[dose_data.combined@meta.data$GBC_values > 0,])[1]
  print(annot_freq)
  nExp <- round(ncol(data.filt) * 0.08 * annot_freq)  # expect 8.0% doublets for 10k recovery * annotation freq in total set for homotypic
  data.filt <- doubletFinder_v3(data.filt, pN = 0.25, pK = optimal_pk, nExp = nExp, PCs = 1:15)
  
  # name of the DF prediction can change, so extract the correct column name.
  DF.name = colnames(data.filt@meta.data)[grepl("DF.classification", colnames(data.filt@meta.data))]
  cowplot::plot_grid(ncol = 2, DimPlot(data.filt, group.by = "orig.ident") + NoAxes(),
                     DimPlot(data.filt, group.by = DF.name) + NoAxes())
  saveRDS(data.filt@meta.data,paste0("Dose_GBCvalues_",i,"_Doublets.rds"))
}

##############################################################################
### Visualize Doublets against fully processed data (Requires 1x Scripts) ####
##############################################################################

Dose_GBC_1_Doublets <- readRDS("Dose_GBCvalues_1_Doublets.rds")
Dose_GBC_2_Doublets <- readRDS("Dose_GBCvalues_2_Doublets.rds")
Dose_GBC_3_Doublets <- readRDS("Dose_GBCvalues_3_Doublets.rds")
Dose_GBC_4_Doublets <- readRDS("Dose_GBCvalues_4_Doublets.rds")
Dose_GBC_8_Doublets <- readRDS("Dose_GBCvalues_8_Doublets.rds")
Dose_GBC_9_Doublets <- readRDS("Dose_GBCvalues_9_Doublets.rds")
Dose_GBC_12_Doublets <- readRDS("Dose_GBCvalues_12_Doublets.rds")
Dose_GBC_13_Doublets <- readRDS("Dose_GBCvalues_13_Doublets.rds")

dose_meta <- readRDS("Dose_MetaData.rds")
UMAP_values <- readRDS("Dose_UMAP_values.rds")

# Add doublets to met data
dose_meta$singlet <- "unlabeled"
dose_meta[row.names(Dose_GBC_1_Doublets),]$singlet <- Dose_GBC_1_Doublets$DF.classifications_0.25_0.05_16
dose_meta[row.names(Dose_GBC_2_Doublets),]$singlet <- Dose_GBC_2_Doublets$DF.classifications_0.25_0.19_18
dose_meta[row.names(Dose_GBC_3_Doublets),]$singlet <- Dose_GBC_3_Doublets$DF.classifications_0.25_0.01_9
dose_meta[row.names(Dose_GBC_4_Doublets),]$singlet <- Dose_GBC_4_Doublets$DF.classifications_0.25_0.1_9
dose_meta[row.names(Dose_GBC_8_Doublets),]$singlet <- Dose_GBC_8_Doublets$DF.classifications_0.25_0.02_12
dose_meta[row.names(Dose_GBC_9_Doublets),]$singlet <- Dose_GBC_9_Doublets$DF.classifications_0.25_0.13_9
dose_meta[row.names(Dose_GBC_12_Doublets),]$singlet <- Dose_GBC_12_Doublets$DF.classifications_0.25_0.27_7
dose_meta[row.names(Dose_GBC_13_Doublets),]$singlet <- Dose_GBC_13_Doublets$DF.classifications_0.25_0.05_10
table(dose_meta$singlet)

# Plot Doublets on the UMAP
UMAP_annot <- merge(UMAP_values,dose_meta,by=0)
row.names(UMAP_annot) <- UMAP_annot$Row.names
UMAP_annot$GBC_pM_factor <- factor(UMAP_annot$GBC_pM)
UMAP_annot$Time_factor <- factor(UMAP_annot$Time)

ggplot(UMAP_annot, aes(x = UMAP_1, y = UMAP_2, colour=singlet)) + geom_point(size=0.5)
UMAP_annot$clusters_and_doublets <- as.character(UMAP_annot$seurat_clusters)
UMAP_annot[UMAP_annot$singlet=='Doublet',]$clusters_and_doublets <- "Doublet"

manual_colors <- hue_pal()(7)
manual_colors <- c(manual_colors,"#000000")
ggplot(UMAP_annot, aes(x = UMAP_1, y = UMAP_2, colour=clusters_and_doublets)) + geom_point(size=0.5) +
  scale_colour_manual(values=manual_colors) + theme_bw()

### Load Integrated Data ###

integated_meta <- readRDS("Integrated_MetaData.rds")
UMAP_values <- readRDS("Integrated_UMAP_values.rds")

dose_samples <- readRDS("DoseSamples.rds")
time_samples <- readRDS("TimeSamples.rds")

integated_meta <- integated_meta[dose_samples,]

# Add Doublet data to meta data
integated_meta$singlet <- "unlabeled"
integated_meta[row.names(Dose_GBC_1_Doublets),]$singlet <- Dose_GBC_1_Doublets$DF.classifications_0.25_0.05_16
integated_meta[row.names(Dose_GBC_2_Doublets),]$singlet <- Dose_GBC_2_Doublets$DF.classifications_0.25_0.19_18
integated_meta[row.names(Dose_GBC_3_Doublets),]$singlet <- Dose_GBC_3_Doublets$DF.classifications_0.25_0.01_9
integated_meta[row.names(Dose_GBC_4_Doublets),]$singlet <- Dose_GBC_4_Doublets$DF.classifications_0.25_0.1_9
integated_meta[row.names(Dose_GBC_8_Doublets),]$singlet <- Dose_GBC_8_Doublets$DF.classifications_0.25_0.02_12
integated_meta[row.names(Dose_GBC_9_Doublets),]$singlet <- Dose_GBC_9_Doublets$DF.classifications_0.25_0.13_9
integated_meta[row.names(Dose_GBC_12_Doublets),]$singlet <- Dose_GBC_12_Doublets$DF.classifications_0.25_0.27_7
integated_meta[row.names(Dose_GBC_13_Doublets),]$singlet <- Dose_GBC_13_Doublets$DF.classifications_0.25_0.05_10

# Plot Doublets on the UMAP
UMAP_annot <- merge(UMAP_values,integated_meta,by=0)
row.names(UMAP_annot) <- UMAP_annot$Row.names
UMAP_annot$GBC_pM_factor <- factor(UMAP_annot$GBC_pM)
UMAP_annot$Time_factor <- factor(UMAP_annot$Time)

ggplot(UMAP_annot, aes(x = UMAP_1, y = UMAP_2, colour=singlet)) + geom_point(size=0.5)
UMAP_annot$clusters_and_doublets <- as.character(UMAP_annot$seurat_clusters)
UMAP_annot[UMAP_annot$singlet=='Doublet',]$clusters_and_doublets <- "Doublet"

manual_colors <- hue_pal()(15)
manual_colors <- c(manual_colors,"#000000")
ggplot(UMAP_annot, aes(x = UMAP_1, y = UMAP_2, colour=clusters_and_doublets)) + geom_point(size=0.5) +
  scale_colour_manual(values=manual_colors) + geom_point(size=0.25) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size = 10),axis.title=element_text(size=11),aspect.ratio=1,legend.position = "bottom")
