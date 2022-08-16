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

###########################
### Load Processed Data ###
###########################

integated_meta <- readRDS("Integrated_MetaData.rds")
UMAP_values <- readRDS("Integrated_UMAP_values.rds")

dose_samples <- readRDS("DoseSamples.rds")
time_samples <- readRDS("TimeSamples.rds")

time_doublets <- readRDS("time_doublets.rds")
time_wo_doublets <- setdiff(time_samples,time_doublets)
time_samples <- time_wo_doublets

Integrated_E_mat_nsprcomp <- readRDS("Integrated_E_mat_nsprcomp.rds")
Integrated_M_mat_nsprcomp <- readRDS("Integrated_M_mat_nsprcomp.rds")

####################
### Bubble Plots ###
####################

### Bubble Plots ###

GetBoundarySamples <- function(Scores,row,column,bonds_df){
  
  # Select samples based on score
  if (row == 1){ # If first values, get less than
    samples <- names(Scores[Scores < bonds_df[row,column]])
  }
  else if (row == dim(bonds_df)[1]+1){ # if beyond the last value, > last value
    samples <- names(Scores[Scores >= bonds_df[row-1,column]])
  }
  else { # else between current value and previous value
    tmp <- Scores[Scores >= bonds_df[row-1,column]]
    samples <- names(tmp[tmp< bonds_df[row,column]])
  }
  return(samples)
}

EMCount_Heatmap <- function(Escores,Mscores,bonds_df){
  #plot(overlapping_Escores$PC1, overlapping_Mscores$PC1,pch=16,cex=0.5)
  plot(c(min(Escores),max(Escores)),c(min(Mscores),max(Mscores)))
  manual_colors <- hue_pal()(16)
  
  out_matrix = matrix(, nrow = 16, ncol = 5)
  index = 1
  for(i in 1:(dim(bonds_df)[1]+1)){
    for(j in 1:(dim(bonds_df)[1]+1)){
      E_boundary_samples <- GetBoundarySamples(Escores,i,1,bonds_df)
      M_boundary_sampels <- GetBoundarySamples(Mscores,j,2,bonds_df)
      intersect_samples = intersect(E_boundary_samples, M_boundary_sampels)
      
      # Basic Counts
      #print(i)
      #print(j)
      #print(length(intersect_samples))
      
      # Inspect boundaries by ploting
      points(Escores[intersect_samples],Mscores[intersect_samples],col=manual_colors[i+(j-1)*4],cex=0.25)
      
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
      
      # This is the old code where I made 'i' columns and 'j' rows
      # for quick visulization via matrix but it drives me crazy
      # try to rectify this with X,Y coordinates
      #odds_mat[j,i] = test$estimate
      #pv_mat[j,i] = test$p.value
      #counts_mat[j,i] = cluster_time_samples + cluster_con_samples
      
      # New code to ouput a longformat with all of the data
      out_matrix[index,] <- c(i,j,test$estimate,test$p.value,cluster_time_samples + cluster_con_samples)
      index = index + 1
    }
  }
  colnames(out_matrix) <- c("i","j","Odds","pv","Count")
  
  return(out_matrix)
}

### Bubble Plot Sample Sets ###

# Define our sample sets
Con_Samples <- dose_samples
Time_Samples <- time_samples
ConAll_samples <- rownames(integated_meta[integated_meta$Batch == "C",])
length(Time_Samples)
length(Con_Samples)
length(ConAll_samples)

# Get all samples with a label
Labeled_Samples <- c(Time_Samples,Con_Samples)
length(Labeled_Samples)

# Set PC rank, from previous code for looking at different combo
m = 1
n = 1

E_scores <- Integrated_E_mat_nsprcomp$x[,m]
M_scores <- Integrated_M_mat_nsprcomp$x[,n]

Escore_IQR_med <- quantile(E_scores[Labeled_Samples],probs=c(0.25,0.50,0.75))
Mscore_IQR_med <- quantile(M_scores[Labeled_Samples],probs=c(0.25,0.50,0.75))
Bounds_IQR_med <- as.data.frame(cbind(Escore_IQR_med,Mscore_IQR_med))
print(Bounds_IQR_med)

### Bubble Plot EM-scores ###
results <- EMCount_Heatmap(E_scores,M_scores,Bounds_IQR_med)

# Make data frame and define Odds and X/Y coordinates
odd_mat_long <- as.data.frame(results)
odd_mat_long$Odds <- log2(odd_mat_long$Odds)
odd_mat_long$X <- odd_mat_long$i
odd_mat_long$Y <- odd_mat_long$j

# Make Bubble Plot
bubble_p <- ggplot(odd_mat_long, aes(x = X, y = Y,size=Count,fill=Odds)) + 
  geom_point(alpha=1.0, shape=21, color="black") +
  scale_size(range = c(5, 15)) +
  scale_fill_gradient2(low='blue',high='red',mid='white',midpoint=0,limits=c(-3,3),oob=squish) +
  theme_bw() + coord_fixed(ratio = 1) + ylim(c(0.75,4.25)) + xlim(c(0.75,4.25))
bubble_p

# Check sample location
check_p <- ggplot(odd_mat_long, aes(x = X, y = Y,size=Count,fill=Odds,label=Odds)) + 
  geom_point(alpha=1.0, shape=21, color="black") +
  scale_size(range = c(5, 15)) +
  scale_fill_gradient2(low='blue',high='red',mid='white',midpoint=0,limits=c(-3,3),oob=squish) +
  theme_bw() + coord_fixed(ratio = 1) + ylim(c(0.75,4.25)) + xlim(c(0.75,4.25)) +geom_text(hjust=0, vjust=0,size=2)
check_p

#########################
### Quartile Clusters ###
#########################

# Get samples for each of the quartile clusters
quart_clusters <- rep(NA, 16)
index = 1  
for(i in 1:(dim(Bounds_IQR_med)[1]+1)){
  for(j in 1:(dim(Bounds_IQR_med)[1]+1)){
    # Get samples in E and M quartiles
    E_boundary_samples <- GetBoundarySamples(E_scores,i,1,Bounds_IQR_med)
    M_boundary_sampels <- GetBoundarySamples(M_scores,j,2,Bounds_IQR_med)
    # Get intersection of samples
    intersect_samples = intersect(E_boundary_samples, M_boundary_sampels)
    # Filter by labeled samples
    filter_samples <- intersect(intersect_samples,union(Time_Samples,Con_Samples))
    quart_clusters[index] <- list(filter_samples)
    index = index + 1
  }
}

### Check uniqueness and completeness
check_mat <- matrix(,nrow=length(Labeled_Samples),ncol=16)
n= 0
for (sample in Labeled_Samples){
  n=n+1
  for (i in 1:16){
    if(sample %in% quart_clusters[[i]]){
      check_mat[n,i] = 1
    }
    else{
      check_mat[n,i] = 0
    }
  }
}
sum(check_mat)               # Check that sum = # of sampes
summary(rowSums(check_mat))  # Check that row sums are all 1

# Save Results for Reuse
saveRDS(Bounds_IQR_med,"Bounds_IQR_med.rds")
saveRDS(odd_mat_long,"BubbleOddsMat.rds")
saveRDS(quart_clusters,"EMQuartileClusters.rds")