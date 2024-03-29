# Folder
setwd("/Users/nicholaspanchy/Documents/Work_UTK/DoseTime_FromScratch2/CoxModels_Updated/")

#Name: dataset_code.R
#Author: Andrew Willems <awillems@vols.utk.edu>
#Purpose: Concordance index (C-index) analysis for Nick's clusters of genes
#Load needed packages----
library(survival)
library(survminer)
library(tidyverse)
library(Seurat)

filter_TCGA_files = c('ACC_df_fpkm_update.RData',
                      'BLCA_df_fpkm_update.RData',
                      'BRCA_df_fpkm_update.RData',
                      'CESC_df_fpkm_update.RData',
                      'CHOL_df_fpkm_update.RData',
                      'COAD_df_fpkm_update.RData',
                      'ESCA_df_fpkm_update.RData',
                      'GBM_df_fpkm_update.RData',   
                      'HNSC_df_fpkm_update.RData',  
                      'KIRC_df_fpkm_update.RData',
                      'KIRP_df_fpkm_update.RData',
                      'LAML_df_fpkm_update.RData', 
                      'LGG_df_fpkm_update.RData',
                      'LIHC_df_fpkm_update.RData',
                      'LUAD_df_fpkm_update.RData',
                      'LUSC_df_fpkm_update.RData',
                      'MESO_df_fpkm_update.RData',
                      'OV_df_fpkm_update.RData',
                      'PAAD_df_fpkm_update.RData',
                      'READ_df_fpkm_update.RData',
                      'SARC_df_fpkm_update.RData',
                      'SKCM_df_fpkm_update.RData',
                      'STAD_df_fpkm_update.RData',
                      'THCA_df_fpkm_update.RData',
                      'UCEC_df_fpkm_update.RData',
                      'UCS_df_fpkm_update.RData',
                      'UVM_df_fpkm_update.RData')

#Loading the cox_model_fitter() function from its file----
source("Code/cox_model.R")
source("Code/tcga_replicateFilter.R")

#Loading Nick's data and removing any data frames that have no genes----
clust_int <- readRDS(file = "Data/Clusters_IntegrationFirst_Quart_3rdQuartFC_Alt.rds")

# Update gene map
library(HGNChelper)
source("../NonOverlappingReplacement.R")
currentmamp <- getCurrentHumanMap()

# Get all genes
all_genes <- c()
for (df in clust_int){
  current_df <- df
  all_genes <- c(all_genes,current_df$gene)
}
all_genes <- unique(all_genes)

# Store variables
result_lists <- vector(mode = "list", length = 50)
active_genes <- vector(mode = "list", length = 50)
gene_lists <- vector(mode = "list", length = 50)

# Make > 1 for more than one replicate
for (n in 1:1){
  results_mat = matrix(, nrow = length(filter_TCGA_files), ncol = 21)
  activegenes_mat = matrix(, nrow = length(filter_TCGA_files), ncol = 19)
  activegenes_list <- vector(mode = "list",length=length(filter_TCGA_files)*15)
  i_list <- vector(mode = "list",length=length(filter_TCGA_files)*15)
  j_list <- vector(mode = "list",length=length(filter_TCGA_files)*15)
  list_index=1
  
  #for (i in 1:3){
  for (i in 1:length(filter_TCGA_files)){
    print(filter_TCGA_files[i])
    
    #Loading the cancer specific data frame----
    load(paste("Data/",filter_TCGA_files[i],sep=""),
         verbose = TRUE)
    #print(dim(cox_df[1]))
    
    # Fix genes
    check_cox <- checkGeneSymbols(colnames(cox_df),unmapped.as.na = FALSE,map=currentmamp)
    check_cox_NoOverlap <- NonOverlappingReplacement(check_cox)
    colnames(cox_df) <- check_cox_NoOverlap$Suggested.Symbol
    
    # Filter Samples
    dim(cox_df)
    cox_df <- cox_df[ , !(names(cox_df) %in% c("ajcc.n","ajcc.m"))]
    filter_samples <- tcga_replicateFilter(row.names(cox_df),"RNA")
    cox_df <- cox_df[filter_samples,]
    cox_df <- cox_df[cox_df$days.to.last.follow.up > 0,]
    cox_df <- cox_df[cox_df$sample.type != 'Solid Tissue Normal',]
    dim(cox_df)
    
    # Convert to TPM
    expression_data <- subset(cox_df, select=c(TSPAN6:AC007389.5))
    expression_data_TPM <- t(apply(expression_data,1, function(x) 1000000*(x/sum(x))))
    expression_data_TPM_LogP1 <- log1p(expression_data_TPM)
    #check
    #summary(rowSums(expression_data))
    #summary(rowSums(expression_data_TPM))
    cox_df[, colnames(expression_data_TPM_LogP1)] <- expression_data_TPM_LogP1
    
    #print(table(cox_df$vital.status))
    
    # Save sample number and size
    results_mat[i,1] = i
    results_mat[i,2] = dim(cox_df)[1]
    
    activegenes_mat[i,1] = i
    activegenes_mat[i,2] = dim(cox_df)[1]
    
    # Baseline model using all genes, we don't really use this
    cox_model <- cox_model_fitter(my.seed = as.numeric(n), cox.df = cox_df,
                                  gene.num = length(all_genes),
                                  cox.predictors = all_genes,
                                  tumor.stage = FALSE,
                                  tumor.n = FALSE,
                                  tumor.m = FALSE,
                                  regular.cox = FALSE,
                                  save.regular.cox.genes = FALSE,
                                  # Output from previous version of code, not used here so dump everything to the same file
                                  # and don't worry about
                                  my.filename = paste0("Data/Outputs/tmpAll_clust_int_nod_read_222022.csv")
    )
    #max(cox_model$CV$cvm)
    
    # Save Performance and # of active genes
    results_mat[i,3] = max(cox_model$CV$cvm)
    activegenes_mat[i,3] = length(cox_model$`Active Genes`)
    
    #For loop to get all of the gene names from each data frame and submit them
    #to the cox_model_fitter() function to get concordance index----
    counter <- 1
    outputs <- list()
    my_cindicies <- c()
    my_activegenes <- c()
    
    j = 1
    for (df in clust_int){
      current_df <- df
      if(dim(current_df)[1] != 0){
        
        # Run the model for each cluster
        cox_model <- cox_model_fitter(my.seed = as.numeric(n), cox.df = cox_df,
                                      gene.num = length(rownames(current_df)),
                                      cox.predictors = current_df$gene,
                                      tumor.stage = FALSE,
                                      tumor.n = FALSE,
                                      tumor.m = FALSE,
                                      regular.cox = FALSE,
                                      save.regular.cox.genes = FALSE,
                                      # Output from previous version of code, not used here so dump everything to the same file
                                      # and don't worry about
                                      my.filename = paste0("Data/Outputs/",counter,"_clust_integrated_all_read.csv"))
        
        outputs[[counter]] <- cox_model
        perc_done <- paste0(round((100*counter)/length(clust_int)), "% done with this dataset")
        print(perc_done)
        counter <- counter + 1
        
        #Storing all of the c-index values in a vector that we can use later
        c_finder <-cox_model$CV$index[1]
        current_c <- cox_model$CV$cvm[c_finder]
        current_c <- round(current_c, digits = 4)
        my_cindicies <- c(my_cindicies, current_c)
        my_activegenes <- c(my_activegenes,length(cox_model$`Active Genes`))
        
        # Save actual actives genes + sample and cluster index
        activegenes_list[list_index] <- list(cox_model$`Active Genes`)
        i_list[list_index] <- i
        j_list[list_index] <- j
        list_index=list_index+1
        
      }else{
        print("This input doesn't have any genes and the cox model is skipping to the next data frame")
        perc_done <- paste0(round((100*counter)/length(clust_int)), "% done with this dataset")
        print(perc_done)
        counter <- counter + 1
        
      }
      
      #Getting the top performing data frame 
      top_cindex <-max(my_cindicies)
      top_index <- which(my_cindicies==top_cindex)
      j=j+1
    }
    
    # Save results from the loop for c-index
    results_mat[i,4] = my_cindicies[1]
    results_mat[i,5] = my_cindicies[2]
    results_mat[i,6] = my_cindicies[3]
    results_mat[i,7] = my_cindicies[4]
    results_mat[i,8] = my_cindicies[5]
    results_mat[i,9] = my_cindicies[6]
    results_mat[i,10] = my_cindicies[7]
    results_mat[i,11] = my_cindicies[8]
    results_mat[i,12] = my_cindicies[9]
    results_mat[i,13] = my_cindicies[10]
    results_mat[i,14] = my_cindicies[11]
    results_mat[i,15] = my_cindicies[12]
    results_mat[i,16] = my_cindicies[13]
    results_mat[i,17] = my_cindicies[14]
    results_mat[i,18] = my_cindicies[15]
    results_mat[i,19] = my_cindicies[16]
    
    # Save results from the loop for active gene count
    activegenes_mat[i,4] = my_activegenes[1]
    activegenes_mat[i,5] = my_activegenes[2]
    activegenes_mat[i,6] = my_activegenes[3]
    activegenes_mat[i,7] = my_activegenes[4]
    activegenes_mat[i,8] = my_activegenes[5]
    activegenes_mat[i,9] = my_activegenes[6]
    activegenes_mat[i,10] = my_activegenes[7]
    activegenes_mat[i,11] = my_activegenes[8]
    activegenes_mat[i,12] = my_activegenes[9]
    activegenes_mat[i,13] = my_activegenes[10]
    activegenes_mat[i,14] = my_activegenes[11]
    activegenes_mat[i,15] = my_activegenes[12]
    activegenes_mat[i,16] = my_activegenes[13]
    activegenes_mat[i,17] = my_activegenes[14]
    activegenes_mat[i,18] = my_activegenes[15]
    activegenes_mat[i,19] = my_activegenes[16]
    
    # Record patient classes
    results_mat[i,20] = table(cox_df$vital.status)[[1]]
    results_mat[i,21] = table(cox_df$vital.status)[[2]]
    
  }
  
  # Store results in list, only use for > 1 replicate
  result_lists[[n]] <- results_mat
  active_genes[[n]] <- activegenes_mat
}

results_mat
activegenes_mat

save(results_mat, file = "MultiCancerResultsTable_QuartALL_LogP1_Alt.RData")
save(activegenes_mat, file = "MultiCancerGeneCountTable_QuartAL_LogP1_AltL.RData")

save(activegenes_list, file = "MultiCancerGeneLists_QuartALL_LogP1_Alt.RData")
save(i_list, file = "MultiCancerGeneilist_QuartAL_LogP1L_Alt.RData")
save(j_list, file = "MultiCancerjlist_QuartALL_LogP1_Alt.RData")

#load("MultiCancerGeneLists_QuartALL_Alt.RData")
#load("MultiCancerGeneilist_QuartALL_Alt.RData")
#load("MultiCancerjlist_QuartALL_Alt.RData")

### Process Results into Figures ###

library(ggplot2)
library(scales)

names = c('ACC',
                      'BLCA',
                      'BRCA',
                      'CESC',
                      'CHOL',
                      'COAD',
                      'ESCA',
                      'GBM',   
                      'HNSC',  
                      'KIRC',
                      'KIRP',
                      'LAML', 
                      'LGG',
                      'LIHC',
                      'LUAD',
                      'LUSC',
                      'MESO',
                      'OV',
                      'PAAD',
                      'READ',
                      'SARC',
                      'SKCM',
                      'STAD',
                      'THCA',
                      'UCEC',
                      'UCS',
                      'UVM')

load("MultiCancerResultsTable_QuartALL_LogP1_Alt.RData")
load("MultiCancerGeneCountTable_QuartAL_LogP1_AltL.RData")

result_df <- as.data.frame(results_mat)
activegenes_mat_df <- as.data.frame(activegenes_mat)
result_df$Names <- names
activegenes_mat_df$Names <- names

results_mat_quad <- matrix(, nrow = 16, ncol =5)
for (i in 4:19){
  cluster_num <- i-3
  x_position <- ceiling(cluster_num/4)
  y_position <- cluster_num - 4*(x_position-1)
  cluster_cindex <- mean(result_df[,i])
  cluster_cindex_GT10 <- mean(result_df[result_df$V21 > 10,i])
  ActiveGenes <- mean(activegenes_mat[,i])
  
  results_mat_quad[cluster_num,1] <- x_position
  results_mat_quad[cluster_num,2] <- y_position
  results_mat_quad[cluster_num,3] <- cluster_cindex
  results_mat_quad[cluster_num,4] <- cluster_cindex_GT10
  results_mat_quad[cluster_num,5] <- ActiveGenes
}
results_mat_quad_df <- as.data.frame(results_mat_quad)

#ggplot() + geom_point(data = results_mat_quad_df, aes(x = V1, y = V2,fill=V3), color='black',pch=21, size=8) +
#  scale_fill_gradient2(low='red',high='blue',mid='white',midpoint=0.6,limits=c(0.5,0.7),oob=squish) +
#  coord_fixed(ratio = 1)

ggplot() + geom_point(data = results_mat_quad_df, aes(x = V1, y = V2,fill=V4,size=V5), color='black',pch=21) +
  scale_fill_gradient2(low='black',high='red',mid='white',midpoint=0.645,limits=c(0.57,0.70),oob=squish) +
  coord_fixed(ratio = 1) + scale_size(range = c(3, 12)) + theme_bw() + ylim(c(0.75,4.25)) + xlim(c(0.75,4.25))

######################
### Scaled C-index ###
######################

### Data ###
result_df_scale <- result_df
result_df_scale[c(4:19)] <- t(scale(t(result_df[c(4:19)])))
result_df_scale$AvgClusterIndex <- rowMeans(result_df_scale[,c(4:19)])

### Heatmap ###

results_mat_long_score_order <- matrix(, nrow = 448, ncol =10) # for adding averages
index= 1
for (i in 1:27){
  for (j in 4:19){
    results_mat_long_score_order[index,1] <- j
    results_mat_long_score_order[index,2] <- i
    #results_mat_long_score_order[index,3] <- result_df[i,j] # Raw Cindex
    results_mat_long_score_order[index,3] <- result_df_scale[i,j] # Scaled Cindex
    results_mat_long_score_order[index,4] <- result_df_scale$V2[i]
    results_mat_long_score_order[index,5] <- mean(result_df_scale[result_df_scale$V21 > 10,j])
    results_mat_long_score_order[index,6] <- j-3
    results_mat_long_score_order[index,7] <- result_df_scale$Names[i]
    
    x = ceiling((j-3)/4)
    y = (j-3) %% 4
    if (y==0){
      y=4
    }
    results_mat_long_score_order[index,8] <- paste("(",x,',',y,")",sep="")
    
    index = index + 1
  }
}

# Add averages #
for (k in 1:16){
  
  results_mat_long_score_order[index,] <- results_mat_long_score_order[k,]
  
  results_mat_long_score_order[index,3] <- results_mat_long_score_order[index,5]
  results_mat_long_score_order[index,4] <- -1
  results_mat_long_score_order[index,7] <- "Mean"
  
  index = index + 1
}

results_mat_long_score_order_df <- as.data.frame(results_mat_long_score_order)
results_mat_long_score_order_df[,1] <- as.numeric(results_mat_long_score_order_df[,1])
results_mat_long_score_order_df[,2] <- as.numeric(results_mat_long_score_order_df[,2])
results_mat_long_score_order_df[,3] <- as.numeric(results_mat_long_score_order_df[,3])
results_mat_long_score_order_df[,4] <- as.numeric(results_mat_long_score_order_df[,4])
results_mat_long_score_order_df[,5] <- as.numeric(results_mat_long_score_order_df[,5])
results_mat_long_score_order_df[,6] <- as.numeric(results_mat_long_score_order_df[,6])

colnames(results_mat_long_score_order_df) <- c("X","Y","Cindex","Count","AverageScore", "Cluster","Names", "Model", "Score_Rank","Count_Rank")

results_mat_long_score_order_df$AverageScore <- round(results_mat_long_score_order_df$AverageScore,digits=4)
table(results_mat_long_score_order_df$AverageScore)

results_mat_long_score_order_df$Score_Rank <- rank(results_mat_long_score_order_df$AverageScore,ties.method = "min")
results_mat_long_score_order_df$Count_Rank <- rank(results_mat_long_score_order_df$Count,ties.method = "min")

#results_mat_long_score_order_df[results_mat_long_score_order_df$AverageScore == -1.3198,    9] <- 1
#results_mat_long_score_order_df[results_mat_long_score_order_df$AverageScore == -1.1886,  9] <- 2
#results_mat_long_score_order_df[results_mat_long_score_order_df$AverageScore == -0.3205,  9] <- 3
#results_mat_long_score_order_df[results_mat_long_score_order_df$AverageScore == -0.1886,    9] <- 4
#results_mat_long_score_order_df[results_mat_long_score_order_df$AverageScore == -0.1213, 9] <- 5

#results_mat_long_score_order_df[results_mat_long_score_order_df$AverageScore == -0.0657, 9] <- 6
#results_mat_long_score_order_df[results_mat_long_score_order_df$AverageScore == -0.0506, 9] <- 7
#results_mat_long_score_order_df[results_mat_long_score_order_df$AverageScore == 0.0442, 9] <- 8
#results_mat_long_score_order_df[results_mat_long_score_order_df$AverageScore == 0.0655, 9] <- 9
#results_mat_long_score_order_df[results_mat_long_score_order_df$AverageScore == 0.2608,    9] <- 10

#results_mat_long_score_order_df[results_mat_long_score_order_df$AverageScore == 0.3014, 9] <- 11
#results_mat_long_score_order_df[results_mat_long_score_order_df$AverageScore == 0.3782 ,9] <- 12
#results_mat_long_score_order_df$AverageScore == 0.461 ,9] <- 13
#results_mat_long_score_order_df[results_mat_long_score_order_df$AverageScore == 0.5469 ,9] <- 14
#results_mat_long_score_order_df[results_mat_long_score_order_df$AverageScore == 0.5605 ,9] <- 15

#results_mat_long_score_order_df[results_mat_long_score_order_df$AverageScore == 0.6366 ,9] <- 16

#table(results_mat_long_score_order_df$Count)
#results_mat_long_score_order_df[results_mat_long_score_order_df$Count == -1,10] <- 0
#results_mat_long_score_order_df[results_mat_long_score_order_df$Count == 36,10] <- 1
#results_mat_long_score_order_df[results_mat_long_score_order_df$Count == 56,10] <- 2
#results_mat_long_score_order_df[results_mat_long_score_order_df$Count == 77,10] <- 4

#results_mat_long_score_order_df[results_mat_long_score_order_df$Count == 79,10] <- 5
#results_mat_long_score_order_df[results_mat_long_score_order_df$Count == 85,10] <- 6

#results_mat_long_score_order_df[results_mat_long_score_order_df$Count == 140,10] <- 7
#results_mat_long_score_order_df[results_mat_long_score_order_df$Count == 142,10] <- 8

#results_mat_long_score_order_df[results_mat_long_score_order_df$Count == 156,10] <- 9
#results_mat_long_score_order_df[results_mat_long_score_order_df$Count == 166,10] <- 10
#results_mat_long_score_order_df[results_mat_long_score_order_df$Count == 175,10] <- 11
#results_mat_long_score_order_df[results_mat_long_score_order_df$Count == 259,10] <- 12

#results_mat_long_score_order_df[results_mat_long_score_order_df$Count == 263,10] <- 13
#results_mat_long_score_order_df[results_mat_long_score_order_df$Count == 306,10] <- 14
#results_mat_long_score_order_df[results_mat_long_score_order_df$Count == 344,10] <- 15
#results_mat_long_score_order_df[results_mat_long_score_order_df$Count == 348,10] <- 16
#results_mat_long_score_order_df[results_mat_long_score_order_df$Count == 377,10] <- 17

#results_mat_long_score_order_df[results_mat_long_score_order_df$Count == 409,10] <- 18
#results_mat_long_score_order_df[results_mat_long_score_order_df$Count == 426,10] <- 19
#results_mat_long_score_order_df[results_mat_long_score_order_df$Count == 432,10] <- 20
#results_mat_long_score_order_df[results_mat_long_score_order_df$Count == 466,10] <- 21
#results_mat_long_score_order_df[results_mat_long_score_order_df$Count == 492,10] <- 22


#results_mat_long_score_order_df[results_mat_long_score_order_df$Count == 508,10] <- 24
#results_mat_long_score_order_df[results_mat_long_score_order_df$Count == 518,10] <- 25
#results_mat_long_score_order_df[results_mat_long_score_order_df$Count == 526,10] <- 26
#results_mat_long_score_order_df[results_mat_long_score_order_df$Count == 536,10] <- 27

#results_mat_long_score_order_df[results_mat_long_score_order_df$Count == 551,10] <- 28
#results_mat_long_score_order_df[results_mat_long_score_order_df$Count == 1096,10] <- 29

#results_mat_long_score_order_df[,9] <- as.numeric(results_mat_long_score_order_df[,9])
#results_mat_long_score_order_df[,10] <- as.numeric(results_mat_long_score_order_df[,10])

results_mat_long_score_order_df$Names <- factor(results_mat_long_score_order_df$Names,levels=unique(results_mat_long_score_order_df$Names[order(results_mat_long_score_order_df$Count_Rank)]))
results_mat_long_score_order_df$Model <- factor(results_mat_long_score_order_df$Model,levels=unique(results_mat_long_score_order_df$Model[order(results_mat_long_score_order_df$Score_Rank)]))

# Drop <= 10 Positive Samples
#results_mat_long_score_order_df <- results_mat_long_score_order_df[!(results_mat_long_score_order_df$Names %in% c("PRAD","KICH")),]

# Scaled Cindex
ggplot(results_mat_long_score_order_df, aes(x = Model , y = Names)) + 
  geom_tile(aes(fill  = Cindex))  +
  scale_fill_gradient2(low='black',high='red',mid='white',midpoint=0.0,limits=c(-3,3),oob=squish) +
  theme_bw() + ylab("Data Set (ordered by size)") + xlab("Model (ordered by components)") + theme(axis.text.x = element_text(angle = 90))
