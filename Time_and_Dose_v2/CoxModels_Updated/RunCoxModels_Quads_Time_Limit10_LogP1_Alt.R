# Folder
setwd("/Users/nicholaspanchy/Documents/Work_UTK/DoseTime_FullCorr_FromScratch/CoxModels_Updated/")

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
clust_int <- readRDS(file = "Data/Clusters_IntegrationFirst_Quart_3rdQuarFC_TimeOnlyAlt.rds")

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

save(results_mat, file = "MultiCancerResultsTable_QuartTimeLogP1_Alt.RData")
save(activegenes_mat, file = "MultiCancerGeneCountTable_QuartTime_LogP1_AltL.RData")

save(activegenes_list, file = "MultiCancerGeneLists_QuartTime_LogP1_Alt.RData")
save(i_list, file = "MultiCancerGeneilist_QuartTime_LogP1L_Alt.RData")
save(j_list, file = "MultiCancerjlist_QuartTime_LogP1_Alt.RData")


### Process ###

library(ggplot2)
library(scales)

load("MultiCancerResultsTable_QuartTimeLogP1_Alt.RData")
load("MultiCancerGeneCountTable_QuartTime_LogP1_AltL.RData")

result_df <- as.data.frame(results_mat)

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
summary(results_mat_quad_df$V4)

ggplot() + geom_point(data = results_mat_quad_df, aes(x = V1, y = V2,fill=V4,size=V5), color='black',pch=21) +
  scale_fill_gradient2(low='black',high='red',mid='white',midpoint=0.645,limits=c(0.57,0.70),oob=squish) +
  coord_fixed(ratio = 1) + scale_size(range = c(3, 12)) + theme_bw() + ylim(c(0.75,4.25)) + xlim(c(0.75,4.25))
