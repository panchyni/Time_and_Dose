#STAD_dataset.R
#Purpose: TCGA-STAD analysis

### Setup Environemt ###

# Folder
setwd("E:/2021/Work_06_17_2021_TGFBTimeSeries/CoxModels_Updated/")

# Install Packages
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#
#BiocManager::install("TCGAbiolinks")

#Loading needed packages----
library(ggplot2)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(tidyverse)

#Loading the merged data frame----
coad_query <- GDCquery(project       = "TCGA-STAD",
                       data.category = "Transcriptome Profiling",
                       data.type     = "Gene Expression Quantification",
                       workflow.type = "HTSeq - FPKM")

#Downloading the files
GDCdownload(query           = coad_query,
            method          = "api",
            files.per.chunk = 10,
            directory       = "Data/")

#Making the SummarizedExperiment object
coad_data_se <- GDCprepare(coad_query, summarizedExperiment = TRUE,
                           directory = "Data/")
coad_data_df <- as.data.frame(colData(coad_data_se))
coad_data_df$vital_status <- factor(coad_data_df$vital_status,
                                    levels = c("Alive", "Dead"),
                                    labels = c(0,1))

coad_data_df$vital_status <- as.numeric(as.character(coad_data_df$vital_status))
bulk_rna_df <- coad_data_se@assays@data@listData[["HTSeq - FPKM"]]
colnames(bulk_rna_df) <- coad_data_se@colData@rownames
rownames(bulk_rna_df) <- coad_data_se@rowRanges@elementMetadata@listData[["external_gene_name"]]
bulk_rna_df <- t(bulk_rna_df)
bulk_rna_df <- as.data.frame(bulk_rna_df)
bulk_rownames <- rownames(bulk_rna_df)
bulk_rna_df$barcode <- bulk_rownames
bulk_rna_df_unique <- subset(bulk_rna_df,select = unique(colnames(bulk_rna_df)))

coad_data_df_unique <- subset(coad_data_df,select = unique(colnames(coad_data_df)))
merged_df <- merge(bulk_rna_df_unique, coad_data_df_unique, by = 'barcode')
rownames(merged_df) <- merged_df$barcode
merged_df <- merged_df[,2:length(colnames(merged_df))]
merged_df$days_to_last_follow_up <- ifelse(merged_df$vital_status==1, merged_df$days_to_death, merged_df$days_to_last_follow_up)
merged_df <- filter(merged_df, days_to_last_follow_up != "NA")
cox_time <- merged_df$days_to_last_follow_up

cox_event <- merged_df$vital_status
cox_tumor <- merged_df$ajcc_pathologic_stage
cox_tumor_n <- merged_df$ajcc_pathologic_n
cox_tumor_m <- merged_df$ajcc_pathologic_m
cox_gender <- merged_df$gender
cox_eth <- merged_df$ethnicity
cox_race <- merged_df$race
cox_type <- merged_df$definition
cox_df <- subset(merged_df, select=c(TSPAN6:AC007389.5))
cox_df$days.to.last.follow.up <- cox_time
cox_df$vital.status <- cox_event
cox_df$tumor.stage <- cox_tumor
cox_df$ajcc.m <- cox_tumor_m
cox_df$ajcc.n <- cox_tumor_n
cox_df$race <- cox_race
cox_df$ethnicity <- cox_eth
cox_df$gender <- cox_gender
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="A", replacement="")
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="B", replacement="")
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern="C", replacement="")
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "Stage iv", replacement = 4)
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "Stage iii", replacement = 3)
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "Stage ii", replacement = 2)
cox_df$tumor.stage <- gsub(cox_df$tumor.stage, pattern = "Stage i", replacement = 1)
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="a", replacement="")
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="b", replacement="")
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="c", replacement="")
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="N0", replacement=0)
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="N1", replacement=1)
cox_df$ajcc.n <- gsub(cox_df$ajcc.n, pattern="N2", replacement=2)
cox_df$sample.type <- cox_type
cox_df <- filter(cox_df, !tumor.stage=="not reported")
#cox_df <- cox_df[complete.cases(cox_df[, "ajcc.m"]), ]
cox_df$days.to.last.follow.up <- ifelse(cox_df$days.to.last.follow.up < 1, 1,cox_df$days.to.last.follow.up)
save(cox_df, file = "Data/TCGA-STAD/STAD_df_fpkm_update.RData")