# Time_and_Dose
R code for Time and Dose Project

Instruction
  - At the top of each script is setwd command, change the file path to the location of the Time_and_Dose folder
  - Scripts should be run in order from 0 to 7 although processing single-cell data through Seurat is data intensive and some tasks
  are best run on a server/data machine if available.
  - After running script 6, move the results to CoxModels_Updated and run the scripts there
    1. The get_dataset scripts will pull data from TCGA
    2. The RunCoxModel scripts will generate the cox models
  - Note: Putative doublet checks may fail if run before the doublet check scripts, but this is supplemental and can be ignored

Notes:

No raw sequencing data is included because of file size limitations
  - scRNA data can be obtained from GEO (this will be addressed in the paper)
  - script for obtaining TCGA data are included (FPKM data may be deprecated so might have to be calculated from raw counts)
  - Intermediate files are also excluded because most are larger than the 25Mb limit but can be generated from the raw data

Time_and_Dose.zip contains the first version of the code which has fewer corrections for confunding variables
