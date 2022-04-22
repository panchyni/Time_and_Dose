# Time_and_Dose
R code for Time and Dose Project

Instructions (brief)
  - At the top of each script is setwd command, change the file path to the location of the Time_and_Dose folder
  - Scripts should be run in order from 0 to 4 (2 and 2a can be run together)
  - Afte running script 4, move the results to CoxModels_Updated and run the scripts there
    1. The get_dataset scripts will pull data from TCGA
    2. The RunCoxModel scripts will generate the cox models

Notes:

No raw sequencing data is included because of file size limitations
  - scRNA data can be obtained from GEO (this will be addressed in the paper)
  - script for obtaining TCGA data are included
  
Panels for figures 1 to 3 will be generated automatically, but 4 and 5 need to be manually
saved to preseve spacing between points (otherwise there is alot of dead space)
