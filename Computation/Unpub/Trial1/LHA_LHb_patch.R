library(dplyr)
library(Seurat)
library(patchwork)
library(readxl)

#Load data from xlsx
LHA_LHb_patch.data <- read_xlsx('~/Google\ Drive/My\ Drive/Meletis/Unpub_RNA-seq/LHA-LHb_Patch-Seq/TPM_FINAL_169.xlsx', col_names = FALSE)

#Create new df for Ephys_Cathegory
LHA_LHb_patch.ephys <- LHA_LHb_patch.data[2,2:length(LHA_LHb_patch.data)]
row.names(LHA_LHb_patch.ephys) <- LHA_LHb_patch.data[2,1]

#Data clean up and add row & col names
LHA_LHb_patch.data_clean <- LHA_LHb_patch.data[5:dim(LHA_LHb_patch.data)[1],]
row.names(LHA_LHb_patch.data_clean) <- make.names(LHA_LHb_patch.data_clean$...1, unique = TRUE)
LHA_LHb_patch.data_clean$...1 <- NULL
colnames(LHA_LHb_patch.data_clean) <- LHA_LHb_patch.data[4,2:length(LHA_LHb_patch.data)]


