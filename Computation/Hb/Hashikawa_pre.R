library(dplyr)
library(patchwork)
library(ggplot2)
library(scales)
library(Seurat)
library(SeuratObject)
library(DoubletDecon)

# Get the list of folders
folders <- list.files("/Users/pax/Google\ Drive/My\ Drive/Meletis/Data_reproduce/GSE137478/GSE137478_RAW", full.names = T)
proj.names <- list.files("/Users/pax/Google\ Drive/My\ Drive/Meletis/Data_reproduce/GSE137478/GSE137478_RAW")
# Read 10X data in the folders into a list
data.list <- lapply(folders, Read10X)
names(data.list) <- paste0(substr(folders, nchar(folders)-9, nchar(folders)), ".data")

# Create Seurat object from the list of data
se.list1 <- list()
for (i in seq_along(data.list)) {
  se.list1[[i]] <- CreateSeuratObject(data.list[[i]], project = substr(folders, nchar(folders)-9, nchar(folders)), min.cells = 3, min.features = 200)
}
names(se.list) <- paste0(substr(folders, nchar(folders)-9, nchar(folders)), ".se")

location='./DoubletDecon/'
for (i in seq_along(se.list)) {
  se.list[[i]][["source"]] <- "Hashikawa"
  se.list[[i]][["percent.mt"]] <- PercentageFeatureSet(se.list[[i]], pattern = "^Mt.")
  se.list[[i]] <- subset(se.list[[i]], subset = nCount_RNA > 700 & nCount_RNA < 15000 & percent.mt < 20)
  se.list[[i]] <- NormalizeData(se.list[[i]])
  se.list[[i]] <- FindVariableFeatures(se.list[[i]], selection.method = "vst", nfeatures = 2000)
  se.list[[i]] <- ScaleData(se.list[[i]], verbose = FALSE)
  se.list[[i]] <- RunPCA(se.list[[i]], npcs = 30, verbose = FALSE)
  se.list[[i]] <- RunUMAP(se.list[[i]], reduction = "pca", dims = 1:30)
  se.list[[i]] <- FindNeighbors(se.list[[i]], reduction = "pca", dims = 1:30)
  se.list[[i]] <- FindClusters(se.list[[i]], resolution = 0.5)
  
  # Run DoubletDecon workdflow
  newFiles <- Improved_Seurat_Pre_Process(se.list[[i]], num_genes=50, write_files=FALSE)
  filename <- paste0("", proj.names[i])
  write.table(newFiles$newExpressionFile, paste0(location, filename, "_expression"), sep="\t")
  write.table(newFiles$newFullExpressionFile, paste0(location, filename, "_fullExpression"), sep="\t")
  write.table(newFiles$newGroupsFile, paste0(location, filename , "_groups"), sep="\t", col.names = F)
  results <- Main_Doublet_Decon(rawDataFile=newFiles$newExpressionFile, 
                             groupsFile=newFiles$newGroupsFile, 
                             filename=filename, 
                             location=location,
                             fullDataFile=NULL, 
                             removeCC=FALSE, 
                             species="hsa", 
                             rhop=1.1, 
                             write=TRUE, 
                             PMF=TRUE, 
                             useFull=FALSE, 
                             heatmap=FALSE,
                             centroids=TRUE,
                             num_doubs=100, 
                             only50=FALSE,
                             min_uniq=4,
                             nCores=-1)
  
  # Rename cell ID of results$DRS_doublet_table to match with se.combined
  rownames(results$DRS_doublet_table) <- gsub("\\.","-",rownames(results$DRS_doublet_table))
  
  # Make a list of cell IDs That are doublet
  is.doublet <- rownames(results$DRS_doublet_table)[results$DRS_doublet_table$isADoublet == "TRUE"]
  
  # Add doublet to metadata
  se.list[[i]][['doublet']] <- 0
  se.list[[i]]$doublet[which(names(se.list[[i]]$nCount_RNA)%in% is.doublet)] <- 1
  
  # Subset to remove doublet
  se.list[[i]] <- subset(se.list[[i]], subset = doublet == 1, invert = T)
}

# Select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = se.list, nfeatures = 2000)

# Find integration anchors
se.anchors <- FindIntegrationAnchors(object.list = se.list, anchor.features = features)

# Creates an 'integrated' data assay
se.combined <- IntegrateData(anchorset = se.anchors)

# Add source metadata
se.combined[["source"]] <- "Hashikawa"

# Save RDS
saveRDS(se.combined, "./RDS_obj/Hashikawa_combined.RDS")

DefaultAssay(se.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
se.combined <- ScaleData(se.combined, verbose = FALSE)
se.combined <- RunPCA(se.combined, npcs = 30, verbose = FALSE)
se.combined <- RunUMAP(se.combined, reduction = "pca", dims = 1:30)
se.combined <- FindNeighbors(se.combined, reduction = "pca", dims = 1:30)
se.combined <- FindClusters(se.combined, resolution = 0.5)
DimPlot(se.combined, reduction = 'umap')

# Save RDS
saveRDS(se.combined, "./RDS_obj/Hashikawa_combined.RDS")

