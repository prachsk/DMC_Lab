library(Seurat)
library(dplyr)
library(Matrix)
library(abind)
library(purrr)

# Read files into list
folder <- list.files("/Users/pax/Google\ Drive/My\ Drive/Meletis/Data_reproduce/GSE130597/GSE130597", pattern="*.gz", full.names=TRUE)
data.list <- lapply(folder, read.table, sep="\t", header=TRUE, row.names=1)
names(data.list) <- paste0(substr(folder, 125, nchar(folder)-7), ".data")
data.list <- data.list[order(names(data.list))]

# Create Seurat Object
dge.se.list <- list()
proj_name <- c()
for (i in seq_along(data.list)){ 
  proj_name <- c(proj_name, substr(names(data.list[i]), 1, nchar(names(data.list[i]))-5))
  dge.se.list[[i]] <- CreateSeuratObject(data.list[[i]], project = proj_name[i], min.cells = 10)}
names(dge.se.list) <- substr(names(data.list), 1, 8)

# Creat a list for objects to merge and merge all the objects in the list
list.to.merge <- list()
for (i in seq(2, length(dge.se.list))) {
  list.to.merge <- append(list.to.merge, dge.se.list[[i]])
}
dge.merge <- merge(x = dge.se.list[[1]], y = list.to.merge, add.cell.ids = c(proj_name), project = "control.highfat.counts")

# Calculate percent mt
dge.merge[["percent.mt"]] <- PercentageFeatureSet(dge.merge, pattern = "^mt-")

# Create metadata for control and high-fat group
group.list <- list()

for (i in 1:length(dge.se.list)){
  if (i %% 2 == 1) 
    group.list[[i]] <- array(rep("control",length(names(dge.se.list[[i]]$nCount_RNA))))
  if (i %% 2 != 1) 
    group.list[[i]] <- array(rep("highfat",length(names(dge.se.list[[i]]$nCount_RNA))))
}

names(group.list) <- paste0(substr(folder[i], 125, nchar(folder)-7), ".group")

for (i in seq_along(dge.se.list)){
  dimnames(group.list[[i]]) <- list(names(dge.se.list[[i]]$nCount_RNA))
  }      
group <- do.call("abind", group.list)

# Add group metadata to the object
dge.merge <- AddMetaData(dge.merge, metadata = group, col.name = "group")

# Create batch metadata
batch.list <- list()

batch.A <- list("1" = c(1,2), "2" = c(3,4), "3" = c(5,6), "4" = c(7,8), "5" = c(9,10), "6" = c(11,12), "7" = c(13,14))

for (i in 1:length(dge.se.list)){
  for(j in 1:length(batch.A)){
    if (i %in% batch.A[[j]]) 
      batch.list[[i]] <- array(rep(paste0(names(batch.A[j])),length(names(dge.se.list[[i]]$nCount_RNA))))}}

names(batch.list) <- paste0(substr(folder[i], 125, nchar(folder)-7), ".batch.A")

for (i in seq_along(dge.se.list)){
  dimnames(batch.list[[i]]) <- list(names(dge.se.list[[i]]$nCount_RNA))}       

batch.A <- do.call("abind", batch.list)

# Add batch metadata to the object
dge.merge <- AddMetaData(dge.merge, batch.A, "batch.A")

# Visualize QC metrics as a violin plot
VlnPlot(dge.merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Visualize feature-feature relationships with FeatureScatter.
plot1 <- FeatureScatter(dge.merge, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(dge.merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Sub set for quality cells
dge.merge <- subset(dge.merge, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA < 15000 & percent.mt > 0.1 & percent.mt < 10)

# Normalize the data
dge.merge <- NormalizeData(dge.merge)

# Find variable features
dge.merge <- FindVariableFeatures(dge.merge, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(dge.merge), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(dge.merge)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scaling data
all.genes <- rownames(dge.merge)
dge.merge <- ScaleData(dge.merge, features = all.genes)

# Run PCA
dge.merge <- RunPCA(dge.merge, npcs = 120, features = VariableFeatures(object = dge.merge))
VizDimLoadings(dge.merge, dims = 1:2, reduction = "pca")
ElbowPlot(dge.merge, ndims = 50)

# Find Cluster
dge.merge <- FindNeighbors(dge.merge, dims = 1:100)
dge.merge <- FindClusters(dge.merge, resolution = 0.5)

# Run TSNE
dge.merge <- RunTSNE(dge.merge, dims = 1:100, check_duplicate = FALSE)
DimPlot(dge.merge, reduction = "tsne")

#Save RDS
saveRDS(dge.merge, "../Result_RDS/dge.merge.RDS")

