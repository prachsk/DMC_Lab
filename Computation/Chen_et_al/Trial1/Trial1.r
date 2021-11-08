library(dplyr)
library(Seurat)
library(patchwork)

#Read in count matrix from .txt file
cell.data <- read.delim(file = '/Users/pax/Google\ Drive/My\ Drive/Meletis/Data/GSE87544/GSE87544_Merged_17samples_14437cells_count.txt.gz', header = TRUE, sep = '\t')
cell.genes <- cell.data$Gene
cell.data['Gene'] <- NULL
dim(cell.data)

#Create Seurat obj from the count matrix
cell.se <- CreateSeuratObject(counts = cell.data, project = 'LHA', row.names = cell.genes)

VlnPlot(cell.se, features = c('nFeature_RNA', 'nCount_RNA'), ncol = 2)

#Select cells with at least 2000 genes
cell.se <- subset(cell.se, subset = nFeature_RNA > 2000)

#Find high variable genes
cell.se <- FindVariableFeatures(cell.se)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(cell.se), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(cell.se)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Scaling the data
all.genes <- rownames(cell.se)
cell.se <- ScaleData(cell.se, features = all.genes)

#Run PCA
cell.se <- RunPCA(cell.se, npcs = 50, features = VariableFeatures(object = cell.se))
VizDimLoadings(cell.se, dims = 1:2, reduction = "pca")
ElbowPlot(cell.se, ndims = 25)

#Find cluster
cell.se <- FindNeighbors(cell.se, dims = 1:17)
cell.se <- FindClusters(cell.se, resolution = 0.1)

#Run tsne
cell.se <- RunTSNE(cell.se, dims = 1:17)
DimPlot(cell.se, reduction = "tsne")

#Save Data
saveRDS(cell.se, '/Users/pax/Google\ Drive/My Drive/Meletis/Computation/Chen_et_al/Trial1/R_obj/Trial1_merge.rds')




#Iter1
#Subset the largest clustr
iter1.se <- subset(cell.se, idents = 0)

#Find high variable genes
iter1.se <- FindVariableFeatures(iter1.se)

#Scaling the data
all.genes_c0 <- rownames(iter1.se)
iter1.se <- ScaleData(iter1.se, features = all.genes_c0)

#Run PCA
iter1.se <- RunPCA(iter1.se, npcs = 50, features = VariableFeatures(object = iter1.se))
VizDimLoadings(iter1.se, dims = 1:2, reduction = "pca")
ElbowPlot(iter1.se, ndims = 25)

#Find cluster
iter1.se <- FindNeighbors(iter1.se, dims = 1:15)
iter1.se <- FindClusters(iter1.se, resolution = 0.5)

#Run tsne
iter1.se <- RunTSNE(iter1.se, dims = 1:15)
DimPlot(iter1.se, reduction = "tsne")



#Iter2
#Subset the largest cluster
iter2.se <- subset(iter1.se, idents = 0)

#Find high variable genes
iter2.se <- FindVariableFeatures(iter2.se)

#Scaling the data
all.genes_c0 <- rownames(iter2.se)
iter2.se <- ScaleData(iter2.se, features = all.genes_c0)

#Run PCA
iter2.se <- RunPCA(iter2.se, npcs = 50, features = VariableFeatures(object = iter2.se))
VizDimLoadings(iter2.se, dims = 1:2, reduction = "pca")
ElbowPlot(iter2.se, ndims = 20)

#Find cluster
iter2.se <- FindNeighbors(iter2.se, dims = 1:15)
iter2.se <- FindClusters(iter2.se, resolution = 0.5)

#Run tsne
iter2.se <- RunTSNE(iter2.se, dims = 1:15)
DimPlot(iter2.se, reduction = "tsne")




#Iter3
#Subset the largest cluster
iter3.se <- subset(iter2.se, idents = 0)

#Find high variable genes
iter3.se <- FindVariableFeatures(iter3.se)

#Scaling the data
all.genes_c0 <- rownames(iter3.se)
iter3.se <- ScaleData(iter3.se, features = all.genes_c0)

#Run PCA
iter3.se <- RunPCA(iter3.se, npcs = 50, features = VariableFeatures(object = iter3.se))
VizDimLoadings(iter3.se, dims = 1:2, reduction = "pca")
ElbowPlot(iter3.se, ndims = 20)

#Find cluster
iter3.se <- FindNeighbors(iter3.se, dims = 1:15)
iter3.se <- FindClusters(iter3.se, resolution = 1)

#Run tsne
iter3.se <- RunTSNE(iter3.se, dims = 1:15)
DimPlot(iter3.se, reduction = "tsne")