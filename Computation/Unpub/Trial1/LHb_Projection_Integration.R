library(dplyr)
library(Seurat)
library(patchwork)

#Read data from .txt and cleanup
LHb_lskra.data <- read.table('/Users/pax/Google\ Drive/My\ Drive/Meletis/Unpub_RNA-seq/LHb_Projection_data/LHb_Iskra.txt', header = TRUE)
LHb_lskra.genes <- LHb_lskra.data$GeneID
LHb_lskra.data['GeneID'] <- NULL
LHb_lskra.se <- CreateSeuratObject(LHb_lskra.data, project = 'LHb_lskra', row.names = LHb_lskra.genes)

LH_lskra.data <- read.table('/Users/pax/Google\ Drive/My\ Drive/Meletis/Unpub_RNA-seq/LHb_Projection_data/Iskra_LH.txt', header = TRUE)
LH_lskra.genes <- rownames(LH_lskra.data)
row.names(LH_lskra.data) <- NULL
LH_lskra.se <- CreateSeuratObject(LH_lskra.data, project = 'LH_lskra', row.names = LH_lskra.genes)

#Subset LHb_lskra.se
VlnPlot(LHb_lskra.se, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
plot1 <- FeatureScatter(LHb_lskra.se, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = 'orig.ident')
plot1
LHb_lskra.se <- subset(LHb_lskra.se, subset = nFeature_RNA > 1250)

#Subset LH_lskra.se
VlnPlot(LH_lskra.se, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
plot1 <- FeatureScatter(LH_lskra.se, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = 'orig.ident')
plot1
LH_lskra.se <- subset(LH_lskra.se, subset = nFeature_RNA > 2500)

#Split LHb_lskra.se and LH_lskra.se into the list of 2 obj each by its orig.ident
LHb.list <- SplitObject(LHb_lskra.se, split.by = "orig.ident")
LH.list <- SplitObject(LH_lskra.se, split.by = "orig.ident")

#Normalize and identify variable features for each dataset in LHb.list and lH.list
LHb.list <- lapply(X = LHb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

LH.list <- lapply(X = LH.list, FUN = function(y) {
  y <- NormalizeData(y)
  y <- FindVariableFeatures(y, selection.method = "vst", nfeatures = 2000)
})

#Select features that are repeatedly variable across datasets in LHb.list and Lh.list for integration
LHb.list_features <- SelectIntegrationFeatures(object.list = LHb.list)
LH.list_features <- SelectIntegrationFeatures(object.list = LH.list)

#Find anchors
LHb.anchors <- FindIntegrationAnchors(object.list = LHb.list, anchor.features = LHb.list_features)
LH.anchors <- FindIntegrationAnchors(object.list = LH.list, anchor.features = LH.list_features)

#This command creates an 'integrated' data assay
LHb.combined <- IntegrateData(anchorset = LHb.anchors)
LH.combined <- IntegrateData(anchorset = LH.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(LHb.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
LHb.combined <- ScaleData(LHb.combined, verbose = FALSE)
LHb.combined <- RunPCA(LHb.combined, npcs = 30, verbose = FALSE)
LHb.combined <- RunUMAP(LHb.combined, reduction = "pca", dims = 1:30)
LHb.combined <- FindNeighbors(LHb.combined, reduction = "pca", dims = 1:30)
LHb.combined <- FindClusters(LHb.combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(LHb.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(LHb.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2


# original unmodified data still resides in the 'RNA' assay
DefaultAssay(LH.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
LH.combined <- ScaleData(LH.combined, verbose = FALSE)
LH.combined <- RunPCA(LH.combined, npcs = 30, verbose = FALSE)
LH.combined <- RunUMAP(LH.combined, reduction = "pca", dims = 1:30)
LH.combined <- FindNeighbors(LH.combined, reduction = "pca", dims = 1:30)
LH.combined <- FindClusters(LH.combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(LH.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(LH.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2















