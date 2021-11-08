library(dplyr)
library(Seurat)
library(patchwork)
library(readxl)

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

# Create list of the 2 se objs
LH_LHb.list <- list(LH_lskra.se, LH_lskra.se)

#Normalize and identify variable features for each dataset independently
LH_LHb.list <- lapply(X = LH_LHb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

#Select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = LH_LHb.list)

LH_LHb.anchors <- FindIntegrationAnchors(object.list = LH_LHb.list, anchor.features = features)

#Creates an 'integrated' data assay
LH_LHb.combined <- IntegrateData(anchorset = LH_LHb.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
DefaultAssay(LH_LHb.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
LH_LHb.combined <- ScaleData(LH_LHb.combined, verbose = FALSE)
LH_LHb.combined <- RunPCA(LH_LHb.combined, npcs = 30, verbose = FALSE)
LH_LHb.combined <- RunUMAP(LH_LHb.combined, reduction = "pca", dims = 1:30)
LH_LHb.combined <- FindNeighbors(LH_LHb.combined, reduction = "pca", dims = 1:30)
LH_LHb.combined <- FindClusters(LH_LHb.combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(LH_LHb.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(LH_LHb.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
DimPlot(LH_LHb.combined, reduction = "umap", split.by = "orig.ident")

#Performing differential expression after integration, switch back to the original data
DefaultAssay(LH_LHb.combined) <- "RNA"
c1.markers <- FindConservedMarkers(LH_LHb.combined, ident.1 = 1, grouping.var = "orig.ident", verbose = FALSE)
head(c1.markers)

rownames(head(c1.markers))

FeaturePlot(LH_LHb.combined, features = c(rownames(head(c1.markers))), min.cutoff = "q9")

saveRDS(LH_LHb.combined, './LH_LHb_integ.rds')


