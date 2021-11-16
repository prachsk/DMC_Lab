library(dplyr)
library(Seurat)
library(patchwork)
library(readxl)

# Read LHA data from .tab and cleanup
lh.data <- read.table('~/Google\ Drive/My\ Drive/Meletis/Unpub_RNA-seq/LHA_sc_RNA-seq/lh_rpkms.tab', header = TRUE)
lh.nm <- make.names(lh.data$gene, unique = TRUE)
lh.data["gene"] <- NULL

# Create Seurat object
lh.se <- CreateSeuratObject(lh.data, project = "LHA", row.names = lh.nm)
lh.se

# Add source to metadata
lh.se[["source"]] <- rep(c("Unpub"), times = ncol(lh.se))

# Mt genes QC
lh.se[["percent.mt"]] <- PercentageFeatureSet(lh.se, pattern = "^Mt.")

# Visualize QC metrics as a violin plot
VlnPlot(lh.se, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(lh.se, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1

# Subset the data
lh.se <- subset(lh.se, subset = nFeature_RNA > 2000 & percent.mt < 1.0)

# Normalize the data
lh.se <- NormalizeData(lh.se)

# Identify high variable genes
lh.se <- FindVariableFeatures(lh.se, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(lh.se), 10)
plot1 <- VariableFeaturePlot(lh.se)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scaling the data
all.genes <- rownames(lh.se)
lh.se <- ScaleData(lh.se, features = all.genes)

# Run PCA
lh.se <- RunPCA(lh.se, npcs = 50,features = VariableFeatures(object = lh.se))
VizDimLoadings(lh.se, dims = 1:2, reduction = "pca")
ElbowPlot(lh.se, ndims = 25)

# Find clusters
lh.se <- FindNeighbors(lh.se, dims = 1:17)
lh.se <- FindClusters(lh.se, resolution = 0.5)

# Run UMAP
lh.se <- RunUMAP(lh.se, dims = 1:17)
DimPlot(lh.se, reduction = "umap")

# Save RDS
saveRDS(lh.se, '../RDS_obj/LHA_scRNA.RDS')

# Finding all markers within the clusters
lh.markers <- FindAllMarkers(lh.se, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Get the top 2 markers of each cluster
top2.genes <- lh.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# Get the top 1 markers of each cluster
top1.genes <- lh.markers %>%
  group_by(cluster) %>%
  slice_max(n = 1, order_by = avg_log2FC)

# Violin plot and feature plot of the top 1&2 markers
VlnPlot(lh.se, features = top2.genes$gene)
VlnPlot(lh.se, features = top1.genes$gene)
FeaturePlot(lh.se, features = top1.genes$gene)

# Making heatmap
lh.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(lh.se, features = top10$gene) + NoLegend()


