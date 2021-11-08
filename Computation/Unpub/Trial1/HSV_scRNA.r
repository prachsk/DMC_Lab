library(dplyr)
library(Seurat)
library(patchwork)
library(readxl)

#Read data from .tab and cleanup
hsv.data <- read.table('~/Google\ Drive/My\ Drive/Meletis/Unpub_RNA-seq/LHA_HSV_sc_RNA-seq/hsv_rpkms.tab', header = TRUE)
hsv.genes <- make.names(hsv.data$gene, unique = TRUE)
hsv.data["gene" ] <- NULL

#Create Seurat Object
hsv.se <- CreateSeuratObject(hsv.data, project = "HSV_sc_RNA-seq", row.names = hsv.genes)
hsv.se

# Visualize QC metrics as a violin plot
VlnPlot(hsv.se, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
plot1 <- FeatureScatter(hsv.se, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1

#Subset the data
hsv.se <- subset(hsv.se, subset = nFeature_RNA > 1500)

#Normalize the data
hsv.se <- NormalizeData(hsv.se)

#Identify high variable genes
hsv.se <- FindVariableFeatures(hsv.se, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(hsv.se), 10)
plot1 <- VariableFeaturePlot(hsv.se)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Scaling the data
all.genes <- rownames(hsv.se)
hsv.se <- ScaleData(hsv.se, features = all.genes)

#Run PCA
hsv.se <- RunPCA(hsv.se, npcs = 50,features = VariableFeatures(object = hsv.se))
VizDimLoadings(hsv.se, dims = 1:2, reduction = "pca")
ElbowPlot(hsv.se, ndims = 25)

#Find clusters
hsv.se <- FindNeighbors(hsv.se, dims = 1:15)
hsv.se <- FindClusters(hsv.se, resolution = 1.2)

#Run UMAP
hsv.se <- RunUMAP(hsv.se, dims = 1:15)
DimPlot(hsv.se, reduction = "umap")

#Save
saveRDS(hsv.se, './HSV_scRNA.rds')

#Finding all markers within the clusters
hsv.markers <- FindAllMarkers(hsv.se, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#Get the top 2 markers of each cluster
top2.genes <- hsv.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

#Get the top 1 markers of each cluster
top1.genes <- hsv.markers %>%
  group_by(cluster) %>%
  slice_max(n = 1, order_by = avg_log2FC)

#Violin plot and feature plot of the top 1&2 markers
VlnPlot(hsv.se, features = top2.genes$gene)
VlnPlot(hsv.se, features = top1.genes$gene)
FeaturePlot(hsv.se, features = top1.genes$gene)

#Making heatmap
hsv.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(hsv.se, features = top10$gene) + NoLegend()














