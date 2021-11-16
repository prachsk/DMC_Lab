library(dplyr)
library(Seurat)
library(patchwork)
library(readxl)

# Read data from .txt and cleanup
LHb_lskra.data <- read.table('/Users/pax/Google\ Drive/My\ Drive/Meletis/Unpub_RNA-seq/LHb_Projection_data/LHb_Iskra.txt', header = TRUE)
LHb_lskra.genes <- LHb_lskra.data$GeneID
LHb_lskra.data['GeneID'] <- NULL
LHb_lskra.se <- CreateSeuratObject(LHb_lskra.data, project = 'LHb_lskra', row.names = LHb_lskra.genes)

# Add source to metadata
LHb_lskra.se[["source"]] <- rep(c("Unpub"), times = ncol(LHb_lskra.se))

# copy origi.ident to area
LHb_lskra.se[["area"]] <- LHb_lskra.se$orig.ident

# change orig.ident to LHb_lskra
LHb_lskra.se$orig.ident <- rep(c("LHb_lskra"), times = ncol(LHb_lskra.se))

# mt qc and Subset LHb_lskra.se
LHb_lskra.se[["percent.mt"]] <- PercentageFeatureSet(LHb_lskra.se, pattern = "^Mt.")
VlnPlot(LHb_lskra.se, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(LHb_lskra.se, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = 'orig.ident')
plot1
LHb_lskra.se <- subset(LHb_lskra.se, subset = nFeature_RNA > 1250 & nFeature_RNA < 1.0)

# Normalize the data
LHb_lskra.se <- NormalizeData(LHb_lskra.se)

# Identify high variable genes
LHb_lskra.se <- FindVariableFeatures(LHb_lskra.se, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(LHb_lskra.se), 10)
plot1 <- VariableFeaturePlot(LHb_lskra.se)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scaling the data
all.genes <- rownames(LHb_lskra.se)
LHb_lskra.se <- ScaleData(LHb_lskra.se, features = all.genes)

# Run PCA
LHb_lskra.se <- RunPCA(LHb_lskra.se, npcs = 50,features = VariableFeatures(object = LHb_lskra.se))
VizDimLoadings(LHb_lskra.se, dims = 1:2, reduction = "pca")
ElbowPlot(LHb_lskra.se, ndims = 25)

# Find clusters
LHb_lskra.se <- FindNeighbors(LHb_lskra.se, dims = 1:20)
LHb_lskra.se <- FindClusters(LHb_lskra.se, resolution = 0.5)

# Run UMAP
LHb_lskra.se <- RunUMAP(LHb_lskra.se, dims = 1:20)
DimPlot(LHb_lskra.se, reduction = "umap")

# Save RDS
saveRDS(LHb_lskra.se, '../RDS_obj/LHb_lskra.RDS')







