library(dplyr)
library(Seurat)
library(patchwork)
library(readxl)

# Load data from xlsx
LHA_LHb_patch.data <- read_xlsx('~/Google\ Drive/My\ Drive/Meletis/Unpub_RNA-seq/LHA-LHb_Patch-Seq/TPM_FINAL_169.xlsx', col_names = FALSE)

# Create df for Lib_info_Sample_Name
LHA_LHb_patch.lib <- LHA_LHb_patch.data[1,2:length(LHA_LHb_patch.data)]
row.names(LHA_LHb_patch.lib) <- LHA_LHb_patch.data[1,1]
LHA_LHb_patch.lib <- t(LHA_LHb_patch.lib)
LHA_LHb_patch.lib <- as.data.frame(LHA_LHb_patch.lib)

# Create new df for Ephys_Cathegory
LHA_LHb_patch.ephys <- LHA_LHb_patch.data[2,2:length(LHA_LHb_patch.data)]
row.names(LHA_LHb_patch.ephys) <- LHA_LHb_patch.data[2,1]
LHA_LHb_patch.ephys <- t(LHA_LHb_patch.ephys)
LHA_LHb_patch.ephys <- as.data.frame(LHA_LHb_patch.ephys)

# Create df for gene names
LHA_LHb_patch.genes <- LHA_LHb_patch.data$...1[5:dim(LHA_LHb_patch.data)[1]]
LHA_LHb_patch.genes <- make.names(LHA_LHb_patch.genes, unique = TRUE)

# Create new df for Animal_Cell ID
LHA_LHb_patch.AC_ID <- LHA_LHb_patch.data[3,2:length(LHA_LHb_patch.data)]

# Data clean up and add row & col names
LHA_LHb_patch.data_clean <- LHA_LHb_patch.data[5:dim(LHA_LHb_patch.data)[1],2:dim(LHA_LHb_patch.data)[2]]
rownames(LHA_LHb_patch.data_clean) <- LHA_LHb_patch.genes
colnames(LHA_LHb_patch.data_clean) <- LHA_LHb_patch.AC_ID

# Create Seurat obj
LHA_LHb_patch.se <- CreateSeuratObject(LHA_LHb_patch.data_clean, project = "LHA_LHb_patch")

# Create Metadata for Lib_info_Sample_Name and add metadata to the object
row.names(LHA_LHb_patch.lib) <- rownames(LHA_LHb_patch.se[[]])
LHA_LHb_patch.se <- AddMetaData(LHA_LHb_patch.se, LHA_LHb_patch.lib, col.name = "Lib_info_Sample_Name")

# Create Metadata for Ephys_Cathegory and add metadata to the object
row.names(LHA_LHb_patch.ephys) <- rownames(LHA_LHb_patch.se[[]])
LHA_LHb_patch.se <- AddMetaData(LHA_LHb_patch.se, LHA_LHb_patch.ephys, col.name = "Ephys_Cathegory")

# Add data source to metadata
LHA_LHb_patch.se[["source"]] <- rep(c("Unpub"), times = ncol(LHA_LHb_patch.se))

# Transform the transcript count into TPM
#LHA_LHb_patch.se$nCount_RNA <- (LHA_LHb_patch.se$nCount_RNA)/1e6

# Find percent of mt genes
LHA_LHb_patch.se[["percent.mt"]] <- PercentageFeatureSet(LHA_LHb_patch.se, pattern = "^mt.")

# Visualize QC metrics as a violin plot
VlnPlot(LHA_LHb_patch.se, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Visualize feature-feature relationships with FeatureScatter plot
plot1 <- FeatureScatter(LHA_LHb_patch.se, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(LHA_LHb_patch.se, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Subset the data
LHA_LHb_patch.se <- subset(LHA_LHb_patch.se, subset = percent.mt <= 20)

# Normalizing the data
LHA_LHb_patch.se <- NormalizeData(LHA_LHb_patch.se)

# Find variable features
LHA_LHb_patch.se <- FindVariableFeatures(LHA_LHb_patch.se, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(LHA_LHb_patch.se), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(LHA_LHb_patch.se)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scaling the data
all.genes <- rownames(LHA_LHb_patch.se)
LHA_LHb_patch.se <- ScaleData(LHA_LHb_patch.se, features = all.genes)

# Run PCA
LHA_LHb_patch.se <- RunPCA(LHA_LHb_patch.se, npcs = 50, features = VariableFeatures(object = LHA_LHb_patch.se))
VizDimLoadings(LHA_LHb_patch.se, dims = 1:2, reduction = "pca")
ElbowPlot(LHA_LHb_patch.se, ndims = 20)

# Clustering the cells
LHA_LHb_patch.se <- FindNeighbors(LHA_LHb_patch.se, dims = 1:10)
LHA_LHb_patch.se <- FindClusters(LHA_LHb_patch.se, resolution = 0.5)

# Run t-SNE
LHA_LHb_patch.se <- RunTSNE(LHA_LHb_patch.se, dims = 1:10)
DimPlot(LHA_LHb_patch.se, reduction = "tsne")

# Save obj to RDS
saveRDS(LHA_LHb_patch.se, "../RDS_obj/LHA_LHb_patch.RDS")

# Find markers for every cluster compared to all remaining cells, report only the positive ones
LHA_LHb_patch.markers <- FindAllMarkers(LHA_LHb_patch.se, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top2.markers <- LHA_LHb_patch.markers %>% 
  group_by(cluster) %>% 
  slice_max(n = 2, order_by = avg_log2FC)

# Visualize the expression of the top 2 features with VlnPlot and FeaturePlot
VlnPlot(LHA_LHb_patch.se, features = c(top2.markers$gene))
FeaturePlot(LHA_LHb_patch.se, features = c(top2.markers$gene)) + DimPlot(LHA_LHb_patch.se, reduction = "tsne")

# Make a heatmap
LHA_LHb_patch.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(LHA_LHb_patch.se, features = top10$gene) + NoLegend()


