library(dplyr)
library(Seurat)
library(patchwork)

#Read count matrix of AJ17001 with Read10X
AJ17001.data <- Read10X(data.dir = "~/Google\ Drive/My\ Drive/Meletis/Data_reproduce/GSE125065/GSE125065_RAW/GSM3562050_AJ17001/")
AJ17002.data <- Read10X(data.dir = "~/Google\ Drive/My\ Drive/Meletis/Data_reproduce/GSE125065/GSE125065_RAW/GSM3562051_AJ17002/")

#Create SEurat obj
AJ17001.se <- CreateSeuratObject(AJ17001.data, project = "Female", min.cells = 5, min.features = 200)
AJ17002.se <- CreateSeuratObject(AJ17002.data, project = "Male", min.cells = 5, min.features = 200)

#Add sex to the  obj
AJ17001.se[['sex']] <- "F"
AJ17002.se[['sex']] <- "M"

#Merge Seurat obj
LHA.combined <- merge(AJ17001.se, y = AJ17002.se, add.cell.ids = c("F", "M"), project = "LHA_combined")
LHA.combined

#Look for mitochodrial genes in the count matrix
LHA.combined[["percent.mt"]] <- PercentageFeatureSet(LHA.combined, pattern = "^mt-")

#Visualize QC metrics as a violin plot
VlnPlot(LHA.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Visualize feature-feature relationships with FeatureScatter
plot1 <- FeatureScatter(LHA.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(LHA.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Subset the obj
LHA.combined <- subset(LHA.combined, subset = nCount_RNA > 500 & nCount_RNA < 40000 & nFeature_RNA > 600 & nFeature_RNA < 6500  & percent.mt < 40)

#Normalized the gene expression by the total number of transcripts detected in each cell and multiplied by the median transcript count
LHA.combined$nFeature_RNA <- (LHA.combined$nFeature_RNA/sum(LHA.combined$nCount_RNA))*median(LHA.combined$nCount_RNA)

#Log transform
LHA.combined <- NormalizeData(LHA.combined)

#Find top 1000 variable genes
LHA.combined <- FindVariableFeatures(LHA.combined, selection.method = "vst", nfeatures = 1000)

#Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(LHA.combined), 10)

#Plot variable features with and without labels
plot1 <- VariableFeaturePlot(LHA.combined)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Scale data
all.genes <- rownames(LHA.combined)
LHA.combined <- ScaleData(LHA.combined, features = all.genes)

#Run PCA
LHA.combined <- RunPCA(LHA.combined, npcs = 50, features = VariableFeatures(object = LHA.combined))
VizDimLoadings(LHA.combined, dims = 1:2, reduction = "pca")
ElbowPlot(LHA.combined)

#Find Clusters
LHA.combined <- FindNeighbors(LHA.combined, dims = 1:15)
LHA.combined <- FindClusters(LHA.combined, resolution = 0.5)

#Run and visualize t-SNE
LHA.combined <- RunTSNE(LHA.combined, dims = 1:15)
DimPlot(LHA.combined, reduction = "tsne")

#Save the LHA.combined obj
saveRDS(LHA.combined, file = "./Mickelsen_merge.rds")

#Find markers in every cluster
LHA.combined.markers <- FindAllMarkers(LHA.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#Find top1 markers in each cluster
top1.markers <- LHA.combined.markers %>%
  group_by(cluster) %>%
  slice_max(n = 1, order_by = avg_log2FC)

#Visualize top1 markers in each cluster in violin plot
VlnPlot(LHA.combined, features = c(top1.markers$gene))

#Find the median of neuron markers
neuron.markers <- AverageExpression(LHA.combined, features = c("Snap25", "Syp", "Tubb3", "Elavl2"))

#Create a list of cell names that are neurons
neuron.list <- WhichCells(LHA.combined, idents = c(0,1,3,4,5,7,9,13,15))

#Add cell_type to metadata. Every cell is classified as non-neurons at the beginning
LHA.combined[["cell_type"]] <- "Non-neurons"

#Classify cell with names in neuron.list to be neurons
LHA.combined$cell_type[which(names(LHA.combined$nCount_RNA) %in% neuron.list)] <- "Neurons"

#Visualize cell type in t-SNE with DimPlot
DimPlot(LHA.combined, reduction = "tsne", group.by = "cell_type")

#Save the LHA.combined obj
saveRDS(LHA.combined, file = "./Mickelsen_merge.rds")


