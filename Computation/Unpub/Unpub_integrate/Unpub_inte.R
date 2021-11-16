library(Seurat)
library(SeuratData)
library(patchwork)


# Read RDS of unpub data into a list
RDS.list <- list()
RDS <- list.files("/Users/pax/Google\ Drive/My\ Drive/Meletis/Computation/Unpub/RDS_obj", pattern = "*.RDS")

for (i in seq_along(RDS)){
 RDS.list[i] <- readRDS(paste0("/Users/pax/Google\ Drive/My\ Drive/Meletis/Computation/Unpub/RDS_obj/", RDS[i]))
}

names(RDS.list) <- substr(RDS, 1, nchar(RDS)-4)

# Identify variable features for each dataset independently
RDS.list <- lapply(X = RDS.list, FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


# Select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = RDS.list)

# Find integration anchors
unpub.anchors <- FindIntegrationAnchors(object.list = RDS.list, anchor.features = features)

# Creates an 'integrated' data assay
unpub.combined <- IntegrateData(anchorset = unpub.anchors)

# Perform integration analysis
DefaultAssay(unpub.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
unpub.combined <- ScaleData(unpub.combined, verbose = FALSE)
unpub.combined <- RunPCA(unpub.combined, npcs = 50, verbose = FALSE)
unpub.combined <- RunUMAP(unpub.combined, reduction = "pca", dims = 1:25)
unpub.combined <- FindNeighbors(unpub.combined, reduction = "pca", dims = 1:25)
unpub.combined <- FindClusters(unpub.combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(unpub.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(unpub.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

# Save RDS
saveRDS(unpub.combined, "./RDS_obj/Unpub_integration.RDS")
