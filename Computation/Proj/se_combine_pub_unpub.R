library(dplyr)
library(Seurat)
library(patchwork)

# Read RDS of pub and unpub data
chen <- readRDS('/Users/pax/Google\ Drive/My\ Drive/Meletis/Computation/Chen_et_al/RDS_obj/Chen_merge.RDS')
mickelsen <- readRDS('/Users/pax/Google\ Drive/My\ Drive/Meletis/Computation/Mickelsen_et_al/RDS_obj/Mickelsen_merge.RDS')
rossi.2019 <- readRDS('/Users/pax/Google\ Drive/My\ Drive/Meletis/Computation/Rossi_et_al_2019/Result_RDS/dge.merge.RDS')
rossi.2021 <- readRDS('/Users/pax/Google\ Drive/My\ Drive/Meletis/Computation/Rossi_et_al_2021/RDS_obj/MLB017_processed.RDS')
unpub <- readRDS('/Users/pax/Google\ Drive/My\ Drive/Meletis/Computation/Unpub/Unpub_integrate/RDS_obj/Unpub_integration.RDS')

# Create list of multiple Seurat object
proj.list <- list(unpub, chen, mickelsen, rossi.2019, rossi.2021)

# Select RNA as default assay, normalize, and identify variable features for each dataset independently
proj.list <- lapply(X = proj.list, FUN = function(x) {
  DefaultAssay(x) <- "RNA"
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = proj.list)

# Find integration anchors
proj.anchors <- FindIntegrationAnchors(object.list = proj.list, anchor.features = features)

# Creates an 'integrated' data assay
proj.combined <- IntegrateData(anchorset = proj.anchors)

# Select integrated as default assay
DefaultAssay(proj.combined) <- "integrated"

# Save RDS
saveRDS(proj.combined, "./RDS_obj/Proj_combined.RDS")

# Load RDS
proj.combined <- readRDS("./RDS_obj/Proj_combined.RDS")

# Run the standard workflow for visualization and clustering
proj.combined <- ScaleData(proj.combined, verbose = FALSE)
proj.combined <- RunPCA(proj.combined, npcs = 200, verbose = FALSE)
proj.combined <- RunUMAP(proj.combined, reduction = "pca", dims = 1:70)
proj.combined <- FindNeighbors(proj.combined, reduction = "pca", dims = 1:70)
proj.combined <- FindClusters(proj.combined, resolution = 0.5)

DimPlot(proj.combined, reduction = "umap", group.by = "stim")








