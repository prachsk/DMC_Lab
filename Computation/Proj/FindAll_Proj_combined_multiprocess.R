library(dplyr)
library(Seurat)
library(patchwork)
library(future)

plan("multicore", workers = 4)
options(future.globals.maxSize = 7000 * 1024^2)

# Read RDS
proj.combined <- readRDS("./RDS_obj/Proj_combined.RDS")

# Run the standard workflow for visualization and clustering
proj.combined <- ScaleData(proj.combined, verbose = FALSE)
proj.combined <- RunPCA(proj.combined, npcs = 200, verbose = FALSE)
proj.combined <- RunUMAP(proj.combined, reduction = "pca", dims = 1:70, n.neighbors = 40, n.epochs = 500, min.dist = 0.5, spread = 0.4, negative.sample.rate = 10)
proj.combined <- FindNeighbors(proj.combined, reduction = "pca", dims = 1:70)
proj.combined <- FindClusters(proj.combined, resolution = 0.5)

# Set Default assay to RNA
DefaultAssay(proj.combined) <- "RNA"

# Find cluster markers
proj.combined.markers <- FindAllMarkers(proj.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

saveRDS(proj.combined.markers, "./RDS_obj/Proj_combined_markers.RDS" )
saveRDS(proj.combined, "./RDS_obj/Proj_combined_FindAll.RDS")