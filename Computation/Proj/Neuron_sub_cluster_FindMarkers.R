library(dplyr)
library(Seurat)
library(patchwork)
library(SeuratObject)
library(future)

plan("multicore", workers = 4)
options(future.globals.maxSize = 7000 * 1024^2)

# Read RDS of neuron.subset and neuron.subset.marker
neuron.subset <- readRDS("./RDS_obj/Proj_combined_neuron_sub.RDS")
neuron.subset.marker <- readRDS("./RDS_obj/Proj_combined_neuron_sub_markers.RDS")

# Get top 2 makers of each cluster
neuron.top2.markers <- neuron.subset.marker %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# Remove cluster 0 and 2 as they have high expression of mt genes. Thus, they are most likely to be low quality cells
neuron.subset <- subset(neuron.subset, idents = c(0, 2), invert = TRUE)

# Re-cluster the neuron.subset after removing cluster 0 and 2
neuron.subset <- FindVariableFeatures(neuron.subset)
DefaultAssay(neuron.subset) <- "integrated"
neuron.subset <- ScaleData(neuron.subset)
neuron.subset <- RunPCA(neuron.subset, npcs = 200)
neuron.subset <- RunUMAP(neuron.subset, reduction = "pca", dims = 1:70, n.neighbors = 40, n.epochs = 500, min.dist = 0.5, spread = 0.4, negative.sample.rate = 10)
neuron.subset <- FindNeighbors(neuron.subset, reduction = "pca", dims = 1:70)
neuron.subset <- FindClusters(neuron.subset, resolution = 0.5)

# Select RNA as default assay for neuron.subset to find markers in each cluster
DefaultAssay(neuron.subset) <- "RNA"
neuron.subset.marker <- FindAllMarkers(neuron.subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Find the median of neuron type (GABAergic or GLUT) markers
neuron.type.median <- as.data.frame(AverageExpression(neuron.subset, features = c("Slc32a1", "Slc17a6"))$RNA)
GABA.ident <- as.integer(colnames(neuron.type.median)[neuron.type.median["Slc32a1",] > neuron.type.median["Slc17a6",]])
GLUT.ident <- as.integer(colnames(neuron.type.median)[neuron.type.median["Slc32a1",] < neuron.type.median["Slc17a6",]])

# Create lists of cell names that are GABA and GLUT
GABA.list <- WhichCells(neuron.subset, idents = c(GABA.ident))
GLUT.list <- WhichCells(neuron.subset, idents = c(GLUT.ident))

# Add neuron_type to metadata
neuron.subset[["neuron_type"]] <- "None"

# Classify cell with names in GABA.list to be GABA
neuron.subset$neuron_type[which(names(neuron.subset$nCount_RNA) %in% GABA.list)] <- "GABA"

# Classify cell with names in GLUT.list to be GLUT
neuron.subset$neuron_type[which(names(neuron.subset$nCount_RNA) %in% GLUT.list)] <- "GLUT"

DimPlot(neuron.subset, reduction = "umap", group.by = "neuron_type") + DimPlot(neuron.subset, reduction = "umap", label = T) + DimPlot(neuron.subset, reduction = "umap", group.by = "source")

# Rename GABA idents
neuron.subset <- RenameIdents(object = neuron.subset, 
                              "0" = "GABA1",
                              "1" = "GABA2",
                              "2" = "GABA3",
                              "3" = "GABA4",
                              "4" = "GABA5",
                              "11" = "GABA6",
                              "12" = "GABA7",
                              "15" = "GABA8",
                              "16" = "GABA9",
                              "17" = "GABA10", 
                              "20" = "GABA11", 
                              "22" = "GABA12",
                              "26" = "GABA13",
                              "29" = "GABA14",
                              "30" = "GABA15",
                              "32" = "GABA16",
                              "33" = "GABA17",
                              "35" = "GABA18",
                              "36" = "GABA19",
                              "37" = "GABA20")

# Rename GLUT idents
neuron.subset <- RenameIdents(object = neuron.subset, 
                              "5" = "GLUT1",
                              "6" = "GLUT2",
                              "7" = "GLUT3",
                              "8" = "GLUT4",
                              "9" = "GLUT5",
                              "10" = "GLUT6",
                              "13" = "GLUT7",
                              "14" = "GLUT8",
                              "18" = "GLUT9",
                              "19" = "GLUT10", 
                              "21" = "GLUT11", 
                              "23" = "GLUT12",
                              "24" = "GLUT13",
                              "25" = "GLUT14",
                              "27" = "GLUT15",
                              "28" = "GLUT16",
                              "31" = "GLUT17",
                              "34" = "GLUT18")

# Save RDS of named neuron.subset
saveRDS(neuron.subset, "./RDS_obj/Proj_combined_neuron_named.RDS")