library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)
library(future)

plan("multiprocess", workers = 4)
options(future.globals.maxSize = 7000 * 1024^2)

# Read RDS of pub and unpub data
chen <- readRDS('/Users/pax/Google\ Drive/My\ Drive/Meletis/Computation/Chen_et_al/RDS_obj/Chen_merge.RDS')
mickelsen <- readRDS('/Users/pax/Google\ Drive/My\ Drive/Meletis/Computation/Mickelsen_et_al/RDS_obj/Mickelsen_merge.RDS')
rossi.2019 <- readRDS('/Users/pax/Google\ Drive/My\ Drive/Meletis/Computation/Rossi_et_al_2019/Result_RDS/dge.merge.RDS')
rossi.2021 <- readRDS('/Users/pax/Google\ Drive/My\ Drive/Meletis/Computation/Rossi_et_al_2021/RDS_obj/MLB017_processed.RDS')
unpub <- readRDS('/Users/pax/Google\ Drive/My\ Drive/Meletis/Computation/Unpub/Unpub_integrate/RDS_obj/Unpub_integration.RDS')

# Merge all Seurat object into one
proj.merge <- merge(x = unpub, y = c(chen, mickelsen, rossi.2019, rossi.2021))

# Selecting RNA as default assay
DefaultAssay(proj.merge) <- "RNA"

# remove the obj loaded from RDS
rm(chen, mickelsen, rossi.2019, rossi.2021, unpub)

# Run standard Seurat workflow
proj.merge <- NormalizeData(proj.merge)
proj.merge <- FindVariableFeatures(proj.merge, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(proj.merge)
proj.merge <- ScaleData(proj.merge)
proj.merge <- RunPCA(proj.merge, npcs = 200, features = VariableFeatures(object = proj.merge))

# Run harmony
proj.harmony <- proj.merge %>% 
  RunHarmony("neuron_type", plot_convergence = TRUE)

# Downstream analysis of harmony
proj.harmony <- proj.harmony %>% 
  RunUMAP(reduction = "harmony", dims = 1:70, n.neighbors = 40, n.epochs = 500, min.dist = 0.5, spread = 0.4, negative.sample.rate = 10) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:70) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

DimPlot(proj.harmony, reduction = "umap", group.by = "source")

# Find cluster markers
proj.markers <- FindAllMarkers(proj.harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
proj.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# Save RDS of proj.markers and proj.harmoony
saveRDS(proj.markers, "./RDS_obj/Harmony_markers.RDS")
saveRDS(proj.harmony, "./RDS_obj/Harmony_pub_unpub.RDS")

# Read RDS
#proj.markers <- readRDS("./RDS_obj/Harmony_markers_neurons.RDS")
proj.harmony <- readRDS("./RDS_obj/Harmony_pub_unpub.RDS")

# Add Cell type metadata
# Find the median of neuron markers
#neuron.markers <- AverageExpression(proj.harmony, features = c("Snap25", "Syp", "Tubb3", "Elavl2"))

#Create a list of cell names that are neurons
#neuron.list <- WhichCells(proj.merge, idents = c(0,1,3,4,5,7,9,13,15))

#Add cell_type to metadata. Every cell is classified as non-neurons at the beginning
#proj.merge[["cell_type"]] <- "Non-neurons"

#Classify cell with names in neuron.list to be neurons
#proj.merge$cell_type[which(names(proj.merge$nCount_RNA) %in% neuron.list)] <- "Neurons"







