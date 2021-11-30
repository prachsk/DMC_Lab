library(dplyr)
library(patchwork)
library(ggplot2)
library(scales)
library(Seurat)
library(SeuratObject)

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
  #x <- SCTransform(x, vars.to.regress = c("percent.mt"))
})

# Select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = proj.list, nfeatures = 2000)

# Find integration anchors
proj.anchors <- FindIntegrationAnchors(object.list = proj.list, anchor.features = features)

# Creates an 'integrated' data assay
proj.combined <- IntegrateData(anchorset = proj.anchors)

# Save RDS
saveRDS(proj.combined, "./RDS_obj/Proj_combined.RDS")

# Read RDS
proj.combined <- readRDS("./RDS_obj/Proj_combined.RDS")

# Select integrated as default assay
DefaultAssay(proj.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
proj.combined <- ScaleData(proj.combined)
proj.combined <- RunPCA(proj.combined, npcs = 200)
proj.combined <- RunUMAP(proj.combined, reduction = "pca", dims = 1:70, n.neighbors = 40, n.epochs = 500, min.dist = 0.5, spread = 0.4, negative.sample.rate = 10)
proj.combined <- FindNeighbors(proj.combined, reduction = "pca", dims = 1:70)
proj.combined <- FindClusters(proj.combined, resolution = 0.5)

# Visualize with DimPlot
DimPlot(proj.combined, reduction = "umap")

# Set Default assay to RNA
DefaultAssay(proj.combined) <- "RNA"

# Find cluster markers
proj.combined.markers <- FindAllMarkers(proj.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Save RDS
saveRDS("./RDS_obj/Proj_combined_markers.RDS")
saveRDS("./RDS_obj/Proj_combined_FindAll.RDS")

# Read RDS
proj.combined.markers <- readRDS("./RDS_obj/Proj_combined_markers.RDS")
proj.combined <- readRDS("./RDS_obj/Proj_combined_FindAll.RDS")

# Get top 2 makers of each cluster
top2.markers <- proj.combined.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# Identify neurons, non-neurons, and none
# Find the median of neuron markers
neuron.markers <- as.data.frame(AverageExpression(proj.combined, features = c("Snap25", "Syp", "Tubb3", "Elavl2"))$RNA)
neuron.idents <- as.integer(colnames(neuron.markers)[colSums(neuron.markers)>15])
non_neuron.idents <- as.integer(colnames(neuron.markers)[colSums(neuron.markers)<5])

# Create lists of cell names that are neurons and non-neurons
neuron.list <- WhichCells(proj.combined, idents = c(neuron.idents))
non_neuron.list <- WhichCells(proj.combined, idents = c(non_neuron.idents))

# Add cell_type to metadata. Every cell is classified as non-neurons at the beginning
proj.combined[["cell_type"]] <- "None"

# Classify cell with names in neuron.list to be neurons
proj.combined$cell_type[which(names(proj.combined$nCount_RNA) %in% neuron.list)] <- "Neurons"

# Classify cell with names in mom_neuron.list to be non-neurons
proj.combined$cell_type[which(names(proj.combined$nCount_RNA) %in% non_neuron.list)] <- "Non-neurons"

# Visualize UMAP with DimPlot
DimPlot(proj.combined, reduction = "umap", group.by = "cell_type") + DimPlot(proj.combined, reduction = "umap", label = TRUE) + DimPlot(proj.combined, reduction = "umap", group.by = "source")

# Suspect that cells with "None" in cell_type are immune cell
# Use VlnPlot to visualize key markers of cluster 5 (almost all "None" cells are in cluster 5)
VlnPlot(proj.combined, features = c("Pdgfra", "C1ql1"))

# Rename "None" to "Immune cells"
proj.combined$cell_type[proj.combined$cell_type == "None"] <- "Immune cells"

# Remove cell with cell_type of Immune cells and visualize wiht DimPlot
proj.combined.subset <- subset(proj.combined, subset = cell_type == "Immune cells", invert = TRUE)
DimPlot(proj.combined, reduction = "umap", group.by = "cell_type") + DimPlot(proj.combined, reduction = "umap", label = T) + DimPlot(proj.combined.subset, reduction = "umap", group.by = "source")


# Subset for only neuron cells
neuron.subset <- subset(proj.combined.subset, subset = cell_type == "Neurons")

# Save neuron.subset RDS
saveRDS(neuron.subset, "./RDS_obj/Proj_combined_neuron_sub.RDS")

# Read RDS of neuron.subset
neuron.subset <- readRDS("./RDS_obj/Proj_combined_neuron_sub.RDS")

# Seurat analysis pipeline on neuron.subset
neuron.subset <- FindVariableFeatures(neuron.subset)
DefaultAssay(neuron.subset) <- "integrated"
neuron.subset <- ScaleData(neuron.subset)
neuron.subset <- RunPCA(neuron.subset, npcs = 200)
neuron.subset <- RunUMAP(neuron.subset, reduction = "pca", dims = 1:70, n.neighbors = 40, n.epochs = 500, min.dist = 0.5, spread = 0.4, negative.sample.rate = 10)
neuron.subset <- FindNeighbors(neuron.subset, reduction = "pca", dims = 1:70)
neuron.subset <- FindClusters(neuron.subset, resolution = 0.5)

# Visualize UMAP with DimPlot
DimPlot(neuron.subset, reduction = "umap", group.by = "source") + DimPlot(neuron.subset, reduction = "umap")

# Save RDS of neuron.subset
saveRDS(neuron.subset, "./RDS_obj/Proj_combined_neuron_sub.RDS")

# Read RDS of neuron.subset
neuron.subset <- readRDS("./RDS_obj/Proj_combined_neuron_sub.RDS")

# Select RNA as default assay for neuron.subset to find markers in each cluster
DefaultAssay(neuron.subset) <- "RNA"
neuron.subset.marker <- FindAllMarkers(neuron.subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Save RDS of neuron.subset
saveRDS(neuron.subset, "./RDS_obj/Proj_combined_neuron_sub.RDS")
saveRDS(neuron.subset.marker, "./RDS_obj/Proj_combined_neuron_sub_markers.RDS")

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

DimPlot(neuron.subset, reduction = "umap", group.by = "neuron_type") + DimPlot(neuron.subset, reduction = "umap", label = T) + DimPlot(neuron.subset, reduction = "umap", group.by = "orig.ident")

# Save RDS of neuron.subset
saveRDS(neuron.subset, "./RDS_obj/Proj_combined_neuron_sub.RDS")

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

# Select RNA as default assay for neuron.subset to find markers in each cluster
DefaultAssay(neuron.subset) <- "RNA"
neuron.subset.marker <- FindAllMarkers(neuron.subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Get top 2 makers of each cluster
neuron.top2.markers <- neuron.subset.marker %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# Save RDS of named neuron.subset
saveRDS(neuron.subset, "./RDS_obj/Proj_combined_neuron_named.RDS")

# Read RDS
neuron.subset <- readRDS("./RDS_obj/Proj_combined_neuron_named.RDS")

Hb.marker.median <- as.data.frame(AverageExpression(neuron.subset, features = c("Tac2", "Pcdh10"))$RNA)
Hb.ident <- colnames(Hb.marker.median)[colSums(Hb.marker.median) > 5]
Hb.list <- WhichCells(neuron.subset, idents = c(Hb.ident))

LHA.marker.median <- as.data.frame(AverageExpression(neuron.subset, features = c("Pax6","Pdyn","Hcrt"))$RNA)
LHA.ident <- colnames(LHA.marker.median)[colSums(LHA.marker.median) > 10]
LHA.list <- WhichCells(neuron.subset, idents = c(LHA.ident))

# FeaturePlot of key expression marker of LHA projecting to LHb (Pax6) and VTA (Pdyn, Hcrt)
FeaturePlot(neuron.subset, features = c("Slc32a1","Gad1","Gad2","Slc17a6")) + DimPlot(neuron.subset, reduction = "umap", group.by = "neuron_type")

# Subset GLUT cells and run Seurat workflow
GLUT.sub <- subset(neuron.subset, subset = neuron_type == "GLUT")
GLUT.sub <- FindVariableFeatures(GLUT.sub)
DefaultAssay(GLUT.sub) <- "integrated"
GLUT.sub <- ScaleData(GLUT.sub)
GLUT.sub <- RunPCA(GLUT.sub, npcs = 200)
GLUT.sub <- RunUMAP(GLUT.sub, reduction = "pca", dims = 1:50)
GLUT.sub <- FindNeighbors(GLUT.sub, reduction = "pca", dims = 1:50)
GLUT.sub <- FindClusters(GLUT.sub, resolution = 0.19)
DimPlot(GLUT.sub, reduction = 'umap')

DefaultAssay(GLUT.sub) <- "RNA"
GLUT.sub.markers <- FindAllMarkers(GLUT.sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Get top 2 makers of each cluster
GLUT.sub.top2.markers <- GLUT.sub.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# Remove cluster 13 as it is low quality cells
GLUT.sub <- subset(GLUT.sub, idents = 13, invert = TRUE)

# Re-clustering GLUT.sub with Seurat workflow
GLUT.sub <- FindVariableFeatures(GLUT.sub)
DefaultAssay(GLUT.sub) <- "integrated"
GLUT.sub <- ScaleData(GLUT.sub)
GLUT.sub <- RunPCA(GLUT.sub, npcs = 200)
GLUT.sub <- RunUMAP(GLUT.sub, reduction = "pca", dims = 1:50)
GLUT.sub <- FindNeighbors(GLUT.sub, reduction = "pca", dims = 1:50)
GLUT.sub <- FindClusters(GLUT.sub, resolution = 0.19)
DimPlot(GLUT.sub, reduction = 'umap')

# Find marker of each cluster in GLUT.sub
DefaultAssay(GLUT.sub) <- "RNA"
GLUT.sub.markers <- FindAllMarkers(GLUT.sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Get top 2 makers of each cluster
GLUT.sub.top2.markers <- GLUT.sub.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# Save Old identity
GLUT.sub[["old.ident"]] <- Idents(GLUT.sub)

GLUT.sub <- RenameIdents(GLUT.sub, 
                         "0" = "Tac1 + Trh",
                         "1" = "Nts + Tcf7l2",
                         "2" = "Tac1 + Pitx2",
                         "3" = "Hcrt + Scg2",
                         "4" = "Penk + Otp",
                         "5" = "Meis2 + Rsrp1",
                         "6" = "Dlk1 + Asb4",
                         "7" = "Bc1 + Pax6",
                         "8" = "Pmch + Cartpt",
                         "9" = "Synpr + Hpcal1",
                         "10" = "Foxb1 + Pitx2",
                         "11" = "Apoe + Atp1a2",
                         "12" = "Trh + Ghrh",
                         "13" = "Tcf4 + Calca",
                         "14" = "Nppc + Nrn1",
                         "15" = "Tcf4 + Grp",
                         "16" = "Avp + Oxt",
                         "17" = "Ucn3 + Trh"
                         )

# Save current active ident as subtypes
GLUT.sub[["subtypes"]] <- Idents(GLUT.sub)

# Save RDS GLUT.sub
saveRDS(GLUT.sub, "./RDS_obj/GLUT_sub_named.RDS")

# Read GLUT.sub RDS
GLUT.sub <- readRDS("./RDS_obj/GLUT_sub_named.RDS")

# DimPlot of GLUT subtypes embedded in UMAP
DimPlot(GLUT.sub, reduction = 'umap', label = T) + labs(title = "UMAP of GLUT Subtypes") + NoLegend()

# # Visualize markers of GLUT clusters with stack VlnPlot
VlnPlot(GLUT.sub , features = c(unique(GLUT.sub.top2.markers$gene)), 
        stack = T, group.by = "old.ident", split.by = "subtypes", flip = T) + 
  NoLegend() + labs(title = "Markers of GLUT clusters")

# Subset GABA cells and run Seurat workflow
GABA.sub <- subset(neuron.subset, subset = neuron_type == "GABA")
GABA.sub <- FindVariableFeatures(GABA.sub)
DefaultAssay(GABA.sub) <- "integrated"
GABA.sub <- ScaleData(GABA.sub)
GABA.sub <- RunPCA(GABA.sub, npcs = 200)
GABA.sub <- RunUMAP(GABA.sub, reduction = "pca", dims = 1:50)
GABA.sub <- FindNeighbors(GABA.sub, reduction = "pca", dims = 1:50)
GABA.sub <- FindClusters(GABA.sub, resolution = 0.25)
DimPlot(GABA.sub, reduction = 'umap')

DefaultAssay(GABA.sub) <- "RNA"
GABA.sub.markers <- FindAllMarkers(GABA.sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Get top 2 makers of each cluster
GABA.sub.top2.markers <- GABA.sub.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# Remove cluster 14 as it is low quality cells
GABA.sub <- subset(GABA.sub, idents = 14, invert = TRUE)

# # Re-clustering GABA.sub with Seurat workflow
GABA.sub <- FindVariableFeatures(GABA.sub)
DefaultAssay(GABA.sub) <- "integrated"
GABA.sub <- ScaleData(GABA.sub)
GABA.sub <- RunPCA(GABA.sub, npcs = 200)
GABA.sub <- RunUMAP(GABA.sub, reduction = "pca", dims = 1:50)
GABA.sub <- FindNeighbors(GABA.sub, reduction = "pca", dims = 1:50)
GABA.sub <- FindClusters(GABA.sub, resolution = 0.25)
DimPlot(GABA.sub, reduction = 'umap')

# Find marker of each cluster in GABA.sub
DefaultAssay(GABA.sub) <- "RNA"
GABA.sub.markers <- FindAllMarkers(GABA.sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Get top 2 makers of each cluster
GABA.sub.top2.markers <- GABA.sub.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# Save Old identity
GABA.sub[["old.ident"]] <- Idents(GABA.sub)

GABA.sub <- RenameIdents(GABA.sub,
                         "0" = "Tac2 + A030009H04Rik",
                         "1" = "Nts + Cartpt",
                         "2" = "Tac1 + Meis2",
                         "3" = "Ddc + Sncg",
                         "4" = "Gal + Lmo3",
                         "5" = "Sst + Meis2",
                         "6" = "Tcf4 + Synpr",
                         "7" = "Gm42418 + 2900097C17Rik",
                         "8" = "Nnat + C1ql2",
                         "9" = "Tcf7l2 + Calb1",
                         "10" = "Ebf1 + Nkx2-4",
                         "11" = "Plp1 + Ptgds",
                         "12" = "Sst + Arpp21",
                         "13" = "Fam81a + Rprm",
                         "14" = "Pmch + Olig1",
                         "15" = "Cck + Vtn",
                         "16" = "Cst3 + C1qa",
                         "17" = "Slc5a7 + Acly",
                         "18" = "Calb2 + Synpr",
                         "19" = "Hdc + Slc18a2"
                         )

# Save current active ident as subtypes
GABA.sub[["subtypes"]] <- Idents(GABA.sub)

# Save RDS GABA.sub
saveRDS(GABA.sub, "./RDS_obj/GABA_sub_named.RDS")

# Read GABA.sub RDS
GABA.sub <- readRDS("./RDS_obj/GABA_sub_named.RDS")

# DimPlot of GABA subtypes embedded in UMAP
DimPlot(GABA.sub, reduction = 'umap', label = T) + labs(title = "UMAP of GABA Subtypes") + NoLegend()

# Visualize markers of GABA clusters with stack VlnPlot
VlnPlot(GABA.sub , features = c(unique(GABA.sub.top2.markers$gene)),
        stack = T, group.by = "old.ident", split.by = "subtypes", flip = T) + 
  NoLegend() + labs(title = "Markers of GABA clusters")

GLUTGABA.merge <- merge(GLUT.sub, y = GABA.sub,)




#FeaturePlot(neuron.subset, features = c("Slc32a1","Slc17a6"))

#GLUT.features <-c('Pmch','Nrgn','Zic1','Tac1','Ebf3','Hcrt','Gpr101','Trh','Synpr','Grp','Col27a1','Syt2','Otp','Sst','Cartpt','Gda',
                  #'Pitx2','Pdyn','Tcf4','Onecut2','Cbln2','Cck','Cbln1','Meis2','Nkx2-1')

#VlnPlot(neuron.subset, features = c(unique(neuron.top2.markers$gene)), stack = T) + RotatedAxis()

#DotPlot(neuron.subset, features = c(unique(neuron.top2.markers$gene))) + RotatedAxis()


