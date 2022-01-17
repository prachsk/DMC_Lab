library(dplyr)
library(patchwork)
library(ggplot2)
library(scales)
library(Seurat)
library(SeuratObject)
library(dittoSeq)
library(Nebulosa)

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

# Assign new Sources
proj.combined[["Source"]] <- proj.combined$source
unpub.LHA <- rownames(proj.combined@meta.data)[proj.combined$orig.ident == "LHA" & proj.combined$source == "Unpub"]
unpub.patch <- rownames(proj.combined@meta.data)[proj.combined$orig.ident == "LHA_LHb_patch" & proj.combined$source == "Unpub"]
unpub.LHb <- rownames(proj.combined@meta.data)[proj.combined$orig.ident == "LHb_lskra" & proj.combined$source == "Unpub"]
unpub.HSV <- rownames(proj.combined@meta.data)[proj.combined$orig.ident == "HSV" & proj.combined$source == "Unpub"]

# Classify cell according to the unpub areas
proj.combined$Source[which(names(proj.combined$nCount_RNA) %in% unpub.LHA)] <- "Unpub_LHA"
proj.combined$Source[which(names(proj.combined$nCount_RNA) %in% unpub.patch)] <- "Unpub_LHA_LHb_patch"
proj.combined$Source[which(names(proj.combined$nCount_RNA) %in% unpub.LHb)] <- "Unpub_LHb"
proj.combined$Source[which(names(proj.combined$nCount_RNA) %in% unpub.HSV)] <- "Unpub_HSV"

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

# FeaturePlot of "Pdgfra" & "C1ql1"
FeaturePlot(proj.combined, features = c("Pdgfra", "C1ql1"))

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

DimPlot(neuron.subset, reduction = "umap", group.by = "neuron_type") + DimPlot(neuron.subset, reduction = "umap", label = T) + DimPlot(neuron.subset, reduction = "umap", group.by = "Source")

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

neuron.subset[["subtypes"]] <- Idents(neuron.subset)

# UMAP plots of cell from different sources in neuron.subset
p1 <- DimPlot(neuron.subset, cells.highlight = rownames(neuron.subset@meta.data)[neuron.subset$Source == "Unpub_LHA"], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Unpub_LHA") + NoLegend()
p2 <- DimPlot(neuron.subset, cells.highlight = rownames(neuron.subset@meta.data)[neuron.subset$Source == "Unpub_LHb"], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Unpub_LHb") + NoLegend()
p3 <- DimPlot(neuron.subset, cells.highlight = rownames(neuron.subset@meta.data)[neuron.subset$Source == "Unpub_HSV"], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Unpub_HSV") + NoLegend()
p4 <- DimPlot(neuron.subset, cells.highlight = rownames(neuron.subset@meta.data)[neuron.subset$Source == "Unpub_LHA_LHb_patch"], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Unpub_LHA_LHb_patch") + NoLegend()
p5 <- DimPlot(neuron.subset, cells.highlight = rownames(neuron.subset@meta.data)[neuron.subset$Source == "Chen"], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Chen") + NoLegend()
p6 <- DimPlot(neuron.subset, cells.highlight = rownames(neuron.subset@meta.data)[neuron.subset$Source == "Mickelsen"], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Mickelsen") + NoLegend()
p7 <- DimPlot(neuron.subset, cells.highlight = rownames(neuron.subset@meta.data)[neuron.subset$Source == "Rossi_2019"], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Rossi_2019") + NoLegend()
p8 <- DimPlot(neuron.subset, cells.highlight = rownames(neuron.subset@meta.data)[neuron.subset$Source == "Rossi_2021"], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Rossi_2021") + NoLegend()
p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + (DimPlot(neuron.subset, reduction = "umap", label = T) + NoLegend())

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
                         "0" = "GLUT1",
                         "1" = "GLUT2",
                         "2" = "GLUT3",
                         "3" = "GLUT4",
                         "4" = "GLUT5",
                         "5" = "GLUT6",
                         "6" = "GLUT7",
                         "7" = "GLUT8",
                         "8" = "GLUT9",
                         "9" = "GLUT10",
                         "10" = "GLUT11",
                         "11" = "GLUT12",
                         "12" = "GLUT13",
                         "13" = "GLUT14",
                         "14" = "GLUT15",
                         "15" = "GLUT16",
                         "16" = "GLUT17",
                         "17" = "GLUT18"
                         )

# Save current active ident as subtypes
GLUT.sub[["subtypes"]] <- Idents(GLUT.sub)

GLUT.sub <- RenameIdents(GLUT.sub, 
                         "GLUT1" = "Tac1 + Trh",
                         "GLUT2" = "Nts + Tcf7l2",
                         "GLUT3" = "Tac1 + Pitx2",
                         "GLUT4" = "Hcrt + Scg2",
                         "GLUT5" = "Penk + Otp",
                         "GLUT6" = "Meis2 + Rsrp1",
                         "GLUT7" = "Dlk1 + Asb4",
                         "GLUT8" = "Bc1 + Pax6",
                         "GLUT9" = "Pmch + Cartpt",
                         "GLUT10" = "Synpr + Hpcal1",
                         "GLUT11" = "Foxb1 + Pitx2",
                         "GLUT12" = "Apoe + Atp1a2",
                         "GLUT13" = "Trh + Ghrh",
                         "GLUT14" = "Tcf4 + Calca",
                         "GLUT15" = "Nppc + Nrn1",
                         "GLUT16" = "Tcf4 + Grp",
                         "GLUT17" = "Avp + Oxt",
                         "GLUT18" = "Ucn3 + Trh"
                         )

# Save RDS GLUT.sub
saveRDS(GLUT.sub, "./RDS_obj/GLUT_sub_named.RDS")

# Read GLUT.sub RDS
GLUT.sub <- readRDS("./RDS_obj/GLUT_sub_named.RDS")

# DimPlot of GLUT subtypes embedded in UMAP
DimPlot(GLUT.sub, reduction = 'umap', label = T) + labs(title = "UMAP of GLUT Subtypes") + NoLegend()
DimPlot(GLUT.sub, reduction = 'umap', group.by = "subtypes" , label = T) + labs(title = "UMAP of GLUT Subtype Numbers") + NoLegend()

# DimPlot of GLUT subtypes embedded in UMAP group by source
DimPlot(GLUT.sub, reduction = 'umap', group.by = "Source") + labs(title = "UMAP of GLUT Subtypes")

DimPlot(GLUT.sub, reduction = 'umap', group.by = "Ephys_Cathegory") + labs(title = "UMAP of Electro Physiological Property of GLUT Subtypes")
dittoBarPlot(object = GLUT.sub, var = "Ephys_Cathegory", scale = "count", group.by = "subtypes")

p1 <- DimPlot(GLUT.sub, cells.highlight = rownames(GLUT.sub@meta.data)[GLUT.sub$Ephys_Cathegory == 1], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Patch-seq with 1 EC") + NoLegend()
p2 <- DimPlot(GLUT.sub, cells.highlight = rownames(GLUT.sub@meta.data)[GLUT.sub$Ephys_Cathegory == 2], cols.highlight = "blue", cols = "gray") + labs(title = "Cells from Patch-seq with 2 EC") + NoLegend()
p3 <- DimPlot(GLUT.sub, cells.highlight = rownames(GLUT.sub@meta.data)[GLUT.sub$Ephys_Cathegory == 3], cols.highlight = "dark green", cols = "gray") + labs(title = "Cells from Patch-seq with 3 EC") + NoLegend()
p4 <- DimPlot(GLUT.sub, cells.highlight = rownames(GLUT.sub@meta.data)[GLUT.sub$Ephys_Cathegory == 4], cols.highlight = "yellow", cols = "gray") + labs(title = "Cells from Patch-seq with 4 EC") + NoLegend()
p5 <- DimPlot(GLUT.sub, cells.highlight = rownames(GLUT.sub@meta.data)[GLUT.sub$Ephys_Cathegory == 5], cols.highlight = "purple", cols = "gray") + labs(title = "Cells from Patch-seq with 5 EC") + NoLegend()
p6 <- DimPlot(GLUT.sub, cells.highlight = rownames(GLUT.sub@meta.data)[GLUT.sub$Ephys_Cathegory == 6], cols.highlight = "orange", cols = "gray") + labs(title = "Cells from Patch-seq with 6 EC") + NoLegend()
p7 <- DimPlot(GLUT.sub, cells.highlight = rownames(GLUT.sub@meta.data)[GLUT.sub$Ephys_Cathegory == 7], cols.highlight = "black", cols = "gray") + labs(title = "Cells from Patch-seq with 7 EC") + NoLegend()
p1 + p2 + p3 + p4 + p5 + p6 + p7 + (DimPlot(GLUT.sub, reduction = "umap", group.by = "subtypes", label = T) + NoLegend())

# # Visualize markers of GLUT clusters with stack VlnPlot
VlnPlot(GLUT.sub , features = c(unique(GLUT.sub.top2.markers$gene)), 
        stack = T, group.by = "old.ident", split.by = "subtypes", flip = T) + 
  NoLegend() + labs(title = "Markers of GLUT clusters")

# UMAP plots of cell from different sources in GLUT.se
p1 <- DimPlot(GLUT.sub, cells.highlight = rownames(GLUT.sub@meta.data)[GLUT.sub$Source == "Unpub_LHA"], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Unpub_LHA in GLUT") + NoLegend()
p2 <- DimPlot(GLUT.sub, cells.highlight = rownames(GLUT.sub@meta.data)[GLUT.sub$Source == "Unpub_LHb"], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Unpub_LHb in GLUT") + NoLegend()
p3 <- DimPlot(GLUT.sub, cells.highlight = rownames(GLUT.sub@meta.data)[GLUT.sub$Source == "Unpub_HSV"], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Unpub_HSV in GLUT") + NoLegend()
p4 <- DimPlot(GLUT.sub, cells.highlight = rownames(GLUT.sub@meta.data)[GLUT.sub$Source == "Unpub_LHA_LHb_patch"], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Unpub_LHA_LHb_patch in GLUT") + NoLegend()
p5 <- DimPlot(GLUT.sub, cells.highlight = rownames(GLUT.sub@meta.data)[GLUT.sub$Source == "Chen"], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Chen in GLUT") + NoLegend()
p6 <- DimPlot(GLUT.sub, cells.highlight = rownames(GLUT.sub@meta.data)[GLUT.sub$Source == "Mickelsen"], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Mickelsen in GLUT") + NoLegend()
p7 <- DimPlot(GLUT.sub, cells.highlight = rownames(GLUT.sub@meta.data)[GLUT.sub$Source == "Rossi_2019"], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Rossi_2019 in GLUT") + NoLegend()
p8 <- DimPlot(GLUT.sub, cells.highlight = rownames(GLUT.sub@meta.data)[GLUT.sub$Source == "Rossi_2021"], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Rossi_2021 in GLUT") + NoLegend()
p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + (DimPlot(GLUT.sub, reduction = "umap", group.by = "subtypes", label = T) + NoLegend())

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
                         "0" = "GABA1",
                         "1" = "GABA2",
                         "2" = "GABA3",
                         "3" = "GABA4",
                         "4" = "GABA5",
                         "5" = "GABA6",
                         "6" = "GABA7",
                         "7" = "GABA8",
                         "8" = "GABA9",
                         "9" = "GABA10",
                         "10" = "GABA11",
                         "11" = "GABA12",
                         "12" = "GABA13",
                         "13" = "GABA14",
                         "14" = "GABA15",
                         "15" = "GABA16",
                         "16" = "GABA17",
                         "17" = "GABA18",
                         "18" = "GABA19",
                         "19" = "GABA20"
                         )

# Save current active ident as subtypes
GABA.sub[["subtypes"]] <- Idents(GABA.sub)

GABA.sub <- RenameIdents(GABA.sub,
                         "GABA1" = "Tac2 + A030009H04Rik",
                         "GABA2" = "Nts + Cartpt",
                         "GABA3" = "Tac1 + Meis2",
                         "GABA4" = "Ddc + Sncg",
                         "GABA5" = "Gal + Lmo3",
                         "GABA6" = "Sst + Meis2",
                         "GABA7" = "Tcf4 + Synpr",
                         "GABA8" = "Gm42418 + 2900097C17Rik",
                         "GABA9" = "Nnat + C1ql2",
                         "GABA10" = "Tcf7l2 + Calb1",
                         "GABA11" = "Ebf1 + Nkx2-4",
                         "GABA12" = "Plp1 + Ptgds",
                         "GABA13" = "Sst + Arpp21",
                         "GABA14" = "Fam81a + Rprm",
                         "GABA15" = "Pmch + Olig1",
                         "GABA16" = "Cck + Vtn",
                         "GABA17" = "Cst3 + C1qa",
                         "GABA18" = "Slc5a7 + Acly",
                         "GABA19" = "Calb2 + Synpr",
                         "GABA20" = "Hdc + Slc18a2"
                         )

# Save RDS GABA.sub
saveRDS(GABA.sub, "./RDS_obj/GABA_sub_named.RDS")

# Read GABA.sub RDS
GABA.sub <- readRDS("./RDS_obj/GABA_sub_named.RDS")

# DimPlot of GABA subtypes embedded in UMAP
DimPlot(GABA.sub, reduction = 'umap', label = T) + labs(title = "UMAP of GABA Subtypes") + NoLegend()
DimPlot(GABA.sub, reduction = 'umap', group.by = "subtypes" , label = T) + labs(title = "UMAP of GABA Subtype Numbers") + NoLegend()

# DimPlot of GABA subtypes embedded in UMAP group by source
DimPlot(GABA.sub, reduction = 'umap', group.by = "Source") + labs(title = "UMAP of GABA Subtypes")

FeaturePlot(GABA.sub, features = c("Slc32a1", "Slc17a6"))
FeaturePlot(GLUT.sub, features = c("Slc32a1", "Slc17a6"))

# Visualize markers of GABA clusters with stack VlnPlot
VlnPlot(GABA.sub , features = c(unique(GABA.sub.top2.markers$gene)),
        stack = T, group.by = "old.ident", split.by = "subtypes", flip = T) + 
  NoLegend() + labs(title = "Markers of GABA clusters")

# # UMAP plots of cell from different sources in GABA.se
p1 <- DimPlot(GABA.sub, cells.highlight = rownames(GABA.sub@meta.data)[GABA.sub$Source == "Unpub_LHA"], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Unpub_LHA in GABA") + NoLegend()
p2 <- DimPlot(GABA.sub, cells.highlight = rownames(GABA.sub@meta.data)[GABA.sub$Source == "Unpub_LHb"], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Unpub_LHb in GABA") + NoLegend()
p3 <- DimPlot(GABA.sub, cells.highlight = rownames(GABA.sub@meta.data)[GABA.sub$Source == "Unpub_HSV"], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Unpub_HSV in GABA") + NoLegend()
p4 <- DimPlot(GABA.sub, cells.highlight = rownames(GABA.sub@meta.data)[GABA.sub$Source == "Unpub_LHA_LHb_patch"], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Unpub_LHA_LHb_patch in GABA") + NoLegend()
p5 <- DimPlot(GABA.sub, cells.highlight = rownames(GABA.sub@meta.data)[GABA.sub$Source == "Chen"], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Chen in GABA") + NoLegend()
p6 <- DimPlot(GABA.sub, cells.highlight = rownames(GABA.sub@meta.data)[GABA.sub$Source == "Mickelsen"], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Mickelsen in GABA") + NoLegend()
p7 <- DimPlot(GABA.sub, cells.highlight = rownames(GABA.sub@meta.data)[GABA.sub$Source == "Rossi_2019"], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Rossi_2019 in GABA") + NoLegend()
p8 <- DimPlot(GABA.sub, cells.highlight = rownames(GABA.sub@meta.data)[GABA.sub$Source == "Rossi_2021"], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Rossi_2021 in GABA") + NoLegend()
p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + (DimPlot(GABA.sub, reduction = "umap", group.by = "subtypes", label = T) + NoLegend())

# Plots of GLUT fave markers
FeaturePlot(neuron.subset, features = c('Trp73', 'Glp1r', 'Gpr101', 'Samd3', 'Nov', 'Ndnf', 'Chat', 'Slc5a7'))
FeaturePlot(GLUT.sub, features = c('Trp73', 'Glp1r', 'Gpr101', 'Samd3', 'Nov', 'Ndnf', 'Chat', 'Slc5a7'))
FeaturePlot(GABA.sub, features = c('Trp73', 'Glp1r', 'Gpr101', 'Samd3', 'Nov', 'Ndnf', 'Chat', 'Slc5a7'))
VlnPlot(neuron.subset, features = c('Trp73', 'Glp1r', 'Gpr101', 'Samd3', 'Nov', 'Ndnf', 'Chat', 'Slc5a7'), stack = T, flip = T) + NoLegend()
VlnPlot(GLUT.sub, features = c('Trp73', 'Glp1r', 'Gpr101', 'Samd3', 'Nov', 'Ndnf', 'Chat', 'Slc5a7'), split.by = "subtypes", stack = T, flip = T) + NoLegend()
VlnPlot(GABA.sub, features = c('Trp73', 'Glp1r', 'Gpr101', 'Samd3', 'Nov', 'Ndnf', 'Chat', 'Slc5a7'), split.by = "subtypes", stack = T, flip = T) + NoLegend()

# Plots of GABA fave markers
FeaturePlot(neuron.subset, features = c('Prlr', 'Esr1', 'Pex5l', 'Pvalb', 'Gpr149', 'Calcr', 'Npy', 'Kcnab1'))
FeaturePlot(GLUT.sub, features = c('Prlr', 'Esr1', 'Pex5l', 'Pvalb', 'Gpr149', 'Calcr', 'Npy', 'Kcnab1'))
FeaturePlot(GABA.sub, features = c('Prlr', 'Esr1', 'Pex5l', 'Pvalb', 'Gpr149', 'Calcr', 'Npy', 'Kcnab1'))
VlnPlot(neuron.subset, features = c('Prlr', 'Esr1', 'Pex5l', 'Pvalb', 'Gpr149', 'Calcr', 'Npy', 'Kcnab1'), stack = T, flip = T) + NoLegend()
VlnPlot(GLUT.sub, features = c('Prlr', 'Esr1', 'Pex5l', 'Pvalb', 'Gpr149', 'Calcr', 'Npy', 'Kcnab1'), split.by = "subtypes", stack = T, flip = T) + NoLegend()
VlnPlot(GABA.sub, features = c('Prlr', 'Esr1', 'Pex5l', 'Pvalb', 'Gpr149', 'Calcr', 'Npy', 'Kcnab1'), split.by = "subtypes", stack = T, flip = T) + NoLegend()

# Plots of mk from Cell paper
VlnPlot(neuron.subset, features = c('Map1b', 'Slc17a6', 'Slc32a1', 'Gad1', 'Pmch', 
                               'Cartpt', 'Hcrt', 'Trh', 'Sst', 'Gpr83', 'Gpr101',
                               'Meis2', 'Otp', 'Calb1', 'Calb2', 'Nrgn', 'Tac1',
                               'Tac2', 'Bdnf', 'Synpr', 'Nts', 'Col25a1', 'Gal', 'Th'),
        stack = T, fill.by = 'ident', flip = T) + labs(title = "Expression of Markers from Wang et al., in Neuron Subclass") + NoLegend()

FeaturePlot(neuron.subset, features = c('Map1b', 'Slc17a6', 'Slc32a1', 'Gad1', 'Pmch', 
                                   'Cartpt', 'Hcrt', 'Trh', 'Sst', 'Gpr83', 'Gpr101',
                                   'Meis2', 'Otp', 'Calb1', 'Calb2', 'Nrgn', 'Tac1',
                                   'Tac2', 'Bdnf', 'Synpr', 'Nts', 'Col25a1', 'Gal', 'Th'))

DotPlot(neuron.subset, features = c('Map1b', 'Slc17a6', 'Slc32a1', 'Gad1', 'Pmch', 
                                    'Cartpt', 'Hcrt', 'Trh', 'Sst', 'Gpr83', 'Gpr101',
                                    'Meis2', 'Otp', 'Calb1', 'Calb2', 'Nrgn', 'Tac1',
                                    'Tac2', 'Bdnf', 'Synpr', 'Nts', 'Col25a1', 'Gal', 'Th'),
        split.by = "Source", cols = "RdBu") + RotatedAxis() + labs(title = "DotPlot of Markers from Wang et al., in Neuron Subclass by Source")

VlnPlot(GLUT.sub, features = c('Map1b', 'Slc17a6', 'Slc32a1', 'Gad1', 'Pmch', 
                               'Cartpt', 'Hcrt', 'Trh', 'Sst', 'Gpr83', 'Gpr101',
                               'Meis2', 'Otp', 'Calb1', 'Calb2', 'Nrgn', 'Tac1',
                               'Tac2', 'Bdnf', 'Synpr', 'Nts', 'Col25a1', 'Gal', 'Th'),
        stack = T, fill.by = 'ident', flip = T) + labs(title = "Expression of Markers from Wang et al., in GLUT Subtypes") + NoLegend()

FeaturePlot(GLUT.sub, features = c('Map1b', 'Slc17a6', 'Slc32a1', 'Gad1', 'Pmch', 
                               'Cartpt', 'Hcrt', 'Trh', 'Sst', 'Gpr83', 'Gpr101',
                               'Meis2', 'Otp', 'Calb1', 'Calb2', 'Nrgn', 'Tac1',
                               'Tac2', 'Bdnf', 'Synpr', 'Nts', 'Col25a1', 'Gal', 'Th'))


VlnPlot(GABA.sub, features = c('Map1b', 'Slc17a6', 'Slc32a1', 'Gad1', 'Pmch', 
                               'Cartpt', 'Hcrt', 'Trh', 'Sst', 'Gpr83', 'Gpr101',
                               'Meis2', 'Otp', 'Calb1', 'Calb2', 'Nrgn', 'Tac1',
                               'Tac2', 'Bdnf', 'Synpr', 'Nts', 'Col25a1', 'Gal', 'Th'),
        stack = T, fill.by = 'ident', flip = T) + labs(title = "Expression of Markers from Wang et al., in GABA Subtypes") + NoLegend()

FeaturePlot(GABA.sub, features = c('Map1b', 'Slc17a6', 'Slc32a1', 'Gad1', 'Pmch', 
                                   'Cartpt', 'Hcrt', 'Trh', 'Sst', 'Gpr83', 'Gpr101',
                                   'Meis2', 'Otp', 'Calb1', 'Calb2', 'Nrgn', 'Tac1',
                                   'Tac2', 'Bdnf', 'Synpr', 'Nts', 'Col25a1', 'Gal', 'Th'))


# Bar plot of cell count from source in each cluster
dittoBarPlot(object = neuron.subset, var = "Source", scale = "count", group.by = "subtypes") + labs(title = "Bar Plot of Absolute Cell Count by Source in All Neuronal Cells")
dittoBarPlot(object = neuron.subset, var = "Source", scale = "percent", group.by = "subtypes") + labs(title = "Bar Plot of Percent Count by Source in All Neuronal Cells")
dittoBarPlot(object = GLUT.sub, var = "Source", scale = "count", group.by = "subtypes") + labs(title = "Bar Plot of Absolute Cell Count by Source in GLUT")
dittoBarPlot(object = GLUT.sub, var = "Source", scale = "percent", group.by = "subtypes") + labs(title = "Bar Plot of Percent Count by Source in GLUT")
dittoBarPlot(object = GABA.sub, var = "Source", scale = "count", group.by = "subtypes") + labs(title = "Bar Plot of Absolute Cell Count by Source in GABA")
dittoBarPlot(object = GABA.sub, var = "Source", scale = "percent", group.by = "subtypes") + labs(title = "Bar Plot of Percent Count by Source in GABA")

# DotPlot of neuron.subset with key markers and split by Source
DotPlot(neuron.subset, features = c(unique(neuron.top2.markers$gene)), split.by = "Source", cols = "RdBu") + RotatedAxis()
DotPlot(neuron.subset, features = c('Trp73', 'Glp1r', 'Gpr101', 'Samd3', 'Nov', 'Ndnf', 'Chat', 'Slc5a7', 'Prlr', 'Esr1', 'Pex5l', 'Pvalb', 'Gpr149', 'Calcr', 'Npy', 'Kcnab1')) + RotatedAxis()
DotPlot(neuron.subset, features = c('Trp73', 'Glp1r', 'Gpr101', 'Samd3', 'Nov', 'Ndnf', 'Chat', 'Slc5a7', 'Prlr', 'Esr1', 'Pex5l', 'Pvalb', 'Gpr149', 'Calcr', 'Npy', 'Kcnab1'), split.by = "Source", cols = "RdBu") + RotatedAxis()


# Co-expression in GLUT
plot_density(GLUT.sub, c('Slc5a7', 'Ndnf'), joint = T, combine = F)
plot_density(GLUT.sub, c('Samd3', 'Nov'), joint = T, combine = F)
plot_density(GLUT.sub, c('Gpr149', 'Trp73'), joint = T, combine = F)
plot_density(GLUT.sub, c('Gpr149', 'Samd3'), joint = T, combine = F)
plot_density(GLUT.sub, c('Gpr101', 'Pvalb'), joint = T, combine = F)
FeaturePlot(GLUT.sub, c("Chrna1", "Chrnb2"))

# Co-expression in GABA
plot_density(GABA.sub, c('Slc5a7', 'Chat'), joint = T, combine = F)
plot_density(GABA.sub, c('Gpr101', 'Calcr'), joint = T, combine = F)
FeaturePlot(GABA.sub, c("Chrna1", "Chrnb2"))

# FeaturePlot of the top 2 markers for each cluster
FeaturePlot(GLUT.sub, features = c(unique(GLUT.sub.top2.markers$gene)))
FeaturePlot(GABA.sub, features = c(unique(GABA.sub.top2.markers$gene)))







