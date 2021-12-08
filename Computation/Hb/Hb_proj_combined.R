library(dplyr)
library(patchwork)
library(ggplot2)
library(scales)
library(Seurat)
library(SeuratObject)

# Read RDS of pub and unpub data
unpub <- readRDS('/Users/pax/Google\ Drive/My\ Drive/Meletis/Computation/Unpub/Unpub_integrate/RDS_obj/Unpub_integration.RDS')
hashikawa <- readRDS('/Users/pax/Google\ Drive/My\ Drive/Meletis/Computation/Hb/RDS_obj/Hashikawa_combined.RDS')
wallace <- readRDS('/Users/pax/Google\ Drive/My\ Drive/Meletis/Computation/Hb/RDS_obj/hab_batch1.rds')

# Update Seurat obj of wallace from v2 to v3 and add source metadata
wallace <- UpdateSeuratObject(object = wallace)
wallace[["source"]] <- "Wallace"
wallace[["percent.mt"]] <- wallace$percent.mito
wallace[["percent.mito"]] <- NULL

# Create list of multiple Seurat object
Hb.list <- list(unpub, hashikawa, wallace)

# Select RNA as default assay, normalize, and identify variable features for each dataset independently
Hb.list <- lapply(X = Hb.list, FUN = function(x) {
  DefaultAssay(x) <- "RNA"
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = Hb.list, nfeatures = 2000)

# Find integration anchors
Hb.anchors <- FindIntegrationAnchors(object.list = Hb.list, anchor.features = features)

# Creates an 'integrated' data assay
Hb.combined <- IntegrateData(anchorset = Hb.anchors)

# Assign new Sources
Hb.combined[["Source"]] <- Hb.combined$source
unpub.LHA <- rownames(Hb.combined@meta.data)[Hb.combined$orig.ident == "LHA" & Hb.combined$source == "Unpub"]
unpub.patch <- rownames(Hb.combined@meta.data)[Hb.combined$orig.ident == "LHA_LHb_patch" & Hb.combined$source == "Unpub"]
unpub.LHb <- rownames(Hb.combined@meta.data)[Hb.combined$orig.ident == "LHb_lskra" & Hb.combined$source == "Unpub"]
unpub.HSV <- rownames(Hb.combined@meta.data)[Hb.combined$orig.ident == "HSV" & Hb.combined$source == "Unpub"]

# Classify cell according to the unpub areas
Hb.combined$Source[which(names(Hb.combined$nCount_RNA) %in% unpub.LHA)] <- "Unpub_LHA"
Hb.combined$Source[which(names(Hb.combined$nCount_RNA) %in% unpub.patch)] <- "Unpub_LHA_LHb_patch"
Hb.combined$Source[which(names(Hb.combined$nCount_RNA) %in% unpub.LHb)] <- "Unpub_LHb"
Hb.combined$Source[which(names(Hb.combined$nCount_RNA) %in% unpub.HSV)] <- "Unpub_HSV"

# Save RDS of Hb.combined
saveRDS(Hb.combined, "./RDS_obj/Hb_proj.RDS")

# Read RDS of Hb.combined
Hb.combined <- readRDS("./RDS_obj/Hb_proj.RDS")

# Select integrated as default assay
DefaultAssay(Hb.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
Hb.combined <- ScaleData(Hb.combined)
Hb.combined <- RunPCA(Hb.combined, npcs = 100)
Hb.combined <- RunUMAP(Hb.combined, reduction = "pca", dims = 1:70)
Hb.combined <- FindNeighbors(Hb.combined, reduction = "pca", dims = 1:70)
Hb.combined <- FindClusters(Hb.combined, resolution = 0.5)

DimPlot(Hb.combined, reduction = 'umap') + DimPlot(Hb.combined, reduction = 'umap', group.by = 'Source')

# Set Default assay to RNA
DefaultAssay(Hb.combined) <- "RNA"

# Find cluster markers
Hb.combined.markers <- FindAllMarkers(Hb.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Save RDS of Hb.combined and Hb.combined.markers
saveRDS(Hb.combined, "./RDS_obj/Hb_proj.RDS")
saveRDS(Hb.combined.markers, "./RDS_obj/Hb_proj_mk.RDS")

# Read RDS of Hb.combined and Hb.combined.markers
Hb.combined <- readRDS("./RDS_obj/Hb_proj.RDS")
Hb.combined.markers <- readRDS("./RDS_obj/Hb_proj_mk.RDS")

# Get top 2 makers of each cluster
top2.markers <- Hb.combined.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# Identify neurons, non-neurons, and none
# Find the median of neuron markers
neuron.markers <- as.data.frame(AverageExpression(Hb.combined, features = c("Tac2", "Slc17a7", "Slc17a6", "Snap25", "Gap43"))$RNA)
neuron.idents <- as.integer(colnames(neuron.markers)[colSums(neuron.markers)>15])
non_neuron.idents <- as.integer(colnames(neuron.markers)[colSums(neuron.markers)<5.5])

# Create lists of cell names that are neurons and non-neurons
neuron.list <- WhichCells(Hb.combined, idents = c(neuron.idents))
non_neuron.list <- WhichCells(Hb.combined, idents = c(non_neuron.idents))

# Add cell_type to metadata. Every cell is classified as non-neurons at the beginning
Hb.combined[["cell_type"]] <- "None"

# Classify cell with names in neuron.list to be neurons
Hb.combined$cell_type[which(names(Hb.combined$nCount_RNA) %in% neuron.list)] <- "Neurons"

# Classify cell with names in mom_neuron.list to be non-neurons
Hb.combined$cell_type[which(names(Hb.combined$nCount_RNA) %in% non_neuron.list)] <- "Non-neurons"

DimPlot(Hb.combined, reduction = 'umap', label = T) + DimPlot(Hb.combined, reduction = 'umap', group.by = 'cell_type')

# Remove cell_type of None
Hb.combined <- subset(Hb.combined, subset = cell_type == "None", invert = T)

# Save RDS of Hb.combined
saveRDS(Hb.combined, "./RDS_obj/Hb_proj_cell_type.RDS")

# Read RDS of Hb.combined
Hb.combined <- readRDS("./RDS_obj/Hb_proj_cell_type.RDS")

# FeaturePlot of cluster 6 and 16
FeaturePlot(Hb.combined, features = c("Ccl4", "Ccl7", "Pdgfra", "Cspg5"))

# Select neuron cells
neuron.sub <- subset(Hb.combined, subset = cell_type == "Neurons")

# Seurat analysis pipeline on neuron.subset
neuron.sub <- FindVariableFeatures(neuron.sub)
DefaultAssay(neuron.sub) <- "integrated"
neuron.sub <- ScaleData(neuron.sub, vars.to.regress = "percent.mt")
neuron.sub <- RunPCA(neuron.sub, npcs = 100)
neuron.sub <- RunUMAP(neuron.sub, reduction = "pca", dims = 1:70)
neuron.sub <- FindNeighbors(neuron.sub, reduction = "pca", dims = 1:70)
neuron.sub <- FindClusters(neuron.sub, resolution = 0.5)

DimPlot(neuron.sub, reduction = 'umap')

# Find cluster markers
DefaultAssay(neuron.sub) <- "RNA"
neuron.sub.markers <- FindAllMarkers(neuron.sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Get top 2 makers of each cluster
neuron.sub.top2.markers <- neuron.sub.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# Find the median of Hb area (MHb or Lhb) markers
Hb.area.median <- as.data.frame(AverageExpression(neuron.sub, features = c("Tac2", "Pcdh10"))$RNA)
MHb.ident <- as.integer(colnames(Hb.area.median)[Hb.area.median["Tac2",] > Hb.area.median["Pcdh10",]])
LHb.ident <- as.integer(colnames(Hb.area.median)[Hb.area.median["Tac2",] < Hb.area.median["Pcdh10",]])

# Create lists of cell names that are MHb and LHb
MHb.list <- WhichCells(neuron.sub, idents = c(MHb.ident))
LHb.list <- WhichCells(neuron.sub, idents = c(LHb.ident))

# Add Hb_area to metadata
neuron.sub[["Hb_area"]] <- "None"

# Classify cell with names in GABA.list to be GABA
neuron.sub$Hb_area[which(names(neuron.sub$nCount_RNA) %in% MHb.list)] <- "MHb"

# Classify cell with names in GLUT.list to be GLUT
neuron.sub$Hb_area[which(names(neuron.sub$nCount_RNA) %in% LHb.list)] <- "LHb"

# Save RDS
saveRDS(neuron.sub, "./RDS_obj/neuron_Hb_area_proj.RDS")

# Read RDS
neuron.sub <- readRDS("./RDS_obj/neuron_Hb_area_proj.RDS")

DimPlot(neuron.sub, reduction = "umap", group.by = "Hb_area") + DimPlot(neuron.sub, reduction = "umap", label = T) + DimPlot(neuron.sub, reduction = "umap", group.by = "Source")

FeaturePlot(neuron.sub, features = c("Tac2", "Pcdh10"))

# Save Old identity
neuron.sub[["old.ident"]] <- Idents(neuron.sub)

# Rename MHb idents
neuron.sub <- RenameIdents(object = neuron.sub, 
                              "0" = "MHb1",
                              "1" = "MHb2",
                              "2" = "MHb3",
                              "3" = "MHb4",
                              "10" = "MHb5",
                              "13" = "MHb6")

# Rename LHb idents
neuron.sub <- RenameIdents(object = neuron.sub, 
                              "4" = "LHb1",
                              "5" = "LHb2",
                              "6" = "LHb3",
                              "7" = "LHb4",
                              "8" = "LHb5",
                              "9" = "LHb6",
                              "11" = "LHb7",
                              "12" = "LHb8",
                              "14" = "LHb9")

DimPlot(neuron.sub, reduction = "umap", label = T) + labs(title = "UMAP of Hb area") + NoLegend() 
DimPlot(neuron.sub, reduction = "umap", group.by = "Source")

# Find cluster markers
DefaultAssay(neuron.sub) <- "RNA"
neuron.sub.markers <- FindAllMarkers(neuron.sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Do heat map
neuron.sub.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DefaultAssay(neuron.sub) <- "integrated"
DoHeatmap(neuron.sub, features = top10$gene) + NoLegend()

# Save RDS
saveRDS(neuron.sub, "./RDS_obj/neuron_Hb_area_named.RDS")

# Read RDS
neuron.sub <- readRDS("./RDS_obj/neuron_Hb_area_named.RDS")
DotPlot(neuron.sub, features = c('Agt', 'Cadps2', 'Chat', 'Chrna3', 'Chrna4', 'Chrnb3', 'Cplx1', 
                                 'Gabbr1', 'Napb', 'Nppa', 'Slc17a7', 'Slc18a3',
                                 'Slc5a7', 'Snca', 'Sncg', 'Stxbp5l', 'Sv2b', 'Syt6',
                                 'Syt9', 'Ywhaz', 'Abat', 'Atp2a2', 'Atp2b2', 'Baiap3',
                                 'Cacna1b', 'Cadps', 'Camk2a', 'Cd47', 'Clu', 'Htr2c',
                                 'Kcnc4', 'Lin7a', 'Nat8l', 'Nlgn1', 'Nrxn3', 'Ntrk2',
                                 'Pak1', 'Pde1b', 'Rph3a', 'Sirpa', 'Slc17a6', 'Slc1a2',
                                 'Slc6a1', 'Slc6a7', 'Sv2a', 'Syn2', 'Syp', 'Syt11', 'Vamp1'),
        split.by = "Source", cols = "RdBu") + RotatedAxis() + labs(title = "Neurotransmitter Related Genes Split")

# DotPlot of key markers 
DotPlot(neuron.sub, features = c('Gabbr1', 'Gabbr2', 'Gria3', 'Grid2', 'Grin2b', 'Grm5', 'Gabra1', 'Gabrb1', 'Gabrg2', 'Gabrg3')) + RotatedAxis()
DotPlot(neuron.sub, features = c('Gabbr1', 'Gabbr2', 'Gria3', 'Grid2', 'Grin2b', 'Grm5', 'Gabra1', 'Gabrb1', 'Gabrg2', 'Gabrg3'), split.by = "Source", cols = "RdBu") + RotatedAxis()

# Subset to select cells from LHb
LHb.se <- subset(neuron.sub, subset = Hb_area == "LHb")

# Remove eGFP
counts <- GetAssayData(LHb.se, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c('eGFP'))),]
LHb.se <- subset(LHb.se, features = rownames(counts))

# Run Seurat analysis pipeline on LHb.se
LHb.se <- FindVariableFeatures(LHb.se)
DefaultAssay(LHb.se) <- "integrated"
LHb.se <- ScaleData(LHb.se, vars.to.regress = "percent.mt")
LHb.se <- RunPCA(LHb.se, npcs = 100)
LHb.se <- RunUMAP(LHb.se, reduction = "pca", dims = 1:70)
LHb.se <- FindNeighbors(LHb.se, reduction = "pca", dims = 1:70)
LHb.se <- FindClusters(LHb.se, resolution = 0.33)
DimPlot(LHb.se, reduction = "umap", label = T)
DimPlot(LHb.se, reduction = "umap", group.by = "Source")

# Find cluster markers
DefaultAssay(LHb.se) <- "RNA"
LHb.se.markers <- FindAllMarkers(LHb.se, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Get top 2 makers of each cluster
LHb.se.top2.markers <- LHb.se.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# Save Old identity
LHb.se[["old.ident"]] <- Idents(LHb.se)

# Rename LHb subclass idents
LHb.se <- RenameIdents(object = LHb.se, 
                           "0" = "Ndfip1+Avp",
                           "1" = "Pcp4+Nwd2",
                           "2" = "Cbln1+Rph3a",
                           "3" = "Vps8+Ahi1",
                           "4" = "Sst+Gpr151",
                           "5" = "Ret+Ppargc1a",
                           "6" = "Nptxr+Syn2",
                           "7" = "Ntng1+C1ql3",
                           "8" = "Ctxn3+Foxp1",
                           "9" = "Tac2+Scube1",
                           "10" = "Nrgn+Ppp3ca")

DimPlot(LHb.se, reduction = "umap", label = T) + labs(title = "UMAP of LHb Subclass") + NoLegend() 

# FeaturePlot of high Slc17a6 and low Slc32a1. Thus most are GLUT
FeaturePlot(LHb.se, features = c('Slc17a6', 'Slc32a1'))

# Subset to select cells from MHb
MHb.se <- subset(neuron.sub, subset = Hb_area == "MHb")

# Remove eGFP
counts <- GetAssayData(MHb.se, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c('eGFP'))),]
MHb.se <- subset(MHb.se, features = rownames(counts))

# Run Seurat analysis pipeline on MHb.se
MHb.se <- FindVariableFeatures(MHb.se)
DefaultAssay(MHb.se) <- "integrated"
MHb.se <- ScaleData(MHb.se, vars.to.regress = "percent.mt")
MHb.se <- RunPCA(MHb.se, npcs = 100)
MHb.se <- RunUMAP(MHb.se, reduction = "pca", dims = 1:70)
MHb.se <- FindNeighbors(MHb.se, reduction = "pca", dims = 1:70)
MHb.se <- FindClusters(MHb.se, resolution = 0.5)
DimPlot(MHb.se, reduction = "umap", label = T)
DimPlot(MHb.se, reduction = "umap", group.by = "Source")

# Find cluster markers
DefaultAssay(MHb.se) <- "RNA"
MHb.se.markers <- FindAllMarkers(MHb.se, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Get top 2 makers of each cluster
MHb.se.top2.markers <- MHb.se.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# Save Old identity
MHb.se[["old.ident"]] <- Idents(MHb.se)

# Rename MHb subclass idents
MHb.se <- RenameIdents(object = MHb.se, 
                            "0" = "Nrgn+Lrrc55",
                            "1" = "Kcnmb4+Syt15",
                            "2" = "Cck+Adcyap1",
                            "3" = "Synpr+Prkcb",
                            "4" = "Pam+Actb",
                            "5" = "Car2+Pcp4l1",
                            "6" = "Rn7s1+Gm37376",
                            "7" = "Ptgds+X6330403A02Rik",
                            "8" = "Olig1+Pdgfra")

DimPlot(MHb.se, reduction = "umap", label = T) + labs(title = "UMAP of MHb Subclass") + NoLegend()


FeaturePlot(MHb.se, features = "Meis2")
FeaturePlot(MHb.se, features = c('Slc32a1', 'Slc17a6', 'Slc17a7', 'Tac2'))
FeaturePlot(MHb.se, features = c('Slc18a3', 'Chat'))
DotPlot(MHb.se, features = c('Slc17a6', 'Slc17a7', 'Tac2'), split.by = "Source", cols = "RdBu") + RotatedAxis()


FeaturePlot(LHb.se, features = c('Slc32a1', 'Slc17a6'))
FeaturePlot(LHb.se, features = c('Chrm3', 'Vgf', 'Gpr151', 'Peg10', 'Cartpt', 'Sst')) + DimPlot(LHb.se, reduction = 'umap', label = T) + NoLegend()











