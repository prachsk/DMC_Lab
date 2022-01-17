library(dplyr)
library(patchwork)
library(ggplot2)
library(scales)
library(Seurat)
library(SeuratObject)
library(dittoSeq)
library(Nebulosa)

# Read RDS of pub and unpub data
LHb.rds <- readRDS('//Users/pax/Google\ Drive/My\ Drive/Meletis/Computation/Unpub/RDS_obj/LHb_lskra.RDS')
hashikawa <- readRDS('/Users/pax/Google\ Drive/My\ Drive/Meletis/Computation/Hb/RDS_obj/Hashikawa_combined.RDS')
wallace <- readRDS('/Users/pax/Google\ Drive/My\ Drive/Meletis/Computation/Hb/RDS_obj/hab_batch1.rds')

# Update Seurat obj of wallace from v2 to v3 and add source metadata
wallace <- UpdateSeuratObject(object = wallace)
wallace[["source"]] <- "Wallace"
wallace[["percent.mt"]] <- wallace$percent.mito
wallace[["percent.mito"]] <- NULL

# Create list of multiple Seurat object
Hb.list <- list(LHb.rds, hashikawa, wallace)

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
Sert.LHb <- rownames(Hb.combined@meta.data)[Hb.combined$area == "Sert" & Hb.combined$source == "Unpub"]
TH.LHb <- rownames(Hb.combined@meta.data)[Hb.combined$area == "TH" & Hb.combined$source == "Unpub"]

# Classify cell according to the unpub areas
Hb.combined$Source[which(names(Hb.combined$nCount_RNA) %in% Sert.LHb)] <- "Unpub_LHb_Sert"
Hb.combined$Source[which(names(Hb.combined$nCount_RNA) %in% TH.LHb)] <- "Unpub_LHb_TH"

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

# Seurat analysis pipeline on neuron.sub
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

# Remove cluster 7 as it is low quality
neuron.sub <- subset(neuron.sub, idents = 7, invert = T)

# Re-run analysis pipeline on neuron.sub
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
                              "6" = "MHb5",
                              "10" = "MHb6",
                              "11" = "MHb7",
                              "13" = "MHb8",
                              "14" = "MHb9")

# Rename LHb idents
neuron.sub <- RenameIdents(object = neuron.sub, 
                              "4" = "LHb1",
                              "5" = "LHb2",
                              "7" = "LHb3",
                              "8" = "LHb4",
                              "9" = "LHb5",
                              "12" = "LHb6")

DimPlot(neuron.sub, reduction = "umap", label = T) + labs(title = "UMAP of Hb area") + NoLegend() 
DimPlot(neuron.sub, reduction = "umap", group.by = "Source")


# UMAP plots of cell from different sources in neuron.sub
p1 <- DimPlot(neuron.sub, cells.highlight = rownames(neuron.sub@meta.data)[neuron.sub$Source == "Unpub_LHb_Sert"], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Unpub_LHb_Sert") + NoLegend()
p2 <- DimPlot(neuron.sub, cells.highlight = rownames(neuron.sub@meta.data)[neuron.sub$Source == "Unpub_LHb_TH"], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Unpub_LHb_TH") + NoLegend()
p3 <- DimPlot(neuron.sub, cells.highlight = rownames(neuron.sub@meta.data)[neuron.sub$Source == "Hashikawa"], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Hashikawa") + NoLegend()
p4 <- DimPlot(neuron.sub, cells.highlight = rownames(neuron.sub@meta.data)[neuron.sub$Source == "Wallace"], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Wallace") + NoLegend()

p1 + p2 + p3 + p4 + (DimPlot(neuron.sub, reduction = "umap", label = T) + NoLegend())

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

# Rename LHb subtypes idents
LHb.se <- RenameIdents(object = LHb.se, 
                       "0" = "LHb1",
                       "1" = "LHb2",
                       "2" = "LHb3",
                       "3" = "LHb4",
                       "4" = "LHb5",
                       "5" = "LHb6",
                       "6" = "LHb7",
                       "7" = "LHb8")

LHb.se[["subtypes"]] <- Idents(LHb.se)

# Rename LHb subclass idents
LHb.se <- RenameIdents(object = LHb.se, 
                       "LHb1" = "Cbln1 + Arpp21",
                       "LHb2" = "Pcp4 + Nwd2",
                       "LHb3" = "Nrgn + Ahi1",
                       "LHb4" = "Sst + Fstl1",
                       "LHb5" = "Nptxr + Necab1",
                       "LHb6" = "Cplx1 + Nr2f2",
                       "LHb7" = "Ntng1 + C1ql3",
                       "LHb8" = "Ctxn3 + Cbln4")

# Save RDS
saveRDS(LHb.se, "./RDS_obj/LHb_sub.RDS")

# Read RDS
LHb.se <- readRDS("./RDS_obj/LHb_sub.RDS")

DimPlot(LHb.se, reduction = "umap", label = T) + labs(title = "UMAP of LHb Subclass") + NoLegend() 

# UMAP plots of cell from different sources in LHb.se
p1 <- DimPlot(LHb.se, cells.highlight = rownames(LHb.se@meta.data)[LHb.se$Source == "Unpub_LHb_Sert"], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Unpub_LHb_Sert") + NoLegend()
p2 <- DimPlot(LHb.se, cells.highlight = rownames(LHb.se@meta.data)[LHb.se$Source == "Unpub_LHb_TH"], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Unpub_LHb_TH") + NoLegend()
p3 <- DimPlot(LHb.se, cells.highlight = rownames(LHb.se@meta.data)[LHb.se$Source == "Hashikawa"], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Hashikawa") + NoLegend()
p4 <- DimPlot(LHb.se, cells.highlight = rownames(LHb.se@meta.data)[LHb.se$Source == "Wallace"], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Wallace") + NoLegend()

p1 + p2 + p3 + p4 + (DimPlot(LHb.se, reduction = "umap", group.by = "subtypes", label = T) + NoLegend())

FeaturePlot(LHb.se, features = c(unique(LHb.se.top2.markers$gene)))

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

# Remove cluster 4 as they are low quality cells
MHb.se <- subset(MHb.se, idents = 4, invert = TRUE)

# Re-runun Seurat analysis pipeline on MHb.se
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

# Rename MHb subtypes idents
MHb.se <- RenameIdents(object = MHb.se, 
                            "0" = "MHb1",
                            "1" = "MHb2",
                            "2" = "MHb3",
                            "3" = "MHb4",
                            "4" = "MHb5",
                            "5" = "MHb6",
                            "6" = "MHb7",
                            "7" = "MHb8",
                            "8" = "MHb9",
                            "9" = "MHb10",
                            "10" = "MHb11")

# Save Old identity
MHb.se[["subtypes"]] <- Idents(MHb.se)

# Rename MHb subtypes idents
MHb.se <- RenameIdents(object = MHb.se, 
                       "MHb1" = "Nrgn + Ctnnbip1",
                       "MHb2" = "Synpr + Rapgef4",
                       "MHb3" = "Kcnmb4 + Nnat",
                       "MHb4" = "Cck + Adcyap1",
                       "MHb5" = "Car2 + Pcp4l1",
                       "MHb6" = "Fosb + Fos",
                       "MHb7" = "Plp1 + Atp1a2",
                       "MHb8" = "X6330403A02Rik + Ptgds",
                       "MHb9" = "Apoe + Dusp1",
                       "MHb10" = "ERCC-00130 + ERCC-00074",
                       "MHb11" = "Olig1 + Gpr17")

# Save RDS
saveRDS(MHb.se, "./RDS_obj/MHb_sub.RDS")

# Read RDS
MHb.se <- readRDS("./RDS_obj/MHb_sub.RDS")

DimPlot(MHb.se, reduction = "umap", group.by = "subtypes", label = T) + labs(title = "UMAP of MHb Subclass") + NoLegend()

# UMAP plots of cell from different sources in MHb.se
p1 <- DimPlot(MHb.se, cells.highlight = rownames(MHb.se@meta.data)[MHb.se$Source == "Unpub_LHb_Sert"], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Unpub_LHb_Sert") + NoLegend()
p2 <- DimPlot(MHb.se, cells.highlight = rownames(MHb.se@meta.data)[MHb.se$Source == "Unpub_LHb_TH"], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Unpub_LHb_TH") + NoLegend()
p3 <- DimPlot(MHb.se, cells.highlight = rownames(MHb.se@meta.data)[MHb.se$Source == "Hashikawa"], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Hashikawa") + NoLegend()
p4 <- DimPlot(MHb.se, cells.highlight = rownames(MHb.se@meta.data)[MHb.se$Source == "Wallace"], cols.highlight = "red", cols = "gray") + labs(title = "Cells from Wallace") + NoLegend()

p1 + p2 + p3 + p4 + (DimPlot(MHb.se, reduction = "umap", group.by = "subtypes", label = T) + NoLegend())

FeaturePlot(MHb.se, features = c(unique(MHb.se.top2.markers$gene)))



FeaturePlot(MHb.se, features = "Meis2")
FeaturePlot(LHb.se, features = c('Slc32a1', 'Slc17a6', 'Slc17a7', 'Tac2'))
FeaturePlot(MHb.se, features = c('Slc18a3', 'Chat'))
DotPlot(MHb.se, features = c('Slc17a6', 'Slc17a7', 'Tac2'), split.by = "Source", cols = "RdBu") + RotatedAxis()


FeaturePlot(LHb.se, features = c('Slc32a1', 'Slc17a6'))
FeaturePlot(LHb.se, features = c('Chrm3', 'Vgf', 'Gpr151', 'Peg10', 'Cartpt', 'Sst')) + DimPlot(LHb.se, reduction = 'umap', label = T) + NoLegend()

# Bar plot of cell count from source in each cluster
dittoBarPlot(object = LHb.se, var = "Source", scale = "count", group.by = "subtypes") + labs(title = "Bar Plot of Absolute Cell Count by Source in LHb Areas")
dittoBarPlot(object = LHb.se, var = "Source", scale = "percent", group.by = "subtypes") + labs(title = "Bar Plot of Percent Count by Source in LHb Areas")
dittoBarPlot(object = MHb.se, var = "Source", scale = "count", group.by = "subtypes") + labs(title = "Bar Plot of Absolute Cell Count by Source in MHb Areas")
dittoBarPlot(object = MHb.se, var = "Source", scale = "percent", group.by = "subtypes") + labs(title = "Bar Plot of Percent Count by Source in MHb Areas")







