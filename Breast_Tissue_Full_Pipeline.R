# Breast_Tissue_Full_Pipeline.R
# Created by Andrew Sanford
# sanfordac2@vcu.edu
# 
# Full pipeline: SpotSweeper QC → Downstream Analysis
# This script contains save points -- uncomment load points for crashes and resume
#
# Required packages:
# BiocManager::install(c("SpotSweeper", "SpatialExperiment", "SingleCellExperiment",
#                        "ggspavis", "scuttle", "scater", "scran", "SingleR", "celldex"))
# install.packages(c("Seurat", "ggplot2", "patchwork", "dplyr", "jsonlite",
#                    "SeuratObject", "presto", "pheatmap", "igraph",
#                    "networkD3", "openxlsx", "RColorBrewer", "plyr"))


###///////////////////###
### LOAD ALL PACKAGES ###
###///////////////////###

library(Seurat)
library(SeuratObject)
library(SpotSweeper)
library(SpatialExperiment)
library(SingleCellExperiment)
library(ggplot2)
library(patchwork)
library(dplyr)
library(jsonlite)
library(ggspavis)
library(scuttle)
library(scater)
library(scran)
library(igraph)
library(presto)
library(pheatmap)
library(SingleR)
library(celldex)
library(RColorBrewer)
library(networkD3)
library(openxlsx)


###////////////////////////////////////////////////////###
### SECTION 1: LOADING DATA WITH SEURAT                ###
###////////////////////////////////////////////////////###

# Set working directory to where your 8 µm (or other bins) data is located
setwd("C:/Users/Andrew/Desktop/binned_outputs/square_008um")

# Load Seurat object from the spatial dataset
seurat_obj <- Load10X_Spatial(
  data.dir = "C:\\Users\\Andrew\\Desktop\\binned_outputs\\square_008um",
  filename = "filtered_feature_bc_matrix.h5",
  image = "tissue_lowres_image.png",
  assay = "Spatial"
)

# Attach spatial image manually if not already embedded
# seurat_obj[["image"]] <- Read10X_Image(
#   image.dir = "spatial",
#   assay = "Spatial"
# )

# Calculate mitochondrial percentage
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Verify image attachment when loaded
names(seurat_obj@images)   # should list image names e.g. slice1
head(seurat_obj@meta.data) # check metadata: percent.mt, nCount_Spatial, etc.

# Check counts dimension
dim(GetAssayData(seurat_obj, assay = "Spatial", slot = "counts"))

# Check data matrix dimension (normalized)
dim(GetAssayData(seurat_obj, assay = "Spatial", slot = "data"))

# Save Seurat load
saveRDS(seurat_obj, "seurat_initialLoad.rds")

# LOAD POINT -- resume here after a crash:
# setwd("C:/Users/Andrew/Desktop/binned_outputs/square_008um")
# seurat_obj <- readRDS("seurat_initialLoad.rds")


###/////////////////////////////////###
### SECTION 2: SPE DEVELOPMENT      ###
###/////////////////////////////////###

# Extract required components from Seurat
counts <- GetAssayData(seurat_obj, assay = "Spatial", slot = "counts")
coldata <- seurat_obj@meta.data
coords  <- GetTissueCoordinates(seurat_obj)

# Convert coordinates to numeric matrix
spatial_coords <- as.matrix(coords[, c("x", "y")])
rownames(spatial_coords) <- colnames(counts)

# Build SpatialExperiment
spe <- SpatialExperiment(
  assays       = list(counts = counts),
  colData      = coldata,
  spatialCoords = spatial_coords
)

# Rename QC fields for SpotSweeper compatibility
colData(spe)$sum                    <- seurat_obj$nCount_Spatial
colData(spe)$detected               <- seurat_obj$nFeature_Spatial
colData(spe)$subsets_Mito_percent   <- seurat_obj$percent.mt

# Consistency check
stopifnot(all(rownames(spatialCoords(spe)) == colnames(counts)))


###///////////////////////////////////////////###
### SECTION 3: ATTACH IMAGE TO SPE            ###
###///////////////////////////////////////////###

# Ensure spe was created correctly (should return a data frame with 0 rows if image missing)
imgData(spe)

# Use the actual path to the image file
image_path <- "C:/Users/Andrew/Desktop/binned_outputs/square_008um/spatial/tissue_lowres_image.png"

# Read scalefactors
scale_json    <- fromJSON("C:/Users/Andrew/Desktop/binned_outputs/square_008um/spatial/scalefactors_json.json")
scale_factors <- scale_json$tissue_lowres_scalef

# Confirm it's numeric
str(scale_factors)

# Add image to SpatialExperiment
sampleID <- colData(spe)$orig.ident[1]

spe <- addImg(
  spe,
  sample_id   = sampleID,
  image_id    = "lowres_image",
  imageSource = image_path,
  scaleFactor = scale_factors
)

# Verify attachment
imgData(spe)
colnames(seurat_obj@meta.data)

# Save SPE prior to running localOutliers
saveRDS(spe, "spe_preQC.rds")

# LOAD POINT:
# spe <- readRDS("spe_preQC.rds")


###///////////////////////////////////////////###
### SECTION 4: RUNNING SPOTSWEEPER QC         ###
###///////////////////////////////////////////###

# Run localOutliers for mito, sum, detected -- this will take awhile
run_sweep <- function(spe_obj, n) {
  spe_obj <- localOutliers(spe_obj, metric = "subsets_Mito_percent", direction = "higher", log = FALSE, n_neighbors = n, samples = "sample_id")
  spe_obj <- localOutliers(spe_obj, metric = "sum",                   direction = "lower",  log = TRUE,  n_neighbors = n, samples = "sample_id")
  spe_obj <- localOutliers(spe_obj, metric = "detected",              direction = "lower",  log = TRUE,  n_neighbors = n, samples = "sample_id")
  spe_obj$local_outliers <- spe_obj$subsets_Mito_percent_outliers |
    spe_obj$sum_outliers |
    spe_obj$detected_outliers
  return(spe_obj)
}

spe <- run_sweep(spe, 12) # adjust n_neighbors if needed

# Verify output
colnames(colData(spe))
length(spe$subsets_Mito_percent_outliers)
length(spe$sum_outliers)
length(spe$detected_outliers)

# QC summary
cat(
  "Total spots:",   ncol(spe), "\n",
  "Retained spots:", sum(!spe$local_outliers), "\n",
  "Removed spots:",  sum(spe$local_outliers), "\n",
  "Percent removed:", round(mean(spe$local_outliers) * 100, 2), "%\n"
)

# Save SPE with outlier flags before subsetting
saveRDS(spe, "spe_all_outliers.rds")

# LOAD POINT:
# spe <- readRDS("spe_all_outliers.rds")


###///////////////////////////////////////////###
### SECTION 5: QC VISUALIZATION (Optional)   ###
###///////////////////////////////////////////###

# Create factor column for plotting (keep original logical for subsetting)
spe$local_outliers_plot <- factor(
  ifelse(spe$local_outliers, "outlier", "retained"),
  levels = c("retained", "outlier")
)

# Check coordinates and spot names line up
# NOTE: only valid if spe still has the full unsubsetted spots (i.e. not resumed from spe_postQC load point)
stopifnot(all(rownames(spatialCoords(spe)) == colnames(seurat_obj)))

# Extract spatial coordinates and build plotting dataframe
coords_plot <- spatialCoords(spe)
df_qc <- data.frame(
  x      = coords_plot[, "x"],
  y      = coords_plot[, "y"],
  status = spe$local_outliers_plot
)

# Plot retained vs outlier spots
ggplot(df_qc, aes(x = y, y = -x)) +
  geom_point(aes(color = status), size = 0.1) +
  scale_color_manual(values = c("grey80", "#D7301F")) + # retained = gray, outliers = red
  coord_fixed() +
  theme_void() +
  labs(title = "QC of Data Using SpotSweeper", color = "QC Status")


###///////////////////////////////////////////###
### SECTION 6: SUBSET & SAVE CLEAN SPE/SEURAT###
###///////////////////////////////////////////###

# Extract barcodes to keep (non-outlier spots)
spots_to_keep <- colnames(spe)[spe$local_outliers == FALSE]

# Subset original Seurat object
seurat_clean <- subset(seurat_obj, cells = spots_to_keep)

# Checks for proper subsetting
length(colnames(seurat_clean)) # should equal length(spots_to_keep)
length(spots_to_keep)
names(seurat_clean@images)     # should show e.g. "slice1"

# Save clean SPE (for SingleR later) and clean Seurat
saveRDS(spe[, spots_to_keep], "spe_postQC.rds")
saveRDS(seurat_clean, "seurat_subPreProcessing.rds")

# LOAD POINT:
# seurat_clean <- readRDS("seurat_subPreProcessing.rds")


###///////////////////////////////////////////////////////////###
### SECTION 7: NORMALIZATION & DIMENSIONALITY REDUCTION       ###
###///////////////////////////////////////////////////////////###

# Normalize and pre-process
seurat_clean <- NormalizeData(seurat_clean, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_clean <- FindVariableFeatures(seurat_clean, selection.method = "vst", nfeatures = 2000)
seurat_clean <- ScaleData(seurat_clean)
seurat_clean <- RunPCA(seurat_clean)

# Neighborhood graph, clustering, and UMAP
seurat_clean <- FindNeighbors(seurat_clean, dims = 1:15)
seurat_clean <- FindClusters(seurat_clean, resolution = 0.5)
seurat_clean <- RunUMAP(seurat_clean, dims = 1:15)

# Plot clusters in UMAP space
DimPlot(seurat_clean, reduction = "umap", label = TRUE) +
  ggtitle("Seurat Clusters Post SpotSweeper")

# Save fully processed Seurat object
saveRDS(seurat_clean, "seuratPostQC_subAndNormed.rds")

# LOAD POINT:
# seurat_clean <- readRDS("seuratPostQC_subAndNormed.rds")


###///////////////////////////////////////////###
### SECTION 8: FIND CLUSTER MARKERS           ###
###///////////////////////////////////////////###

# Confirm active assay name (will be "Spatial" when coming from seurat_obj subset)
names(seurat_clean@assays) # confirm

# Find all markers
markers <- FindAllMarkers(
  seurat_clean,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25
)

# Save markers -- recommended: downstream sections (9, 13) depend on this object
write.csv(markers, "cluster_markers.csv", row.names = FALSE)

# LOAD POINT (skip FindAllMarkers and load from disk instead):
# markers <- read.csv("cluster_markers.csv")

# Top 10 markers per cluster
top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  pull(gene)


###///////////////////////////////////////////###
### SECTION 9: HEATMAP BUILDING               ###
###///////////////////////////////////////////###

# Scale only top10 features to reduce memory load
active_assay <- DefaultAssay(seurat_clean)
seurat_clean <- ScaleData(seurat_clean, features = unique(top10), assay = active_assay)

# Extract normalized expression for top10 genes (assay-name agnostic)
active_assay <- DefaultAssay(seurat_clean)
expr_mat <- as.matrix(GetAssayData(seurat_clean, assay = active_assay, slot = "data")[unique(top10), ])

# Extract cluster identities
clusters_ids <- Idents(seurat_clean)

# Average expression per cluster
expr_df <- data.frame(cluster = clusters_ids, t(expr_mat))
avg_expr_df <- expr_df %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean))

avg_expr_mat <- t(as.matrix(avg_expr_df[, -1]))
colnames(avg_expr_mat) <- avg_expr_df$cluster

# Z-score scale by gene (row)
avg_expr_z <- t(scale(t(avg_expr_mat)))

# Plot heatmap (first 100 genes)
pheatmap(
  avg_expr_z[1:min(100, nrow(avg_expr_z)), ],
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 6
)

# --- Downsample heatmap with top 5 markers ---
top5 <- markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup()

# Downsample spots per cluster for tractable heatmap
set.seed(123)
spots_subset <- subset(
  seurat_clean,
  cells = unlist(
    lapply(split(Cells(seurat_clean), Idents(seurat_clean)), function(cells) {
      sample(cells, min(length(cells), 100))
    })
  )
)

DoHeatmap(spots_subset, features = unique(top10), group.by = "seurat_clusters")


###///////////////////////////////////////////###
### SECTION 10: SPATIAL CLUSTER PLOTS        ###
###///////////////////////////////////////////###

# All clusters on slide
SpatialDimPlot(
  seurat_clean,
  group.by      = "seurat_clusters",
  pt.size.factor = 1,
  image.alpha   = 0.0,
  alpha         = c(0.9, 0.9)
) + theme(legend.position = "right")

# Top 3 clusters by spot count
top3_clusters <- names(sort(table(seurat_clean$seurat_clusters), decreasing = TRUE))[1:3]
top3_seurat   <- subset(seurat_clean, idents = top3_clusters)

SpatialDimPlot(top3_seurat, label = TRUE, label.size = 5, image.alpha = 0) +
  ggtitle("Spatial Plot of Top 3 Clusters")


###///////////////////////////////////////////###
### SECTION 11: SPOTS PER CLUSTER CHART       ###
###///////////////////////////////////////////###

cluster_counts <- table(Idents(seurat_clean))
df_counts <- as.data.frame(cluster_counts)
colnames(df_counts) <- c("Cluster", "Count")
df_counts <- df_counts %>%
  mutate(Percentage = Count / sum(Count) * 100)

ggplot(df_counts, aes(x = Percentage, y = reorder(Cluster, Percentage))) +
  geom_col(fill = "grey80", color = "black") +
  geom_text(aes(label = sprintf("%.1f%% (n = %s)", Percentage, Count)),
            hjust = -0.1, size = 3) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Percentage of Spots", y = "Cluster") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )


###///////////////////////////////////////////###
### SECTION 12: SINGLER ANNOTATION            ###
###///////////////////////////////////////////###

# Load post-QC SPE
spe_clean <- readRDS("spe_postQC.rds")

# Compute size factors and log-normalize
sf <- colSums(counts(spe_clean))
sf[sf == 0] <- 1
sizeFactors(spe_clean) <- sf
spe_clean <- logNormCounts(spe_clean)
assayNames(spe_clean) # should show "counts" and "logcounts"

# Load reference (broad human atlas; swap for BlueprintEncodeData() if immune-heavy)
ref <- HumanPrimaryCellAtlasData()

# Run SingleR at spot level
singler_res <- SingleR(
  test            = spe_clean,
  ref             = ref,
  labels          = ref$label.main,
  assay.type.test = "logcounts"
)

# Match barcodes and bring labels into Seurat
stopifnot(all(colnames(seurat_clean) %in% colnames(spe_clean)))

singler_labels <- singler_res$labels
names(singler_labels) <- colnames(spe_clean)
seurat_clean$SingleR_label <- singler_labels[colnames(seurat_clean)]
table(seurat_clean$SingleR_label, useNA = "always")

# Compute majority-vote consensus label per cluster
majority_label <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  tbl <- table(x)
  names(tbl)[which.max(tbl)]
}

df_singler <- data.frame(
  cluster = as.character(Idents(seurat_clean)),
  label   = as.character(seurat_clean$SingleR_label)
)

consensus <- df_singler %>%
  group_by(cluster) %>%
  summarise(consensus_label = majority_label(label), .groups = "drop")

print(consensus) # sanity check

# Map consensus back to each spot
map             <- setNames(as.character(consensus$consensus_label), consensus$cluster)
clusters_vec    <- as.character(Idents(seurat_clean))
consensus_labels <- map[clusters_vec]
seurat_clean$SingleR_consensus <- unname(consensus_labels)

table(seurat_clean$SingleR_consensus, clusters_vec, useNA = "always")
table(seurat_clean$SingleR_consensus, useNA = "always")

# Use only the first image (avoids duplicate image issues)
seurat_clean@images <- seurat_clean@images[1]

# SingleR Visualizations

# Generate palette — colorRampPalette fallback handles >12 unique labels safely
labels_sr  <- unique(seurat_clean$SingleR_consensus)
base_colors <- brewer.pal(n = 12, name = "Set3")
my_colors  <- setNames(
  colorRampPalette(base_colors)(length(labels_sr)),
  labels_sr
)

DimPlot(seurat_clean, reduction = "umap", group.by = "SingleR_consensus",
        label = TRUE, cols = my_colors) +
  ggtitle("SingleR Consensus Annotations (per cluster)")

SpatialDimPlot(seurat_clean, group.by = "SingleR_consensus",
               image.alpha = 0, cols = my_colors) +
  ggtitle("Spatial SingleR Consensus Annotations")

# Save annotated Seurat
# saveRDS(seurat_clean, "seurat_with_SingleR_consensus.rds")


###///////////////////////////////////////////###
### SECTION 13: GPT / MANUAL ANNOTATIONS     ###
###///////////////////////////////////////////###

# Export top 5 markers per cluster for manual/GPT annotation
top5_export <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%
  ungroup() %>%
  select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj)

write.xlsx(top5_export, file = "top5_genes_per_cluster.xlsx", rowNames = FALSE)

# Define human-readable cluster labels (update after annotation)
cluster_ids <- c(
  "0"  = "(0) Fibroblasts",
  "1"  = "(1) Basal epithelial",
  "2"  = "(2) Luminal epithelial",
  "3"  = "(3) Tumor-like epithelial / stress",
  "4"  = "(4) Enterocyte-like / metabolic",
  "5"  = "(5) T cells & macrophages",
  "6"  = "(6) Plasma B cells",
  "7"  = "(7) Endothelial cells",
  "8"  = "(8) Cycling / proliferating",
  "9"  = "(9) Mitochondrial / oxidative stress",
  "10" = "(10) Immediate early response",
  "11" = "(11) Noise / low-quality"
)

# High-contrast color palette
cluster_colors <- c(
  "(0) Fibroblasts"                       = "#FF5500",
  "(1) Basal epithelial"                  = "#377EB8",
  "(2) Luminal epithelial"                = "#4DAF4A",
  "(3) Tumor-like epithelial / stress"    = "#E41A1C",
  "(4) Enterocyte-like / metabolic"       = "#984EA3",
  "(5) T cells & macrophages"             = "#F781BF",
  "(6) Plasma B cells"                    = "#A65628",
  "(7) Endothelial cells"                 = "#FFFF33",
  "(8) Cycling / proliferating"           = "#1B9E77",
  "(9) Mitochondrial / oxidative stress"  = "#E7298A",
  "(10) Immediate early response"         = "#66C2FF",
  "(11) Noise / low-quality"              = "#808080"
)

# Map cluster IDs to metadata
seurat_clean$celltype <- plyr::mapvalues(
  x    = as.character(Idents(seurat_clean)),
  from = names(cluster_ids),
  to   = as.vector(cluster_ids)
)

# Restrict palette to observed labels
observed_labels     <- unique(seurat_clean$celltype)
cluster_colors_plot <- cluster_colors[observed_labels]

# Spatial plot with manual annotations (image already trimmed to first in Section 12)
SpatialDimPlot(
  seurat_clean,
  group.by       = "celltype",
  label          = FALSE,
  label.size     = 2.5,
  cols           = cluster_colors_plot,
  repel          = TRUE,
  image.alpha    = 0
) + ggtitle("Cluster 0–11 Plot With GPT Annotations")


###///////////////////////////////////////////###
### SECTION 14: SANKEY DIAGRAM               ###
###///////////////////////////////////////////###

# Seurat cluster → cell type Sankey
meta_df <- as.data.frame(seurat_clean@meta.data)

links_sk <- meta_df %>%
  mutate(seurat_clusters = as.character(seurat_clusters)) %>%
  group_by(seurat_clusters, celltype) %>%
  summarise(value = n(), .groups = "drop")

nodes_sk <- data.frame(
  name  = c(unique(links_sk$seurat_clusters), unique(links_sk$celltype)),
  group = c(
    rep("Seurat",   length(unique(links_sk$seurat_clusters))),
    rep("Celltype", length(unique(links_sk$celltype)))
  )
)

links_sk$IDsource <- match(links_sk$seurat_clusters, nodes_sk$name) - 1
links_sk$IDtarget <- match(links_sk$celltype,        nodes_sk$name) - 1

sankeyNetwork(
  Links      = links_sk,
  Nodes      = nodes_sk,
  Source     = "IDsource",
  Target     = "IDtarget",
  Value      = "value",
  NodeID     = "name",
  NodeGroup  = "group",
  sinksRight = TRUE,
  nodeWidth  = 30,
  fontSize   = 12
)


