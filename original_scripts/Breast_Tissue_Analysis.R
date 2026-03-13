# Breast_Tissue_Analysis.R
# Created by Andrew Sanford on 07/20/2025

# Libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SingleCellExperiment)
library(SpotSweeper)
library(SpatialExperiment)
library(ggspavis)
library(scuttle)
library(scater)
library(scran)
library(igraph)
library(SeuratObject)
library(presto)
library(pheatmap)
library(SingleR)
library(celldex)

# Set working directory
setwd("C:\\Users\\Andrew\\Desktop\\binned_outputs\\square_008um")

### PRE-PROCESSING ###
----------------------

# Load data and attach image with Seurat

seurat_obj <- Load10X_Spatial(
  data.dir = "C:\\Users\\Andrew\\Desktop\\binned_outputs\\square_008um",
  filename = "filtered_feature_bc_matrix.h5",
  image = "tissue_lowres_image.png"
)

# Load post-QC spotsweeper data -- SPE FORMAT
spe_clean <- readRDS("spe_postQC.rds")

# Convert from SpatialExperiment to Seurat
seurat_clean <- as.Seurat(spe_clean, counts = "counts", data = NULL, assay = "Spatial")

if ("data" %in% colnames(imgData(spe))) {
  seurat_obj[["image"]] <- imgData(spe)$data[[1]]
}

# Rerun normalization
seurat_clean <- NormalizeData(seurat_clean, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_clean <- FindVariableFeatures(seurat_clean, selection.method = "vst", nfeatures = 2000)
seurat_clean <- ScaleData(seurat_clean)
seurat_clean <- RunPCA(seurat_clean)

# Run neighborhood graph and clustering
seurat_clean <- FindNeighbors(seurat_clean, dims = 1:15)
seurat_clean <- FindClusters(seurat_clean, resolution = 0.5)

# Run UMAP for visualization
seurat_clean <- RunUMAP(seurat_clean, dims = 1:15)

# Plot clusters in UMAP space
DimPlot(seurat_clean, reduction = "umap", label = TRUE) + ggtitle("Seurat Clusters Post SpotSweeper")

# Save the clean seurat with image attached
#saveRDS(seurat_clean, file = "seurat_clean_with_image.rds")

# Load the clean seurat
seurat_clean <- readRDS(file = "seuratPostQC_subAndNormed")

# Find markers in clusters
markers <- FindAllMarkers(
  seurat_clean,
  only.pos = TRUE,       # only keep positive markers
  min.pct = 0.25,        # gene expressed in at least 25% of spots in that cluster
  logfc.threshold = 0.25 # at least log2FC of 0.25
)

# Save markers for downstream interpretation
#write.csv(markers, "cluster_markers.csv", row.names = FALSE)

# Load markers
markers <- read.csv("cluster_markers.csv")

# Take top 10 markers per cluster and extract gene names
top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  pull(gene)



### HEATMAP BUILDING ###
------------------------
  
# Scale only these features to reduce memory load
seurat_clean <- ScaleData(seurat_clean, features = top10)

# # Plot heatmap -- way too intensive
# DoHeatmap(seurat_clean, features = top10, raster = TRUE) + NoLegend()

# Rename assay type as it was lost in SPE -> Seurat conversion

# Rename the lone assay to "RNA"
# Rename the existing (unnamed) assay to "RNA"
seurat_clean@assays <- list(RNA = seurat_clean@assays[[1]])

# Confirm
names(seurat_clean@assays)


# Extract normalized expression matrix for top10 genes
expr_mat <- as.matrix(seurat_clean@assays$RNA@data[unique(top10), ])

# Extract cluster identities (assumes Idents are set)
clusters <- Idents(seurat_clean)

# Make a data frame: each row = spot with cluster, columns = genes
expr_df <- data.frame(cluster = clusters, t(expr_mat))

# Average expression per cluster
avg_expr_df <- expr_df %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean))

# Convert back to matrix: rows = genes, columns = clusters
avg_expr_mat <- t(as.matrix(avg_expr_df[,-1]))
colnames(avg_expr_mat) <- avg_expr_df$cluster

# Scale by row (genes) if desired
avg_expr_z <- t(scale(t(avg_expr_mat)))

# Plot heatmap
pheatmap(avg_expr_z[1:100, ], # first 40 genes
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_row = 6)


# Filter for strong markers and take top 5 per cluster
top5 <- markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup()

# Scale data for these top marker genes
seurat_clean <- ScaleData(seurat_clean, assay = "Spatial", features = top5$gene)

# Plot heatmap with larger text for readability
DoHeatmap(seurat_obj_subset, assay = "Spatial", features = top5$gene, size = 2.5) +
  theme(axis.text = element_text(size = 10))

# Downsample spots per cluster for plotting
set.seed(123)
spots_subset <- seurat_clean %>%
  subset(
    cells = unlist(
      lapply(split(Cells(seurat_clean), Idents(seurat_clean)), function(cells) {
        sample(cells, min(length(cells), 100))
      })
    )
  )

DoHeatmap(spots_subset, features = top10, group.by = "seurat_clusters")



### PLOTTING CLUSTERS ON SLIDE ###
----------------------------------
  
SpatialDimPlot(
  seurat_clean,
  group.by = "seurat_clusters",
  pt.size.factor = 1,        # smaller points for HD density
  image.alpha = 0.0,         # make histology invisible
  alpha = c(0.9, 0.9)        # control spot transparency
) +
  theme(legend.position = "right")

## top 3 clusters plotted on slide

top3 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 3) %>%
  pull(gene) %>%
  unique()

# Get top 3 clusters by spot count
top3_clusters <- names(sort(table(seurat_clean$seurat_clusters), decreasing = TRUE))[1:3]

# Subset Seurat object to top 3 clusters
top3_seurat <- subset(seurat_clean, idents = top3_clusters)

# Plot spatial clusters on image
SpatialDimPlot(top3_seurat, label = TRUE, label.size = 5, image.alpha = 0) + 
  ggtitle("Spatial plot of Top 3 Clusters")



### SPOTS PER CLUSTER HORIZONTAL CHART ###
------------------------------------------

# Convert table to df with explicit colnames
# Count spots per cluster
cluster_counts <- table(Idents(seurat_clean))

# Convert table to df with explicit colnames
df <- as.data.frame(cluster_counts)
colnames(df) <- c("Cluster", "Count")

# Add percentages
df <- df %>%
  mutate(Percentage = Count / sum(Count) * 100)

# Horizontal bar chart
ggplot(df, aes(x = Percentage, y = reorder(Cluster, Percentage))) +
  geom_col(fill = "grey80", color = "black") +
  geom_text(aes(label = sprintf("%.1f%% (n = %s)", Percentage, Count)),
            hjust = -0.1, size = 3) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Percentage of Spots", y = "Cluster") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())



### ANNOTATING CLUSTERS ###
---------------------------
  
# SingleR on SpatialExperiment
# Install once if needed:
# BiocManager::install(c("SingleR","celldex","scuttle"))

library(SingleR)
library(celldex)
library(scuttle)

# Load SPE
spe_clean <- readRDS("spe_postQC.rds")

# Compute size factors
sf <- colSums(counts(spe_clean))

# Replace zeros with 1 (or small epsilon)
sf[sf == 0] <- 1
sizeFactors(spe_clean) <- sf

# Now log-normalize
spe_clean <- logNormCounts(spe_clean)

# Confirm assays
assayNames(spe_clean)
# Should show: "counts" "logcounts"

# Load a reference (pick one)
# For human data, two common choices:
ref <- HumanPrimaryCellAtlasData()      # broad, good starting point
# ref <- BlueprintEncodeData()          # immune-heavy, more granular for immune contexts

# Run SingleR at the spot level
# Make sure genes overlap; SingleR handles intersect internally
singler_res <- SingleR(
  test = spe_clean,
  ref  = ref,
  labels = ref$label.main,   # or ref$label.fine for finer labels
  assay.type.test = "logcounts"
)

# Bring labels into Seurat
# Load your Seurat object that you cluster/UMAP’d already
seurat_clean <- readRDS("seuratPostQC_subAndNormed")

# Check spot/barcode names should match between objects
stopifnot(all(colnames(seurat_clean) %in% colnames(spe_clean)))

# Assign SingleR labels
singler_labels <- singler_res$labels
names(singler_labels) <- colnames(spe_clean)   # keep same spot/barcode names

# assign to Seurat (order by colnames)
seurat_clean$SingleR_label <- singler_labels[colnames(seurat_clean)]

# check labels
table(seurat_clean$SingleR_label, useNA = "always")


# Build a dataframe of cluster + SingleR label
df <- data.frame(
  cluster = as.character(Idents(seurat_clean)),
  label   = as.character(seurat_clean$SingleR_label)
)


# Majority vote function
majority_label <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  tbl <- table(x)
  names(tbl)[which.max(tbl)]
}


# Compute consensus label per cluster
consensus <- df %>%
  group_by(cluster) %>%
  summarise(consensus_label = majority_label(label), .groups = "drop")

print(consensus)  # sanity check


# Map consensus back to each cell/spot
map <- setNames(as.character(consensus$consensus_label), consensus$cluster)

# convert clusters to characters
clusters <- as.character(Idents(seurat_clean))

# lookup consensus labels
consensus_labels <- map[clusters]

# assign safely to meta.data
seurat_clean$SingleR_consensus <- unname(consensus_labels)

# check per cluster
table(seurat_clean$SingleR_consensus, clusters, useNA = "always")

# check overall counts
table(seurat_clean$SingleR_consensus, useNA = "always")

# Visualize annotations using cluster-level consensus labels
# use only one image (from beta version of spe w/ two same images)
seurat_clean@images <- seurat_clean@images[1]


# UMAP with cluster-level consensus labels
DimPlot(seurat_clean, reduction = "umap", group.by = "SingleR_consensus", label = TRUE) +
  ggtitle("SingleR consensus annotations (per cluster)")

# Spatial view
SpatialDimPlot(seurat_clean, group.by = "SingleR_consensus", image.alpha = 0) +
  ggtitle("Spatial SingleR consensus annotations")

# Save annotated Seurat for later
#saveRDS(seurat_clean, "seurat_with_SingleR_consensus.rds")



# make it pretty
library(RColorBrewer)

# Get unique consensus labels
labels <- unique(seurat_clean$SingleR_consensus)

# Generate a palette with enough distinct colors
my_colors <- setNames(brewer.pal(n = max(3, length(labels)), name = "Set3")[1:length(labels)], labels)

# UMAP with cluster-level consensus labels and custom colors
DimPlot(seurat_clean, reduction = "umap", group.by = "SingleR_consensus", label = TRUE, cols = my_colors) +
  ggtitle("SingleR consensus annotations (per cluster)")

# Spatial view with custom colors
SpatialDimPlot(seurat_clean, group.by = "SingleR_consensus", image.alpha = 0, cols = my_colors) +
  ggtitle("Spatial SingleR consensus annotations")



### GPT Annotations ###

library(dplyr)
library(openxlsx)  # for Excel output

# Filter for strong markers and take top 5 per cluster
top5 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%
  ungroup()

# Optional: keep only relevant columns
top5_export <- top5 %>%
  select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj)

# Write to Excel
write.xlsx(top5_export, file = "top5_genes_per_cluster.xlsx", rowNames = FALSE)


# Keep only the first image
seurat_clean@images <- seurat_clean@images[1]

# Define human-readable cluster labels
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

# Define high-contrast colors
cluster_colors <- c(
  "(0) Fibroblasts" = "#FF5500",
  "(1) Basal epithelial" = "#377EB8",
  "(2) Luminal epithelial" = "#4DAF4A",
  "(3) Tumor-like epithelial / stress" = "#E41A1C",
  "(4) Enterocyte-like / metabolic" = "#984EA3",
  "(5) T cells & macrophages" = "#F781BF",
  "(6) Plasma B cells" = "#A65628",
  "(7) Endothelial cells" = "#FFFF33",
  "(8) Cycling / proliferating" = "#1B9E77",
  "(9) Mitochondrial / oxidative stress" = "#E7298A",
  "(10) Immediate early response" = "#66C2FF",
  "(11) Noise / low-quality" = "#808080"
)

# Map cluster IDs to metadata column safely
library(plyr)
seurat_clean$celltype <- plyr::mapvalues(
  x = as.character(Idents(seurat_clean)),
  from = names(cluster_ids),
  to   = as.vector(cluster_ids)
)

# Keep only colors for labels actually present
observed_labels <- unique(seurat_clean$celltype)
cluster_colors_plot <- cluster_colors[observed_labels]

# Plot
SpatialDimPlot(
  seurat_clean,
  group.by = "celltype",
  label = FALSE,
  label.size = 2.5,
  cols = cluster_colors_plot,
  repel = TRUE,
  image.alpha = 0
) + ggtitle("Cluster 0-11 Plot With GPT Annotations")


### ----- Sankey Seurat ----- ###
library(dplyr)
library(networkD3)
library(dplyr)
library(networkD3)

# Convert meta.data to standard data frame
meta_df <- as.data.frame(seurat_clean@meta.data)

# Build links table
links <- meta_df %>%
  mutate(seurat_clusters = as.character(seurat_clusters)) %>%
  group_by(seurat_clusters, celltype) %>%
  summarise(value = n(), .groups = "drop")

# Build nodes table
nodes <- data.frame(
  name = c(unique(links$seurat_clusters), unique(links$celltype)),
  group = c(rep("Seurat", length(unique(links$seurat_clusters))),
            rep("Celltype", length(unique(links$celltype))))
)

# Map labels to indices (zero-based)
links$IDsource <- match(links$seurat_clusters, nodes$name) - 1
links$IDtarget <- match(links$celltype, nodes$name) - 1

# Sankey plot
sankeyNetwork(
  Links = links, Nodes = nodes,
  Source = "IDsource", Target = "IDtarget",
  Value = "value", NodeID = "name",
  NodeGroup = "group",
  sinksRight = TRUE,
  nodeWidth = 30, fontSize = 12
)









