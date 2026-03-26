# Application of SpotSweeper Dynamic QC in a Spatial Transcriptomic Pipeline of Visium HD Data

**Author:** Andrew Sanford, sanfordac2@vcu.edu/sanford.andrew.520@gmail.com  
**Institution:** Virginia Commonwealth University  
**Data:** 10x Visium HD -- Breast Tissue, 8 µm bins

---

## Overview

This repository contains a full end-to-end R pipeline for spatial transcriptomics analysis of human breast tissue using 10x Visium HD data. The pipeline covers everything from raw data loading and quality control through clustering, cell type annotation, and visualization.

The pipeline is divided into two major phases:

1. **QC Phase** -- Spatially-aware quality control using [SpotSweeper](https://github.com/MicTott/SpotSweeper) to flag and remove low-quality spots based on local neighborhood statistics rather than global thresholds.
2. **Analysis Phase** -- Normalization, dimensionality reduction, clustering, marker detection, cell type annotation via SingleR, and a suite of visualizations.

---

## Pipeline Sections

| Section | Description |
|---|---|
| 1 | Load 10x Visium HD data with Seurat |
| 2 | Build SpatialExperiment (SPE) object |
| 3 | Attach tissue image to SPE |
| 4 | Run SpotSweeper QC (`localOutliers`) |
| 5 | Visualize QC outliers on tissue |
| 6 | Subset clean spots, save SPE + Seurat |
| 7 | Normalization, PCA, UMAP, clustering |
| 8 | Find cluster marker genes |
| 9 | Heatmap of top marker expression |
| 10 | Spatial cluster plots on tissue image |
| 11 | Spots-per-cluster bar chart |
| 12 | SingleR automated cell type annotation |
| 13 | Manual / GPT-assisted cluster annotation |
| 14 | Sankey diagram: clusters -> cell types |

---

## Requirements

**R version:** 4.3+

Install Bioconductor packages:
```r
BiocManager::install(c(
  "SpotSweeper", "SpatialExperiment", "SingleCellExperiment",
  "ggspavis", "scuttle", "scater", "scran", "SingleR", "celldex"
))
```

Install CRAN packages:
```r
install.packages(c(
  "Seurat", "ggplot2", "patchwork", "dplyr", "jsonlite",
  "SeuratObject", "presto", "pheatmap", "igraph",
  "networkD3", "openxlsx", "RColorBrewer"
))
```

---

## Input Data

Expected directory structure from a 10x Visium HD output:

```
square_008um/
├── filtered_feature_bc_matrix.h5
└── spatial/
    ├── tissue_lowres_image.png
    ├── tissue_positions.csv
    └── scalefactors_json.json
```

Update the `setwd()` and file path variables at the top of the script to point to your local data directory before running.

---

## Save / Load Points

The script includes `saveRDS` checkpoints throughout so you can resume after a crash without re-running expensive steps. Each save point has a corresponding commented-out `readRDS` load line directly below it. Key checkpoints:

| File | Description |
|---|---|
| `seurat_initialLoad.rds` | Raw Seurat object after loading |
| `spe_preQC.rds` | SPE before SpotSweeper |
| `spe_all_outliers.rds` | SPE with outlier flags |
| `spe_postQC.rds` | Clean SPE (retained spots only) |
| `seurat_subPreProcessing.rds` | Clean Seurat before normalization |
| `seuratPostQC_subAndNormed.rds` | Fully processed Seurat |
| `cluster_markers.csv` | FindAllMarkers output |
| `top5_genes_per_cluster.xlsx` | Top 5 markers per cluster for annotation |

> `.rds` and `.h5` files are excluded from this repo via `.gitignore` due to file size. Only the pipeline script is tracked.

---

## Cell Type Annotations

Automated annotation is performed with **SingleR** using the Human Primary Cell Atlas as a reference. Cluster-level consensus labels are assigned by majority vote across spots in each cluster.

Manual annotations in Section 13 reflect GPT-assisted interpretation of top marker genes and are specific to this dataset — update `cluster_ids` and `cluster_colors` as needed for your own data.

---

## Notes

- This pipeline was developed for **8 µm bin** Visium HD data. SpotSweeper neighbor count (`n_neighbors = 12`) may need tuning for other bin sizes or tissue types.
- `dplyr` must be loaded before any package that masks `summarise` (e.g. `plyr`). The script uses explicit `plyr::mapvalues()` calls to avoid conflicts.
- The `colorRampPalette` fallback in Section 12 handles datasets with more than 12 unique SingleR labels gracefully.

- The main pipeline was merged with tools and the original scripts are in the subfolder.

## Special Thanks

A very special thanks to [Dr. Katarzyna Tyc](https://github.com/katarzynatyc) for mentoring me and supporting me during the creation of this pipeline.

A thank you to my collaborators on the T.Y.C. Genomics Team
- Francisco Meersohm
- Maxwell Moberg
- Sameen 

Thank you to my wife Baeleigh for keeping me sane during debugging <3

