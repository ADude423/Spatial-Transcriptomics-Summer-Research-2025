# SpotSweeper_Cleaning.R
# Created by Andrew Sanford on 06/2025
# sanfordac2@vcu.edu

# This script contains save points -- uncomment load points for crashes and resume

# Required packages
# BiocManager::install("SpotSweeper")
# BiocManager::install("SpatialExperiment")
# install.packages("Seurat")
# install.packages("ggplot2")

# Load Packages
library(Seurat)
library(SpotSweeper)
library(SpatialExperiment)
library(ggplot2)
library(jsonlite)


###////////////////////////////////////////////////////###
### LOADING DATA WITH SEURAT -- skip if already loaded ###
###////////////////////////////////////////////////////###

# Set working directory to where your 8 µm (or other bins) data is located
setwd("C:/Users/Andrew/Desktop/binned_outputs/square_008um")

# Load Seurat object from the spatial dataset
seurat_obj <- Load10X_Spatial(
  data.dir = ".",
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
names(seurat_obj@images)  # should list image names slice1
head(seurat_obj@meta.data)     # Check metadata present (percent.mt, nCount_Spatial, etc)

# Check counts dimension
dim(GetAssayData(seurat_obj, assay = "Spatial", slot = "counts"))

# Check data matrix dimension (normalized)
dim(GetAssayData(seurat_obj, assay = "Spatial", slot = "data"))

# Save Seurat load
saveRDS(seurat_obj, "seurat_initialLoad.rds")

# Load Seurat load -- will only work if working directory is where the save is --I saved in the same folder I was working in

# Set working directory to where your 8 µm (or other bins) data is located
# setwd("C:/Users/Andrew/Desktop/binned_outputs/square_008um")
# seurat_obj <- readRDS("seurat_initialLoad.rds")


###/////////////////////////////////###
### SPE DEVELOPMENT FOR SPOTSWEEPER ###
###/////////////////////////////////###

# --- Extract required components from Seurat ---
counts <- GetAssayData(seurat_obj, assay = "Spatial", slot = "counts")
coldata <- seurat_obj@meta.data
coords <- GetTissueCoordinates(seurat_obj)

# --- Convert coordinates to numeric matrix ---
spatial_coords <- as.matrix(coords[, c("x", "y")])
rownames(spatial_coords) <- colnames(counts)

# --- Build SpatialExperiment ---
spe <- SpatialExperiment(
  assays = list(counts = counts),
  colData = coldata,
  spatialCoords = spatial_coords
)

# --- Rename QC fields for SpotSweeper compatibility ---
colData(spe)$sum <- seurat_obj$nCount_Spatial
colData(spe)$detected <- seurat_obj$nFeature_Spatial
colData(spe)$subsets_Mito_percent <- seurat_obj$percent.mt

# --- Optional: consistency check ---
stopifnot(all(rownames(spatialCoords(spe)) == colnames(counts)))


###///////////////////////////////////////////###
### RUNNING SPOTSWEEPER ON SPATIAL EXPERIMENT ###
###///////////////////////////////////////////###

# Ensure spe (SpatialExperiment) was created correctly
imgData(spe) # if returned dataFrame with 0 rows and 1 column, image is missing

# Use the actual path to the image file
image_path <- "C:/Users/Andrew/Desktop/binned_outputs/square_008um/spatial/tissue_lowres_image.png"

# Get scaleFactor from image
img <- Read10X_Image(
  image.dir = "C:/Users/Andrew/Desktop/binned_outputs/square_008um/spatial",
  image.name = "tissue_lowres_image.png"
)

# Manually load scalefactors_json.json
scale_json <- fromJSON("C:/Users/Andrew/Desktop/binned_outputs/square_008um/spatial/scalefactors_json.json")

# Extract lowres scaleFactor
scale_factors <- scale_json$tissue_lowres_scalef

# Confirm it's numeric
str(scale_factors)

# Add image to SpatialExperiment
# - sample_id: arbitrary identifier for this sample, must match the 'samples' column if used in SpotSweeper
# - image_id: internal name for the image, used to reference this image within the SPE object
sampleID <- colData(spe)$orig.ident[1]  # just take the first sample if one

spe <- addImg(
  spe,
  sample_id = sampleID,
  image_id = "lowres_image",    # arbitrary but unique label for this image
  imageSource = image_path,
  scaleFactor = scale_factors
)

# Verify attachment
imgData(spe)

# Verify naming
colnames(seurat_obj@meta.data)

# Should return: [1] "orig.ident"          "nCount_Spatial"      "nFeature_Spatial"    "percent.mt"         
# [5] "Spatial_snn_res.0.5" "seurat_clusters"  

# Save SPE prior to running localOutliers
saveRDS(spe, "spe_preQC.rds")

# Load SPE before running localOutliers -- see first save spot for help
# spe <- readRDS("spe_preQC.rds")

# run localOutliers for mito, sum, detected -- this will take awhile
run_sweep <- function(spe_obj, n){
  spe_obj <- localOutliers(spe_obj, metric = "subsets_Mito_percent", direction = "higher", log = FALSE, n_neighbors = n, samples = "sample_id")
  spe_obj <- localOutliers(spe_obj, metric = "sum", direction = "lower", log = TRUE, n_neighbors = n, samples = "sample_id")
  spe_obj <- localOutliers(spe_obj, metric = "detected", direction = "lower", log = TRUE, n_neighbors = n, samples = "sample_id")
  spe_obj$local_outliers <- spe_obj$subsets_Mito_percent_outliers | spe_obj$sum_outliers | spe_obj$detected_outliers
  return(spe_obj)
}

spe <- run_sweep(spe, 12) # Change neighbors if needed

# Check data
colnames(colData(spe))
colnames(spe)

# Verify lengths
length(spe$subsets_Mito_percent_outliers)
length(spe$sum_outliers)
length(spe$detected_outliers)

# Summary of QC
cat(
  "Total spots:", ncol(spe), "\n",
  "Retained spots:", sum(!spe$local_outliers), "\n",
  "Removed spots:", sum(spe$local_outliers), "\n",
  "Percent removed:", round(mean(spe$local_outliers) * 100, 2), "%"
)

#save spe prior to subsetting in case of memory issues
saveRDS(spe, "spe_all_outliers")

#load the spe
#spe <- readRDS("spe_all_outliers")

# Extract barcodes/spots to keep (retained spots)
spots_to_keep <- colnames(spe)[spe$local_outliers == FALSE]  

# Subset original Seurat object to keep only retained spots
seurat_clean <- subset(seurat_obj, cells = spots_to_keep)

# Checks for proper subsetting
length(colnames(seurat_clean))  # should be length of spots_to_keep
length(spots_to_keep)


# Confirm image is still attached
names(seurat_clean@images)  # should show images like "slice1"
seurat_clean@images$image   # check the image object exists

# Save prior to preprocessing
saveRDS(seurat_clean, "seurat_subPreProcessing")

# Load pre processed and subsetted seurat
# seurat_clean <- readRDS("seurat_subPreProcessing")

# Data normalization and pre-processing
seurat_clean <- RunPCA(seurat_clean)
seurat_clean <- FindNeighbors(seurat_clean, dims = 1:15)
seurat_clean <- FindClusters(seurat_clean, resolution = 0.5)
seurat_clean <- RunUMAP(seurat_clean, dims = 1:15)

# save final product for 12 neighbors
saveRDS(seurat_clean, "seuratPostQC_subAndNormed")


###////////////////////////////////////////////###
### VISUALIZATION OF FAILED QC DATA (Optional) ###
###////////////////////////////////////////////###

# Keep the original logical local_outliers column for subsetting
# Create a new factor column for plotting
spe$local_outliers_plot <- factor(
  ifelse(spe$local_outliers, "outlier", "retained"),
  levels = c("retained", "outlier")
)

# Sanity check: coordinates and spot names line up
stopifnot(all(rownames(spatialCoords(spe)) == colnames(seurat_obj)))

# Extract spatial coordinates
coords <- spatialCoords(spe)

# Build a plotting dataframe
df <- data.frame(
  x = coords[, "x"],
  y = coords[, "y"],
  status = spe$local_outliers_plot
)

# Plot retained vs outlier spots
ggplot(df, aes(x = y, y = -x)) +  # rotate clockwise if needed
  geom_point(aes(color = status), size = 0.1) +
  scale_color_manual(
    values = c("grey80", "#D7301F")  # retained = gray, outliers = red
  ) +
  coord_fixed() +
  theme_void() +
  labs(title = "QC of Data Using SpotSweeper", color = "QC Status")

