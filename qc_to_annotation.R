BiocManager::install("SingleCellExperiment")
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)

#load 10x data
counts.data <- Read10X(data.dir = "data/")

#create a Seurat object
counts <- CreateSeuratObject(counts = counts.data,
                           min.features = 100)
#see meta.data
head(counts@meta.data) 

#calculate mitochondria gene percentage and add it to the meta data
counts[["percent.mt"]] <- PercentageFeatureSet(counts, pattern = "^MT-")

# Add number of genes per UMI for each cell to metadata
counts[["log10GenesPerUMI"]] <- log10(counts$nFeature_RNA) / log10(counts$nCount_RNA)

# make QC based on visualizations:
metadata <- counts@meta.data
CalculateBarcodeInflections(counts)

# UMI (counts) per cell

# Filter out low quality reads using selected thresholds (which will change)
counts_filtered <- subset(counts, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) 

# Create .RData object to load at any time
save(counts_filtered, file="data/counts_filtered.RData")

## Normalization, Scaling, and PCA
# normalization with logarithm
counts_normalized <- NormalizeData(counts_filtered, normalization.method = "LogNormalize", scale.factor = 10000)

counts_normalized <- FindVariableFeatures(counts_normalized, selection.method = "vst", nfeatures = 2000)
# FindVariableFeatures = find unusual features on the mean variability plot

# Scale data
counts_scaled <- ScaleData(counts_normalized, features = rownames(counts_normalized))

# Run PCA
counts_PCA <- RunPCA(counts_scaled, features = VariableFeatures(object = counts_scaled))

########

# Cell cycle scoring

# counts_cellcycle <- CellCycleScoring(counts_PCA, g2m.features=g2m_genes, s.features=s_genes)

# scTransform
counts_transformed<- SCTransform(counts_PCA, vars.to.regress = c("percent.mt"))

#Integration for samples/ datasets
#integ_features <- SelectIntegrationFeatures(object.list = counts_transformed, nfeatures = 3000) 
# PCAPlot(counts_PCA)

## Get no. of PC

# get pca from elbow plot
ElbowPlot(object = counts_transformed, 
          ndims = 40)

# Determine percent of variation associated with each PC
pct <- counts_transformed[["pca"]]@stdev / sum(counts_transformed[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

co1

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2

# Minimum of the two calculation
pcs <- min(co1, co2)
pcs

## Cluster the cell



