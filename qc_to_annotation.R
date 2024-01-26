BiocManager::install("SingleCellExperiment")
BiocManager::install("HGNChelper")
install.packages("writexl")
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(writexl)
library(HGNChelper)

#load 10x data
counts.data <- Read10X(data.dir = "data/")

#create a Seurat object
counts <- CreateSeuratObject(counts = counts.data,  min.cells = 3, min.features = 200)
#see meta.data
#head(counts@meta.data) 

#calculate mitochondria gene percentage and add it to the meta data
counts[["percent.mt"]] <- PercentageFeatureSet(counts, pattern = "^MT-")

# Add number of genes per UMI for each cell to metadata
counts[["log10GenesPerUMI"]] <- log10(counts$nFeature_RNA) / log10(counts$nCount_RNA)

## make QC based on visualizations:
#metadata <- counts@meta.data
#CalculateBarcodeInflections(counts)

#Violin plot of the features in meta data
VlnPlot(counts, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Get your IQR (Interquartile range) and lower/upper quartile using:
percent.mt = c(counts[["percent.mt"]])
percent.mt <- unlist(percent.mt) # unlist to make it atomic, contains only one type of element

nfeature = c(counts[["nFeature_RNA"]])
nfeature <- unlist(nfeature)

#1 Get values of Q1, Q3, and IQR
q1.mt = quantile(percent.mt, 0.25)
q3.mt = quantile(percent.mt, 0.75)

q1.feature = quantile(nfeature, 0.25)
q3.feature = quantile(nfeature, 0.75)

#2 get IQR (= q3-q1)
iqr.mt = IQR(percent.mt) 
iqr.feature = IQR(nfeature)

#3 get threshold values for outliers
min_mt <- q1.mt - 1.5 * iqr.mt
max_mt <- q3.mt + 1.5 * iqr.mt

min_feature <- q1.feature - 1.5 * iqr.feature
max_feature <- q3.feature + 1.5 * iqr.feature

#Correlation (scatter) plot of the counts and features / percentage of mito
#library(patchwork)
#plot1 <- FeatureScatter(counts, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(counts, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2

# Filter out low quality reads using selected thresholds (which will change)
counts <- subset(counts, subset = nFeature_RNA > min_feature & nFeature_RNA < 5000 & percent.mt < max_mt) 

# Create .RData object to load at any time
# save(counts_filtered, file="data/counts_filtered.RData")

## Normalization, Scaling, and PCA
# normalization with logarithm
counts <- NormalizeData(counts, normalization.method = "LogNormalize", scale.factor = 10000)

counts <- FindVariableFeatures(counts, selection.method = "vst", nfeatures = 2000)
# FindVariableFeatures = find unusual features on the mean variability plot

# Scale data
counts <- ScaleData(counts, features = rownames(counts))

# Run PCA
counts <- RunPCA(counts, features = VariableFeatures(object = counts))

########

# Cell cycle scoring

# counts_cellcycle <- CellCycleScoring(counts_PCA, g2m.features=g2m_genes, s.features=s_genes)

# scTransform
#counts_transformed<- SCTransform(counts, vars.to.regress = c("percent.mt"))

#Integration for samples/ datasets
#integ_features <- SelectIntegrationFeatures(object.list = counts_transformed, nfeatures = 3000) 
# PCAPlot(counts_PCA)

## Get no. of PC

# get pca from elbow plot
#ElbowPlot(object = counts, ndims = 40)

# Determine percent of variation associated with each PC
#pct <- counts[["pca"]]@stdev / sum(counts[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
#cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
#co1 <- which(cumu > 90 & pct < 5)[1]

#co1

# Determine the difference between variation of PC and subsequent PC
#co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
#co2

# Minimum of the two calculation
#pcs <- min(co1, co2)
#pcs

## Cluster the cell
# DefaultAssay(counts_transformed) <- "integrated"
counts <- FindNeighbors(counts, dims = 1:10) #dims = dimensions
counts <- FindClusters(counts, resolution = 0.8)
counts <- RunUMAP(counts, dims = 1:10)
DimPlot(counts, reduction = "umap")

## Annotation

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")

# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# load database and select desired tissue
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Immune system"

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

# get cell type score
es.max = sctype_score(scRNAseqData = counts[["RNA"]]$scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
print(es.max[,1:3])

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(counts@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(counts@meta.data[counts@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(counts@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  


# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

#match cell types with clusters
counts@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  counts@meta.data$customclassif[counts@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
p <- DimPlot(counts, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif',) 
p <- p + ggtitle("Immune Cells")
print(p)

View(sctype_scores$cluster)
