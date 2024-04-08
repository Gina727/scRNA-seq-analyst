library(plumber)

#* @apiTitle scRNA-seq api

  library(Seurat)
  library(HGNChelper)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(openxlsx)
  # load gene set preparation function
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  # load cell type annotation function
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
  db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"

uploaded_file <- NULL 

library(uuid)
key <- NULL

#* @get /userkey
userkey <- function(key, res) {
  key <<- key
  if (is.null(key)) {
    res$status <- 400  # Bad Request
    return(list(error = "The api key is missing."))
  }
  res$setHeader("X-API-KEY", key)
}

#* @param f:file
#* @parser multi
#* @parser rds
#* @post /upload_files
#* @serializer print
# add header 
function(f, res) {
  file_name <- paste0(key, '.RDS')
  uploaded_file <<- f[[1]]
  f
  SaveSeuratRds(uploaded_file, file = file_name)
  res$setHeader("X-API-KEY", key)
}


#* show quality control
#* @post /qcplot
#* @serializer png
# find local folder 
qcplot <- function(req, res){
  seurat_obj <- readRDS(paste0(key, ".RDS"))
  if(is.null(seurat_obj)){
    stop("No Seurat object uploaded", call. = FALSE)
  }
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  file_name <- paste0("qcplot-", key, '.RDS')
  SaveSeuratRds(seurat_obj, file = file_name)

  d = VlnPlot(seurat_obj, features = c("nFeature_RNA", "percent.mt"), ncol = 2, pt.size = 0.1)
  print(d)

  res$setHeader("X-API-KEY", key)
}

#* quality control
#* @serializer print
#* @get /qc
qc <- function(min.features, max.features, max.mtpercent){
  seurat_obj <- readRDS(paste0("qcplot-", key, ".RDS"))
  min.features <- as.numeric(min.features)
  max.features <- as.numeric(max.features)
  max.mtpercent <- as.numeric(max.mtpercent)
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > min.features & nFeature_RNA < max.features & percent.mt < max.mtpercent)
  
  file_name <- paste0("qc-", key, '.RDS')
  SaveSeuratRds(seurat_obj, file = file_name)
}

#* normalization and run PCA reduction
#* @serializer png
#* @post /norm_pca
norm_pca <- function(scaling_factor, num_hvgs, norm_method, hvg_method, res){

  seurat_obj <- readRDS(paste0("qc-", key, ".RDS"))
  scaling_factor <- as.numeric(scaling_factor)
  num_hvgs <- as.numeric(num_hvgs)

  seurat_obj <- NormalizeData(seurat_obj, normalization.method = norm_method, scale.factor = scaling_factor)
    
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = hvg_method, nfeatures = num_hvgs)
    
  # Scale data
  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))

  # Run PCA
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

  d = ElbowPlot(object = seurat_obj, ndims = 40)
  print(d)

  file_name <- paste0("norm_pca-", key, '.RDS')
  SaveSeuratRds(seurat_obj, file = file_name)
  res$setHeader("X-API-KEY", key)
}

#* clustering_umap
#* @post /clustering_umap
clustering_umap <- function(dim, resolution, res){
  
  seurat_obj <- readRDS(paste0("norm_pca-", key, ".RDS"))

  dim <- as.numeric(dim)
  resolution <- as.numeric(resolution)

  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:dim)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:dim)

  file_name <- paste0("clustering_umap-", key, '.RDS')
  SaveSeuratRds(seurat_obj, file = file_name)

  res$setHeader("X-API-KEY", key)
}

#* clustering_tsne
#* @serializer png
#* @post /clustering_tsne
clustering_tsne <- function(dims, resolution, res){
  seurat_obj <- readRDS(paste0("norm_pca-", key, ".RDS"))

  dims <- as.numeric(dims)
  resolution <- as.numeric(resolution)

  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:dims)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
  seurat_obj <- RunTSNE(seurat_obj, dims = 1:dims)

  file_name <- paste0("clustering_tsne-", key, '.RDS')
  SaveSeuratRds(seurat_obj, file = file_name)

  t = TSNEPlot(seurat_obj)
  print(t)

  res$setHeader("X-API-KEY", key)
}

#* annotation_sctype_umap
#* @serializer png
#* @post /annotation_sctype_umap
annotation_sctype_umap <- function(tissue, res){

  seurat_obj <- readRDS(paste0("clustering_umap-", key, ".RDS"))

  # prepare gene sets
    gs_list = gene_sets_prepare(db_, tissue)
    es.max = sctype_score(scRNAseqData = seurat_obj[["RNA"]]$scale.data, scaled = TRUE, 
                          gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
  # merge by cluster
    cL_resutls = do.call("rbind", lapply(unique(seurat_obj@meta.data$seurat_clusters), function(cl){
      es.max.cl = sort(rowSums(es.max[ ,rownames(seurat_obj@meta.data[seurat_obj@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
      head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_obj@meta.data$seurat_clusters==cl)), 10)
    }))
    sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

    # set low-confident (low ScType score) clusters to "unknown"
    sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"

    #match cell types with clusters
    seurat_obj@meta.data$customclassif = ""
    for(j in unique(sctype_scores$cluster)){
      cl_type = sctype_scores[sctype_scores$cluster==j,]; 
      seurat_obj@meta.data$customclassif[seurat_obj@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
    }

    file_name <- paste0("annotation_umap-", key, '.RDS')
    SaveSeuratRds(seurat_obj, file = file_name)

    d <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')+ggtitle(tissue)
    print(d)

    res$setHeader("X-API-KEY", key)
}

#* annotation_sctype_tsne
#* @serializer png
#* @post /annotation_sctype_tsne
annotation_sctype_tsne <- function(tissue, res){

  seurat_obj <- readRDS(paste0("clustering_tsne-", key, ".RDS"))

  # prepare gene sets
    gs_list = gene_sets_prepare(db_, tissue)
    es.max = sctype_score(scRNAseqData = seurat_obj[["RNA"]]$scale.data, scaled = TRUE, 
                          gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
    
    # merge by cluster
    cL_resutls = do.call("rbind", lapply(unique(seurat_obj@meta.data$seurat_clusters), function(cl){
      es.max.cl = sort(rowSums(es.max[ ,rownames(seurat_obj@meta.data[seurat_obj@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
      head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_obj@meta.data$seurat_clusters==cl)), 10)
    }))
    sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
    
    # set low-confident (low ScType score) clusters to "unknown"
    sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
    
    #match cell types with clusters
    seurat_obj@meta.data$customclassif = ""
    for(j in unique(sctype_scores$cluster)){
      cl_type = sctype_scores[sctype_scores$cluster==j,]; 
      seurat_obj@meta.data$customclassif[seurat_obj@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
    }
    file_name <- paste0("annotation_tsne", key, '.RDS')
    SaveSeuratRds(seurat_obj, file = file_name)

    b = TSNEPlot(seurat_obj, group.by = 'customclassif')+ggtitle(tissue)
    print(b)

    res$setHeader("X-API-KEY", key)
}

# download
