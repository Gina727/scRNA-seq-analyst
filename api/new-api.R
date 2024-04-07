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
user_key <- NULL

api_key_list <- c()

#* @get /userkey
userkey <- function() {
  user_key <<- UUIDgenerate()
  api_key_list <<- c(api_key_list, user_key)
  return(user_key)
}

#* @param f:file
#* @parser multi
#* @parser rds
#* @post /upload_files
#* @serializer print
# add header 
upload_files <- function(f, HTTP_API_KEY) {
  uploaded_file <<- f[[1]]
  res$setHeader("HTTP_API_KEY", HTTP_API_KEY)
}

# * check api key
# * @post /check
# check <- function(api_key){
  # Extract the API key from request headers
  # if (is.null(api_key) || !(api_key %in% api_key_list)) {
  #   res$status <- 401  # Unauthorized status code
  #   list(error = "Invalid API key")
  # }

#* show quality control
#* @post /qcplot
#* @serializer png
# find local folder 
qcplot <- function(){
  seurat_obj <- uploaded_file
  if(is.null(seurat_obj)){
    stop("No Seurat object uploaded", call. = FALSE)
  }
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  d = VlnPlot(seurat_obj, features = c("nFeature_RNA", "percent.mt"), ncol = 2, pt.size = 0.1)
  print(d)

  uploaded_file <<- seurat_obj
}

#* quality control
#* @serializer print
#* @get /qc
qc <- function(min.features, max.features, max.mtpercent, res){
  seurat_obj <- uploaded_file
  min.features <- as.numeric(min.features)
  max.features <- as.numeric(max.features)
  max.mtpercent <- as.numeric(max.mtpercent)
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > min.features & nFeature_RNA < max.features & percent.mt < max.mtpercent)
  
  uploaded_file <<- seurat_obj

  max.features
}

#* normalization and run PCA reduction
#* @serializer png
#* @post /norm_pca
norm_pca <- function(scaling_factor, num_hvgs, norm_method, hvg_method, res){
  seurat_obj <- uploaded_file

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

  uploaded_file <<- seurat_obj
}

#* clustering_umap
#* @post /clustering_umap
clustering_umap <- function(dim, resolution){
  seurat_obj <- uploaded_file

  dim <- as.numeric(dim)
  resolution <- as.numeric(resolution)

  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:dim)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:dim)

  uploaded_file <<- seurat_obj
  SaveSeuratRds(seurat_obj, file='clustering_umap.RDS')
}

#* clustering_tsne
#* @serializer png
#* @post /clustering_tsne
clustering_tsne <- function(dims, resolution){
  seurat_obj <- uploaded_file

  dims <- as.numeric(dims)
  resolution <- as.numeric(resolution)

  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:dims)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
  seurat_obj <- RunTSNE(seurat_obj, dims = 1:dims)
  
  t = TSNEPlot(seurat_obj)
  print(t)
  uploaded_file <<- seurat_obj

}

#* annotation_sctype_umap
#* @serializer png
#* @post /annotation_sctype_umap
annotation_sctype_umap <- function(tissue){

  seurat_obj <- uploaded_file

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
    
    d <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')+ggtitle(tissue)
    print(d)

    uploaded_file <<- seurat_obj
}

#* annotation_sctype_tsne
#* @serializer png
#* @post /annotation_sctype_tsne
annotation_sctype_tsne <- function(tissue){

  seurat_obj <- uploaded_file

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
    
    b = TSNEPlot(seurat_obj, group.by = 'customclassif')+ggtitle(tissue)
    print(b)

    uploaded_file <<- seurat_obj
}

# download
