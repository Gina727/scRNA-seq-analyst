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
destfile <- NULL
library(uuid)
key <- NULL

#* @get /check_files
#* This api is used to check the files in the directory and return the latest stage
check_files <- function(req, res){
    key <- req$HTTP_KEY
    files <- list.files("./")
    stages <- c("annotation_umap", "clustering_umap", "norm_pca", "qc", "qcplot")
    counts <- rep(0, length(stages))
    names(counts) <- stages

    for (file_name in files){
        if (grepl(key, file_name)){
            for (stage in stages){
                if (grepl(stage, file_name)){
                    counts[stage] <- counts[stage] + 1
                }
            }
        }
    }

    if (any(counts > 0)){
        print(names(counts)[which.max(counts)])
    } else {
        print("upload")
    }
}

#* @param link:str
#* @post /user_url_download
#* This api is used to download the file uploaded by the client to the same directory of this R script on the server.
function(link, sys.user_id, req, res){
  sys.user_id <<- sys.user_id
  destfile <<- paste0("./", sys.user_id, ".RDS")
  download.file(link, file.path("rds", destfile))
  destfile
}

#* @post /qcplot
#* @serializer unboxedJSON
#* This api is used to analyze the percentage of mitochondria gene, distribution of number of gene and cells, and plot a violin graph for quality control.
qcplot <- function(req, res){
  seurat_obj <- readRDS(file.path("rds", paste0(sys.user_id, ".RDS")))
  if(is.null(seurat_obj)){
    stop("No Seurat object uploaded", call. = FALSE)
  }
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  file_name <- paste0(sys.user_id, "-qcplot", '.RDS')
  SaveSeuratRds(seurat_obj, file = file.path("rds", file_name))

  d <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "percent.mt"), ncol = 2, pt.size = 0.1)
  graph_name = paste0(sys.user_id,"-vlnplot",".png")
  ggsave(filename = file.path("images", graph_name), plot = d, width = 10, height = 10, dpi = 300)

  list(success = true, message = (content = paste0("http://scrna.m2mda.com/images/", graph_name), type = "image", status = true))
}
# Handler for loading picture, response resource type what
# Define response type
# Human check
# Test ai 
# How to show pic? Render session

#* @get /qc
#* This api is used for quality control, selecting desired range of number of genes and cells, and a maximum threshold for mitochondria genes
qc <- function(min.features, max.features, max.mtpercent){
  seurat_obj <- readRDS(file.path("rds", paste0(sys.user_id, "-qcplot", ".RDS")))
  min.features <- as.numeric(min.features)
  max.features <- as.numeric(max.features)
  max.mtpercent <- as.numeric(max.mtpercent)
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > min.features & nFeature_RNA < max.features & percent.mt < max.mtpercent)
  
  file_name <- paste0(sys.user_id, "-qc", '.RDS')
  SaveSeuratRds(seurat_obj, file = file.path("rds", file_name))
}

#* @serializer unboxedJSON
#* @post /norm_pca
#* This api is used for normalization and pca reduction so that the expression matrix can be normalized and find top variable genes. PCA reduction is then used to reduce multidimensional matrix to a 2D graph.
norm_pca <- function(scaling_factor, num_hvgs, norm_method, hvg_method, res){

  seurat_obj <- readRDS(file.path("rds", paste0(sys.user_id, "-qc", ".RDS")))
  scaling_factor <- as.numeric(scaling_factor)
  num_hvgs <- as.numeric(num_hvgs)

  seurat_obj <- NormalizeData(seurat_obj, normalization.method = norm_method, scale.factor = scaling_factor)
    
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = hvg_method, nfeatures = num_hvgs)
    
  # Scale data
  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))

  # Run PCA
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

  d <- ElbowPlot(object = seurat_obj, ndims = 40)
  graph_name = paste0(sys.user_id, "-elbowplot",".png")
  ggsave(filename = file.path("images", graph_name), plot = d, width = 10, height = 10, dpi = 300)

  file_name <- paste0(sys.user_id,"-norm_pca", '.RDS')
  SaveSeuratRds(seurat_obj, file = file.path("rds", file_name))

  list(type = "image", content = paste0("http://scrna.m2mda.com/images/", graph_name))
}

#* @post /clustering_umap
#* This api is used for clustering (with UMAP visualization) the Seurat object with client's parameter input, dim , which is chosen from the previous elbow plot, and resolution, the extent of how many clusters will be generated (high resolution: more clusters)
clustering_umap <- function(dim, resolution, res){
  
  seurat_obj <- readRDS(file.path("rds", paste0(sys.user_id,"-norm_pca",".RDS")))

  dim <- as.numeric(dim)
  resolution <- as.numeric(resolution)

  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:dim)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:dim)

  file_name <- paste0(sys.user_id, "-clustering_umap", '.RDS')
  SaveSeuratRds(seurat_obj, file = file.path("rds", file_name))
}

#* @post /clustering_tsne
#* This api is used for clustering (with t-SNE visualization) the Seurat object with client's parameter input, dim , which is chosen from the previous elbow plot, and resolution, the extent of how many clusters will be generated (high resolution: more clusters)
clustering_tsne <- function(dims, resolution, res){
  seurat_obj <- readRDS(file.path("rds", paste0(sys.user_id, "-norm_pca", ".RDS")))

  dims <- as.numeric(dims)
  resolution <- as.numeric(resolution)

  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:dims)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
  seurat_obj <- RunTSNE(seurat_obj, dims = 1:dims)

  file_name <- paste0(sys.user_id, "-clustering_tsne", '.RDS')
  SaveSeuratRds(seurat_obj, file = file.path("rds", file_name))
}

#* annotation_sctype_umap
#* @serializer unboxedJSON
#* @post /annotation_sctype_umap
#* This api is used for annotating the clustered object with scType algorithm, calculating the sctype score with positive and negative marker gene sets to classify the cell type. UMAP is used for visualization.
annotation_sctype_umap <- function(tissue, res){

  seurat_obj <- readRDS(file.path("rds", paste0(sys.user_id, "-clustering_umap", ".RDS")))

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

    file_name <- paste0(sys.user_id, "-annotation_umap", '.RDS')
    SaveSeuratRds(seurat_obj, file = file.path("rds", file_name))

    d <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')+ggtitle(tissue)
    graph_name = paste0(sys.user_id,"-annotation_umap",".png")
    ggsave(filename = file.path("images", graph_name), plot = d, width = 10, height = 10, dpi = 300)

    list(type = "image", content = paste0("http://scrna.m2mda.com/images/", graph_name))
}

#* annotation_sctype_tsne
#* @serializer unboxedJSON
#* @post /annotation_sctype_tsne
#* This api is used for annotating the clustered object with scType algorithm, calculating the sctype score with positive and negative marker gene sets to classify the cell type. t-SNE is used for visualization.
annotation_sctype_tsne <- function(tissue, res){

  seurat_obj <- readRDS(file.path("rds", paste0(sys.user_id, "-clustering_tsne", ".RDS")))

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
    file_name <- paste0(sys.user_id, "-annotation_tsne", '.RDS')
    SaveSeuratRds(seurat_obj, file = file.path("rds", file_name))

    b <- TSNEPlot(seurat_obj, group.by = 'customclassif')+ggtitle(tissue)
    graph_name = paste0(sys.user_id,"-annotation_tsne",".png")
    ggsave(filename = file.path("images", graph_name), plot = d, width = 10, height = 10, dpi = 300)

    list(type = "image", content = paste0("http://scrna.m2mda.com/images/", graph_name))
}

# download
