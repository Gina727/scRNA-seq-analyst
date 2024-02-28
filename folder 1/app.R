options(shiny.maxRequestSize=30*1024^2)
library(rsconnect)
library(shiny)
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(writexl)
library(HGNChelper)

###
normmethodlist <- c("LogNormalize")
hvglist <- c("vst", "mvp", "disp")
tissuelist <- c("Immune system")
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
###

ui <- navbarPage("scRNA-seq Analyst",
               tabPanel(
                 "Lobby",
                 sidebarLayout(
                     navlistPanel(
                       id = "panels",
                       "Input and Preprocessing",
                       tabPanel("1. Upload",
                                value = "1",
                                p("Please upload .RDS file. To convert other files into .RDS, see 'How to use page'."),
                                fileInput("file", h3("File input")),
                                actionButton("to2", "Process file")
                                ),
                       tabPanel("2. Quality Control",
                                value = "2",
                                plotOutput("qc_plot"),
                                sliderInput("nFeature_RNA", "Num of Genes",
                                                     min = 0, max = 10000, value = c(200, 5000)),
                                sliderInput("percent.mt", "Max. Percentage of mito-genes",
                                                     min = 0, max = 100, value = 5),
                                actionButton("to3", "Filter")
                                ),
                       tabPanel("3. Normalization and Reduce Dimension",
                                value = "3",
                                p("Choose scaling factor, normalizartion method and desired highly variable gene number."),
                                  selectInput("normmethod", "Normalization method", normmethodlist),
                                  br(),
                                  numericInput("scalingfactor", "Scaling factor", value = 10000),
                                  br(),
                                  selectInput("hvgmethod", "Selection method", hvglist),
                                  mainPanel(
                                    actionButton("explainhvg", "Explain methods"),
                                    verbatimTextOutput("explainhvg"),
                                  ),
                                  br(),
                                  numericInput("HVGs", "HVGs", value = 2000),
                                  actionButton("to4", "Scale and Run PCA"),
                                ),
                       "Clustering and Annotation",
                       tabPanel("4. Clustering",
                                value = "4",
                                plotOutput("elbow_plot"),
            
                                numericInput("dim", "Number of PCA", value = 10),
                                br(),
                                numericInput("resolution", "Resoultion", value = 0.8),
                                br(),
                                actionButton("umap", "Run UMAP reduction"),
                                actionButton("tsne", "Run t-SNE reduction"),
                                ),
                       tabPanel("5. Annotation",
                                value = "5",
                                column(4, plotOutput("umapcluster")),
                                column(4, plotOutput("tsnecluster")),
                                selectInput("tissue","Tissue name", tissuelist),
                                actionButton("sctype_umap", "scType annotation using UMAP reduction"),
                                actionButton("sctype_tsne", "scType annotation using t-SNE reduction")
                                ),
                       "Output",
                       tabPanel("Output",
                                value = "6",
                                plotOutput("output"))
                   ),
                   mainPanel(),
                 ),
                 ),
               tabPanel("How to use",
                        titlePanel("Single Cell RNA sequence analyst"),
                        mainPanel(
                          h2("How to use?"),
                          p("Upload your .rds file or choose a given dataset to try."),
                          p("Before you start, create a seurat object of your data and export them as rds file."),
                          p("You may follow the example below to get your rds file."),
                          p("1. Open a new R script."),
                          p("2. Put all 10x files (barcodes.tsv, genes.tsv and matrix.mtx) into a folder inside the same folder as your R script."),
                          p("3. Install and load packages SingleCellExperiment and Seurat using"),
                          code("install.packages () and library( )"),
                          p("4. Run the following code"),
                          code("counts.data <- Read10X(data.dir = 'path_to_the_file_saving_10x_files')"),
                          br(),
                          code("counts <- CreateSeuratObject(counts = counts.data,  min.cells = 3, min.features = 200)"),
                          br(),
                          code("SaveSeuratRds(counts, file='counts.RDS')"),
                          
                          h2("Workflow"),
                          p("This is the structure of scRNA-seq analysing."),
                          img(src = "scRNA-seq.png", height = 600, width = 500)
                        ),),
  )
server <- shinyServer(function(input, output, session) {
  seurat_obj_react <- reactiveVal()
  
  # section 1 upload
  observeEvent(input$to2, {
    # Read the file as a data frame
    seurat_obj <- LoadSeuratRds(input$file$datapath)
    # Store the seurat object in the reactive value
    # Perform some QC and normalization steps
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
    
    seurat_obj_react(seurat_obj)
    
    # Next panel
    updateNavlistPanel(
      session,
      inputId = "panels",
      selected = paste0("2")
    )
    })
  
  # Render the QC plot
  output$qc_plot <- renderPlot({
    # Get the seurat object from the reactive value
    qc <- seurat_obj_react()
    # Plot the number of features and percent mitochondrial content per cell
    VlnPlot(qc, features = c("nFeature_RNA", "percent.mt"), ncol = 2, pt.size = 0.1)
  })
  
  # section 2 qc
  observeEvent(input$to3, {
    seurat_obj <- seurat_obj_react()
    seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > min(input$nFeature_RNA) & nFeature_RNA < max(input$nFeature_RNA) & percent.mt < max(input$percent.mt))
    seurat_obj_react(seurat_obj)
    
    # Next panel
    updateNavlistPanel(
      session,
      inputId = "panels",
      selected = paste0("3")
    )
  })
  
  # section 3 norm & HVG
  observeEvent(input$explainhvg, {
    output$explainhvg <- renderText({
      paste("vst: First, fits a line to the relationship of log(variance) and log(mean) using local polynomial regression (loess). Then standardizes the feature values using the observed mean and expected variance (given by the fitted line). Feature variance is then calculated on the standardized values after clipping to a maximum (see clip.max parameter)." , "mean.var.plot (mvp): First, uses a function to calculate average expression (mean.function) and dispersion (dispersion.function) for each feature. Next, divides features into num.bin (deafult 20) bins based on their average expression, and calculates z-scores for dispersion within each bin. The purpose of this is to identify variable features while controlling for the strong relationship between variability and average expression.", "dispersion (disp): selects the genes with the highest dispersion values", sep="\n")
    })
  })

  observeEvent(input$to4, {
    seurat_obj <- seurat_obj_react()
    
    # normalization with logarithm
    seurat_obj <- NormalizeData(seurat_obj, normalization.method = input$normmethod, scale.factor = input$scalingfactor)
    
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = input$hvgmethod, nfeatures = input$HVGs)
    
    # Scale data
    seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
    
    # Run PCA
    seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
    
    seurat_obj_react(seurat_obj)
    
    # Plot Elbow plot to choose PCA
    output$elbow_plot <- renderPlot({
      elbow <- seurat_obj_react()
      # Plot the number of features and percent mitochondrial content per cell
      ElbowPlot(object = elbow, ndims = 40)
    })
    
    # next panel
    updateNavlistPanel(
      session,
      inputId = "panels",
      selected = paste0("4")
    )
  })
  
  observeEvent(input$umap, {
    seurat_obj <- seurat_obj_react()
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:input$dim) #dims = dimensions
    seurat_obj <- FindClusters(seurat_obj, resolution = input$resolution)
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:input$dim)
    seurat_obj_react(seurat_obj)
    
    output$umapcluster <- renderPlot({
      DimPlot(seurat_obj, reduction = "umap")
    })
    
    updateTabsetPanel(
      session,
      inputId = "panels",
      selected = paste0("5")
    )
  })
  
  observeEvent(input$tsne, {
    seurat_obj <- seurat_obj_react()
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:input$dim) #dims = dimensions
    seurat_obj <- FindClusters(seurat_obj, resolution = input$resolution)
    seurat_obj <- RunTSNE(seurat_obj, dims = 1:input$dim)
    seurat_obj_react(seurat_obj)
    
    output$tsnecluster <- renderPlot({
      TSNEPlot(seurat_obj)
    })
    updateTabsetPanel(
      session,
      inputId = "panels",
      selected = paste0("5")
    )
  })
  
  observeEvent(input$sctype_umap, {
    seurat_obj <- seurat_obj_react()
    tissue = input$tissue
    
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
    p <- DimPlot(counts, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif',) 
    p <- p + ggtitle(input$tissue)
    
    output$sctype <- renderPlot({
      print(p)
    })
    seurat_obj_react(seurat_obj)
    
    updateTabsetPanel(
      session,
      inputId = "panels",
      selected = paste0("6")
    )
  })
  })

shinyApp(ui = ui, server = server)
