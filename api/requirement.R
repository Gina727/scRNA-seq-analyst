options(repos = c(CRAN = "https://cran.rstudio.com/"))

# 创建一个包含所有需要安装的包的向量
packages <- c("rsconnect", "shiny", "SingleCellExperiment", "Seurat",
              "tidyverse", "Matrix", "scales", "cowplot", "RCurl",
              "writexl", "HGNChelper", "openxlsx", "plumber")

# 使用lapply()遍历包列表，安装尚未安装的包
lapply(packages, function(pkg) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        install.packages(pkg)
    }
})

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("SingleCellExperiment", "Seurat"))