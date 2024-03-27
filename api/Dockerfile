FROM rstudio/plumber

# 安装依赖的系统库
RUN apt-get update && \
    apt-get install -y libglpk40 libcurl4-openssl-dev libssl-dev libxml2-dev && \
    rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('devtools', 'remotes'))"
RUN R -e "remotes::install_version('SeuratObject', version = '5.0.1')"
RUN R -e "install.packages(c('SingleCellExperiment', 'Seurat', 'HGNChelper'), dependencies=TRUE)"

COPY / / 

EXPOSE 8000

CMD ["app-api.R"]

