FROM rstudio/plumber

RUN apt-get update && \
    apt-get install -y libglpk40 libcurl4-openssl-dev libssl-dev libxml2-dev && \
    rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('devtools', 'remotes'))"

RUN R -e "remotes::install_version('SeuratObject', version = '5.0.1')"

RUN R -e "install.packages(c('SingleCellExperiment', 'Seurat', 'HGNChelper','dplyr','ggplot2'), dependencies=TRUE)"
RUN R -e "install.packages(c('openxlsx','SeuratDisk','SeuratData'))"
COPY / / 

EXPOSE 8000

CMD ["new_api.R"]

