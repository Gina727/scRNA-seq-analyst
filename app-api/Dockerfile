FROM rstudio/plumber

# 安装依赖的系统库
RUN apt-get update && \
    apt-get install -y libglpk40 libcurl4-openssl-dev libssl-dev libxml2-dev && \
    rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('devtools', 'remotes'))"

# 使用 devtools 或 remotes 安装特定版本的 R 包
RUN R -e "remotes::install_version('SeuratObject', version = '5.0.1')"

# 继续安装其他包
RUN R -e "install.packages(c('SingleCellExperiment', 'Seurat', 'HGNChelper','dplyr','ggplot2','uuid'), dependencies=TRUE)"
RUN R -e "install.packages(c('openxlsx'))"
# 复制你的应用文件到镜像中
COPY / / 

# 暴露端口
EXPOSE 8000

CMD ["new_api.R"]

