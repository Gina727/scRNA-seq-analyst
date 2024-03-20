FROM rstudio/plumber

RUN R -e "install.packages('plumber')"
RUN R -e "install.packages('SingleCellExperiment')"
RUN R -e "install.packages('Seurat')"
RUN R -e "install.packages('HGNChelper')"
RUN apt-get update

COPY / /

EXPOSE 80

ENTRYPOINT ["Rscript", "app-api.R"]