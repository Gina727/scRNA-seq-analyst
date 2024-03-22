FROM rstudio/plumber

RUN R -e "install.packages('plumber')"
RUN R -e "install.packages('SingleCellExperiment')"
RUN R -e "install.packages('Seurat')"
RUN R -e "install.packages('SeuratObject')"
RUN R -e "install.packages('HGNChelper')"
RUN apt-get update
RUN apt-get install -y libglpk40

CMD Rscript /app-api/requirement.R

COPY / /

EXPOSE 80

ENTRYPOINT ["Rscript", "app-api.R"]
