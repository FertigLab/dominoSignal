# use --platform=linux/amd64 to avoid 'no match for platform in the manifest' on M1
FROM rocker/tidyverse:4

COPY . /dominoSignal
WORKDIR /dominoSignal

RUN sudo apt-get update -y && \
    apt-get upgrade -y && \
    apt-get install libhdf5-dev build-essential patch libglpk-dev libxml2-dev gfortran -y

#packages below didn't install with devtools::install_deps, needed BiocManager;
RUN Rscript -e 'BiocManager::install(c("biomaRt", "ComplexHeatmap", "S4Arrays", "SingleCellExperiment", "SummarizedExperiment"),ask=FALSE)'

#install all other dependencies
RUN Rscript -e 'devtools::install_deps(".", dependencies=TRUE)'

#need to restart R sometimes https://github.com/r-lib/devtools/issues/2395
RUN Rscript -e 'devtools::install(".", dependencies=TRUE)'

