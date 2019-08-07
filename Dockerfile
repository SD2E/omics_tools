#FROM sd2e/python3:ubuntu18
FROM r-base:3.6.1
# ENV PATH "/opt/bin/:$PATH"
# ADD config.yml /config.yml
# ADD src /opt/src

#RUN sed -i -re 's/([a-z]{2}\.)?archive.ubuntu.com|security.ubuntu.com/old-releases.ubuntu.com/g' /etc/apt/sources.list
# Install R
# RUN export DEBIAN_FRONTEND=noninteractive \
#     && apt-get update \
#     && apt-get install -y apt-transport-https software-properties-common \
#     && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
#     && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' \
#     && apt-get update \
#     && apt-get -y install libssl-dev \
#                           libcurl4-openssl-dev \
#                           libc6 \
#                           r-base-core \
#                           r-base \
#                           r-base-dev \
#                           libxml2-dev \
#                           wget



# Install nsmblR and other R packages
#RUN Rscript -e 'install.packages("devtools", type="source", repos = "http://cran.us.r-project.org")'
#RUN Rscript -e 'install.packages("devtools")'
#RUN Rscript -e 'library("devtools"); install_github("RcppCore/RcppEigen")'
RUN Rscript -e 'install.packages("RcppEigen")'
RUN Rscript -e 'install.packages("BiocManager"); BiocManager::install("edgeR")'
RUN Rscript -e 'BiocManager::install("clusterProfiler")'
#RUN Rscript -e 'Sys.setenv(R_LIBS_USER="/usr/local/lib/R/site-library"); BiocManager::install("clusterProfiler", dependencies=TRUE, type = "source")'
RUN apt -y autoremove && \
    apt -y clean
#ERROR: dependency 'enrichplot' is not available for package
#RUN Rscript -e 'BiocManager::install("clusterProfiler", dependencies=TRUE, type = "source", checkBuilt = TRUE)
#RUN Rscript -e 'install.packages(c("igraph", "optparse"))'
#RUN Rscript -e 'install.packages(c("corrr","dplyr", "ggplot2","magrittr","Matrix","readr","tibble","votesys"))'
#RUN Rscript -e 'library(devtools); devtools::install_github("joshuaurrutia/nsmblR")'

# Installing package(s) 'clusterProfiler'
#
# R[write to console]: also installing the dependencies 'sys', 'backports', 'ellipsis', 'zeallot', 'formatR', 'askpass', 'farver', 'bit', 'prettyunits', 'vctrs', 'lambda.r', 'futile.options', 'curl', 'mime', 'openssl', 'hms', 'triebeard', 'tweenr', 'polyclip', 'RcppEigen', 'colorspace', 'utf8', 'bit64', 'blob', 'memoise', 'pkgconfig', 'BH', 'plogr', 'futile.logger', 'snow', 'data.table', 'fastmatch', 'stringr', 'httr', 'jsonlite', 'progress', 'urltools', 'xml2', 'gridGraphics', 'ggforce', 'ggrepel', 'viridis', 'labeling', 'munsell', 'R6', 'cli', 'crayon', 'fansi', 'pillar', 'assertthat', 'BiocGenerics', 'Biobase', 'IRanges', 'DBI', 'RSQLite', 'S4Vectors', 'BiocParallel', 'DO.db', 'fgsea', 'reshape2', 'cowplot', 'europepmc', 'ggplotify', 'ggraph', 'ggridges', 'gridExtra', 'igraph', 'purrr', 'RColorBrewer', 'UpSetR', 'digest', 'gtable', 'lazyeval', 'rlang', 'scales', 'tibble', 'viridisLite', 'withr', 'dplyr', 'glue', 'stringi', 'tidyselect', 'AnnotationDbi', 'DOSE', 'enrichplot', 'ggplot2', 'GO.db', 'GOSemSim', 'magrittr', 'plyr', 'qvalue', 'rvcheck', 'tidyr'
# RUN apt-get -y install wget && \
#     apt autoremove && \
#     apt clean

ADD omics_tools /install/omics_tools
RUN cd /install/omics_tools && \
          mkdir data && cd data && \
          wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz && \
          cd ../ && \
          python3 -c 'from annotate_GO_KEGG import check_data_files; check_data_files()'
#/usr/local/lib/python3.6/dist-packages/omics_tools/data/Bacteria.gene_info
COPY README.md setup.py requirements.txt /install/
RUN cd /install \
  && pip3 install .

# Add test data
ADD tests /tests
#ADD src /src
#RUN chmod +x /src/*
