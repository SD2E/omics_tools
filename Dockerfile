FROM sd2e/python3:ubuntu18
#FROM r-base:3.6.1
#RUN sed -i -re 's/([a-z]{2}\.)?archive.ubuntu.com|security.ubuntu.com/old-releases.ubuntu.com/g' /etc/apt/sources.list \
# Install R
RUN export DEBIAN_FRONTEND=noninteractive \
     && apt-get update \
     && apt-get install -y apt-transport-https software-properties-common \
     && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
     && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' \
     && apt-get update \
     && apt-get -y install libssl-dev \
                           libcurl4-openssl-dev \
                           libc6 \
                           r-base-core \
                           r-base \
                           r-base-dev \
                           libxml2-dev \
                           wget \
			                     rsync



# Install nsmblR and other R packages
RUN Rscript -e 'install.packages("RcppEigen")'
RUN Rscript -e 'install.packages("BiocManager"); BiocManager::install("edgeR")'
RUN Rscript -e 'BiocManager::install("clusterProfiler")'

ADD omics_tools /install/omics_tools
RUN cd /install/omics_tools && \
          mkdir data && cd data && \
          wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz

#/usr/local/lib/python3.6/dist-packages/omics_tools/data/Bacteria.gene_info
COPY README.md setup.py requirements.txt /install/
RUN cd /install \
  && pip3 install .


    # brittle to the location of they python install, there's definitely
    # a better way to do this, should eventually by a paramater or seperate app
RUN rsync -rltv /install/omics_tools/data /usr/local/lib/python3.6/dist-packages/omics_tools/ && \
    cd /usr/local/lib/python3.6/dist-packages/omics_tools/ && \
    python3 -c 'from annotate_GO_KEGG import check_data_files; check_data_files()' && \
    gunzip /usr/local/lib/python3.6/dist-packages/omics_tools/data/gene2go.gz

RUN pip3 install fisher
# Add test data
ADD tests /tests
#ADD src /src
#RUN chmod +x /src/*
