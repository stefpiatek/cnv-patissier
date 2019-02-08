FROM r-base:3.5.1

RUN apt-get update \
    && apt-get install -y --no-install-recommends libcurl4-openssl-dev libxml2-dev libssl-dev \
    && rm -rf /var/lib/apt/lists/*

RUN R -e 'install.packages(c("BiocManager", "devtools", "optparse")) ; \
    BiocManager::install("GenomicAlignments", version = "3.8") ; \
    devtools::install_version("ExomeDepth", version = "1.1.10", repos = "http://cran.ma.imperial.ac.uk")'

WORKDIR /mnt/cnv-caller-resources/exome-depth/
