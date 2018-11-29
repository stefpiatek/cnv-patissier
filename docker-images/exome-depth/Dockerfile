FROM r-base:3.5.1

COPY . /usr/scripts/
WORKDIR /usr/scripts/

RUN apt-get update \ 
    && apt-get install -y --no-install-recommends libcurl4-openssl-dev libxml2-dev

RUN R -e 'install.packages(c("bitops_1.0-6.tar.gz", "RCurl_1.95-4.11.tar.gz", "matrixStats_0.54.0.tar.gz", "BH_1.66.0-1.tar.gz", "formatR_1.5.tar.gz", "lambda.r_1.2.3.tar.gz", "getopt_1.20.2.tar.gz", "optparse_1.6.0.tar.gz", "futile.options_1.0.1.tar.gz", "futile.logger_1.4.3.tar.gz", "zlibbioc_1.28.0.tar.gz", "aod_1.3.tar.gz", "snow_0.4-3.tar.gz", "GenomeInfoDbData_1.2.0.tar.gz", "BiocGenerics_0.28.0.tar.gz", "Biobase_2.42.0.tar.gz", "S4Vectors_0.20.0.tar.gz", "IRanges_2.16.0.tar.gz", "BiocParallel_1.16.0.tar.gz", "XVector_0.22.0.tar.gz", "Biostrings_2.50.0.tar.gz", "VGAM_1.0-6.tar.gz", "GenomeInfoDb_1.18.0.tar.gz", "GenomicRanges_1.34.0.tar.gz", "DelayedArray_0.8.0.tar.gz", "SummarizedExperiment_1.12.0.tar.gz", "Rsamtools_1.34.0.tar.gz",  "GenomicAlignments_1.18.0.tar.gz", "ExomeDepth_1.1.10.tar.gz"), repos = NULL, type = "source")'

WORKDIR /mnt/cnv-caller-resources/exome-depth/
