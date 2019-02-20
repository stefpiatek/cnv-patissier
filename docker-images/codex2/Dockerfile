FROM r-base:3.5.1

RUN apt-get update \
    && apt-get install -y --no-install-recommends libcurl4-openssl-dev libxml2-dev libssl-dev \
    && rm -rf /var/lib/apt/lists/*


RUN R -e "install.packages(c('optparse', 'BiocManager', 'devtools')); \
    BiocManager::install('BSgenome', version = '3.8'); \
    BiocManager::install('BSgenome.Hsapiens.1000genomes.hs37d5', version = '3.8'); \
    BiocManager::install('BSgenome.Hsapiens.UCSC.hg19', version = '3.8'); \
    BiocManager::install('CODEX', version = '3.8'); \
    devtools::install_github('yuchaojiang/CODEX2/package', ref='26e796cf222d1104a1274ee7a3f8c4af7bf37c62');" 


CMD ["bash"]
