FROM r-base:3.5.1

RUN apt-get update \ 
    && apt-get install -y --no-install-recommends libcurl4-openssl-dev libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

RUN R -e 'install.packages(c("optparse", "ExomeDepth"))'

WORKDIR /mnt/cnv-caller-resources/exome-depth/
