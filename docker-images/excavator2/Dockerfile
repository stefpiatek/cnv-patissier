FROM r-base:3.5.1

RUN R -e "install.packages('Hmisc', repos = 'http://cran.uk.r-project.org')"

RUN apt-get update \
    && apt-get install -y --no-install-recommends  bedtools\
    && rm -rf /var/lib/apt/lists/*

# Fix to get libpng12.so.0 installed
RUN wget http://fr.archive.ubuntu.com/ubuntu/pool/main/libp/libpng/libpng12-0_1.2.54-1ubuntu1_amd64.deb \
    && dpkg -x libpng12-0_1.2.54-1ubuntu1_amd64.deb /root/libpng12 \
    && ln -s /root/libpng12/lib/x86_64-linux-gnu/libpng12.so.0 /usr/lib/libpng12.so.0 \
    && ln -s /root/libpng12/lib/x86_64-linux-gnu/libpng12.so.0.54.0 /usr/lib/libpng12.so.0.54.0 \
    && ldconfig

RUN wget -q https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 \
    && tar xjvf samtools-1.9.tar.bz2 \
    && cd /samtools-1.9/ \
    && ./configure && make \
    && mv /samtools-1.9/samtools /bin/ \
    && rm -rf /samtools-1.9 

RUN cd /usr \
    && wget --no-check-certificate https://sourceforge.net/projects/excavator2tool/files/EXCAVATOR2_Package_v1.1.2.tgz \
    && tar zxvf EXCAVATOR2_Package_v1.1.2.tgz \
    && R CMD SHLIB /usr/EXCAVATOR2_Package_v1.1.2/lib/F77/F4R.f \
    && R CMD SHLIB /usr/EXCAVATOR2_Package_v1.1.2/lib/F77/FastJointSLMLibraryI.f \
    && wget --no-check-certificate -O /usr/EXCAVATOR2_Package_v1.1.2/data/ucsc.hg19.bw "https://www.dropbox.com/s/djltemd4pb1051x/ucsc.hg19.bw?dl=1" \
    && wget  --no-check-certificate -O /usr/EXCAVATOR2_Package_v1.1.2/data/GCA_000001405.15_GRCh38.bw "https://www.dropbox.com/s/soshiqrayggu2bb/GCA_000001405.15_GRCh38.bw?dl=1"

ENV PATH $PATH:/usr/EXCAVATOR2_Package_v1.1.2

WORKDIR /usr/EXCAVATOR2_Package_v1.1.2

CMD ["bash"]
