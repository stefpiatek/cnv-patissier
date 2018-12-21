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

COPY . /usr/

RUN R CMD SHLIB /usr/EXCAVATOR2_Package_v1.1.2/lib/F77/F4R.f; \
    R CMD SHLIB /usr/EXCAVATOR2_Package_v1.1.2/lib/F77/FastJointSLMLibraryI.f;

ENV PATH $PATH:/usr/EXCAVATOR2_Package_v1.1.2

WORKDIR /usr/EXCAVATOR2_Package_v1.1.2

CMD ["bash"]
