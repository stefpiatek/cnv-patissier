FROM ubuntu:18.04

RUN apt-get update && apt-get install -y \
    bzip2 \
    g++ \
    libbz2-dev \
    liblzma-dev \
    make \
    ncurses-dev \
    wget \
    zlib1g-dev

WORKDIR /tmp

RUN wget -q https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 \
    && tar xjvf samtools-1.9.tar.bz2 \
    && cd /tmp/samtools-1.9/ \
    && ./configure && make \
    && mv /tmp/samtools-1.9/samtools /bin/ \
    && rm -rf /tmp//samtools-1.9

CMD ["bash"]
