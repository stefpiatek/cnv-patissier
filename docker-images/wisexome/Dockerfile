FROM continuumio/miniconda2:4.5.4

WORKDIR /usr/

RUN apt-get update && apt-get install -y gcc \
    && pip install --upgrade pip \
    && pip install \
        matplotlib==2.0.1 \
        numpy==1.13.0 \
        pysam==0.10.0 \
        scipy==0.19.0 \
        xlwt==1.3.0 \
    && git clone https://github.com/VUmcCGP/wisexome.git

WORKDIR /usr/wisexome
COPY . /usr/wisexome/reference

CMD ["bash"]
