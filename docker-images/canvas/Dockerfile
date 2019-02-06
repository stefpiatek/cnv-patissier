FROM mono:5.14

ARG DEBIAN_FRONTEND=noninteractive

RUN echo "deb http://ftp.de.debian.org/debian testing main" >> /etc/apt/sources.list \
    && apt-get update \
    && apt-get install -y --no-install-recommends wget unzip \
    && apt-get -t testing install -y --no-install-recommends python3.7 \
    && rm -rf /var/lib/apt/lists/*

RUN wget --no-check-certificate https://github.com/Illumina/canvas/releases/download/v1.11.0/Canvas-1.11.0_x64.zip \
    && unzip Canvas-1.11.0_x64.zip -d /usr/canvas/ \
    && rm Canvas-1.11.0_x64.zip

WORKDIR /usr/canvas

