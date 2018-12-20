FROM mono:5.14
COPY canvas-1.11.0 /usr/canvas/

RUN echo "deb http://ftp.de.debian.org/debian testing main" >> /etc/apt/sources.list \
    && apt-get update \
    && apt-get -t testing install -y --no-install-recommends python3.6 \
    && rm -rf /var/lib/apt/lists/*


WORKDIR /usr/canvas

