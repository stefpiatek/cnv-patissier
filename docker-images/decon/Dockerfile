FROM r-base:3.1.2

WORKDIR /usr/scripts/

RUN R -e "install.packages(c('packrat'))" \
    && wget --no-check-certificate -q https://zenodo.org/record/1254945/files/RahmanTeam/DECoN-v1.0.2.zip \
    && unzip DECoN-v1.0.2.zip \
    && cd RahmanTeam-DECoN-8d740c5/Linux \
    && Rscript sessionInfo.R --bootstrap-packrat \
    && cp packrat/packrat_source/.Rprofile .

WORKDIR /usr/scripts/RahmanTeam-DECoN-8d740c5/Linux

CMD ["bash"]
