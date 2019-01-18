# later versions have errors in dependency installation for CopywriteR
FROM r-base:3.2.5  

RUN R -e "install.packages(c('optparse', 'dplyr', 'tidyr', 'stringr')); source('http://bioconductor.org/biocLite.R'); biocLite('CopywriteR');"

CMD ["bash"]
