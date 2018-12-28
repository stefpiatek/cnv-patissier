# later versions have errors in dependency installation for CopywriteR
FROM r-base:3.2.5  

RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('CopywriteR'); install.packages('optparse')"

CMD ["bash"]
