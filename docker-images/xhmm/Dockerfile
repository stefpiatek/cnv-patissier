FROM broadinstitute/gatk3:3.8-1

RUN wget --no-check-certificate https://psychgen.u.hpc.mssm.edu/plinkseq_downloads/plinkseq-x86_64-0.10.zip \
    && unzip plinkseq-x86_64-0.10.zip -d /usr/plinkseq-0.10 \
    && rm plinkseq-x86_64-0.10.zip \
    && wget --no-check-certificate https://bitbucket.org/statgen/xhmm/get/cc14e528d909.zip \
    && unzip cc14e528d909.zip -d /usr/xhmm_temp \
    && mv /usr/xhmm_temp/statgen-xhmm-cc14e528d909 /usr/xhmm \
    && rm cc14e528d909.zip

RUN make -C /usr/xhmm
ENV PATH $PATH:/usr/xhmm:/usr/plinkseq-0.10
