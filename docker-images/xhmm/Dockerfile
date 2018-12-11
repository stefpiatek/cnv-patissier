FROM broadinstitute/gatk3:3.8-1

COPY xhmm /usr/xhmm
COPY plinkseq-0.10 /usr/plinkseq-0.10
RUN make -C /usr/xhmm
ENV PATH $PATH:/usr/xhmm:/usr/plinkseq-0.10
