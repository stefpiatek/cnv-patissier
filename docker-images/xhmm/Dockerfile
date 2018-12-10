FROM broadinstitute/gatk3:3.8-1

COPY xhmm /usr/xhmm
RUN make -C /usr/xhmm
ENV PATH $PATH:/usr/xhmm