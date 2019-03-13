FROM java:8-jdk

ENV CLASSPATH=/usr/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar:/usr/savvycnv

# install GATK
RUN wget -O /usr/gatk.tar.bz2 "https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.8-1-0-gf15c1c3ef" \
    && tar -xjf /usr/gatk.tar.bz2 -C /usr/ 

# install savvy cnv
RUN apt-get update \
    && apt-get install -y git \
    && rm -rf /var/lib/apt/lists/* \
    && git clone http://github.com/rdemolgen/SavvySuite.git /usr/savvycnv \
    && cd /usr/savvycnv \
    && git checkout f996a83 \
    && javac *.java



