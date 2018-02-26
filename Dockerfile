## bamdiff blastmapper bwa hisat hisat-build STAR samtools 
## makeblastdb samtools magicblast fastq-dump

FROM ubuntu:latest
RUN rm /bin/sh && ln -s /bin/bash /bin/sh
MAINTAINER Steve Tsang <mylagimail2004@yahoo.com>

RUN apt-get update && apt-get install --yes \
 build-essential \
 gcc-multilib \
 apt-utils \
 zlib1g-dev \
 vim-common \
 wget \
 libncurses5-dev \
 autotools-dev \
 autoconf \
 git \
 perl \
 libbz2-dev \
 liblzma-dev \
 apt-utils \
 libz-dev \
 ncurses-dev \
 zlib1g-dev \
 libcurl3 

WORKDIR /opt
RUN git clone https://github.com/samtools/htslib.git
WORKDIR /opt/htslib
RUN autoheader
RUN autoconf
RUN ./configure
RUN make
RUN make install
ENV PATH "$PATH:/opt/htslib/"

WORKDIR /opt
RUN git clone https://github.com/samtools/samtools.git
WORKDIR /opt/samtools
RUN autoheader
RUN autoconf -Wno-syntax
RUN ./configure    # Optional, needed for choosing optional functionality
RUN make
RUN make install
ENV PATH "$PATH:/opt/samtools/"

#WORKDIR /opt
#RUN git clone git://github.com/samtools/bcftools.git
#WORKDIR /opt/bcftools
#RUN autoheader
#RUN autoconf
#RUN ./configure
#RUN make
#RUN make install
#ENV PATH "$PATH:/opt/bcftools/"

WORKDIR /opt
RUN git clone https://github.com/lh3/bwa.git
WORKDIR /opt/bwa
RUN make
#RUN cp bwa /usr/local/bin
ENV PATH "$PATH:/opt/bwa/"
#./bwa index ref.fa
#./bwa mem ref.fa read-se.fq.gz | gzip -3 > aln-se.sam.gz
#./bwa mem ref.fa read1.fq read2.fq | gzip -3 > aln-pe.sam.gz

RUN apt-get install -y r-base python

#WORKDIR /opt/
#RUN git clone git://github.com/statgen/bamUtil.git
#WORKDIR /opt/bamUtil
#RUN make cloneLib
#RUN make all

WORKDIR /opt/
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.7.1+-x64-linux.tar.gz
RUN tar xvzf ncbi-blast-2.7.1+-x64-linux.tar.gz
WORKDIR /opt/ncbi-blast-2.7.1+
ENV PATH "$PATH:/opt/ncbi-blast-2.7.1+/"

WORKDIR /opt/
RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
RUN tar xvzf sratoolkit.current-ubuntu64.tar.gz
WORKDIR /opt/sratoolkit.2.9.0-ubuntu64
ENV PATH "$PATH: /opt/sratoolkit.2.9.0-ubuntu64/bin/"
RUN apt-get install -y unzip

WORKDIR /opt/
RUN wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip
RUN unzip hisat2-2.1.0-Linux_x86_64.zip
WORKDIR /opt/hisat2-2.1.0
ENV PATH "$PATH:/opt/hisat2-2.1.0/"

WORKDIR /opt/
RUN git clone https://github.com/alexdobin/STAR.git
WORKDIR /opt/STAR/source
RUN make STAR
ENV PATH "$PATH:/opt/STAR/source/"

WORKDIR /
RUN git clone https://github.com/NCBI-Hackathons/RNA_mapping.git

