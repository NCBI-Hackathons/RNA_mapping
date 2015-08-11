#!/bin/bash

cd /home/ubuntu/rawdata
./doGetRawdata.sh

cd /home/ubuntu/genome
./doGetChr20.sh

cd /home/ubuntu/bwa
./doBWA.sh

cd /home/ubuntu/hisat
./doHISAT.sh

cd /home/ubuntu/star
./doSTAR.sh

cd /home/ubuntu/blastmapper
./doBlastmapper.sh

cd /home/ubuntu/csaw
Rscript runCsaw.R ../blastmapper/blastmapper.sorted.bam ../hisat/hisat.sorted.bam ../star/star.sorted.bam ../bwa/bwa.sorted.bam

cd /home/ubuntu/bamdiff
./bamDiff.py ../csaw/csaw.bam.results.csv ../hisat/hisat.sorted.bam ../star/star.sorted.bam > hisatVSstar
