#!/bin/bash

#cd /home/ubuntu/rawdata
date
echo "start rawdata..."
rawdata/doGetRawdata.sh

#cd /home/ubuntu/genome
date
echo "Getting genome..."
genome/doGetChr20.sh

#cd /home/ubuntu/bwa
date
echo "start bwa..."
bwa/doBWA.sh

#cd /home/ubuntu/hisat
date
echo "start hisat2..."
hisat/doHISAT.sh

#cd /home/ubuntu/star
date
echo "start star..."
star/doSTAR.sh

#cd /home/ubuntu/blastmapper
#./doBlastmapper.sh
date
echo "start magicblast..."
magicblast/doMagicblast.sh

#cd /home/ubuntu/csaw
date
echo "start csaw..."
Rscript csaw/runCsaw.R ../magicblast/magicblast.sorted.bam ../hisat/hisat.sorted.bam ../star/star.sorted.bam ../bwa/bwa.sorted.bam

#cd /home/ubuntu/bamdiff
date
echo "start bamdiff..."
bamdiff/bamDiff.py ../csaw/csaw.bam.results.csv ../hisat/hisat.sorted.bam ../star/star.sorted.bam
date
echo "Done!"
