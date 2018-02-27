#!/bin/bash

#cd /home/ubuntu/rawdata
date
echo "start rawdata..."
cd /RNA_mapping/build/rawdata
./doGetRawdata.sh

#cd /home/ubuntu/genome
date
echo "Getting genome..."
cd /RNA_mapping/build/genome
./doGetChr20.sh

#cd /home/ubuntu/bwa
date
echo "start bwa..."
cd /RNA_mapping/build/bwa
./doBWA.sh

#cd /home/ubuntu/hisat
date
echo "start hisat2..."
cd /RNA_mapping/build/hisat
./doHISAT.sh

#cd /home/ubuntu/star
date
echo "start star..."
cd /RNA_mapping/build/star
./doSTAR.sh

#cd /home/ubuntu/blastmapper
#./doBlastmapper.sh
date
echo "start magicblast..."
cd /RNA_mapping/build/magicblast
./doMagicblast.sh

#cd /home/ubuntu/csaw
date
echo "start csaw..."
cd /RNA_mapping/build/csaw
Rscript runCsaw.R ../magicblast/magicblast.sorted.bam ../hisat/hisat.sorted.bam ../star/star.sorted.bam ../bwa/bwa.sorted.bam

#cd /home/ubuntu/bamdiff
date
echo "start bamdiff..."
cd /RNA_mapping/build/bamdiff
bamDiff.py ../csaw/csaw.bam.results.csv ../hisat/hisat.sorted.bam ../star/star.sorted.bam
date
echo "Done!"
