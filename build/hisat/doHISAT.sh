#!/bin/bash

cd /home/ubuntu/genome
hisat-build chr20.fa chr20.hisat

cd /home/ubuntu/hisat
hisat -p 16 -x ../genome/chr20.hisat -U ../rawdata/SRR1153470_1.fastq -S hisat.sam
samtools view -bSF4 hisat.sam | samtools sort - hisat.sorted
samtools index hisat.sorted.bam

