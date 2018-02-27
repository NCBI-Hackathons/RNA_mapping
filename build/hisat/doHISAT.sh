#!/bin/bash

#cd /home/ubuntu/genome
hisat2-build ../genome/chr20.fa ../genome/chr20.hisat

#cd /home/ubuntu/hisat
hisat2 -p 16 -x ../genome/chr20.hisat -U ../rawdata/SRR1153470_1.fastq -S hisat.sam
samtools view -bSF4 hisat.sam | samtools sort - hisat.sorted
samtools index hisat.sorted.bam

