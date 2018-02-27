#!/bin/bash

#cd /home/ubuntu/genome
bwa index ../genome/chr20.fa

#cd /home/ubuntu/bwa
cd /RNA_mapping/build/bwa
bwa mem -t 16 ../genome/chr20.fa ../rawdata/SRR1153470_1.fastq > bwa.sam
samtools view -bSF4 bwa.sam | samtools sort -o bwa.sorted.bam
samtools index bwa.sorted.bam

