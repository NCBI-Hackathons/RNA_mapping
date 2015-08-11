#!/bin/bash

cd /home/ubuntu/genome
bwa index chr20.fa

cd /home/ubuntu/bwa
bwa mem -t 16 ../genome/chr20.fa ../rawdata/SRR1153470_1.fastq > bwa.sam
samtools view -bSF4 bwa.sam | samtools sort - bwa.sorted
samtools index bwa.sorted.bam

