#!/bin/bash

cd /home/ubuntu/genome
mkdir chr20.STAR
STAR --runMode genomeGenerate --runThreadN 16 --genomeDir chr20.STAR --genomeFastaFiles chr20.fa

cd /home/ubuntu/star
STAR --genomeDir ../genome/chr20.STAR --runThreadN 16 --readFilesIn ../rawdata/SRR1153470_1.fastq --outFileNamePrefix star
samtools view -bSF4 starAligned.out.sam | samtools sort - star.sorted
samtools index star.sorted.bam

