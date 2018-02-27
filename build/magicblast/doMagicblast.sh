#!/bin/bash

###
### use magicblast https://ncbiinsights.ncbi.nlm.nih.gov/2016/10/13/introducing-magic-blast/
###

cd /home/ubuntu/genome
makeblastdb -parse_seqids -dbtype nucl -title chr20 -out chr20.blast -in ../genome/chr20.fa

cd /home/ubuntu/blastmapper
sed -n '1~4s/^@/>/p;2~4p' ../rawdata/SRR1153470_1.fastq > SRR1153470_1.fasta
#blastmapper -num_threads 16 -query SRR1153470_1.fasta -db ../genome/chr20.blast -out blastmapper.sam -outfmt 15
magicblast -sra SRR1153470 -db /home/ubuntu/genome/chr20.blast -num_threads 16 -outfmt sam -out magicblast.sam
samtools view -bSF4 magicblast.sam | samtools sort - magicblast.sorted
samtools index magicblast.sorted.bam

