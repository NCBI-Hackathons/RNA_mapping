#!/bin/bash

cd /home/ubuntu/genome
makeblastdb -parse_seqids -dbtype nucl -title chr20 -out chr20.blast -in ../genome/chr20.fa

cd /home/ubuntu/blastmapper
sed -n '1~4s/^@/>/p;2~4p' ../rawdata/SRR1153470_1.fastq > SRR1153470_1.fasta
blastmapper -num_threads 16 -query SRR1153470_1.fasta -db ../genome/chr20.blast -out blastmapper.sam -outfmt 15
samtools view -bSF4 blastmapper.sam | samtools sort - blastmapper.sorted
samtools index blastmapper.sorted.bam

