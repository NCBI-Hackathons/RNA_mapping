#!/bin/bash

cd /home/ubuntu/blastmapper
rm blastmapper.*
rm SRR1153470_1.fasta

cd /home/ubuntu/bwa
rm bwa.*

cd /home/ubuntu/csaw
rm csaw.bam.results.csv

cd /home/ubuntu/genome
rm -r chr20*.*
rm Log.out

cd /home/ubuntu
rm -r ncbi

cd /home/ubuntu/hisat
rm hisat.*

cd /home/ubuntu/rawdata
rm SRR1153470_?.fastq

cd /home/ubuntu/star
rm star*.*
