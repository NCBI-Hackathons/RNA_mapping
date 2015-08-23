#!/bin/bash

# http://ncbiinsights.ncbi.nlm.nih.gov/2013/12/24/introducing-the-new-human-genome-assembly-grch38/

# http://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.30
# http://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/

wget -O GRCh38.p4.gff.gz ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000001405.30_GRCh38.p4/GCF_000001405.30_GRCh38.p4_genomic.gff.gz
gunzip GRCh38.p4.gff.gz
grep "^NC_000020.11" GRCh38.p4.gff > GRCh38.p4.chr20.gff 

wget -O chr20.fa.gz ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000001405.30_GRCh38.p4/GCF_000001405.30_GRCh38.p4_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/chr20.fna.gz
gunzip chr20.fa.gz

#wget -O chr20.fa.gz ftp.ensembl.org/pub/release-81/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.20.fa.gz
#gunzip chr20.fa.gz

#wget -O GRCh38.Ensembl.81.gtf.gz ftp.ensembl.org/pub/release-81/gtf/homo_sapiens/Homo_sapiens.GRCh38.81.gtf.gz
#gunzip GRCh38.Ensembl.81.gtf.gz

