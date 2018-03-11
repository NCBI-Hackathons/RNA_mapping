#!/bin/bash

# https://sites.stanford.edu/abms/content/giab-reference-materials-and-data
# https://catalog.coriell.org/0/Sections/Search/Sample_Detail.aspx?Ref=GM12878

# http://www.ncbi.nlm.nih.gov/sra/SRX457730
# 100 ng of Poly-A-RNA
# 115359773 reads (2 x 101bp Illumina paired end)

/opt/sratoolkit.2.9.0-ubuntu64/bin/fastq-dump --split-files SRR1153470

