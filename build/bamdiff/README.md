#bamDiff

This takes as input a comma separated value (csv) file and two or more .bam 
files. The csv should be a file output by [runCsaw.R](../csaw/runCsaw.R) (e.g. column 2 seqid, column 3 start, column 4 stop, column -1 p-value) containing regions in which the 
.bam files show statistically different peaks. 

###STEP ONE:
SAM_comparitor calculates some summary statistics for each bam file, including 
the total number of reads contained, the total number of alignments (can be more 
than the number of reads if reads align multiple times), the proportion of 
reads that map uniquely (only once), and the proportion of reads which are 
unmapped. 

###STEP TWO:
For each region identified in the csv, provides text-based tables displaying the 
results for pairwise comparisons of the .bam files. This includes the number of 
reads in the identified region in one .bam file which are unmapped in the other. 
In addition, the number and location of mapping partners in the other .bam are 
analyzed for each mapping read, and the median values are reported. The location 
metric is simply the ration of mapping partners which fall within the identified 
region in the partner .bam file, out of all mapping partners in that file. 

###STEP THREE (yet to be implemented):
Include information about read quality (based on CIGAR strings -- perfect mapping,
short matches, indels, etc). 
Use annotation file (.gff) to look at where these alternate mappings are:
- pseudo genes
- repetitive regions
- alternative exon splicing sites?

###*NOTE:*
You can look at the behavior of reads in your favorite gene in the following way:
- Create a csv file:
	- Include a header (one line, can be anything)
	- Make the second line correspond to your gene, with these fields:
		- index,"chrom name",start,stop

--------------------------------------------------------------------------------

```Usage: bamDiff [OPTIONS] [-n int] [-o outputFile] regions.csv file1.bam file2.bam [file3.bam ...]```


```-h, --help```  
Display options and more info about bamDiff


```-s, --simple```  
Displays only the table reporting the number of 
unmapped reads.

```-S, --summarize```  
Displays only summary statistics for each of the bam
files.

```-v, --verbose```  
Displays summary statistics, tables for unmapped, 
map counts, and location info, as well as most
mapped regions.
DEFAULT VIEW BEHAVIOR: Unmapped tables, most mapped
regions.


```-o, --output```  
__*string*__	Specify rootname for output files. Otherwise, results 
will print to standard out.

```-n, --numRegions```  
__*int*__ 	Allows the user to specify the number of regions in
the csv to be examined.

```-p, --p-value```  
__*float*__	Sets a p-value threshold. All regions in the csv
file will a more significant (lower) p-value will be
examined.


```-a, --annotate```  
__*file*__ 	Not yet functional. In the future will allow you to 
examine the regions most mapped to and extract the
annotations in the gtf/gff file.

To see example usage and output, go to [our wiki tutorial](https://github.com/DCGenomics/RNA_mapping/wiki/7.-Compare-alignments-with-bamDiff).
