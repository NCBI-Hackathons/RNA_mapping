#!/usr/bin/env python
intro = """
################################################################################
#                               SAM_COMPARITOR                                 #
################################################################################
This takes as input a comma separated value (csv) file and two or more .bam 
files. The csv should be a file output by CSAW containing regions in which the 
.bam files show statistically different peaks.

STEP ONE:
SAM_comparitor calculates some summary statistics for each bam file, including 
the total number of reads contained, the total number of alignments (can be more 
than the number of reads if reads align multiple times), the proportion of 
reads that map uniquely (only once), and the proportion of reads which are 
unmapped.

STEP TWO:
For each region identified in the csv, provides text-based tables displaying the 
results for pairwise comparisons of the .bam files. This includes the number of
reads in the identified region in one .bam file which are unmapped in the other.
In addition, the number and location of mapping partners in the other .bam are 
analyzed for each mapping read, and the median values are reported. The location 
metric is simply the ration of mapping partners which fall within the identified
region in the partner .bam file, out of all mapping partners in that file.

STEP THREE (yet to be implemented):
Include information about read quality (based on CIGAR strings -- perfect mapping,
short matches, indels, etc)
Use annotation file (.gff) to look at where these alternate mappings are:
	- pseudo genes
	- repetitive regions
	- alternative exon splicing sites?

** Include command line options to simplify/expand output.
** Allow options for specifying how many regions in the csv to analyze (use 
   head or something).
** Allow to pick reads in the region, and let you see WHERE those reads are 
   mapping in the other file.

################################################################################

*NOTE:
You can look at the behavior of reads in your favorite gene in the following way:
- Create a csv file:
	- Include a header (one line, can be anything)
	- Make the second line correspond to your gene, with these fields:
		- index,"chrom name",start,stop

################################################################################
"""

import sys
import subprocess

if len(sys.argv) < 4:
	print "\nUsage: %s regions.csv file1.bam file2.bam [file3.bam ...]"
	if len(sys.argv) > 1:
		if "-h" in sys.argv[1]:
			print intro
	quit()

inBed = sys.argv[1]
n = len(sys.argv)

text_only = True

bed = open(inBed)
bed.readline()

################################################################################

def median(l):
	List = sorted(l)
	i = len(l)
	if i == 0:
		return -1
	if i%2 == 0:	med = float(List[i//2] + List[i//2-1])/2
	else:			med = List[i//2]
	return med

################################################################################

print "File key:"
print "#\tNAME"
i = 0
for bam in sys.argv[2:]:
	i += 1
	print "%d\t%s" % (i, sys.argv[i+1])

################################################################################

print "\nSummary statistics:"
for i in range(2,n):
	f_name = sys.argv[i]
	subprocess.call("samtools sort -n %s tmp" % (f_name), shell=True)
	subprocess.call("samtools view tmp.bam > tmp1", shell=True)
	subprocess.call('rm tmp.*', shell=True)
	with open("tmp1") as f:
		unmapped = count = primary = not_unique = 0
		lastName = ""

		while True:
			line = f.readline()
			if line == "": break
			elif not line.startswith("@"):
				fields = line.split()
				name = fields[0]
				flag = int(fields[1])
				count += 1
				if flag & 4:       unmapped += 1
				if not flag & 256: primary  += 1
				if name != lastName:
					currPos = f.tell()
					nextName = f.readline()
					if nextName == "": break
					else: nextName = nextName.split()[0]
					if name == nextName:
						not_unique += 1
					f.seek(currPos)
				lastName = name

			unique = primary - not_unique
			pcnt_unique = float(unique)/primary*100
			pcnt_um = float(unmapped)/primary*100

	print f_name
	print "    Total Number of reads:   %d" % (primary)
	print "    Number of alignments:    %d" % (count - unmapped)
	print "    Uniquely mapping reads:  %1.1f %% (%d/%d)" % (pcnt_unique, unique, primary)
	print "    Unmapped reads:          %1.1f %% (%d/%d)\n" % (pcnt_um, unmapped, primary)


################################################################################
# For interesting regions checks where each read maps in the other aligners' outputs
img_count = 0

for line in bed:
	img_count += 1
	mapped = []
	unmapped = []
	reads = []

	fields = line.strip().split(",")
	chrm   = fields[1][1:-1]
	start  = fields[2]
	end    = fields[3]

	#Grab the reads that are mapped to the given regions
	#print "Obtaining reads for region"
	for i, bam in enumerate(sys.argv):
		if i > 1:
			subprocess.call("samtools view -o tmp%d %s %s:%s-%s" % (i-1, bam, chrm, start, end), shell=True)
			#subprocess.call("wc -l tmp%d" % (i-1), shell=True)

	#For each bam file
	for i in range(2,n):
		mapped.append([])
		reads.append([])
		unmapped.append([])
		with open("tmp%d" % (i-1)) as g:
			names = []
			for line2 in g:
				name = line2.split()[0]
				names.append(name)
			for j in range(2,n):
				if i != j:
					#print "Comparing file", i-1, "vs. file", j-1
					mapped[-1].append([])
					reads[-1].append([])
					unmapped[-1].append(0)
					inBam = {}
					inTmp = {}
					subprocess.call("samtools view %s > tmp-1" % sys.argv[j], shell = True)
					with open("tmp-1") as h:
						for l in h:
							name = l.split()[0]
							if name in names:
								if name not in inBam:
									inBam[name] = []
								inBam[name].append(l)
					subprocess.call("cat tmp%d > tmp-1" % (j-1), shell = True)
					with open("tmp-1") as h:
						for l in h:
							name = l.split()[0]
							if name in names:
								if name not in inTmp:
									inTmp[name] = []
								inTmp[name].append(l)

					for name in names:
						um = False
						if name not in inBam:
							unmapped[-1][-1] += 1
							um = True
						else:
							if name in inTmp:
								oc = len(inTmp[name])
							else:
								oc = 0
							lc = len(inBam[name])

							if lc == 1:
								flag = int(inBam[name][0].split()[1])
								if flag & 4:
									unmapped[-1][-1] += 1
									um = True
						if not um:
							#outside = lc - oc
							#inside  = oc
							mapped[-1][-1].append(float(oc)/float(lc))
							reads[-1][-1].append(lc)


	subprocess.call("rm tmp*", shell=True)

	#print mapped
	#print unmapped
	#print reads

	print ""
	print "Region %s:%s-%s" % (chrm, start, end)
	print "="*(16+8*(n-2))
	print " "*16 + "|" + " "*7 + "vs. file"
	print "File\t# Reads\t|"+"\t".join([str(x) for x in range(1,n-1)])
	print "-"*16 + "+" + "-"*(8*(n-2)-1)
	print "Proportion of unmapped reads"
	print "-"*16 + "+" + "-"*(8*(n-2)-1)
	for i in range(n-2):
		l = -1
		print "%d\t%d\t|" % (i+1, unmapped[i][0] + len(reads[i][0])),
		for j in range(n-2):
			if i != j:
				l += 1
				um = unmapped[i][l]
				tot = um + len(reads[i][l])
				if tot == 0:
					tot = 1
				#print "%d/%d\t" % (um, tot),
				print "%1.2f\t" % (float(um)/tot),
			else:
				print "-\t",
		print ""
	print "-"*16 + "+" + "-"*(8*(n-2)-1)
	print "Median # times each read mapped" # include ranges here as well?
	print "-"*16 + "+" + "-"*(8*(n-2)-1)
	for i in range(n-2):
		l = -1
		print "%d\t%d\t|" % (i+1, unmapped[i][0] + len(reads[i][0])),
		for j in range(n-2):
			if i != j:
				l += 1
				print "%1.1f\t" % (median(reads[i][l])),
			else: print "-\t",
		print ""
	print "-"*16 + "+" + "-"*(8*(n-2)-1)
	print "Median proportion times each read mapped inside region" # include ranges here as well?
	print "-"*16 + "+" + "-"*(8*(n-2)-1)
	for i in range(n-2):
		l = -1
		print "%d\t%d\t|" % (i+1, unmapped[i][0] + len(reads[i][0])),
		for j in range(n-2):
			if i != j:
				l += 1
				print "%1.1f\t" % (median(mapped[i][l])),
			else: print "-\t",
		print ""
	print "="*(16+8*(n-2))

#	################################################################################
#	if not text_only:
#		import matplotlib.pyplot as plt
#
#		k = n-2
#
#		plt.figure()
#		for i in range(k):
#			l=-1
#			for j in range(k):
#				if i == j :
#					plt.subplot( k, k, (i)*k + (j+1) )
#					plt.plot(0,0)
#					# Summary stuff goes here
#				elif i != j:
#					l+=1
#					plt.subplot( 2*k, k, 2*i*k + (j+1))
#					#plt.pie([unmapped[i][l], len(reads[i][l])],autopct='%1.1f%%')
#					um = float(unmapped[i][l])
#					tot = um + len(reads[i][l])
#					if tot == 0:
#						tot = 1
#
#					plt.barh(0,1,1,color="blue")
#					plt.barh(0,um/tot,1,color="orange")
#
#					plt.subplot( 2*k, k, 2*i*k + (j+1) + k)
#					if len(mapped[i][l]) != 0:
#						plt.hist(mapped[i][l])
#
#		plt.savefig("comp%d.png" % img_count)
#		plt.close()


