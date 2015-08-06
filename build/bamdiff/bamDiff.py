#!/usr/bin/env python

intro = """
################################################################################
#                                   bamDiff                                    #
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

################################################################################

*NOTE:
You can look at the behavior of reads in your favorite gene in the following way:
- Create a csv file:
    - Include a header (one line, can be anything)
    - Make the second line correspond to your gene, with these fields:
        - index,"chrom name",start,stop

################################################################################
"""

usage = """
\nUsage: %s [OPTIONS] [-n int] [-o outputFile] regions.csv file1.bam file2.bam [file3.bam ...]

-h, --help                  Display options and more info about bamDiff

-s, --simple                Displays only the table reporting the number of 
                            unmapped reads.
-S, --summarize             Displays only summary statistics for each of the bam
                            files.
-v, --verbose               Displays summary statistics, tables for unmapped, 
                            map counts, and location info, as well as most
                            mapped regions.
                            DEFAULT VIEW BEHAVIOR: Unmapped tables, most mapped
                            regions.

-o, --output        string  Specify rootname for output files. Otherwise, results
                            will print to standard out.
-n, --numRegions    int     Allows the user to specify the number of regions in
                            the csv to be examined.
-p, --p-value       float   Sets a p-value threshold. All regions in the csv
                            file will a more significant (lower) p-value will be
                            examined.

-a, --annotate      file    Not yet functional. In the future will allow you to 
                            examine the regions most mapped to and extract the
                            annotations in the gtf/gff file.
"""


import sys
import subprocess
import getopt

################################################################################

summary = False
simple  = False
default = True
verbose = False

f_out   = False

num_regions = 1
pThresh    = 10**10

help       = False
annotation = False

################################################################################

options,remainder = getopt.gnu_getopt(
                                        sys.argv[1:], 
                                        "o:n:p:a:vsSh", 
                                        [
                                            'verbose', 
                                            'simple',
                                            'summarize',
                                            'numRegions=',
                                            'output=',
                                            'help',
                                            'annotate=',
                                            'p-value='
                                            ])

if len(sys.argv) < 4:
    print usage
    print "For more help, enter -h at command line"
    quit()

for opt, arg in options:
    if opt in ('-v', '--verbose'):
        verbose = True
    elif opt in ('-s', '--simple'):
        simple = True
    elif opt in ('-S','--summarize'):
        summary = True
    elif opt in ('-n', '--numRegions'):
        num_regions = int(arg)
    elif opt in ('-o', '--output'):
        f_out = arg
    elif opt in ('-h', '--help'):
        help = True
    elif opt in ('-a', '--annotate'):
        annotation = arg
    elif opt in ('-p', '--p-value'):
        pThresh = float(arg)
        num_regions = 10**100

if help:
    print usage
    print intro
    quit()

if verbose:
    default = simple = summary = False
elif simple:
    default = summary = False
elif summary:
    default = False
elif default:
    pass
else:
    print "Error encountered in specifying output"
    quit()

inBed = remainder[0]
bams  = remainder[1:]
bam_names = [x.split("/")[-1] for x in bams]
k = len(remainder) - 1

################################################################################

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
    if i%2 == 0:    med = float(List[i//2] + List[i//2-1])/2
    else:            med = List[i//2]
    return med

################################################################################
if f_out:
    sys.stdout = open("%s.results" % f_out, "wb")
print "File key:"
print "#\tNAME"
i = 0
for bam in bams:
    i += 1
    print "%d\t%s" % (i, bams[i-1])

################################################################################

if summary or verbose:
    if f_out:
        sys.stdout = open("%s.summary" % f_out, "wb")
    print "\nSummary statistics:"
    for i in range(k):
        f_name = bams[i]
        subprocess.call("samtools sort -n %s tmp" % (f_name), shell=True)
        subprocess.call("samtools view tmp.bam > tmp-1", shell=True)
        subprocess.call('rm tmp.*', shell=True)
        with open("tmp-1") as f:
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

                if not primary == 0:
                    unique = primary - not_unique
                    pcnt_unique = float(unique)/primary*100
                    pcnt_um = float(unmapped)/primary*100
                else:
                    unique = pcnt_unique = pcnt_um = 0

        print f_name
        print "    Total Number of reads:   %d" % (primary)
        print "    Number of alignments:    %d" % (count - unmapped)
        print "    Uniquely mapping reads:  %1.1f %% (%d/%d)" % (pcnt_unique, unique, primary)
        print "    Unmapped reads:          %1.1f %% (%d/%d)\n" % (pcnt_um, unmapped, primary)


################################################################################
# For interesting regions checks where each read maps in the other aligners' outputs
#img_count = 0

csv_count = 0

if f_out:
    sys.stdout = open("%s.results" % f_out, "ab")

while csv_count < num_regions:
    csv_count += 1
    line = bed.readline()
    if line == "":
        break
    #img_count += 1
    mapped       = []
    unmapped     = []
    reads        = []
    superRegions = []

    fields = line.strip().split(",")
    chrm   = fields[1][1:-1]
    start  = fields[2]
    end    = fields[3]
    p      = fields[-1]
    if "e" in p:
    	base, exp = p.split("e")
    	p = float(base)*(10**(float(exp)))
    else:
    	p = float(p)

    if p > pThresh:
        break

    ############################################################
    #WHERE THE MAGIC HAPPENS
    #Grab the reads that are mapped to the given regions
    for i, bam in enumerate(bams):
        subprocess.call("samtools view -o tmp%d %s %s:%s-%s" % (i+1, bam, chrm, start, end), shell=True)
        #subprocess.call("wc -l tmp%d" % (i-1), shell=True)

    #For each bam file
    for i in range(k):
        mapped.append([])
        reads.append([])
        unmapped.append([])
        superRegions.append([])
        with open("tmp%d" % (i+1)) as g:
            names = []
            for line2 in g:
                name = line2.split()[0]
                names.append(name)
            for j in range(k):
                if i != j:
                    #print "Comparing file", i-1, "vs. file", j-1
                    mapped[-1].append([])
                    reads[-1].append([])
                    unmapped[-1].append(0)
                    inBam = {}
                    inTmp = {}
                    subprocess.call("samtools view %s > tmp-1" % bams[j], shell = True)
                    with open("tmp-1") as h:
                        for l in h:
                            name = l.split()[0]
                            if name in names:
                                if name not in inBam:
                                    inBam[name] = []
                                inBam[name].append(l)
                    subprocess.call("cat tmp%d > tmp-1" % (j+1), shell = True)
                    with open("tmp-1") as h:
                        for l in h:
                            name = l.split()[0]
                            if name in names:
                                if name not in inTmp:
                                    inTmp[name] = []
                                inTmp[name].append(l)

                    regions = {}
                    buff = 1000
                    for key in inBam:
                        for l in inBam[key]:
                            rPos   = int(l.split()[3])
                            rChrom = l.split()[2]
                            added = False
                            if rChrom not in regions:
                                regions[rChrom] = []
                            if len(regions[rChrom]) > 0:
                                for i, (rSt, rEnd, c) in enumerate(regions[rChrom]):
                                    nSt = rSt - buff
                                    nEnd = rEnd + buff
                                    if rPos > nSt and rPos < nEnd:
                                        regions[rChrom][i][0] = min(rPos,rSt)
                                        regions[rChrom][i][1] = max(rPos,rEnd)
                                        added = True
                                        regions[rChrom][i][2] += 1
                            if not added:
                                regions[rChrom].append([rPos, rPos+1, 1])

                    for key in regions:
                        regions[key].sort(reverse=True, key=lambda x: x[2])

                    superRegions[-1].append(regions)
                    
                    for name in names:
                        um = False
                        if name not in inBam:
                            unmapped[-1][-1] += 1
                            um = True
                        else:
                            if name in inTmp:
                                inside = len(inTmp[name])
                            else:
                                inside = 0
                            lc = len(inBam[name])

                            if lc == 1:
                                flag = int(inBam[name][0].split()[1])
                                if flag & 4:
                                    unmapped[-1][-1] += 1
                                    um = True
                        if not um:
                            mapped[-1][-1].append(float(inside)/float(lc))
                            reads[-1][-1].append(lc)

    subprocess.call("rm tmp*", shell=True)

    #print mapped
    #print unmapped
    #print reads

    ############################################################
    # REPORT RESULTS
    if not summary:
        print ""
        print "Region %s:%s-%s" % (chrm, start, end)
        print "="*(16+8*(k))
        print " "*16 + "|" + " "*7 + "vs. file"
        print "File\t# Reads\t|"+"\t".join([str(x) for x in range(1,k+1)])
        print "-"*16 + "+" + "-"*(8*(k)-1)
    if simple or default or verbose:
        print "Proportion of unmapped reads"
        print "-"*16 + "+" + "-"*(8*(k)-1)
        for i in range(k):
            l = -1
            print "%d\t%d\t|" % (i+1, unmapped[i][0] + len(reads[i][0])),
            for j in range(k):
                if i != j:
                    l += 1
                    um = unmapped[i][l]
                    tot = um + len(reads[i][l])
                    if tot == 0:
                        tot = 1
                    print "%1.2f\t" % (float(um)/tot),
                else:
                    print "-\t",
            print ""
        print "-"*16 + "+" + "-"*(8*(k)-1)
    if verbose:
        print "Median # times each read mapped" # include ranges here as well?
        print "-"*16 + "+" + "-"*(8*(k)-1)
        for i in range(k):
            l = -1
            print "%d\t%d\t|" % (i+1, unmapped[i][0] + len(reads[i][0])),
            for j in range(k):
                if i != j:
                    l += 1
                    print "%1.1f\t" % (median(reads[i][l])),
                else: print "-\t",
            print ""
        print "-"*16 + "+" + "-"*(8*(k)-1)
        print "Median proportion times each read mapped inside region" # include ranges here as well?
        print "-"*16 + "+" + "-"*(8*(k)-1)
        for i in range(k):
            l = -1
            print "%d\t%d\t|" % (i+1, unmapped[i][0] + len(reads[i][0])),
            for j in range(k):
                if i != j:
                    l += 1
                    print "%1.1f\t" % (median(mapped[i][l])),
                else: print "-\t",
            print ""
    print "="*(16+8*(k))
    if verbose or default:
        print "\nRegions where most identified reads are mapping in second bam file"
        for i in range(k):
            l = -1
            for j in range(k):
                if i != j:
                    l += 1
                    print "Reads in region from file %s mapped to file %s at:" % (bam_names[i] , bam_names[j])
                    for key in superRegions[i][l]:
                        nn = 0
                        for st, stp, ct in superRegions[i][l][key]:
                            nn += 1
                            if nn >= 10:
                                break
                            print "%s\t%d\t%d\t%d" % (key, st, stp, ct)
                    print ""

#    ################################################################################
#    if not text_only:
#        import matplotlib.pyplot as plt
#
#        plt.figure()
#        for i in range(k):
#            l=-1
#            for j in range(k):
#                if i == j :
#                    plt.subplot( k, k, (i)*k + (j+1) )
#                    plt.plot(0,0)
#                    # Summary stuff goes here
#                elif i != j:
#                    l+=1
#                    plt.subplot( 2*k, k, 2*i*k + (j+1))
#                    #plt.pie([unmapped[i][l], len(reads[i][l])],autopct='%1.1f%%')
#                    um = float(unmapped[i][l])
#                    tot = um + len(reads[i][l])
#                    if tot == 0:
#                        tot = 1
#
#                    plt.barh(0,1,1,color="blue")
#                    plt.barh(0,um/tot,1,color="orange")
#
#                    plt.subplot( 2*k, k, 2*i*k + (j+1) + k)
#                    if len(mapped[i][l]) != 0:
#                        plt.hist(mapped[i][l])
#
#        plt.savefig("comp%d.png" % img_count)
#        plt.close()


