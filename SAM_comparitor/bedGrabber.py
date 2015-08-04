#!/usr/bin/env python

import sys
import subprocess

inBed = sys.argv[1]
n = len(sys.argv)

f = open(inBed)

# For interesting regions checks where each read maps in the other aligners' outputs
print "starting analysis"
for line in f:
	mapped = []
	unmapped = []

	fields = line.strip().split()
	chrm   = fields[0]
	start  = fields[1]
	end    = fields[2]

	#Grab the reads that are mapped to the given regions
	print "Obtaining reads for region"
	for i, bam in enumerate(sys.argv):
		if i > 1:
			subprocess.call("samtools view -o tmp%d %s %s:%s-%s" % (i, bam, chrm, start, end), shell=True)
			subprocess.call("wc -l tmp%d" % i, shell=True)

	#For each bam file
	for i in range(2,n-1):
		mapped.append([])
		unmapped.append([])
		with open("tmp%d" % i) as g:
			for j in range(i+1,n):
				g.seek(0)
				if i != j:
					print "Comparing file", i-1, "vs. file", j-1
					mapped[-1].append([])
					unmapped[-1].append(0)
					# For each read within the region of bam file i:
					for line2 in g:
						um = False
						name = line2.split()[0]
						# Grab all locations it maps to in the other bam file
						subprocess.call("samtools view %s | grep %s > tmp%d" % (sys.argv[j], name, n), shell=True)
						# Grab all locations it maps to the other bam file in the same region
						subprocess.call("cat tmp%d | grep %s > tmp%d" % (j, name, n+1), shell=True) 
						with open("tmp%d" % n) as h:
							lc = 0 # The total number of times it was mapped in the other file
							for line3 in h:
								if line3 != "":
									lc += 1
									flag = int(line3.split()[1])
									if flag & 4:
										unmapped[-1][-1] += 1
										um = True
						if lc == 0:
							unmapped[-1][-1] += 1
							um = True

						with open("tmp%d" % (n+1)) as h:
							oc = 0 # The number of times it mapped into the same region
							for line3 in h:
								oc += 1

						outside = lc - oc
						inside  = oc
						if um:
							pass
							#if j == 4:
							#print "Unmapped"
						else:
							mapped[-1][-1].append(float(inside)/float(lc))
							#print mapped[-1][-1]
				 


#subprocess.call("rm tmp*", shell=True)

print mapped

print unmapped

print n
#for bami in mapped:
#	for bamii in bami:
#		print bamii