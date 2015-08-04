#!/usr/bin/env python

import sys
import subprocess

inBed = sys.argv[1]
n = len(sys.argv)

text_only = True

f = open(inBed)
f.readline()

# For interesting regions checks where each read maps in the other aligners' outputs
img_count = 0
print "starting analysis"
for line in f:
	img_count += 1
	mapped = []
	unmapped = []
	reads = []

	fields = line.strip().split(",")
	chrm   = fields[1][1:-1]
	start  = fields[2]
	end    = fields[3]

	#Grab the reads that are mapped to the given regions
	print "Obtaining reads for region"
	for i, bam in enumerate(sys.argv):
		if i > 1:
			subprocess.call("samtools view -o tmp%d %s %s:%s-%s" % (i, bam, chrm, start, end), shell=True)
			subprocess.call("wc -l tmp%d" % i, shell=True)

	#For each bam file
	for i in range(2,n):
		mapped.append([])
		reads.append([])
		unmapped.append([])
		with open("tmp%d" % i) as g:
			for j in range(2,n):
				g.seek(0)
				if i != j:
					print "Comparing file", i-1, "vs. file", j-1
					mapped[-1].append([])
					reads[-1].append([])
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
							#print "Unmapped"
						else:
							mapped[-1][-1].append(float(inside)/float(lc))
							reads[-1][-1].append(lc)
							#print mapped[-1][-1]

	subprocess.call("rm tmp*", shell=True)

	print mapped
	print unmapped
	print reads

	################################################################################
	if not text_only:
		import matplotlib.pyplot as plt

		k = n-2

		plt.figure()
		for i in range(k):
			l=-1
			for j in range(k):
				if i == j :
					plt.subplot( k, k, (i)*k + (j+1) )
					plt.plot(0,0)
					# Summary stuff goes here
				elif i != j:
					l+=1
					plt.subplot( 2*k, k, 2*i*k + (j+1))
					#plt.pie([unmapped[i][l], len(reads[i][l])],autopct='%1.1f%%')
					um = float(unmapped[i][l])
					tot = um + len(reads[i][l])
					if tot == 0:
						tot = 1

					plt.barh(0,1,1,color="blue")
					plt.barh(0,um/tot,1,color="orange")

					plt.subplot( 2*k, k, 2*i*k + (j+1) + k)
					if len(mapped[i][l]) != 0:
						plt.hist(mapped[i][l])

		plt.savefig("comp%d.png" % img_count)
		plt.close()

