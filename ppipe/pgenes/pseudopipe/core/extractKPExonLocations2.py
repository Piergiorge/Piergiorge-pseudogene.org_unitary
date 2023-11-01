#!/usr/bin/env python

import sys, os

feature="CDS"
exons={}
for f in os.listdir(os.getcwd()):
	if f.endswith("gtf"):
		h=open(f)
		for l in h:
			fs=l.split("\t")
			if fs[1]=="protein_coding" and fs[2]==feature:
				if fs[0] not in exons:
					exons[fs[0]]=[]
				exons[fs[0]].append((int(fs[3]),int(fs[4])))

for chr in exons:
	o=open("chr%s_exLocs"%chr, "w")
	el=exons[chr]
	el.sort()
	for exon in el:
		o.write("%s\t%s\t%s\t%s\n" % (feature, chr, exon[0], exon[1]))
	o.close()
