#!/usr/bin/env python

#=================================================================================
# Input: 
#   gff/merged-same-query/chr1.same.query.merged.gff.collapsed.sorted.by.start
#   gff/merged-same-query/chr1.same.query.merged.gff.sorted.by.start
# Read in the collapsed .gff files (containing gene info instead of
# exon info), remove overlaps between merged genes based on expect 
# value and seq identity et al. Matches are sorted by chrom_start.
# Output:
# ./gff/remove-overlap-all-query/chr*.all.query.removed.gff
#================================================================================== [zz]

import os, sys

from openOrFail import *

## labels for a row of data
dataLabels = ['F1', 'F2', 'F3', 'ChromStart', 'ChromEnd', 'Strand',
              'Query', 'GeneId', 'NumOfExons', 'UniqueName',
              'QueryStart', 'QueryEnd', 'QueryLen', 'Frac', 'Expect', 'Score', 'ExonBound']
## define each label as a manifest constant with its position in the list as its value.
for lIndex in xrange(len(dataLabels)): exec '%s = %d' % (dataLabels[lIndex], lIndex)

def fSignum(f):
    if f < 0.0: return -1
    if f > 0.0: return 1
    return 0

def removeCollapsedByExpect(cluster, outFile, geneIdHash):
# given a cluster of overlapping genes, process them in order of
# (lowest expect, hightest score). discard a gene if it has
# significant overlap with a previously selected gene, otherwise add
# it to the list of selected genes.

    # python2.3+ dependency: this works only with a stable sort
    # (potential bug in zz's original, given that perl's sort isn't
    # necessarily stable either)
    cluster.sort(lambda l, r: -fSignum(float(l[Score]) - float(r[Score])))
    cluster.sort(lambda l, r: fSignum(float(l[Expect]) - float(r[Expect])))

    pickedHits = []

    for candidate in cluster:
        cExonBound = candidate[ExonBound]

        for picked in pickedHits:
            if checkOverlap(cExonBound, picked[ExonBound]): break
        else:
            pickedHits.append(candidate)

    for picked in pickedHits:
        geneIdHash[picked[GeneId]] = 1 # keep track of picked genes to filter exons latter.
        outFile.write('\t'.join([str(picked[l]) for l in [F1, F2, F3, ChromStart, ChromEnd, Strand, Query, GeneId, NumOfExons,
                                                     UniqueName, QueryStart, QueryEnd, QueryLen, Frac, Expect, Score, ExonBound]]) + '\n')

def checkOverlap(boundStr1, boundStr2):

  #===============================================================================
  # Given two strings such as "1000..1500 1600..1888", check whether
  # the two strings overlap significantly, Strings are supposedly already sorted
  # return flag: flag = 0 : no overlap ; flag = 1: overlap 
  #=============================================================================== [zz]

    bounds1 = [(int(l), int(u)) for (l, d, u) in [p.split('.') for p in boundStr1.split()]]
    bounds2 = [(int(l), int(u)) for (l, d, u) in [p.split('.') for p in boundStr2.split()]]

    if not bounds1 or not bounds2: return 0
    (b1Lower, b1Upper) = bounds1.pop(0)
    (b2Lower, b2Upper) = bounds2.pop(0)
    
    # loop must terminate because every iteration we either find
    # overlap or eliminate the left-more interval. so even if we never
    # find overlap we will eventually exhaust one or both bounds list
    # and return no overlap. this assumes "supposedly already sorted"
    # comment above is indeed true!
    while 1:
        if b1Upper < b2Lower:
            if not bounds1: return 0
            (b1Lower, b1Upper) = bounds1.pop(0)
            continue

        if b2Upper < b1Lower:
            if not bounds2: return 0
            (b2Lower, b2Upper) = bounds2.pop(0)
            continue

        if b1Lower < b2Lower:
            if b1Upper > b2Upper: return 1
            if b1Upper > (30 + b2Lower - 1): return 1
            if not bounds1: return 0
            (b1Lower, b1Upper) = bounds1.pop(0)
        else:
            if b2Upper > b1Upper: return 1
            if b2Upper > (30 + b1Lower - 1): return 1
            if not bounds2: return 0
            (b2Lower, b2Upper) = bounds2.pop(0)


#main
if len(sys.argv) > 1:
    chrNames = sys.argv[1:]
else:
    chrNames = [str(i) for i in xrange(1,20)] + ['X']

print '...running %s ...' % sys.argv[0]
print '...PART 1 ...'

for chr in chrNames:
    print '...PART 1 %s...' % chr

    inFile = openOrFail('./gff/merged-same-query/%s.same.query.merged.gff.collapsed.sorted.by.start' % chr, 'r')

    try:    os.makedirs('./gff/remove-overlap-all-query/')
    except: pass

    outFile = openOrFail('./gff/remove-overlap-all-query/%s.all.query.removed.collapsed.gff' % chr, 'w')

    # copy header line
    outFile.write(inFile.next())

    # use data from first line to establish initial values.
    f = inFile.next()[:-1].split('\t')
    cluster = [f]
    (f[ChromStart], f[ChromEnd]) = (int(f[ChromStart]), int(f[ChromEnd]))
    lastChromEnd = f[ChromEnd]
    
    # used to track ids of genes that are kept.
    geneIdHash = {}

    for line in inFile:
        f = line[:-1].split('\t')

        # bug in original -- this cutoff does not apply to first line???
        if 1e-4 <= float(f[Expect]): continue

        (f[ChromStart], f[ChromEnd]) = (int(f[ChromStart]), int(f[ChromEnd]))

        if f[ChromStart] <= lastChromEnd:
            # overlap. add gene to current cluster
            cluster.append(f)
            if f[ChromEnd] > lastChromEnd: lastChromEnd = f[ChromEnd]
        else:
            # no overlap. process current cluster and then use this gene to start a new one.
            removeCollapsedByExpect(cluster, outFile, geneIdHash)
            cluster = [f]
            lastChromEnd = f[ChromEnd]

    removeCollapsedByExpect(cluster, outFile, geneIdHash)

    inFile.close()
    outFile.close()

    os.system('sort -g -k 4 %s > %s.sorted.by.start' % (outFile.name, outFile.name))

    # copy header line from the input exon file
    inFileExons = openOrFail('./gff/merged-same-query/%s.same.query.merged.gff' % chr, 'r')
    outFile = openOrFail('./gff/remove-overlap-all-query/%s.all.query.removed.gff' % chr, 'w')
    outFile.write(inFileExons.next())

    for line in inFileExons:
        f = line.split('\t')
        # remove exon counter from gene id, if present
        geneId = '==='.join(f[GeneId].split('===')[:2])

        # only print out those exons whose 'root' gene id is in geneIdHash
        if geneIdHash.has_key(geneId): outFile.write(line)

    inFileExons.close()
    outFile.close()

    os.system('sort -g -k 4 %s > %s.sorted.by.start' % (outFile.name, outFile.name))


