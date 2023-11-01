#!/usr/bin/env python

#===============================================================================
#  Infile:
#        ./gff/ensembl-gene-masked/chr1.ens.gene.masked
#
#  Outfile:
#	 ./gff/remove-overlap-between-same-query/
# Inputs are have been sorted by chrom_start, masked 
# out overlapping with Ensembl gene annotations, this script will remove overlaps
# between matches. 
# However, if two overlapping matches belong to different query, then both are kept.
#================================================================================== [zz]

import os, sys

from openOrFail import *

# labels for a row of data
(Query, ChrId, Ident, Overlap, Positives, Gaps, Strand, QueryStart,
 QueryEnd, QueryLen, MatchStart, MatchEnd, EVal, Score, Original) = xrange(15)


def fSignum(f):
    if f < 0.0: return -1
    if f > 0.0: return 1
    return 0

# reduce a collection of overlapping hits to a set of (nearly)
# non-overlapping hits. this is done by processing the hits from best
# to worse score. a hit is added to the list of picked hits if it does
# not significantly overlap with any already picked hit. there appears
# to be one exception: if the overlap is with a hit from a different
# query and the candidate hit is "sufficiently" good, this overlap is
# allowed and the check for overlaps continues with the next picked
# hit.
#
# the routine writes to the output file the final list of picked hits.
def removeOnlySameQuery(cluster, outFile):
    # bug in original code? perl sort is not necessarily stable
    # python2.3+ dependency: this works only with a stable sort
    cluster.sort(lambda l, r: -fSignum(l[Score] - r[Score]))
    cluster.sort(lambda l, r: fSignum(l[EVal] - r[EVal]))

    # because the keys are sorted, when a candidate hit is picked, it has
    # a higher [evalXscore] ranking than any hit that follows.
    picked = []
    for c in cluster:
        (cStart, cEnd) = (c[MatchStart], c[MatchEnd])

        for p in picked:
            (pStart, pEnd) = (p[MatchStart], p[MatchEnd])

            # check for overlap with this already picked hit. if none,
            # check next picked hit.
            ol  = cEnd <= pEnd and cStart >= pStart            # p completely overlaps c
            ol |= cEnd >= pEnd and cStart <= pStart            # c completely overlaps p
            ol |= cEnd <= pEnd and cEnd   >= (pStart + 30 - 1) # c overlaps p by 30 or more
            ol |= cEnd >= pEnd and cStart <= (pEnd - 30 - 1 )  # p overlaps c by 30 or more
            if not ol: continue

            # there is overlap. if the hits are for the same query,
            # don't add this candidate hit to the picked hit list.
            if p[Query] == c[Query]: break

            # there is overlap, but between different queries. "decide
            # whether to remove the overlapping match based on the
            # value of eval" [zz]
                
            # Unless this is a really good match (either absolutely or
            # relative to the picked scores), discard. (remember,
            # smaller is better for eVal.) this still seems to be
            # broken: consider p[EVal] = 1e-20, then c[EVal] need only
            # be < 1 (not a very demanding eVal).
            if 0.0 == p[EVal]:
                if c[EVal] > 1e-20: break
            else:
                if (c[EVal]/p[EVal]) > 1e20: break
        else:
            # we did not eliminate this candidate, so add it to the picked list.
            picked.append(c)

    for p in picked: outFile.write(p[Original])

# main
if len(sys.argv) > 1:
    chrNames = sys.argv[1:]
else:
    chrNames = [str(i) for i in xrange(1,20)] + ['X']

print 'issues: 1e-20, 1e-3, 30, bugs in eval comparison logic.'

print '...running %s ...' % sys.argv[0]

for chr in chrNames:
    print '...working on %s ...' % chr

    inFile = openOrFail('./gff/ensembl-gene-masked/%s.ens.gene.masked' % chr, 'r')
    print '...open %s...' % inFile.name

    # sorted blastout file [zz]
    try:    os.makedirs('./gff/remove-overlap-same-query/')
    except: pass

    outFile = openOrFail('./gff/remove-overlap-same-query/%s.same.query.removed.gff' % chr, 'w')
    print '...open %s...' % outFile.name

    # skip header line of input.
    inFile.next()
    outFile.write('#query\tchr\tident\toverlap\tpositives\tgaps\tstrand\tquery_start\tquery_end\tquery_len\tmatch_start\tmatch_end\teval\tscore\n')

    # run through the blast hits, collecting overlapping hits into a
    # cluster. call removeOnlySameQuery to process the cluster and
    # write to the outfile non-overlapping hits.
    cluster = []
    for l in inFile:
        f = l[:-1].split() + [l]

        # convert from string to numeric values those quantities that
        # will be used for comparions.
        (f[MatchStart], f[MatchEnd]) = (int(f[MatchStart]), int(f[MatchEnd]))
        (f[EVal], f[Score]) = (float(f[EVal]), float(f[Score]))

        # zz's original code would not have exclude the first entry of
        # a cluster because of a high e value, but that seems more
        # like a glitch than anything else.
        if f[EVal] >= 1e-3: continue

        if not cluster:
            # start first cluster
            cluster = [f]
            lastMatchEnd = f[MatchEnd]
            continue
        
        if f[MatchStart] <= lastMatchEnd:
            cluster.append(f)
            if f[MatchEnd] > lastMatchEnd: lastMatchEnd = f[MatchEnd]
        else:
            # process existing cluster
            removeOnlySameQuery(cluster, outFile)
            # start new cluster
            cluster = [f]
            lastMatchEnd = f[MatchEnd]

    #process last cluster
    removeOnlySameQuery(cluster, outFile)

    inFile.close()
    outFile.close()
