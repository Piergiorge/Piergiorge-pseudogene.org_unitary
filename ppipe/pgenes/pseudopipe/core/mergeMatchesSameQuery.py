#!/usr/bin/env python

#===================================================================
#  Infile:
#     gff/remove-overlap-same-query/chr1.same.query.removed.gff
#  Outfile: 
#  ./gff/merged-same-query/chr*.same.query.merged.gff  
#  ./gff/merged-same-query/chr*.same.query.merged.gff.collapsed
#
#  for matches of the same query, merge adjacent matches or combine
#  them into genes, write out the matches in fasta file, also create
#  new .gff files
#=================================================================== [zz]

import os, sys

from openOrFail import *

# labels for a row of data -- we add last two values here.
dataLabels = ['Query', 'ChrId', 'Ident', 'Overlap', 'Positives', 'Gaps', 'Strand',
              'QueryStart', 'QueryEnd', 'QueryLen', 'ChromStart', 'ChromEnd',
              'Expect', 'Score'] + ['MatchIdBug', 'KeyTag']
# define each label as a manifest constant with its position in the list as its value.
for lIndex in xrange(len(dataLabels)): exec '%s = %d' % (dataLabels[lIndex], lIndex)

# to simplify output comparison, print score using these rules.
def printScore(score):
    s = '%.2f' % score
    if len(s) > 5: s = s[:5] + s[5:].rstrip('0').rstrip('.')
    return s

# Some terminology:
#
#   Query: a protein sequence
#
#   Chromosome: dna sequence
#
#   Match: a blast 'hit' describing an area of similarity between a
# query and a chromosome. A query may give rise to many hits.
#
#   ROS: short for "region of similarity". The actual region of the
# query that blast aligned with the chromsome.  In general something
# less than the full Query sequence.
#
# To some degree, 'hit', 'match' and 'exon' appear to synonymous here.

def checkIfSameGene(prevStart, prevEnd, start, end, strand, queryLen):
# despite the name, this appears to check if two successive hits could
# be successive exons. this check boils down to the degree of overlap
# in the ROS.  returns true if:
#
#   1) previous ROS overlaps and is nearly co-terminant with
#   current ROS and the current ROS is less than or equal to the smaller of
#   20 aas or 20% of the whole current query (not just the ROS).
#   
#   2) current ROS overlaps and extends beyond previous ROS and the
#   length of overlap is less than or equal to the smaller of 20 aas
#   or 20% of the current query as a whole.
#
#   3) There is no overlap at all. 
#
# It is not clear how much cases 1) and 2) really differ.
    flag = 0

    overlapLimit = 20
    if queryLen > 100: overlapLimit = 0.20 * queryLen
    
    if '+' == strand:
        if start >= prevStart \
            and end >= prevEnd - 2 \
            and end <= prevEnd \
            and (end - start + 1) <= overlapLimit: flag = 1
        elif start >= prevStart \
            and start <= prevEnd \
            and end >= prevEnd \
            and (prevEnd - start + 1) <= overlapLimit: flag = 1
        elif start >= prevEnd: flag = 1
    else:
        if prevEnd >= end \
            and start >= prevStart \
            and start <= prevStart + 2 \
            and (end - start + 1 ) <= overlapLimit: flag = 1
        elif  prevEnd >= end \
            and end >= prevStart \
            and start <= prevStart \
            and (end - prevStart + 1) <= overlapLimit: flag = 1
        elif end <= prevStart: flag = 1

    return flag;

#main
print 'need to document: overlap limit (20 or 20%), chrom separation of 10000, gap of 60'

if len(sys.argv) > 1:
    chrNames = sys.argv[1:]
else:
    chrNames = [str(i) for i in xrange(1,20)] + ['X']

print '...running %s ...' % sys.argv[0]

queryGapCutoff = 60

for chr in chrNames:
    print '... %s...' % chr

    inFile = openOrFail('./gff/remove-overlap-same-query/%s.same.query.removed.gff' % chr, 'r')

    # build a table of lists. the first level table is indexed by
    # query and yields a list of matches (blast hits) for that query.
    gffHash = {}
    for line in inFile:
        f = line[:-1].split('\t')
        if 0 == len(f) or '#' == f[0][0]: continue
        f += [-1, ''] # add fields MatchIdBug and KeyTag
        
        # convert some numeric fields from strings to numbers.
        for l in [ChromStart, ChromEnd, QueryStart, QueryEnd, QueryLen]: f[l] = int(f[l])

        try:    gffHash[f[Query]].append(f)
        except: gffHash[f[Query]] = [f]

    inFile.close()

    if not gffHash: continue

    try:    os.makedirs('./gff/merged-same-query/')
    except: pass

    outFile = openOrFail('./gff/merged-same-query/%s.same.query.merged.gff' % chr, 'w')
    outFile.write('\t'.join(['#chr', 'merged', 'merged', 'chrom_start', 'chrom_end', 'strand',
                             'query', 'match_id', 'exonRank', 'num_of_exons', 'unique_name',
                             'query_start', 'query_end', 'query_len', 'expect', 'score'])+'\n')

    outFile2 = openOrFail('./gff/merged-same-query/%s.same.query.merged.gff.collapsed' % chr, 'w')
    outFile2.write('\t'.join(['#chr', 'merged', 'collapsed', 'chrom_start', 'chrom_end', 'strand',
                              'query', 'gene_id', 'num_of_exons', 'unique_name', 'query_start',
                              'query_end', 'query_len', 'frac', 'expect', 'score', 'exon_bound'])+'\n')

    print '... start working on matches of individual query...'

    # process the set of exons (aka matches or hits) resulting from
    # each query to create subsets. each subset contains hits related
    # to one or more putative(?) pgenes.

    # the exons are processed in the order in which they start on the
    # chromosome. there are four possibilities

    # 0) an exon that is for a different strand or is too far apart from
    # the previous exon. In this case, a new subset/pgene is started
    # containing this exon as its first exon, and the previous
    # subset/pgene is closed.

    # 1) Otherwise, if the exon has a significant query overlap with the
    # prior exon, it is taken to mean the protein(query) has two
    # distinct alignments with the chromosome, each corresponding to a
    # different pgene. In this case, we again start a new
    # subset/pgene. This case seems to include some ambiguities that
    # are not addressed (e.g., how should we determine which subset
    # gets related exons?)

    # 2) Otherwise, if the exon is close to the previous, then it is
    # considered a part of the previous and the two exons are merged
    # into one.

    # 3) Otherwise, this exon is a new exon of the current pgene and is
    # added to the subset. It becomes the exon to which the next is
    # compared.

    # Note: pgene membership and exon position is encoded in the
    # geneList structure.
    
    for query in gffHash.keys():
        # sort this query's matches by chrom_start.
        matches = gffHash[query]
        matches.sort(lambda x, y: x[ChromStart] - y[ChromStart])

        currExon = matches.pop(0)
        (cQStart, cQEnd) = (currExon[QueryStart], currExon[QueryEnd])

        geneList = [[]] # a list of pgenes, each of which is a list of exons.
        for newExon in matches:
            # recall we are processing the matches in order of their starting positions on the chromosome.
            distanceOnChrom = newExon[ChromStart] - currExon[ChromEnd] - 1
            (nQStart, nQEnd) = (newExon[QueryStart], newExon[QueryEnd])

            # Need to refine the intron cutoff [zz]
            if newExon[Strand] != currExon[Strand] or distanceOnChrom > 100000:
                # case 0
                currExon[MatchIdBug] = 0 # match id bug does not apply
                geneList[-1].append(currExon) # add current exon to current gene.

                geneList.append([]) # start a new gene.
                currExon = newExon 
                (cQStart, cQEnd) = (currExon[QueryStart], currExon[QueryEnd])
                continue

            # this is really an overlap check
            isSameGene = checkIfSameGene(cQStart, cQEnd, nQStart, nQEnd, newExon[Strand], newExon[QueryLen])

            if not isSameGene:
                # case 1
                currExon[MatchIdBug] = 0    # XXX need to change [zz]
                geneList[-1].append(currExon) # add current exon to current gene.

                geneList.append([]) # start a new gene
                currExon = newExon
                (cQStart, cQEnd) = (currExon[QueryStart], currExon[QueryEnd])
                continue
            else:
                # case 2 or 3.

                # bug? this appears to separation rather than overlap.
                if '+' ==  newExon[Strand]:
                    queryOverlap = nQStart - cQEnd - 1
                else:
                    queryOverlap = cQStart - nQEnd - 1

                ### take into account the missing amino acids in between [zz]
                if queryOverlap < 0: queryOverlap = 0
                queryGapNT = 3 * queryOverlap

                ### if the distance between matches are less than 60 nt, even considering the missing aa [zz]
                if  distanceOnChrom - queryGapNT <= queryGapCutoff:
                    # case 2

                    # is this a bug? --- we know the chroms starts are ordered, but not necessarily the ends.
                    currExon[ChromEnd] = newExon[ChromEnd]

                    # merge data from new into current exon, note new
                    # exon is *not* added to the current gene's list of exons.
                    if nQStart <= cQStart: cQStart = currExon[QueryStart] = nQStart
                    if nQEnd >= cQEnd: cQEnd = currExon[QueryEnd] = nQEnd
                    if float(newExon[Expect]) < float(currExon[Expect]): currExon[Expect] = newExon[Expect]

                    continue
                else:
                    # case 3
                    # bug? normally the '1' exon digit is suppressed in the match id, it wasn't in this case
                    # in zz's original code. we flag this situation by this mysterious '1'.
                    currExon[MatchIdBug] = 1     ### need to change [zz]
                    geneList[-1].append(currExon) # add current exon to current gene.

                    currExon = newExon
                    if cQStart > currExon[QueryStart]: cQStart = currExon[QueryStart]
                    if cQEnd < currExon[QueryEnd]: cQEnd = currExon[QueryEnd]

                    continue

        # don't forget the last exon --- do we need to verify that it is valid?
        currExon[MatchIdBug] = 0
        geneList[-1].append(currExon)

        # generate output for each of the resulting exons.
        for g in xrange(len(geneList)):
            exons = geneList[g]
            numExons = len(exons)
            for e in xrange(numExons):
                exon = exons[e]

            	_id = '%s===%d' % (exon[Query], g+1)
                # check match id bug flag.
                if exon[MatchIdBug] or e > 0: _id += '===%d' % (e+1)

                exon[KeyTag] = '%s===%s===%d' % (exon[Query], chr, exon[ChromStart]) # '===': a crude attempt to avoid ambiguity when working with weird chromosome/sequence names.
                outFile.write('\t'.join([exon[ChrId], 'merged', 'merged'] +
                                    [str(exon[f]) for f in [ChromStart, ChromEnd, Strand, Query]]
                                    + [_id, str(e+1), str(numExons), exon[KeyTag]] +
                                    [str(exon[f]) for f in [QueryStart, QueryEnd, QueryLen]] + [exon[Expect], exon[Score]]) + '\n')

        # print just one line of output for all the exons in a given gene.
        for g in xrange(len(geneList)):
            exons = geneList[g]
            numExons = len(exons)

            exon = exons[0]
            (lastCStart, lastCEnd) = (exon[ChromStart], exon[ChromEnd])
            (lastQStart, lastQEnd) = (exon[QueryStart], exon[QueryEnd])
            (lastExpect, lastScore)= (exon[Expect], float(exon[Score]))
            exonBound = '%d..%d' % (lastCStart, lastCEnd)

            for e in xrange(1, numExons):
                exon = exons[e]
                (cStart, cEnd) = (exon[ChromStart], exon[ChromEnd])
                (qStart, qEnd) = (exon[QueryStart], exon[QueryEnd])

                if cStart < lastCStart:	lastCStart = cStart
                if cEnd > lastCEnd:	lastCEnd = cEnd
                if qStart < lastQStart:	lastQStart = qStart
                if qEnd > lastQEnd:	lastQEnd = qEnd

                if float(exon[Expect]) < float(lastExpect): lastExpect = exon[Expect]
                lastScore += float(exon[Score])

                exonBound += ' %d..%d' % (cStart, cEnd)

            frac = float(lastQEnd - lastQStart + 1)/exon[QueryLen]

            outFile2.write('\t'.join([exon[ChrId], 'merged', 'merged',
                                      str(lastCStart), str(lastCEnd), exon[Strand], exon[Query],
                                      '%s===%d' % (exon[Query], g+1), str(numExons), exon[KeyTag],
                                      str(lastQStart), str(lastQEnd), str(exon[QueryLen]),
                                      '%.2f' % frac, lastExpect, printScore(lastScore), exonBound]) + '\n')

    for f in [outFile, outFile2]:
        f.close()
        os.system('sort -g -k 4 %s > %s.sorted.by.start' % (f.name, f.name))
