#!/usr/bin/env python

#======================================================================================
#  Inputs:
#        gff/remove-overlap-all-query/chr12.all.query.removed.gff.sorted.by.start
#  Outputs:
#        ./gff/extended/chr*.extended.fa
#        ./gff/extended/chr*.extended.gff
#
#  Inputs are merged BLAST out gff file, matches are extended on both sides to match     
#  the length of the query protein plus 30 (60)  nt buffer, write out the matches in fasta   
#  format also create the nre gff file with updates boundaries                        
######################################################################################### [zz]

import mmap, os, sys

from openOrFail import openOrFail

# labels for a row of data, 'Chr' is added here.
dataLabels = ['F1', 'F2', 'F3', 'ChromStart', 'ChromEnd', 'Strand',
              'Query', 'MatchId', 'ExonRank', 'NumOfExons', 'UniqueName',
              'QueryStart', 'QueryEnd', 'QueryLen', 'Expect', 'Score'] + ['Chr']
# define each label as a manifest constant with its position in the list as its value.
for lIndex in xrange(len(dataLabels)): exec '%s = %d' % (dataLabels[lIndex], lIndex)

#main
print 'fix hardcoded chrom lengths. extension pad (buffer) is set to 60 (20 aas?).'

print '...running %s...' % sys.argv[0]

envVars = [('ChromosomeFastaTemplate', 1)]
for ev in envVars:
    (n, req) = ev[:2]
    v = os.getenv(n)
    if req:
        if not v:
            sys.stderr.write('Must set environment variable %s.\n' % n)
            sys.exit(1)
    else:
        v = ev[2]
        sys.stderr.write('Using default value "%s" for environment variable %s.\n' % (v, n))

    exec '%s="%s"' % (n, v)
    
if len(sys.argv) > 1:
    chrNames = sys.argv[1:]
else:
    chrNames = [str(i) for i in xrange(1,20)] + ['X']

# output goes in this directory
try:    os.makedirs('./gff/extended/')
except: pass

for chr in chrNames:
    print '...working on %s...' % chr

    inFile = openOrFail('./gff/remove-overlap-all-query/%s.all.query.removed.gff.sorted.by.start' % chr, 'r')
    inFile.next() # skip header line

    gffHash = {}
    for line in inFile:
        f = line[:-1].split('\t')
        for l in [ChromStart, ChromEnd, ExonRank, NumOfExons, QueryStart, QueryEnd, QueryLen]: f[l] = int(f[l])
        gffHash[f[MatchId]] = f + [chr]

    inFile.close()

    # skip if no entry was read [zz]
    if not gffHash: continue

    # to avoid reading the entire chromosome file, we map it into
    # memory, and then seek to various locations of interest.
    fileName = ChromosomeFastaTemplate % chr
    try:
        chrFd = os.open(fileName, os.O_RDWR)
        chrLen = os.fstat(chrFd)[6]
        chrMap = mmap.mmap(chrFd, chrLen)
    except: assert 0, '%s doen\'t exist or cannot be mmap\'ed\n' % fileName
        
    chrMap.readline() # skip header line.
    firstLine = chrMap.tell()  # where first data line starts
    chrMap.readline() 
    secondLine = chrMap.tell() # where second data line starts
    lLen = secondLine - firstLine - 1 # - 1 because newline is 1 character under unix.
    # adjust len to be nt len.
    chrLen = chrLen - firstLine
    chrLen = chrLen - chrLen/(lLen+1)

    # sort records by order of start location on chromosome.
    sortedKeys = gffHash.keys()
    sortedKeys.sort(lambda l, r: gffHash[l][ChromStart] - gffHash[r][ChromStart])

    buffer = 60
    for id in sortedKeys:
        gff = gffHash[id]
        (query, numOfExons, exonRank) = (gff[Query], gff[NumOfExons], gff[ExonRank])

        geneId = id.split('===')[1] # id format: query.geneId[.exonId] 

        if '+' == gff[Strand]: # forward strand
            if 1 == exonRank:
                # extend left to query start
                leftExtend = gff[QueryStart] - 1
            else:
                # extend left to previous exon
                preId = '%s===%s===%d' % (query, geneId, exonRank - 1)
                leftExtend = gff[QueryStart] - gffHash[preId][QueryEnd] - 1
                if 0 > leftExtend: leftExtend = 0

            if 1 == numOfExons or exonRank == numOfExons:
                # extend right to query end (this seems asymmetric wrt to -1)
                rightExtend = gff[QueryLen] - gff[QueryEnd]
            else:
                # extend right to next exon
                nextId = '%s===%s===%d' % (query, geneId, exonRank + 1)
                rightExtend = gffHash[nextId][QueryStart] - gff[QueryEnd] - 1
                if 0 > rightExtend: rightExtend = 0
                
        else: # reverse strand, exon rank is in chromosome order, but query start/end switch roles.
            if 1 == exonRank:
                # extend left to query end (this seems asymmetric wrt to -1)
                leftExtend = gff[QueryLen] - gff[QueryEnd]
            else:
                # extend left to previous exon
                preId = '%s===%s===%d' % (query, geneId, exonRank - 1)
                leftExtend = gffHash[preId][QueryStart] - gff[QueryEnd] - 1
                if 0 > leftExtend: leftExtend = 0

            if 1 == numOfExons or exonRank == numOfExons:
                # last exon, extend right to query start
                rightExtend = gff[QueryStart]
                # bug?
                if 1 == numOfExons: rightExtend -= 1
            else:
                # extend right to next exon
                nextId = '%s===%s===%d' % (query, geneId, exonRank + 1)
                rightExtend = gff[QueryStart] - gffHash[nextId][QueryEnd] - 1
                if 0 > rightExtend: rightExtend = 0

        # "*3" => convert aa units to nt units.
        gff[ChromStart] -= leftExtend*3 + buffer
        gff[ChromEnd] += rightExtend*3 + buffer

        # keep them in bounds. (why the - 1 ?)
        if 1 > gff[ChromStart]: gff[ChromStart] = 1
        if gff[ChromEnd] >= chrLen: gff[ChromEnd] = chrLen - 1

    print '...writing out %s...' % chr

    
    outFAFile = openOrFail('./gff/extended/%s.extended.fa' % chr, 'w')

    outGFFFile = openOrFail('./gff/extended/%s.extended.gff' % chr, 'w')
    outGFFFile.write('\t'.join(['#chr', 'merged', 'extned', 'chrom_start', 'chrom_end', 'strand', 'query',
                                'match_id', 'exonRank', 'num_of_exons', 'unique_name',
                                'query_start', 'query_end', 'query_len', 'expect', 'score']) + '\n')
    
    for id in sortedKeys:
        gff = gffHash[id]

        outGFFFile.write('\t'.join([gff[Chr], 'merged', 'extd']
                                   + [str(gff[l]) for l in [ChromStart, ChromEnd, Strand, Query]]
                                   + [id]
                                   + [str(gff[l]) for l in [ExonRank, NumOfExons, UniqueName, QueryStart, QueryEnd, QueryLen, Expect, Score]]) + '\n')

        outFAFile.write('  '.join(['>'+gff[UniqueName], gff[Chr], 'merged', 'extd']
                                + [str(gff[l]) for l in [ChromStart, ChromEnd, Strand, Query]]
                                + [id]
                                + [str(gff[l]) for l in [ExonRank, NumOfExons, QueryStart, QueryEnd, QueryLen, Expect, Score]]) + '\n')

        # compute location in memory of chromosome start (relative to
        # start of first data line)
        (lines, offset) = divmod(gff[ChromStart]-1, lLen)
        try:    chrMap.seek(firstLine + lines*(lLen+1) + offset)
        except: print 'pos %s: seeking to %d in a file of %d bytes' %  (gff[ChromStart], firstLine + lines*(lLen+1) + offset, chrMap.size())
            
        seqLen = gff[ChromEnd] - gff[ChromStart] + 1

        # extract sequence from chromosome
        seq = ''
        while len(seq) < seqLen: seq += chrMap.readline()[:-1]
        seq = seq[:seqLen] # trim excess
        
        # now write it out (while lLen might make more sense, 50
        # matches zz's output)
        for x in xrange(0, seqLen, 50): outFAFile.write(seq[x:x+50] + '\n')

    chrMap.close()
    outFAFile.close()
    outGFFFile.close()
