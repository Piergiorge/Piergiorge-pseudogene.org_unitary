#!/usr/bin/env python

##===========================================================
# for each match, print out the 3' downstream nt sequence
# and detect poly-A signal and poly-A tail 
##============================================================ [zz]

import mmap, os, re, string, sys

from openOrFail import openOrFail

POLYAWINDOW = 1000    # length of the nt searched for polya tail [zz]

## labels for a row of gff data plus 'Chr'.
dataLabels = ['F1', 'F2', 'F3', 'ChromStart', 'ChromEnd', 'Strand',
              'Query', 'MatchId', 'GeneId', 'ExonRank', 'NumOfExons', 'UniqueName',
              'QueryStart', 'QueryEnd', 'QueryLen', 'Expect', 'Ident'] +['Chr', 'Line']
## define each label as a manifest constant with its position in the list as its value.
for lIndex in xrange(len(dataLabels)): exec '%s = %d' % (dataLabels[lIndex], lIndex)

signalRE = re.compile(r'AATAAA|ATTAAA|CATAAA|AATATA')

def checkPolya(seq):
  #==========================================
  #  determine the polya class: 1, 2 or 3
  #=========================================== [zz]

    # upcase.
    seq = seq.upper()
   
    maxCount  = 0
    maxPos = '' # if never defined, will print as null (matches original output)
    maxStr = ''

    # zz's code use <= which translate to + 1,
    # but he probably meant <, and the + 1 should be ditched
    # this is O(n^2), but could easily be made O(n)
    for i in xrange(len(seq) - 50 + 1):
        frag = seq[i:i+50]
        count = frag.count('A')
        if count > maxCount: (maxCount, maxPos, maxStr) = (count, i, frag)

    signalPos = '' # if never defined, will print as null (matches original output)
    if maxCount < 30:
        kind = 0
    else:
        mo = signalRE.search(seq[maxPos-1::-1])
        if mo:
            signalPos = mo.span()[0]
            if signalPos < 50: kind = 1
            elif signalPos <= 100: kind = 2
            else: kind = 3
        else:
            kind = 3

    return (kind, maxCount, maxPos, maxStr, signalPos)

# bug in original? zz's seq manipulation module returns something
# for negative lengths (i.e., stop earlier than start). this
# code "reproduces" this bug.
def seqSlice(chrMap, start, stop, firstLine, llen):
    # extract seq for interval [start, stop] from a memory mapped file
    # whose data begins at position firstline, for which each data line contains
    # llen symbols.

    # try to do something reasonable when stop is before start
    if stop < start: (start, stop) = (stop, start)

    if not ((1 <= start <= chrMap.size()) and (1 <= stop <= chrMap.size())):
        print >>sys.stderr, 'Sequence slice (%d, %d) not contained within sequence (1, %d).\n'%(start, stop, chrMap.size())

    # convert to 0-based
    start -= 1
    stop -= 1
    seqLen = (stop - start) + 1
    if seqLen <= 1: print >>sys.stderr, 'Short slice', start, stop, seqLen
    (lines, offset) = divmod(start, llen)
    chrMap.seek(firstLine + lines*(llen+1) + offset)
    seq = ''
    while len(seq) < seqLen: seq += chrMap.readline()[:-1]
    return seq[:seqLen]

#main
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

print '...running %s...' % sys.argv[0]

complementTrans = string.maketrans('acgtACGT', 'tgcaTGCA')

try:    os.makedirs('./polya')
except: pass

for chr in chrNames:
    print '...working on %s...' % chr

    oldGffFile = openOrFail('./step_3/%s.fasta.clean.3.gff' % chr, 'r')
    newFaFile  = openOrFail('./polya/%s.polya.%d.fa' % (chr, POLYAWINDOW), 'w')
    newGffFile = openOrFail('./polya/%s.polya.%d.gff' % (chr, POLYAWINDOW), 'w')

    newGffFile.write('\t'.join(['#chr',        'polya',     'polya',
        'chrom_start', 'chrom_end',   'strand',    'query',
        'match_id',    'gene_id',     'exon_rank', 'num_of_exons',
        'unique_name', 'query_start', 'query_end', 'query_len',
        'expect',      'ident',       'polya']) + '\n')

    #================================================================
    # create a Sequence object and read in chromosome sequence
    # here we use the non-masked chromosomal sequence since the polyA
    # could have been masked out
    #================================================================ [zz]

    # possible bug -- the above comment seems to be at odds with the name of the file.
    fileName = ChromosomeFastaTemplate % chr
    # these files are quite large. use memory mapping to access them.
    try:
        chrFd = os.open(fileName, os.O_RDWR)
        chrLen = os.fstat(chrFd)[6]
        chrMap = mmap.mmap(chrFd, chrLen)
        getByte = chrMap.read_byte
    except:
        sys.stderr.write('%s doen\'t exist or cannot be mmap\'ed\n' % fileName)
        sys.exit(1)
        
    # find the line length (skip the def line)
    chrMap.readline() 
    firstLine = chrMap.tell()
    chrMap.readline()
    secondLine = chrMap.tell()
    llen = secondLine - firstLine - 1 # because newline is 1 character under unix.
    # adjust len to be nt len.
    chrLen = chrLen - firstLine
    chrLen = chrLen - (chrLen + llen)/(llen+1)

    gffHash = {}      # hash table contains the match info [zz]
    sortedIds = []    # array containing the match names, sorted by chrom_start [zz]

    ##===========================================
    ##  the gff file is supposed to be formated,
    ##  use $last_chrom_start to make sure 
    ##============================================ [zz]
    lastChromStart = 0
    for l in oldGffFile:
        if '#' == l[0]: continue

        f = l[:-1].split('\t') + [chr, l]

        ## to make sure the matches are pre-sorted by chrom-start [zz]
        chromStart = int(f[ChromStart])
        if chromStart < lastChromStart: 
            sys.stderr.write('Fatal error, the matches are not sorted by chrom_start:\n' + l)
            sys.exit(1)
        lastChromStart = chromStart
        
        gffHash[f[MatchId]] = f
        ## use the $match_id as hash key [zz]
        sortedIds.append(f[MatchId])

    # bug in original?: this closed a file that didn't exist.
    oldGffFile.close()

    ## determine the poly-A type, see paper  ## [zz]
    numIds = len(sortedIds)
    # bug in original? this loop exceeded bounds of sortedId
    for i in xrange(numIds):
        match = sortedIds[i]
        gff = gffHash[match]

        #XXXX [zz]
        print  chr + '\t' + match

        if '+' == gff[Strand]:
            checkStart = min(int(gff[ChromEnd]) + 1, chrLen)
            checkStop = min(checkStart + POLYAWINDOW - 1, chrLen)

            ## if the next match is less than 1000 bp away [zz]
            if (i + 1) < numIds:
                limit = int(gffHash[sortedIds[i+1]][ChromStart])
                if limit <= checkStop: checkStop = max(limit - 1, 1)
                
            #njc2 # check whether out of bound [zz]
            #njc2 if checkStop >= chrLen: checkStop = chrLen - 1
        
            seq = seqSlice(chrMap, checkStart, checkStop, firstLine, llen)
        else:
            checkStop  = max(int(gff[ChromStart]) - 1, 1)
            checkStart = max(checkStop - POLYAWINDOW + 1, 1)

            ## if the next match is less than 1000 bp away [zz]
            if i:
                limit = int(gffHash[sortedIds[i-1]][ChromEnd])
                if limit >= checkStart: checkStart = min(limit + 1, chrLen)

            #njc2 # check whether out of bound [zz]
            #njc2 if checkStart < 1: checkStart = 1

            seq = seqSlice(chrMap, checkStart, checkStop, firstLine, llen)
            seq = seq[::-1]

        seq = seq.translate(complementTrans)

        (kind, maxCount, maxPos, maxStr, signalPos) = checkPolya(seq)

        newFaFile.write('\t'.join(['>' + match, gff[ChromStart], gff[ChromEnd], gff[Strand],
                                   gff[QueryStart], gff[QueryEnd], gff[QueryLen]] +
                                   [str(n) for n in [checkStart, checkStop, checkStop - checkStart + 1,
                                   maxCount, maxPos, maxStr, signalPos, kind]]) + '\n')
        newFaFile.write(seq + '\n')

        newGffFile.write(gff[Line][:-1].replace('clean3', 'polya') + '\t' + str(kind) + '\n')

    newFaFile.close()
    newGffFile.close()


