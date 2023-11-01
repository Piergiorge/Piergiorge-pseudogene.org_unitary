#!/usr/bin/env python

import operator, os, re, sys

from openOrFail import *

## labels for a row of data
## first set are labels for data from gff input, second set are added fields.
dataLabels = ['F1', 'F2', 'F3', 'ChromStart', 'ChromEnd', 'Strand',
              'Query', 'MatchId', 'GeneID', 'ExonRank', 'NumOfExons', 'UniqueName',
              'QueryStart', 'QueryEnd', 'QueryLen', 'Expect', 'Ident', 'Polya'] + \
             ['Chr', 'Band', 'Frac', 'NumInserts', 'NumDeletes',  'NumShifts',
              'NumStops', 'Disable', 'QuerySeq', 'MatchSeq']

## define each label as a manifest constant with its position in the list as its value.
for lIndex in xrange(len(dataLabels)): exec '%s = %d' % (dataLabels[lIndex], lIndex)

(BandName, BandBegin, BandEnd) = (0, 1, 2)
bandHash = {}

def getBand(chr, matchId, chromStart, chromEnd):

  #==============================================
  # obtain the cytogenic band of the exon/gene 
  #==============================================[zz]

    print '%s\t%d\t%d' % (chr, chromStart, chromEnd)

    chrBandHash = bandHash[chr]
    for x in xrange(len(chrBandHash)):
        (n, b, e) = chrBandHash[x]
        if chromStart >= b and chromEnd <= (e + 10): return n

        # potential bounds bug?
        if chromStart >= b and chromEnd < chrBandHash[x + 1][BandEnd]:
            overlap = e - chromStart + 1

            # rounding arithmetic?
            if overlap >= 0.5 * (chromEnd - chromStart + 1): return  n

            return chrBandHash[x + 1][BandName]

        if chromStart <= e:
            sys.stderr.write('Failed to obtain cytogenic info for %s\n' % matchId)
            sys.exit(1)
    else:
        print 'potential bug, but returning null band just like zz'
        return None

def calIdentity(s1, s2):

    #===========================================================
    #  return indentical number of residues bwtn two sequences
    #============================================================ [zz]
    if not s1 or not s2:
        sys.stderr.write('string empty\n')
        sys.exit(1)
        
    if len(s1) != len(s2):
        sys.stderr.write('the strings have different length\n')
        sys.exit(1)

    return reduce(operator.add, [s1[i] == s2[i] for i in xrange(len(s1))], 0)


alignmentMarksRE = re.compile(r'[/\\Xx*]')
alignmentMarksNoStarRE = re.compile(r'[/\\Xx]')
alignmentMarksNoXRE = re.compile(r'[/\\*]')
manyDashesRE = re.compile(r'-{21,}')

def checkDisable(querySeq, matchSeq):

  #===============================================================================================
  # the idea is to check the n-term and c-term ends of the alignment to exclude
  # the suspicious low-confidence residues, and then check whether the high-confidence
  # sequence conatins any disablement. The high confidence alignment starts with 4 residues
  # with over 2 identity, and also with 5 out of 10 identical residues. Long stretch of insertions 
  # in the match_seq are also considered as low-confidence alignment
  #=============================================================================================== [zz]

    leng = len(matchSeq)

    if not alignmentMarksRE.search(matchSeq): return '0'

    lengTest = 5
    lengTest2 = 10

    # some sort of sliding window would probably be a more efficient approach
    # for the next two loops.
    posNterm = None
    for i in xrange(leng - lengTest2 + 1):
        ident  = calIdentity(querySeq[i:i+lengTest], matchSeq[i:i+lengTest])
        ident2  = calIdentity(querySeq[i:i+lengTest2], matchSeq[i:i+lengTest2])

        if ident >= 2 and ident2 >= 5:
            posNterm = i + 3
            break

    posCterm = None
    for i in xrange(leng - 1, lengTest2 - 2, -1):
        ident  = calIdentity(querySeq[i-lengTest+1:i+1], matchSeq[i-lengTest+1:i+1])
        ident2  = calIdentity(querySeq[i-lengTest2+1:i+1], matchSeq[i-lengTest2+1:i+1])

        # original had ident >= 3
        if ident >= 2 and ident2 >= 5:
            posCterm = i - 3
            break

    if not posNterm or not posCterm: return 'd'

    #
    # I assume perl does not side effect invocee arguments via the following assignments...
    # bounds problem?
    querySeq = querySeq[posNterm:posCterm+1]
    matchSeq = matchSeq[posNterm:posCterm+1]

    posArray = [m.span() for m in manyDashesRE.finditer(querySeq)]

    if posArray:
        start = 0
        tmpSeq = ''
        for (gb, ge) in posArray:
            tmpSeq += matchSeq[start:gb]
            start = ge
        tmpSeq += matchSeq[start:]
        matchSeq = tmpSeq

# tfasty 3.4 change
#    if alignmentMarksNoStarRE.search(matchSeq): return 'D'
    if alignmentMarksNoXRE.search(matchSeq): return 'D'

    return 'd'

#main
if len(sys.argv) > 1:
    chrNames = sys.argv[1:]
else:
    chrNames = [str(i) for i in xrange(1,20)] + ['X']

print '...running %s ...' % sys.argv[0]

## #====================================
## #  read in the cytogeneic band info
## #==================================== [zz]

## bandFile = openOrFail('/home/carriero/bioInformatics/ZhaoleiZhang/ensembl/mouse-030321/mysql/mus_musculus_core_11_3/karyotype.txt.sorted.table', 'r')
## for line in bandFile:
##     f = line.split('\t')
##     try:    bandHash[f[1]].append((f[4], int(f[2]), int(f[3])))
##     except: bandHash[f[1]] = [(f[4], int(f[2]), int(f[3]))]
        
## bandFile.close()

gffHash = {}

try:    os.makedirs('./pexons')
except: pass

for chr in chrNames:
    print '...working on %s...' % chr

    oldGFFFile = openOrFail('./polya/%s.polya.1000.gff' % chr, 'r')
    oldFAFile = openOrFail('./step_3/%s.fasta.clean.3.fa' % chr, 'r')
    newGFFFile = openOrFail('./pexons/%s.pexon.gff' % chr, 'w')
    newFAFile = openOrFail('./pexons/%s.pexon.fa' % chr, 'w');

    oldGFFFile.next()
    sortedKeys = []
    for line in oldGFFFile:
        f = line[:-1].split('\t')
        for l in [ChromStart, ChromEnd, QueryStart, QueryEnd, QueryLen]: f[l] =  int(f[l])
        frac = '%.2f' % ((float(f[QueryEnd] - f[QueryStart] + 1))/f[QueryLen])
#        band = getBand(chr, f[MatchId], f[ChromStart], f[ChromEnd])
	band = 'NoBandData'
        sortedKeys.append(f[MatchId])
        gffHash[f[MatchId]] = f + [chr, chr+band, frac, -1, -1, -1, -1, '', '', '']

    oldGFFFile.close()

    print '...finish reading karyo band...'

    for line in oldFAFile:

        #======================================================================
        #  retrieve number of inserts, delete, frameshits, et al and sequences 
        #====================================================================== [zz]

        try:
            matchId = line.split()[0][1:]
        except:
            sys.stderr.write(('failed to parse following line in %s:\n' % oldFAFile.name) + +line)
            sys.exit(1)

        gffMatch = gffHash[matchId]

        querySeq = oldFAFile.next()[:-1]
        gffMatch[QuerySeq] = querySeq

        matchSeq = oldFAFile.next()[:-1]
        gffMatch[MatchSeq] = matchSeq

        #=========================
        # compute some parameters
        #========================== [zz]
        gffMatch[NumInserts] = querySeq.count('-')
        gffMatch[NumDeletes] = matchSeq.count('-')
        gffMatch[NumShifts] = matchSeq.count('\\') + matchSeq.count('/')
        gffMatch[NumStops] = matchSeq.count('X') + matchSeq.count('x')
# fasta 3.4 change
        gffMatch[NumStops] += matchSeq.count('*')

        #==========  use the calculate identity from FASTA  ========== [zz]
        # meant to comment out non(?) side-effecting cll to cal_identity?

        #========== whether the sequence contains disablement  ========== [zz]
        gffMatch[Disable] = checkDisable(gffMatch[QuerySeq], gffMatch[MatchSeq])
        print matchId + ' ' + gffMatch[Disable]
        
    oldFAFile.close()
    
    ##==========  write out new .gff and new .fa files  ========== [zz]
    newGFFFile.write('\t'.join(['#chr', 'pexon', 'pexon', 'band', 'chrom_start', 'chrom_end', 'strand',
                               'query', 'match_id', 'gene_id', 'exon_rank', 'num_of_exons', 'unique_name',
                               'query_start', 'query_end', 'query_len', 'frac', 'numInserts', 'numDeletes',
                               'numShifts', 'numStops', 'expect', 'ident', 'polya', 'disable']) + '\n')

    for id in sortedKeys:
        gff = gffHash[id]

        newGFFFile.write('\t'.join([gff[Chr], 'exon', 'exon']
                                   + [str(gff[l]) for l in [Band, ChromStart, ChromEnd, Strand, Query]]
                                   + [id]
                                   + [str(gff[l]) for l in [GeneID, ExonRank, NumOfExons, UniqueName,
                                                            QueryStart, QueryEnd, QueryLen, Frac,
                                                            NumInserts, NumDeletes, NumShifts, NumStops,
                                                            Expect, Ident, Polya, Disable]]) + '\n')

        newFAFile.write(' '.join(['>' + id, gff[Chr]]
                                 + [str(gff[l]) for l in [Strand, ChromStart, ChromEnd]]
                                 + [id]
                                 + [str(gff[l]) for l in [GeneID, Query, UniqueName,
                                                          QueryStart, QueryEnd, QueryLen, Frac,
                                                          NumInserts, NumDeletes, NumShifts, NumStops,
                                                          Expect, Ident, Polya, Disable]]) + '\n')
        newFAFile.write(gff[QuerySeq] + '\n')
        newFAFile.write(gff[MatchSeq] + '\n')

    newGFFFile.close()
    newFAFile.close()
