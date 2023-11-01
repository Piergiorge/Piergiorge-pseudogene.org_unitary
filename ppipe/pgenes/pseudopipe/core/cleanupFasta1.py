#!/usr/bin/env python

##===========================================================================================
# The raw parsed outputs from fasta realignment are "cleaned up" in 3 steps
# Step_1: (cleanup_fasta_1.pl, this script)
# 	Input: ./step_3/chr*.fasta.raw.gff (fa)
# 	Output: ./step_3/chr*.fasta.clean.1.gff (fa) 
# Redundant matches are removed, and those matches containing long
# gaps (>60nt) are split.  Dubious residues at both ends of alignments
# are also removed according to criterion. Write out new gff and fa
# files, sorted by the starting position of the matches on chromosome.
#
# Step_2: (cleanup_fasta_2.pl)
# 	Input: ./step_3/chr*.fasta.clean.1.gff (fa)
# 	Output: ./step_3/chr*.fasta.clean.2.gff (fa)
# Overlapping matches are removed according to a procedure similar to
# "remove_overlap.pl"
#
# Step_3: (cleanup_fasta_3.pl)
# 	Input: ./step_3/chr*.fasta.clean.2.gff (fa)
# 	Output: ./step_3/chr*.fasta.clean.3.gff (fa)
# Since some matches are removed or(?) split, the exon and gene
# assignment are messed up.  This script assigns new exon and gene ids
# to matches.
##============================================================================================ [zz]

import operator, re, sys

from openOrFail import openOrFail

# labels for a row of data
# first set are labels for data from extended gff input
dataLabels = ['F1', 'F2', 'F3', 'ChromStart', 'ChromEnd', 'Strand',
              'Query', 'MatchId', 'ExonRank', 'MaxExonRank', 'Match',
              'QueryStart', 'QueryEnd', 'QueryLen', 'Expect', 'Score', 'Ident'] 

# define each label as a manifest constant with its position in the list as its value.
for lIndex in xrange(len(dataLabels)): exec '%s = %d' % (dataLabels[lIndex], lIndex)

# data labels for records output by this script
outDataLabels = ['OutChr', 'OutChromStart', 'OutChromEnd', 'OutStrand',
                 'OutQuery', 'OutGeneName', 'OutKey',
                 'OutQueryStart', 'OutQueryEnd', 'OutQueryLen', 'OutExpect', 'OutScore', 'OutIdent',
                 'OutQSeq', 'OutMSeq']
# define each label as a manifest constant with its position in the list as its value.
for lIndex in xrange(len(outDataLabels)): exec '%s = %d' % (outDataLabels[lIndex], lIndex)

def checkIdent(s1, s2):
# return the number of identical characters in s1 and s2 (see
# calIdentity in makeExons.)

    assert s1 and s2, 'strings empty'
    assert len(s1) == len(s2), 'the strings have different length'

    return reduce(operator.add, [s1[i] == s2[i] for i in xrange(len(s1))], 0)

def _cleanup2(qSeq, mSeq, len1, minIdent1, len2, minIdent2):
# look for a position x such that qSeq and mSeq have at least
# minIdent1 characters in common in their subsequences beginning at x
# of length len1, and at least minIdent2 characters in common in the
# subsequences beginning at x of length len2. it is assumed that len1
# < len2.
#
# return (n, c), the leftmost and rightmost such positions.
#
# we could optimized using a sliding window, but that seems like
# overkill.

    leng = len(qSeq) # is this a bug --- what if m <<< q?

    (posNterm, posCterm) = (None, None)

    for i in xrange(leng - len2): # bug in original? '<= leng' would go past end
        if checkIdent(qSeq[i:i+len1], mSeq[i:i+len1]) >= minIdent1:
            if checkIdent(qSeq[i:i+len2], mSeq[i:i+len2]) >= minIdent2:
                posNterm = i
                break
        
    for i in xrange(leng, len2-1, -1):
        if checkIdent(qSeq[i-len1:i], mSeq[i-len1:i]) >= minIdent1:
            if checkIdent(qSeq[i-len2:i], mSeq[i-len2:i]) >= minIdent2:
                posCterm = i
                break
    
    return (posNterm, posCterm)

def cleanupAlignment(qSeq, mSeq, minGapWidth):
# zz says:
#
# remove long gaps (>20 aa, 60 nt) in the alignment, by breaking into
# fragments. human introns are in general longer than 60 nts. we check
# only the query sequence for gaps (because these gaps indicate elided
# material? --- what would a long gap in the chromosome signify?)
#
# remove dubious residues at ends of the fragments, producing a clean
# (non-dubious?) alignment in which both ends start with a string of
# 10 residues with identity >= 50% and a string of 5 residues with
# identity >= 20%.
#
# return 0-based positions of the cleaned up fragments relative to the
# original. it is assumed the sequences don't contain long gaps (>20
# aa, 60 nt).

    gapRe = re.compile('(-{%d,})' % minGapWidth) # matches minGapWidth or more consecutive '-'s
    start = 0
    pos = [] # will hold intervals defining new (gapless) fragments.
    for m in gapRe.finditer(qSeq):
        (pos1, pos2) = m.span()

        # skip shifts (but what if multiple shifts?, edge conditions?)
        if qSeq[pos1-1] in '\\/': pos1 -= 1
        if qSeq[pos2] in '\\/': pos2 += 1

        pos.append((start, pos1)) # bug: if first gap begins at 0, we end up with (0, 0) as an interval
        start = pos2

    pos.append((start, len(qSeq)))

    # trim dubious regions from the ends of the resulting fragments
    pos2 = []
    for (start, stop) in pos:
        # zz doc bug: min idents are inclusive limits.
        (posN2, posC2) = _cleanup2(qSeq[start:stop], mSeq[start:stop], 5, 2, 10, 5)

        if posN2 == None or posC2 == None:
            print posN2, posC2
            print qSeq
            print mSeq
            continue
        # zz bug: posN2 or posC2 can be undefined, in which case ... ?
        if posN2 < posC2: pos2.append((posN2+start, posC2+start))

    return pos2

def countNts(aaSeq):
# count number of nucleotides in an amino acid sequence by countng
# gaps and shift and then computing nts.

    numGaps = aaSeq.count('-')
    numFwds = aaSeq.count('\\')
    numBcks = aaSeq.count('/')

    return (len(aaSeq) - numGaps - numFwds - numBcks)*3 + numFwds - numBcks

def checkIdentOverlap(s1, s2):
# return the number of identical characters in s1 and s2, as well as
# the number of positions where one or the other has an aa (and so the
# sequences overlapped).

    assert s1 and s2, 'strings empty'
    assert len(s1) == len(s2), 'the strings have different length'

    (ident, overlap) = (0, 0)
    for i in xrange(len(s1)):
        (c1, c2) = (s1[i], s2[i])

        if c1 == c2: ident += 1

        # alphabetic implies aa (what about 'X'?)
        if c1.isalpha() or c2.isalpha(): overlap += 1

    return (ident, overlap)

#main
print 'document lengths and cutoffs for cleanup, gap minimum, need to document identifiers match/matchId/geneName/...'
print '...running %s...' % sys.argv[0]

if len(sys.argv) > 1:
    chrNames = sys.argv[1:]
else:
    chrNames = [str(i) for i in xrange(1,20)] + ['X']

FaMatchRe = re.compile(r'^>(\S+)\s') # matches first token on def line.
for chr in chrNames:
    print '...working on %s...' % chr

    # build a dictionary of fasta match info indexed by match identifier. 
    oldFaFile = openOrFail('./step_3/%s.fasta.raw.fa' % chr, 'r')
    fastaHash = {}
    for l in oldFaFile:
        mo = FaMatchRe.match(l)
        assert mo, 'failed to parse ' + l[:-1]
        fastaHash[mo.group(1)] = (oldFaFile.next()[:-1], oldFaFile.next()[:-1])
    oldFaFile.close()

    oldGffFile = openOrFail('./step_3/%s.fasta.raw.gff' % chr, 'r')

    newHash = {} # collection of clean fragments/exons
    for l in oldGffFile:
        if '#' == l[0]: continue

        f = l[:-1].split('\t')
        (qStart, cStart) = (int(f[QueryStart]), int(f[ChromStart]))

        # retrieve fasta match info for this record.
        (querySeq, matchSeq) = fastaHash[f[Match]]

        # get coordinates of the "clean" fragments within the query
        # sequence using a minimum gap of 20.
        newCoors = cleanupAlignment(querySeq, matchSeq, 20)

        frag = 0
        for (start, stop) in newCoors:
            # adjust start and end positions to account for gaps (what about shifts, only in chromosome?)?
            newQStart = qStart + (start - querySeq[:start].count('-'))
            newQEnd = qStart + (stop - querySeq[:stop-1].count('-')) - 1 # stop-1 reproduces a potential bug in zz's code.

            if '+' == f[Strand]:
                newCStart = cStart + countNts(matchSeq[:start])
                newCEnd   = cStart + countNts(matchSeq[:stop]) - 1
            else:
                newCStart = cStart + countNts(matchSeq[stop:])
                newCEnd   = cStart + countNts(matchSeq[start:]) - 1

            # if we have more than one frag (exon?) add a frag counter to the identifier.
            key = f[Match] ;
            if len(newCoors) > 1: key += '.%d' % frag

            # following script uses gene_name without the exon_id part. [zz]
            geneName = f[MatchId]
            try: geneName = '==='.join(geneName.split('===')[:2])
            except: pass

            qSeq = querySeq[start:stop]
            mSeq = matchSeq[start:stop]
            
            (_ident, _overlap) = checkIdentOverlap(qSeq, mSeq)

            _ident2 = '%.3f' % (float(_ident)/_overlap)

            newHash[key] = [chr, newCStart, newCEnd, f[Strand], f[Query], geneName, key,
                            newQStart, newQEnd, f[QueryLen], f[Expect], f[Score], _ident2,
                            qSeq, mSeq]

            frag += 1

    oldGffFile.close()

    # sort the records in order of chrom start. if more than one
    # record for a given value of chrom start, sort by chrom end
    # (greatest to least).
    #
    # perl may not be using a stable sort, so this may be a bug in the
    # original python from 2.3 uses stable sort.
    keys = newHash.keys()
    keys.sort(lambda l, r: newHash[r][OutChromEnd] - newHash[l][OutChromEnd]) # l,r reversal is intentional.
    keys.sort(lambda l, r: newHash[l][OutChromStart] - newHash[r][OutChromStart]) #

    newGffFile = openOrFail('./step_3/%s.fasta.clean.1.gff' % chr, 'w')
    newFaFile = openOrFail('./step_3/%s.fasta.clean.1.fa' % chr, 'w')

    newGffFile.write('\t'.join(['#chr', 'clean1', 'clean1', 'chrom_start', 'chrom_end', 'strand', 'query',
                        'gene_name', 'unique_name',  'query_start',
                        'query_end', 'query_len', 'expect', 'score', 'ident']) + '\n')

    for key in keys:
        r = newHash[key]
        newGffFile.write('\t'.join([r[OutChr], 'clean1', 'clean1'] + \
                                   [str(r[x]) for x in xrange(OutChromStart,OutQSeq)]) + '\n')
        newFaFile.write('>'+key+'  clean1  clean1  '+ \
                        '  '.join([str(r[x]) for x in xrange(OutChromStart,OutQSeq)]) + '\n')
        newFaFile.write(r[OutQSeq] + '\n')
        newFaFile.write(r[OutMSeq] + '\n')

    newGffFile.close()
    newFaFile.close()
