#!/usr/bin/env python

#==========================================================================================
#  This script does the 2nd step of the "clean-up" procedure 
#  The inputs should have been sorted by chrom_start already
#========================================================================================== [zz]

##===========================================================================================
# The raw parsed outputs from fasta realignment are "cleaned up" in 3 steps
# Step_1: 
# Script: cleanup_fasta_1.pl
# Input: ./step_3/chr*.fasta.raw.gff (fa)
# Output: ./step_3/chr*.fasta.clean.1.gff (fa) 
# Redundant matches are removed, and those matches containing long gaps (>60nt) are split
# dubious residues at both ends of alignments are also removed according to criterion
# write out new gff and fa files, sorted by the starting position on the matches on chromosome
# Step_2:
# Script: cleanup_fasta_2.pl
# Input: ./step_3/chr*.fasta.clean.1.gff (fa)
# Output: ./step_3/chr*.fasta.clean.2.gff (fa)
# Overlapping matches are removed accoring to a procedure similar to "remove_overlap.pl"
# Step_3:
# Script: cleanup_fasta_3.pl
# Input: ./step_3/chr*.fasta.clean.2.gff (fa)
# Output: ./step_3/chr*.fasta.clean.3.gff (fa)
# Since some matches are removed and split, so the exon and gene assignment are messed up.
# This script assign new exon and gene ids to matches
##============================================================================================ [zz]

import os, re, sys

from openOrFail import openOrFail

## labels for a row of data
## first set are labels for data from gff input. second are extra fields.
dataLabels = ['F1', 'F2', 'F3', 'ChromStart', 'ChromEnd', 'Strand',
              'Query', 'GeneName', 'Match', 'QueryStart', 'QueryEnd', 'QueryLen',
              'Expect', 'Score', 'Ident'] + ['QuerySeq','MatchSeq']

## define each label as a manifest constant with its position in the list as its value.
for lIndex in xrange(len(dataLabels)): exec '%s = %d' % (dataLabels[lIndex], lIndex)

def fSignum(f):
    if f < 0.0: return -1
    if f > 0.0: return 1
    return 0

def removeByExpect(gff, outFile):

  #============================================================
  #  simply-mindedly removing matches according to score value
  #============================================================= [zz]

    #==============================================================================
    # we sort the genes by score because they are more useful to separate multiple
    # exon genes or genes after merging
    #============================================================================== [zz]

    # python2.3+ dependency: this works only with a stable sort (potential bug in zz's
    # original)
    candidates = gff.keys()
    candidates.sort(lambda l, r: -fSignum(float(gff[l][Score]) - float(gff[r][Score])))
    candidates.sort(lambda l, r: fSignum(float(gff[l][Expect]) - float(gff[r][Expect])))

    pickedHits = []

    for candidate in candidates:
        (cStart, cEnd) = (gff[candidate][ChromStart], gff[candidate][ChromEnd])
        #==============================================================
        #  decide whether the hit overlaps with an already picked hit
        #  default is no overlap
        #=============================================================== [zz]
        for picked in pickedHits:
            (pStart, pEnd) = (picked[ChromStart], picked[ChromEnd])
            
            if    (cEnd <= pEnd and cStart >= pStart) \
               or (cEnd >= pEnd and cStart <= pStart) \
               or (cEnd <= pEnd and cEnd >= (pStart + 30 - 1)) \
               or (cEnd >= pEnd and cStart <= (pEnd - 30 - 1)): break
        else:
            pickedHits.append(gff[candidate])

    for picked in pickedHits:
        outFile.write('\t'.join([str(picked[l]) for l in [F1, F2, F3, ChromStart, ChromEnd, Strand, Query, GeneName,
                                                     Match, QueryStart, QueryEnd, QueryLen, Expect, Score, Ident]]) + '\n')
#main

print '...running %s...' % sys.argv[0]

if len(sys.argv) > 1:
    chrNames = sys.argv[1:]
else:
    chrNames = [str(i) for i in xrange(1,20)] + ['X']

for chr in chrNames:
    print '...working on %s...' % chr

    oldGffFile = openOrFail('./step_3/%s.fasta.clean.1.gff' % chr, 'r')
    newGffFile = openOrFail('./step_3/%s.fasta.clean.2.gff.tmp' % chr, 'w')

    newGffFile.write('\t'.join([
        '#chr',        'clean2',      'clean2',    'chrom_start',
        'chrom_end',   'strand',      'query',     'gene_name',
        'unique_name', 'query_start', 'query_end', 'querylen',
        'expect',      'score',       'ident']) + '\n')

    #=================================
    #  read in the first line  
    #================================= [zz]
    hashGff =  {}
    oldGffFile.next()

    f = oldGffFile.next()[:-1].split('\t')
    # later comparisons will depend on these fields being integers.
    (f[ChromStart], f[ChromEnd]) = (int(f[ChromStart]), int(f[ChromEnd]))
    f[F2] = f[F3] = 'clean2'
    hashGff[f[Match]] = f
    lastChromEnd = int(f[ChromEnd])
    
    #================================
    #   read rest of the lines  
    #================================= [zz]
    for l in oldGffFile:
        f = l[:-1].split('\t')
        (f[ChromStart], f[ChromEnd]) = (int(f[ChromStart]), int(f[ChromEnd]))

        # cutoff in eval ? [zz]
        # bug in original --- first record is not checked for cutoff...
        if float(f[Expect]) >= 1e-6: continue

        #=================================================================
        # separate cluster of overlapped matches and process them later
        #================================================================= [zz]
        if f[ChromStart] <= lastChromEnd:
            f[F2] = f[F3] = 'clean2'
            hashGff[f[Match]] = f
            if lastChromEnd < f[ChromEnd]: lastChromEnd = f[ChromEnd]
        else:
            # remove overlap, [zz]
            removeByExpect(hashGff, newGffFile)

            #delete the content of the old hash tables
            # put int new data into the hash tables [zz]

            # intential that we don't set F2, F3 to clean2?
            hashGff = { f[Match] : f }
            lastChromEnd = int(f[ChromEnd])

    removeByExpect(hashGff, newGffFile)

    oldGffFile.close()
    newGffFile.close()

    #  sort the new gff by chrom_start [zz]
    os.system('sort  -g -k 4 %s > %s' % (newGffFile.name, newGffFile.name.replace('.tmp', '')))

matchTagRE = re.compile(r'^>(\S+)\s+')

for chr in chrNames:

    #===========================================================
    #  write our NEW_FA, but first need to read in the NEW_GFF
    #============================================================ [zz]

    gffHash = {}
    sortedKeys = []

    newGffFile = openOrFail('./step_3/%s.fasta.clean.2.gff' % chr, 'r')
    newGffFile.next()
    for l in newGffFile:
        f = l[:-1].split('\t')
        gffHash[f[Match]] = f + ['', ''] # place holders for query and match seqs.
        sortedKeys.append(f[Match])
    newGffFile.close()

    oldFaFile = openOrFail('./step_3/%s.fasta.clean.1.fa' % chr, 'r')
    for l in oldFaFile:
        mo = matchTagRE.match(l)
        if not mo:
            sys.stderr.write('Failed to parse following line:\n' + l)
            sys.exit(1)

        match = mo.groups()[0]
        if gffHash.has_key(match):
            gffHash[match][QuerySeq] = oldFaFile.next()
            gffHash[match][MatchSeq] = oldFaFile.next()
        else:
            oldFaFile.next()
            oldFaFile.next()

    oldFaFile.close()

    newFaFile = openOrFail('./step_3/%s.fasta.clean.2.fa' % chr, 'w')
    for key in sortedKeys:
        entry = gffHash[key]
        newFaFile.write('  '.join(['>'+key, 'clean2', 'clean2'] + \
                                  [entry[l] for l in [ChromStart, ChromEnd, Strand, Query, GeneName,
                                                      QueryStart, QueryEnd, QueryLen,
                                                      Expect, Score, Ident]]) + '\n')
        newFaFile.write(entry[QuerySeq])
        newFaFile.write(entry[MatchSeq])
