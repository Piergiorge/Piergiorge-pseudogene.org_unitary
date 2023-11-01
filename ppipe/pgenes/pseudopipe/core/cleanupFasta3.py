#!/usr/bin/env python

import re, sys

from openOrFail import openOrFail

## labels for a row of data
## first set are labels for data from gff input. second are extra fields.
dataLabels = ['F1', 'F2', 'F3', 'ChromStart', 'ChromEnd', 'Strand',
              'Query', 'GeneName', 'UniqueName', 'QueryStart', 'QueryEnd', 'QueryLen',
              'Expect', 'Score', 'Ident'] + ['Chr', 'ExonRank','NumOfExons']

## define each label as a manifest constant with its position in the list as its value.
for lIndex in xrange(len(dataLabels)): exec '%s = %d' % (dataLabels[lIndex], lIndex)

#==========================================================================================
#  This script does the 3rd step of the "clean-up" procedure 
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

#main 
print '...running %s...' % sys.argv[0]

if len(sys.argv) > 1:
    chrNames = sys.argv[1:]
else:
    chrNames = [str(i) for i in xrange(1,20)] + ['X']

nameTagRE = re.compile(r'^>(\S+)\s+')

for chr in chrNames:
    print '...working on %s...' % chr

    oldGffFile = openOrFail('./step_3/%s.fasta.clean.2.gff' % chr, 'r')
    oldFaFile = openOrFail('./step_3/%s.fasta.clean.2.fa' % chr, 'r')
    newGffFile = openOrFail('./step_3/%s.fasta.clean.3.gff' % chr, 'w')

    newGffFile.write('\t'.join([
        '#chr',        'clean3',    'clean3',       'chrom_start',
        'chrom_end',   'strand',    'query',        'match_id',
        'gene_id',     'exon_rank', 'num_of_exons', 'unique_name',
        'query_start', 'query_end', 'query_len',    'expect',
        'ident']) + '\n'),

    newFaFile = openOrFail('./step_3/%s.fasta.clean.3.fa' % chr, 'w')

    #==============================================================
    # read through the gff file to get the number of exons for each 
    # gene 
    #=============================================================
    #  Data structures: 
    #  %fa_hash: contains match and query sequences
    #  %gff_hash: contains info for each exon/match
    #  %gene_structure: contains number of exons for each gene_name
    #  @sorted_names: array containging exon/match names.
    #============================================================== [zz]
    lastChromStart = 0

    geneStructure = {}
    gffHash = {}
    sortedNames = []

    oldGffFile.next()
    for l in oldGffFile:
        f = l[:-1].split('\t') + [chr, 0, 0 ] # initial values for chr, exon rank and number

        ## overly cautious, this is not really necessary [zz]
        chromStart = int(f[ChromStart])
        if chromStart < lastChromStart:
            sys.stderr.write('%d < $%d\nfailed in parsing ' + l)
            sys.exit(1)

        lastChromStart = chromStart

        # determine exon_rank [zz]
        if geneStructure.has_key(f[GeneName]):
            geneStructure[f[GeneName]] += 1
        else:
            geneStructure[f[GeneName]] = 1
        f[ExonRank] = geneStructure[f[GeneName]]

        gffHash[f[UniqueName]] = f
        sortedNames.append(f[UniqueName])

    oldGffFile.close()

    ## read in query and match sequences [zz]
    faHash = {}
    for l in oldFaFile:
        mo = nameTagRE.match(l)
        if not mo:
            sys.stderr.write('Failed to parse following line:\n' + l)
            sys.exit(1)

        name = mo.groups()[0]
        faHash[name] = (oldFaFile.next(), oldFaFile.next())

    oldFaFile.close()
    
    #=============================================
    #  write out the new  gff  and fa
    #================================================= [zz]
    for match in sortedNames:
        geneName = gffHash[match][GeneName]
        gff = gffHash[match]
        numOfExons = geneStructure[geneName]
        gff[NumOfExons] = numOfExons

        matchId = geneName
        if numOfExons > 1: matchId += '.%d' % gff[ExonRank]

        newGffFile.write('\t'.join([gff[Chr], 'clean3', 'clean3',
                                    gff[ChromStart], gff[ChromEnd], gff[Strand], gff[Query],
                                    '%s===%s' % (chr, matchId), '%s===%s' % (chr, geneName),
                                    str(gff[ExonRank]), str(gff[NumOfExons]), match,
                                    gff[QueryStart], gff[QueryEnd], gff[QueryLen],
                                    gff[Expect], gff[Ident]]) + '\n')

        newFaFile.write('  '.join(['>%s===%s' % (chr, matchId),
                                   gff[ChromStart], gff[ChromEnd], gff[Strand], gff[Query],
                                   str(gff[ExonRank]), str(gff[NumOfExons]),'%s===%s' % (chr, match),
                                   gff[QueryStart], gff[QueryEnd], gff[QueryLen],
                                   gff[Expect], gff[Ident]]) + '\n')
        newFaFile.write(''.join(faHash[match]))

    newGffFile.close()
    newFaFile.close()

