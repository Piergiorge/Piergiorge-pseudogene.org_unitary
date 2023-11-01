#!/usr/bin/env python

import cPickle, os, sys

from openOrFail import openOrFail
# this script and pfpo2 will need to be modified a bit to handle tfasty 3.4 output:
#   1) a few attributes (pg_version, pg_name) are multiply defined and so generate warnings
#   2) the bug zz mentions may have been fixed.
#   3) '-' strands appear to have coordinates reversed.
#   4) sy_ vs sx_, ty_ vs tx_ ???

from pfpo2 import *

#================================================================================
#  Infile:
#        ./step_2/chr$chr.fasta.align.log
#        ./gff/extended/chr$chr.extended.gff
#  Output:
#        ./step_3/chr$chr.fasta.raw.gff
#  parse the fasta alignment outputs, matches have been extended and re-aligned
#  with protein, write out new gff and fa files, sorted by the starting position 
#  on the matches on chromosome
#================================================================================ [zz]

# labels for a row of data, 'Chr', 'Identity', 'MatchSeq' and 'QuerySeq' are added here.
dataLabels = ['F1', 'F2', 'F3', 'ChromStart', 'ChromEnd', 'Strand',
              'Query', 'MatchId', 'ExonRank', 'NumOfExons', 'UniqueName',
              'QueryStart', 'QueryEnd', 'QueryLen', 'Expect', 'Score'] + ['Chr', 'Identity', 'MatchSeq', 'QuerySeq']
# define each label as a manifest constant with its position in the list as its value.
for lIndex in xrange(len(dataLabels)): exec '%s = %d' % (dataLabels[lIndex], lIndex)

#main
print '...running %s...' % sys.argv[0]

# get the chromosomes [zz]
if len(sys.argv) > 1:
    chrNames = sys.argv[1:]
else:
    chrNames = [str(i) for i in xrange(1,20)] + ['X']

try:    os.makedirs('./step_3')
except: pass

for chr in chrNames:
   print '...working on %s...' % chr
   # maps from a short, index-based identifier to the longer identifier used by pseudopipe.
   fx2seqId = cPickle.load(open('fastaOut/%s_fx2seqId.p'%chr))


   # data hash table [zz]
   gffHash = {}

   oldGFFFile = openOrFail('./gff/extended/%s.extended.gff' % chr, 'r')
   for line in oldGFFFile:
       line = line.strip()
       if '#' == line[0]: continue;
       f = line.split()
       (f[ChromStart], f[ChromEnd]) = (int(f[ChromStart]), int(f[ChromEnd]))
       gffHash[f[UniqueName]] = f + ['', '', '', '']

   oldGFFFile.close()

   infile = openOrFail('./fastaOut/%s.fasta.align.log' % chr, 'r')
   tFasta = parseTFastyOutput(infile)
   infile.close()
   
   for queryKey in tFasta:
       for (matchId, (proChroPos, header, summaryData, matches, trailer)) in tFasta[queryKey].iteritems():
           (match, mVals, qLab, qVals, qSeq, tLab, tVals, tSeq) = matches[0]
           (tStart, tStop) = (int(tVals['al_start']), int(tVals['al_stop']))

           # the summary line begins with ~45 quoted characters from
           # the input file's def line for this sequence. these
           # characters may include spaces. tfasty then adds four
           # fields: strand, opt, bits, e(1). we want strand, and
           # since we do not know into how many fields the beginning
           # part of summary line will be split, we must count from the
           # end of the field list of white-space separated fields.
           strand = '-+'['[f]' == summaryData[0][-5]]

           matchId = fx2seqId[matchId]
           gff = gffHash[matchId]
           ##============================================================================
           ## there is a bug in fasta log, al_stop - al_start +1 is always 3 nucleotides
	   ## short than the real aligned sequence.  If the match is on forward strand
           ## then match_al_start is unchanged match_al_stop ++3.  If the strand is reverse
           ## then the match_al_start -=2, match_al_stop +=1 
           ##============================================================================ [zz]

           # is this bug still present in fasta 3.4?
           if '+' == strand:
               (gff[ChromStart], gff[ChromEnd]) = (gff[ChromStart] + tStart - 1, gff[ChromStart] + tStop - 1)
           else:
               (gff[ChromStart], gff[ChromEnd]) = (gff[ChromStart] + tStop - 1, gff[ChromStart] + tStart - 1)
 
           gff[Chr] = chr
           gff[Query] = fx2seqId[queryKey].split('===',1)[0]
           gff[QueryStart] =  qVals['al_start']
           gff[QueryEnd] =  qVals['al_stop']
           gff[Strand] = strand

           gff[QueryLen] =  qVals['sq_len']
           gff[QuerySeq] = qSeq
           gff[MatchSeq] = tSeq

           gff[Score] = mVals['sy_score']
           gff[Identity] = mVals['sy_ident']
##           gff[Score] = mVals['sx_score']
##           gff[Identity] = mVals['sx_ident']

   newGFFFile = openOrFail('./step_3/%s.fasta.raw.gff' % chr, 'w')
   newFAFile = openOrFail('./step_3/%s.fasta.raw.fa' % chr, 'w')

   # write out new .gff and new .fa files [zz]
   newGFFFile.write('\t'.join(['#chr', 'extd', 'alned', 'chrom_start', 'chrom_end', 'strand', 'query',
			'match_id', 'exonRank', 'num_of_exons', 'unique_name', 'query_start',
			'query_end', 'query_len', 'expect', 'score', 'ident'])+'\n') 

   # sort the records, first reversly  by chrom_end and then by chrom_start [zz]
   # stable sort bug in original
   gffKeys = gffHash.keys()
   gffKeys.sort(lambda l, r: gffHash[r][ChromEnd] - gffHash[l][ChromEnd]) # reverse via l,r swap
   gffKeys.sort(lambda l, r: gffHash[l][ChromStart] - gffHash[r][ChromStart])

   #================================
   # remove redundant matches
   #================================ [zz]
   prevGffKey = gffKeys.pop(0)
   sortedMatches = []
   for k in gffKeys:
       gff = gffHash[k]
       prevGff = gffHash[prevGffKey]
       if gff[Strand] != prevGff[Strand] \
           or gff[Query] != prevGff[Query] \
           or gff[ChromStart] > prevGff[ChromEnd]:
           
           ##======================================================================
           ## the last_match and the new match have different strand, query et al.
           ##===================================================================== [zz]

           # can two distinct chains of related matches be intertwined?
           sortedMatches.append(prevGffKey)
           prevGffKey = k

       elif gff[Query] == prevGff[Query] \
           and gff[ChromStart] == prevGff[ChromStart] \
           and gff[ChromEnd] > prevGff[ChromEnd]:

           ##======================================================================
           # current match contains last match
           ##====================================================================== [zz]
           print '... %s contains %s, ignored ...' % (k, prevGffKey)

       elif gff[Query] == prevGff[Query] \
           and gff[ChromStart] >=  prevGff[ChromStart] \
           and gff[ChromEnd] <= prevGff[ChromEnd]:

           ##====================================================================
           ## last match contains current match
           #=====================================================================
           print '... %s contains %s, ignored ...' % (prevGffKey, k)

       elif gff[Query] == prevGff[Query] and gff[ChromStart] <= prevGff[ChromEnd]:
           ##====================================================================
           ## this match overlaps with last_match
           ##====================================================================
           print '...%s overlaps with %s, do nothing...' % (k, prevGffKey)
           sortedMatches.append(prevGffKey)
           prevGffKey = k

       else: assert 0, 'logic problem in remove redundant matches'

   sortedMatches.append(prevGffKey)

   for k in sortedMatches:
       gff = gffHash[k]
       (gff[ChromStart], gff[ChromEnd]) = (str(gff[ChromStart]), str(gff[ChromEnd]))

       newGFFFile.write('\t'.join([gff[Chr], 'extd', 'alned'] + \
                                  [gff[f] for f in [ChromStart, ChromEnd, Strand, Query, MatchId, ExonRank, NumOfExons]] + \
                                  [k] + \
                                  [gff[f] for f in [QueryStart, QueryEnd, QueryLen, Expect, Score, Identity]]) + '\n')

       newFAFile.write(' '.join(['>'+k, gff[Chr]] + \
                                [gff[l] for l in [Strand, ChromStart, ChromEnd, Query, MatchId, QueryStart, QueryEnd, QueryLen, Expect, Score]]) + '\n')
	
       newFAFile.write('%s\n%s\n' % (''.join([s[:-1] for s in gff[QuerySeq]]), ''.join([s[:-1] for s in gff[MatchSeq]])))

   newGFFFile.close()
   newFAFile.close()


