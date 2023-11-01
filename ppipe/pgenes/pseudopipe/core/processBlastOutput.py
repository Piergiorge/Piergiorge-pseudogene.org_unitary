#!/usr/bin/env python

# merges three zz processing steps:
#
# 1) reformat.sh
# 
# 2) extract-by-chr.pl
#
# 3) parse-rename-blastout.pl
#
# Notes:
#
# 1) (wip: Renaming now is from dddl_human to [P|Q].... form)
#
# 2) Final sorting is done through external sort command.

##=============================================================
##  read in the TAB-deli BLAST output (from Nick C)
##  reformat it to the "old human" format so that
##  chr_start is always < $chr_end. Also introduce strand.
##  Then split into individual chromosomes
##=============================================================== [zz]

import os, re, sys

from fastaSeqFile import fastaSeqFile
from openOrFail import openOrFail

# first arg is protein query file. second is blast output file

# read in query file and build a dictionary mapping a query to its aa length.
queries = {}
for s in fastaSeqFile(sys.argv[1]):
    # the fasta iterator returns the array of lines for the next
    # sequence.  the first line is the def line, the rest (which
    # include new line characters) hold the sequence.

    # def line looks like '>QueryID safsad fafs ...'
    k = s[0][1:].split()[0]
    if queries.has_key(k):
        print 'duplicate query key %s, ignoring: %s' % (k, s[1])

    queries[k] = str(reduce(lambda x, y: x+y, [len(l)-1 for l in s[1:]], 0))

# now process blast output
ofRe = re.compile(sys.argv[2])

of = [ '%s/%s'%(d, f) for d in sys.argv[3:] for f in os.listdir(d) if ofRe.match(f) ]

## labels for a row of data in blast output
dataLabels = ['Query', 'Chr', 'Ident', 'Overlap', 'Deletes', 'Inserts', 'QueryStart',
              'QueryEnd', 'ChrStart', 'ChrEnd', 'EVal', 'Score']

## define each label as a manifest constant with its position in the list as its value.
for lIndex in xrange(len(dataLabels)): exec '%s = %d' % (dataLabels[lIndex], lIndex)

pm = {'+': 'P', '-': 'M'}

chrFiles = {}
for bf in of:
    for l in open(bf):
        f = l[:-1].split('\t')
        if 12 != len(f):
            print 'warning, ignoring corrupt line: ' + str(f)
            continue

        # throw out weak hits.
        if float(f[EVal]) >= 1e-4: continue

        # assume + strand, but test.
        strand = '+'
        if int(f[ChrEnd]) < int(f[ChrStart]):
            strand = '-'
            (f[ChrStart], f[ChrEnd]) = (f[ChrEnd], f[ChrStart])

        # use following to partition output by (chr X strand).
        k = (f[Chr], pm[strand])
        if not chrFiles.has_key(k):
            chrFiles[k] = openOrFail('./%s_%s_blastHits' % k, 'w')

            ##     k = f[Chr]
            ##     if not chrFiles.has_key(k):
            ##         chrFiles[k] = openOrFail('./chr%s.blastoutPy' % k, 'w')

            # generate a header line.
            chrFiles[k].write('\t'.join(['#query', 'chr', 'ident', 'overlap', 'del', 'insert', 'strand',
                                        'query_start', 'query_end', 'query_len', 'chr_start', 'chr_end', 'eval', 'score'])+'\n')

        chrFiles[k].write('\t'.join(f[:QueryStart]+[strand]+f[QueryStart:ChrStart]
                                    +[queries[f[0]]]+f[ChrStart:])+'\n')

for k in chrFiles:
    chrFiles[k].close()

    # why add 3 to these indices?: we add one to shift to one based
    # indexing used by sort, we add two to account for the fact we
    # insert strand and length columns before ChrStart, ChrEnd, and EVal
    cStartX = ChrStart+3
    EValX = EVal+3
    cEndX = ChrEnd+3

    # sort by chr start, then eval, then chr end (reversed, so longest
    # otherwise equal hit comes first).
    os.system('sort -k %d,%dn -k %d,%dg -k %d,%dnr %s > %s.sorted'
              % (cStartX, cStartX, EValX, EValX, cEndX, cEndX, chrFiles[k].name, chrFiles[k].name))
