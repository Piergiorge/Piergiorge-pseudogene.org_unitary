#!/usr/bin/env python

import os, sys

cwd = os.path.dirname(os.path.realpath(sys.argv[0])) + os.sep

#chrNames = ' '.join([str(i) for i in range(1,20)] + ['X'])
chrNames = ' '.join(sys.argv[1:])

stages = ['filterEnsemblGene.py',
          'removeOverlapSameQuery.py',
          'mergeMatchesSameQuery.py',
          'removeOverlapAllQuery.py',
          'extendMatches.py',
          'fastaAlign.py',
          'parseFastaAlignment.py',
          'cleanupFasta1.py',
          'cleanupFasta2.py',
          'cleanupFasta3.py',
          'checkPolya.py',
          'makeExons.py',
          'makeGenes2.py']

for s in stages:
    print 'running ' + s
    if os.system(cwd + s + ' ' + chrNames): break
else:
    print 'all processing completed.'
    sys.exit(0)

print 'failed during %s stage.' % s
