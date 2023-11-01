#!/usr/bin/env /usr/bin/python2

import sys

def splitFasta(fileName, seqsPerSplit, splitName):
    f = open(fileName)
    seqCount = 0
    fileCount = 0
    of = open(splitName % fileCount, 'w')
    for l in f:
        if '>' == l[0]:
            if seqsPerSplit == seqCount:
                seqCount = 0
                of.close()
                fileCount += 1
                of = open(splitName % fileCount, 'w')
            seqCount += 1
        of.write(l)
    of.close()
    f.close()

splitFasta(sys.argv[1], int(sys.argv[2]), sys.argv[3])
