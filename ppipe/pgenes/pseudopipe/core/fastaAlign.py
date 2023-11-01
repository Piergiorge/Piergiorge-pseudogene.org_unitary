#!/usr/bin/env python

#===============================================================
# Input:
#   ./gff/extended/chr*.extended.fa 
# Output:
#   ./step_2/chr*.fasta.align.log 
#   ./step_2/chr*.fasta.visual.log 
#
#  This script use "tfasta" to align more rigorously the extended
#  chromosomal DNA sequence with query protein
#================================================================ [zz]

import cPickle, os, sys, time

sew = sys.stderr.write

# break up a fasta file into individual (def line, sequence) pairs,
# but don't reformat in anyway. build a dictionary of thses sequences
# indexed by the first element of the def line.
def fasta2dict(fFile):
    (defLine, seq) = ('', '')
    d = {} 
    for l in fFile:
        if '>' == l[0]: # a new def line
            if defLine: # record previous def line and seq
                d[defLine[1:].split()[0]] = (defLine, seq)
            (defLine, seq) = (l, '')
        else: # more of the sequence
            seq += l

    if defLine: d[defLine[1:].split()[0]] = (defLine, seq)

    return d
                
from openOrFail import openOrFail

#main
print '...running %s...' % sys.argv[0]

envVars = [('ProteinQueryFile', 1), ('FastaProgram', 1)]
for ev in envVars:
    (n, req) = ev[:2]
    v = os.getenv(n)
    if req:
        if not v:
            sew('Must set environment variable %s.\n' % n)
            sys.exit(1)
    else:
        v = ev[2]
        sew('Using default value "%s" for environment variable %s.\n' % (v, n))

    exec '%s="%s"' % (n, v)
    
# get the chromosomes [zz]
if len(sys.argv) > 1:
    chrNames = sys.argv[1:]
else:
    chrNames = [str(i) for i in xrange(1,20)] + ['X']

# create a dictionary of query proteins.
queryFile = openOrFail(ProteinQueryFile, 'r')
queryDict = fasta2dict(queryFile)
queryFile.close()

# keep track of blast hit queries that are not in the protein dictionary.
npQueries = []

myPid = os.getpid()
blastMatchTmp = '/tmp/blastmatch_%d.dna.fa'%myPid
queryTmp = '/tmp/query_%d.protein.fa'%myPid

try:    os.makedirs('./fastaOut')
except: pass

# our working identifiers can grow long enough that fasta truncates
# them. we use fx2seqId to map from a short, index-based identifier
# for the chromosome sequence that is used when invoking fasta to
# the longer working identifier.
fastaCount, totalTime = 0, 0.0
for chr in chrNames:
    print '...working on %s...' % chr
    fx2seqId = {}
    
    # set up fasta log files
    fastaLog       = 'fastaOut/%s.fasta.align.log' % chr
    fastaVisualLog = 'fastaOut/%s.fasta.visual.log' % chr
    try: # to delete old log files
        os.unlink(fastaLog)
        os.unlink(fastaVisualLog)
    except: pass

    # create a dictionary of chromosome DNA matchs sequences
    blastMatchFile = openOrFail('./gff/extended/%s.extended.fa' % chr, 'r')
    blastMatchDict = fasta2dict(blastMatchFile)
    blastMatchFile.close()

    for matchId in blastMatchDict.keys():
        # write (extended) chromosome sequence for this blast hit to temp file.
        newId = 'seq_%d'%fastaCount
        bmData = blastMatchDict[matchId]
        fx2seqId[newId] = matchId
        open(blastMatchTmp, 'w').write(('>%s\n'%newId)+''.join(bmData[1:]))

        # extract query from matchId (matchId: queryr_chr_start)
        try:    query = matchId.split('===',1)[0]
        except: assert 0, 'failed to parse %s\n' % matchId

        if queryDict.has_key(query): # write protein query sequence for this hit to temp file.
            open(queryTmp, 'w').write(''.join(queryDict[query]))
        else: # add to list of those not in the protein dictionary and move on.
            npQueries.append(query)
            continue

        # run fasta
        startTime = time.time()
        os.system('''
                  %s -m 10 %s %s >> %s
#                  %s %s %s >> %s
                  rm  -f %s %s
                  ''' % (FastaProgram, queryTmp, blastMatchTmp, fastaLog,
                         FastaProgram, queryTmp, blastMatchTmp, fastaVisualLog,
                         blastMatchTmp, queryTmp))
        totalTime += time.time() - startTime
        fastaCount += 1

    cPickle.dump(fx2seqId, open('fastaOut/%s_fx2seqId.p'%chr, 'w'), -1)

# report queries not found in protein dictionary
print '%d fasta invocations; %f seconds' % (fastaCount, totalTime)

if npQueries:
    sew('unmatched queries:\n')
    sew('\n'.join(npQueries)+'\n')
