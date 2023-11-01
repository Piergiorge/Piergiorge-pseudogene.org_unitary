#!/usr/bin/env python

##=========================================================================
## read in the BLAST parsed ouput for each chr, filter with the masks
## created from EnSembl annotation, remove those hits that overlap
## significantly with the masks [zz]

# Apparently this code (especially locateMatch) assumes that we will
# be handed matches in sorted order (and all matches have the same
# strandedness).

import os, sys

from openOrFail import openOrFail

## The index is passed in and out to exploit the ordering of mask locations
## That is, the next thing we look for will be no earlier than where we have
## already looked. A mask, roughly, corresponds to an exon.
def locateMatch(chrStart, chrEnd, index, maskStart, maskEnd):
    ## find first mask starting beyond chrStart and then back up.
    while maskStart[index] <= chrStart: index += 1
    index -= 1
    (start, startsBetween) = (index, chrStart > maskEnd[index])

    ## find first mask starting beyond chrEnd and then back up.
    while maskStart[index] <= chrEnd: index += 1
    index  -= 1
    (end, endsBetween) = (index, chrEnd > maskEnd[index])

    return (start, startsBetween, end, endsBetween)

# main
print 'need to document overlap parameter (30) and dependency on mask array files.'

# get info for gene masks
envVars = [('BlastoutSortedTemplate', 1), ('ExonMaskFields', 0, '0 1'), ('ExonMaskTemplate', 1)]
for ev in envVars:
    (n, req) = ev[:2]
    v = os.getenv(n)
    if not v:
        if req:
            sys.stderr.write('Must set environment variable %s.\n' % n)
            sys.exit(1)
        else:
            v = ev[2]
            sys.stderr.write('Using default value "%s" for environment variable %s.\n' % (v, n))

    exec '%s="%s"' % (n, v)
    
maskFields = [int(e)for e in ExonMaskFields.split()]
print 'mask fields', maskFields

if len(sys.argv) > 1:
    chrNames = sys.argv[1:]
else:
    chrNames = [str(i) for i in xrange(1,20)] + ['X']

for chr in chrNames:
    # read in the ensembl gene mask file [zz]
    maskFile = openOrFail(ExonMaskTemplate % chr, 'r')

    # arrays to store the start and end of genes (after merging) [zz]
    # initialize to contain two 'virtual' exons at 0 and 1000000000000
    (maskStart, maskEnd) = ([0, 1000000000000], [0, 1000000000000])

    print '... start reading in mask file for %s ...' % chr
    for l in maskFile:
        if '#' == l[0]: continue
        f = l[:-1].split('\t')
        (start, end) = [int(f[x]) for x in maskFields]
        assert start <= end, 'fatal error: mask_start %d > mask_end %d.\n' % (start, end)
        maskStart.insert(-1, start)
        maskEnd.insert(-1, end)

    maskFile.close()
    print '... finish reading in mask file for %s ...' % chr

    inFile = openOrFail(BlastoutSortedTemplate % chr, 'r')

    try:    os.makedirs('./gff/ensembl-gene-masked')
    except: pass
    outFile = openOrFail('./gff/ensembl-gene-masked/%s.ens.gene.masked' % chr, 'w')
    
    # copy the header
    header = inFile.next()
    outFile.write(header)

    # figure out which columns hold data of interest. if naming scheme
    # changes, this will break --- arguably a "good thing".
    colNames = header.split()
    (ChrStart, ChrEnd, Strand) = (colNames.index('chr_start'), colNames.index('chr_end'), colNames.index('strand'))

    ## index of the mask array, declared here so we can keep track
    ## of the position in the mask arrays that we have come to [zz]
    index = 0
    for l in inFile:
        f = l.split()
        (chrStart, chrEnd) = (int(f[ChrStart]), int(f[ChrEnd]))

        # finds the indices of the exons including or immediately to
        # the left of positions chrStart and chrEnd. The 'between'
        # values indicate which condition obtains.
        (start, startsBetween, end, endsBetween) = locateMatch(chrStart, chrEnd, index, maskStart, maskEnd)

        # next locateMatch should start searching beginning with this
        # index.
        index = start

        if startsBetween:
            # the BLAST hit's start lies between two exons with
            # indices 'start' and 'start+1'
            if endsBetween and start == end:
                # the hit ends between the same two exons that stradle
                # its start, so the hit is contained in one non-coding
                # region (no overlap).
                outFile.write(l)

            elif not endsBetween and (start + 1) == end:
                # the hit ends in the right exon of the two that
                # straddle its start. thus there is overlap with just
                # one gene.  keep the hit if the overlap is less than
                # the cutoff 30nt.
                if chrEnd - maskStart[end] + 1 < 30: outFile.write(l)

            elif not endsBetween and (start + 1) < end:
                # the hit ends in an exon at least one exon further to
                # the right, so the hit overlaps at least one complete
                # exon. discard the hit.
                pass

            elif endsBetween and start != end:
                # the hit ends between two exons, and the left of
                # these is not the same as the exon to the left of the
                # its start. thus, it must overlap at least one whole
                # exon. discard the hit.
                pass

            else: assert False, 'match starts between, but then what? (line: "%s")' % l[:-1]
        else:
            # the BLAST hit starts in the exon with index 'start'.
            if not endsBetween and start == end:
                # the hit starts and ends in the same exon, discard [zz]
                pass
            
            elif not endsBetween and start < end:
                # the hit ends in an exon and that exon is to the
                # right of exon in which the hit starts, so there is
                # overlap with portions of at least two
                # exons. (Shouldn't overlap rules apply here?)
                pass
            
            elif endsBetween and start == end:
                # the hit ends between two exons and the left exon of
                # these two is the exon in which the hit starts, thus
                # there is some overlap.  keep if the overlap is less
                # than the cutoff.
                if maskEnd[end] - chrStart + 1 < 30: outFile.write(l)

            elif endsBetween and start < end:
                # the hit ends between two exons and the left exon of
                # these two lies to the right of the exon in which the
                # hit starts, thus there is overlap with at least one
                # while exon. discard the hit.
                pass

            else: assert False, 'match does not start between, but then what? (line: "%s")' % l[:-1]
    outFile.close()
    inFile.close()


