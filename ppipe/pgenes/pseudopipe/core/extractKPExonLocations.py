#!/usr/bin/env python

# we assume we are in a directory containing the mysql files:
# exon.txt.table
# exon_transcript.txt.table
# translation.txt.table
# translation_stable_id.txt.table
# seq_region.txt.table
#
# so fetch from ensembl's mysql directory exon.txt.table.gz exon_transcript.txt.table.gz translation.txt.table.gz translation_stable_id.txt.table.gz seq_region.txt.table.gz
import sys

# ensembl files not longer end with .table, but if you need to work with an older set, you can specify that suffix here.
#old TableSuffix = '.table'
TableSuffix = ''

# first get a list of known peptipes (file name given in argv 1)
peps = [ l.split()[0][1:] for l in open(sys.argv[1]) if '>' == l[0] ]

# build a table that maps pep (aka 'translation') stable ids to internal ids.
tlsid2tlid = {}
bb = [ tlsid2tlid.setdefault(r[1], []).append(r[0]) for l in open('translation_stable_id.txt'+TableSuffix) for r in [l[:-1].split('\t')] ]
probs = [ k for k, v in tlsid2tlid.iteritems() if len(v) > 1 ]
if probs: print 'translation stable to translation internal is one to many!: ' + str(probs)

# build a table that maps translation ids to transcript ids.
tlid2trid = {}
bb = [ tlid2trid.__setitem__(r[0], r[1]) for l in open('translation.txt'+TableSuffix) for r in [l[:-1].split('\t')] ]

# compose these to build a table of known pep transcripts
kptrids = {}
bb = [ kptrids.__setitem__(tlid2trid[tlsid2tlid[p][0]], 1) for p in peps ]

# build a table of exons that are part of a known pep transcript.
# note we don't care about one to many here.
kpex = {}
bb = [ kpex.__setitem__(r[0], 1) for l in open('exon_transcript.txt'+TableSuffix) for r in [l[:-1].split('\t')] if kptrids.has_key(r[1]) ]

# build a list of exon location data for exons in the known pep exon list.
kplocs = [ r+[l] for l in open('exon.txt'+TableSuffix) for r in [l[:-1].split('\t')] if kpex.has_key(r[0]) ]

## # collate by chromosome (well, sequence region id and strand)
## kplocsBySRS = {}
## bb = [ kplocsBySRS.setdefault((r[1], r[4]), []).append(r) for r in kplocs ]

## # sort
## def cmpExLoc(l, r):
##     bl, br, el, er = int(l[2]), int(r[2]), int(l[3]), int(r[3])
##     if bl < br: return -1
##     if bl > br: return 1
##     if el < er: return -1
##     if el > er: return 1
##     return 0

## bb = [ v.sort(cmpExLoc) for v in kplocsBySRS.itervalues() ]

## # dump locations in sorted order by chr and strand
## # build a map from seq region internal id to name  for "coord system" 1 (chromsome-ish) entries in the seq_region table.
## srid2name = {}
## bb = [ srid2name.__setitem__(r[0], r[1]) for l in open('seq_region.txt'+TableSuffix) for r in [l[:-1].split('\t')] ]

## strand = { '1' : 'P', '-1' : 'M' }

## bb = [ open('chr%s_%s_exLocs'%(srid2name[k[0]], strand[k[1]]), 'w').write(''.join([r[-1] for r in kplocsBySRS[k]])) for k in kplocsBySRS ]

# collate by chromosome (well, sequence region id)
kplocsBySRS = {}
bb = [ kplocsBySRS.setdefault(r[1], []).append(r) for r in kplocs ]

# sort
def cmpExLoc(l, r):
    bl, br, el, er = int(l[2]), int(r[2]), int(l[3]), int(r[3])
    if bl < br: return -1
    if bl > br: return 1
    if el < er: return -1
    if el > er: return 1
    return 0

bb = [ v.sort(cmpExLoc) for v in kplocsBySRS.itervalues() ]

# dump locations in sorted order by chr and strand
# build a map from seq region internal id to name  for "coord system" 1 (chromsome-ish) entries in the seq_region table.
srid2name = {}
bb = [ srid2name.__setitem__(r[0], r[1]) for l in open('seq_region.txt'+TableSuffix) for r in [l[:-1].split('\t')] ]

bb = [ open('chr%s_exLocs'%srid2name[k], 'w').write(''.join([r[-1] for r in kplocsBySRS[k]])) for k in kplocsBySRS ]
