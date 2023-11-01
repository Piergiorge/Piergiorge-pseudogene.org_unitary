#!/usr/bin/env python

import os, re, sys

from openOrFail import openOrFail

## labels for a row of data
## first set are labels for data from gff input, second set are added fields.
dataLabels = ['F1', 'F2', 'F3', 'Band', 'ChromStart', 'ChromEnd', 'Strand',
              'Query', 'MatchId', 'GeneId', 'ExonRank', 'NumOfExons', 'UniqueName',
              'QueryStart', 'QueryEnd', 'QueryLen', 'Frac',
              'NumOfIns', 'NumOfDels', 'NumOfShifts', 'NumOfStops',
              'Expect', 'Ident', 'Polya', 'Disable'] + \
             ['Chr', 'QuerySeq', 'MatchSeq']

## the records that will be output contain even more fields.
geneDataLabels = dataLabels + ['ExonBound', 'ExonLen', 'IntronBound', 'IntronLen',
                               'ClassOld', 'ClassNew', 'Comment', 'Comment2', 'ShortID']
## define each label as a manifest constant with its position in the list as its value.
for lIndex in xrange(len(geneDataLabels)): exec '%s = %d' % (geneDataLabels[lIndex], lIndex)

#==========================================================
# collapse pseudo-exons into single gene, and output info
#   The difference btween this script and make-genes.pl is: 
#  need to get the relative protein information such
#  as ribosomal protein, Numt, OR et al
#  These shouldn't be classified as processed pseudogenes
#========================================================== [zz]

# scope issue here
specialHash = {}

## def getBand(chr, geneId, chromStart, chromEnd, refBands):
## ##===============================================
## ## obtain the cytogenic band of the exon/gene
##   #================================================ [zz]

##     # bug in original? zz's version has <= where < is probably what should
##     # be used. more puzzling, he refs the next band too.
##     numBands = len(refBands)
##     for n in xrange(numBands):
##         (name, beg, end) = refBands[n]
##         if chromStart >= beg and chromEnd <= (end + 10):
##             band = name
##             break
##         elif chromStart >= beg and (n+1 < numBands and chromEnd < refBands[n+1][BandEnd]):
##             overlap = end - chromStart + 1;
##             # intended rounding here?
##             if overlap >= 0.5 * (chromEnd - chromStart + 1): band = name
##             else: band = refBands[n+1][BandName]
##             break
##         elif chromStart > end: continue
##         else:
##             sys.stderr.write('%s: %s  Failed to obtain cytogenic info for %s \n' % (chromStart, chromEnd, geneId))
##             sys.exit(1)
##     return chr + band

# note 'class' is a keyword for Python, so we use kind in this routine
def getClass(expCutoff, identCutoff, numOfExons, frac, expect, ident,
             polya, disable, query, queryStart, queryEnd, queryLen):


  #===============================================================
  # decied whether a p-gene is processed, everything is considered 
  # unpreocessed unless proved otherwise. use 70% as a cutoff as 
  # threshold. More conservative way need to consider identity
  # Pseudogenes containing one gap could be a repeat, need to include
  #  as pseudogene in the future
  #=============================================================== [zz]

    ## classes are FP, PSSD, FRAG, GENE-SINGLE GENE-MULT, and DUP [zz]

    #=============================================================
    # use exp=1e-10 and ident= 40% as cutoff for false positive
    # also remove the ORs and Numts and merge RPs to PSSD1 
    #============================================================= [zz]
    if expect > expCutoff or ident < identCutoff: return ('FP', 'FP')

    if 1 == numOfExons:
        if frac > 0.95 and disable == '0' and ident > 0.95: kind = 'GENE-SINGLE'
        elif frac >= 0.7: kind = 'PSSD' + disable + polya
        else: kind = 'FRAG'
    else:
        if (frac > 0.95 or (queryEnd - queryStart + 1) >= (queryLen - 10)) \
            and disable == '0' and ident > 0.95: kind = 'GENE-MULT'
        else: kind = 'DUP'

    kindNew = kind

    #========================================
    #========================================
    ## die  if ( ! exists $special_hash{$query} ) ; [zz]
    try: (comment, comment2) = specialHash[query]
    except: (comment, comment2) = ('', '')

    # should this really be an elif chain? as it stands the logic here is
    # rather subtle.
    if kind.startswith('PSSD') and disable == '0': kindNew = 'PSSD2'

    if kind.startswith('PSSD') and re.search(r'[dD]', disable): kindNew = 'PSSD1'

    if re.search(r'PSSD|GENE-SINGLE', kind) and comment == 'RP': kindNew = 'PSSD1'

    if kind.startswith('PSSD') and comment == 'OR': kindNew = 'OR'

    if kind.startswith('PSSD') and re.search('Numt', comment, re.I): kindNew = 'Numt'

    return (kind, kindNew)

##########################
##	SUBROUTINES	##
########################## [zz]

def printRoutine(gffFile, faFile, geneId, gene, gToTrueGeneBounds=None):
##============================================
  # print single gene record to the output file
##============================================ [zz]
    chr = gene[Chr]

    gffFile.write('\t'.join([chr, 'pgene', 'pgene'] + \
                            [gene[f] for f in [Band, ChromStart, ChromEnd, Strand, Query]] + \
                            [geneId] + \
                            [gene[f] for f in [UniqueName, QueryStart, QueryEnd, QueryLen,
                                               Frac, NumOfIns, NumOfDels, NumOfShifts, NumOfStops,
                                               Expect, Ident, Polya, Disable,
                                               NumOfExons, ExonBound, ExonLen, IntronBound, IntronLen,
                                               ClassOld, ClassNew, ShortID, Comment, Comment2]]) + '\n')
    GeneBound=[]                        
    if(gToTrueGeneBounds is not None):
        GeneBound=[gToTrueGeneBounds[geneId][0], gToTrueGeneBounds[geneId][1]] 
    
    faFile.write('  '.join(['>'+geneId] + \
                           [gene[f] for f in [Band, ChromStart, ChromEnd, Strand,
                                              Query, QueryStart, QueryEnd, QueryLen,
                                              Frac, NumOfIns, NumOfDels, NumOfShifts, NumOfStops,
                                              Expect, Ident, Polya, Disable,
                                              NumOfExons, ExonBound, ExonLen, IntronBound, IntronLen,
                                              ClassOld, ClassNew, ShortID, Comment, Comment2]] + GeneBound) + '\n')
	
                            

    faFile.write(gene[QuerySeq] + '\n')
    faFile.write(gene[MatchSeq] + '\n')

#main
if len(sys.argv) > 1:
    chrNames = sys.argv[1:]
else:
    chrNames = [str(i) for i in xrange(1,20)] + ['X']

## cutoff to determine whether to inclde a sequence as a pseudogene [zz]
EXP_CUTOFF   = 1e-10
IDENT_CUTOFF = 0.40

#=============================================================================
# these queries are contaminated with Alus or other repeats or low complexity 
# need to be removed
#============================================================================ [zz]
badQueryHash = {}

counterHash = {}

## #================================================
## # read in the protein information file
## #================================================ [zz]
## print '... reading special query proteins ...'

## proteinFile = openOrFail('./protein-queries/special.mouse.protein.tbl', 'r')
## proteinFile.next()
## for l in proteinFile:
##     (id, comment, comment2) = l[:-1].split('\t')
##     specialHash[id] = (comment, comment2)
##     counterHash[id] = 0
## proteinFile.close()

## ## read in file containing karyotype band information [zz]
## print '...read cytogenic band information...'

## (BandName, BandBegin, BandEnd) = (0, 1, 2)
## bandHash = {}
## bandFile = openOrFail('/home/carriero/bioInformatics/ZhaoleiZhang/ensembl/mouse-030321/mysql/mus_musculus_core_11_3/karyotype.txt.sorted.table', 'r')
## for l in bandFile:
##     f = l[:-1].split('\t')[0:5] # f1, chr, beg, end, band, ....
##     try: bandHash[f[1]].append((f[4], int(f[2]), int(f[3])))
##     except: bandHash[f[1]] = [(f[4], int(f[2]), int(f[3]))]
## bandFile.close()

try:    os.makedirs('./pgenes')
except: pass

for chr in chrNames:
    print '...working on %s...' % chr

    #==========================================================================
    #  declare data structures
    #==========================================================================
    # %exon_hash : hash table to store info read from pexons.gff
    # @gene_ids :  array of pgene ids
    # %gene_exon_hash : hash table to store exon match_ids for each gene
    ##========================================================================= [zz]

    exonHash = {}
    geneIds = []
    geneExonHash = {}

    gffFile = openOrFail('./pexons/%s.pexon.gff' % chr, 'r')
    faFile = openOrFail('./pexons/%s.pexon.fa' % chr, 'r')

    ##======================================
    #  read in information of the matches 
    ##====================================== [zz]
    for l in gffFile:
        if '#' == l[0]: continue

        f = l[:-1].split('\t') + [chr, '', '']

        #=============================
        ## remove contamination 
        #============================== [zz]

        # not actually used???
        if badQueryHash.has_key(f[Query]): continue

        try: geneExonHash[f[GeneId]].append(f[MatchId])
        except: geneExonHash[f[GeneId]] = [(f[MatchId])]
        if 1 == int(f[ExonRank]): geneIds.append(f[GeneId])

        exonHash[f[MatchId]] = f
    gffFile.close()

    ## read in the query and match sequences [zz]
    for l in faFile:
        if '>' != l[0]:
            sys.stderr.write('failed to parse ' + l)
            sys.exit(1)

        matchId = l.split()[0][1:]
        exonHash[matchId][QuerySeq] = faFile.next()[:-1]
        exonHash[matchId][MatchSeq] = faFile.next()[:-1]
    faFile.close()

    ###
    # my ( $last_match_id, $last_gene_id ); [zz]
    geneHash = {}
    gToTrueGeneBounds = {}
    for geneId in geneIds:
        exonArray = geneExonHash[geneId]

        if 1 == len(exonArray):
            ##=====================================================================
            # in most cases, this is a single exon p'gene, but if there are
            # big chunk of gaps in the middle, then it is considered as multiple 
            # exon candidate 
            ##==================================================================== [zz]

            # [:] creates a shallow copy of the data
            gene = geneHash[geneId] = exonHash[geneId][:] 
            gene[NumOfExons] = '1'
            gene += ['(%s..%s)' % (gene[ChromStart], gene[ChromEnd]),
                     str(int(gene[ChromEnd]) - int(gene[ChromStart]) + 1),
                     '.', '.', '', '', '', '', '']
            gToTrueGeneBounds[geneId] = ['(%s..%s)' % (gene[QueryStart], gene[QueryEnd]),
                                         str(int(gene[QueryEnd]) - int(gene[QueryStart]) + 1)]
            if '-' == gene[Strand]:
                gene[ExonBound] = 'com' + gene[ExonBound]
                gToTrueGeneBounds[geneId][0] = 'com' + gToTrueGeneBounds[geneId][0] 
        else:
            ##=====================================================================
            # multiple exon genes/pseudogenes
            ##==================================================================== [zz]
            numOfExons = len(exonArray)
            ## just to change to a better  handle [zz]
            matchIds = exonArray

            # create an empty record
            gene = geneHash[geneId] = ['' for l in geneDataLabels]
            firstExon = exonHash[matchIds[0]]
            lastExon = exonHash[matchIds[-1]]

            gene[NumOfExons]	= str(numOfExons)
            gene[Strand]	= firstExon[Strand]
            gene[Chr]		= firstExon[Chr]
            gene[UniqueName]	= firstExon[UniqueName]
            gene[ChromStart]	= firstExon[ChromStart]
            gene[ChromEnd]	= lastExon[ChromEnd]
            gene[QueryLen]	= firstExon[QueryLen]
            gene[Query]		= firstExon[Query]

            ##===================================================================
            ## get query start and end on the chromosome
            ## and sort the @match_ids, so the ids are from 5' to 3' on the protein 
            #==================================================================== []
            if '-' == gene[Strand]:
                gene[QueryStart] = lastExon[QueryStart]
                gene[QueryEnd] = firstExon[QueryEnd]
                sortedMatchIds = matchIds[::-1]
            else:
                gene[QueryStart] = firstExon[QueryStart]
                gene[QueryEnd] = lastExon[QueryEnd]
                sortedMatchIds = matchIds

            ## get the query and match sequences [zz]
            for matchId in sortedMatchIds:
                gene[QuerySeq] += exonHash[matchId][QuerySeq] + '  '
                gene[MatchSeq] += exonHash[matchId][MatchSeq] + '  '
            
            ## remove white spaces [zz]

            ## trailing whitespace, that is.
            gene[QuerySeq] = gene[QuerySeq].rstrip()
            gene[MatchSeq] = gene[MatchSeq].rstrip()

            numOfIns = numOfDels = numOfShifts = numOfStops = 0
            exonBound, exonLen, tExonBound, tExonLen = '', '', '', ''
            ## total length of the match, and total identical residues [zz]
            expect   = 1.
            totalLen = totalIdent = 0.
            disable  = '0'
            querySeq = matchSeq = ''

            for id in matchIds:

                #===================================================================
                # here we use @match_ids instead of @sorted_match_ids because we
                # need the match to be sorted on chrom_start, So we can get the gene
                # boundaries et al.
                #==================================================================== [zz]

                exon = exonHash[id]

                # the boundary of exons [zz]
                exonBound += '%s..%s ' % (exon[ChromStart], exon[ChromEnd]) 
                tExonBound += '%s..%s ' % (exon[QueryStart], exon[QueryEnd]) 
                exonLen += '%d ' % (int(exon[ChromEnd]) - int(exon[ChromStart]) + 1)
                tExonLen += '%d ' % (int(exon[QueryEnd]) - int(exon[QueryStart]) + 1)
                numOfIns += int(exon[NumOfIns])
                numOfDels += int(exon[NumOfDels])
                numOfShifts += int(exon[NumOfShifts])
                numOfStops += int(exon[NumOfStops])
                
                if float(exon[Expect]) < expect: expect = float(exon[Expect])

                totalLen += len(exon[MatchSeq])
                totalIdent += len(exon[MatchSeq]) * float(exon[Ident])

                if 'd' == exon[Disable] and 'D' != disable: disable = 'd'
                elif 'D' == exon[Disable]: disable = 'D'


            gene[ExonBound] = '(' + exonBound.rstrip() + ')'
            gene[ExonLen] = '(' + exonLen.rstrip() + ')'
            gToTrueGeneBounds[geneId] = [ '(' + tExonBound.rstrip() + ')', '(' + tExonLen.rstrip() + ')']
            (gene[NumOfIns], gene[NumOfDels]) = (str(numOfIns), str(numOfDels))
            (gene[NumOfShifts], gene[NumOfStops]) = (str(numOfShifts), str(numOfStops))
            gene[Expect] = str(expect)
            gene[Disable] = disable
            gene[Ident] = '%.2f' % (totalIdent/totalLen)
            gene[Frac] = '%.2f' % ((float(gene[QueryEnd]) - float(gene[QueryStart]) + 1)/float(gene[QueryLen]))

            ## to determine polya, we need @sorted_match_ids [zz]
            gene[Polya] = exonHash[sortedMatchIds[-1]][Polya]

            ##============================
            # determine intron boundaries
            ##=========================== [zz]

            # having just put together that nicely formatted string, we take it apart...
            # we should just built and auxillary list

            ## remove the begining of the first exon and the end of the last exon [zz]
            intronArray = [v for v in re.sub(r'\(|\)|\.', ' ', exonBound).split()][1:-1]

            intronBound = ''
            intronLen = ''
            while intronArray:
                start = int(intronArray.pop(0))
                stop = int(intronArray.pop(0))
                intronBound += '%d..%d ' % (start+1, stop-1)
                intronLen += '%d ' % (stop - start - 1)

            gene[IntronBound] = '(' + intronBound.rstrip() + ')'
            gene[IntronLen]   = '(' + intronLen.rstrip() + ')'

            #add "com" if the match is on minus strand [zz]
            if '-' == gene[Strand]:
                gene[ExonBound] = 'com' + gene[ExonBound]
                gene[IntronBound] = 'com' + gene[IntronBound]
                gene[ExonLen] = 'com' + gene[ExonLen]
                gene[IntronLen] = 'com' +  gene[IntronLen]
                gToTrueGeneBounds[geneId][0] = 'com' + gToTrueGeneBounds[geneId][0] 
                gToTrueGeneBounds[geneId][1] = 'com' + gToTrueGeneBounds[geneId][1] 

        ##===========================================================
        # decied whether a p-gene is processed, things are considered
        # unpreocessed unless proved otherwise 
        ##=========================================================== [zz]
        (classOld, classNew) = getClass(EXP_CUTOFF, IDENT_CUTOFF, int(gene[NumOfExons]),
                                        float(gene[Frac]), float(gene[Expect]), float(gene[Ident]),
                                        gene[Polya], gene[Disable], gene[Query],
                                        int(gene[QueryStart]), int(gene[QueryEnd]), int(gene[QueryLen]))

        gene[ClassOld] = classOld
        gene[ClassNew] = classNew
        if specialHash.has_key(gene[Query]):
            (gene[Comment], gene[Comment2]) = specialHash[gene[Query]]

        #get the karyo band name, the genes are already sorted by starting coordinates [zz]
        #        gene[Band] = getBand(chr, gene[GeneId], int(gene[ChromStart]), int(gene[ChromEnd]), bandHash[chr])
        gene[Band] = 'NoBandData'

        #=======================================
        # give it a short name <= 10 characters
        #=======================================
        # don't count on this being unique --- if analysis is run separately for each chrom/seq ids could be re-used.
        if classNew.startswith('PSSD'):
            q = gene[Query]
            c = counterHash.get(q, 1)
            gene[ShortID] =  '%s===%d' % (q, c)
            counterHash[q] = c + 1

    ##===============================
    ##  print out 
    ##=============================== [zz]
    header = '\t'.join([
            '#chr', 'pgene', 'pgene', 'band', 'chrom_start', 'chrom_end', 'strand', 'query', 'gene_id',
            'unique_name', 'query_start', 'query_end', 'query_len', 'frac', 'num_of_ins', 'num_of_dels',
            'num_of_shifts', 'num_of_stops', 'expect', 'ident', 'polya', 'disable', 'num_of_exons',
            'exon_bound', 'exon_len', 'intron_bound', 'intron_len', 'class_old', 'class_new', 'short_ID', 'comment', 'comment2']) + '\n'

    allGffFile = openOrFail('./pgenes/%s.sptrembl.all.gff' % chr, 'w')
    allGffFile.write(header)
    allFaFile = openOrFail('./pgenes/%s.sptrembl.all.fa' % chr, 'w')

    pssd1GffFile = openOrFail('./pgenes/%s.sptrembl.pssd1.gff' % chr, 'w')
    pssd1GffFile.write(header)
    pssd1FaFile = openOrFail('./pgenes/%s.sptrembl.pssd1.fa' % chr, 'w')

    pssd2GffFile = openOrFail('./pgenes/%s.sptrembl.pssd2.gff' % chr, 'w')
    pssd2GffFile.write(header)
    pssd2FaFile = openOrFail('./pgenes/%s.sptrembl.pssd2.fa' % chr, 'w')

    fragGffFile = openOrFail('./pgenes/%s.sptrembl.frag.gff' % chr, 'w')
    fragGffFile.write(header)
    fragFaFile = openOrFail('./pgenes/%s.sptrembl.frag.fa' % chr, 'w')

    dupGffFile = openOrFail('./pgenes/%s.sptrembl.dup.gff' % chr, 'w')
    dupGffFile.write(header)
    dupFaFile = openOrFail('./pgenes/%s.sptrembl.dup.fa' % chr, 'w')

    realGffFile = openOrFail('./pgenes/%s.sptrembl.real.gff' % chr, 'w')
    realGffFile.write(header)
    realFaFile = openOrFail('./pgenes/%s.sptrembl.real.fa' % chr, 'w')

    for geneId in geneIds:
        gene = geneHash[geneId]
        printRoutine(allGffFile, allFaFile, geneId, gene, gToTrueGeneBounds)

        classNew = gene[ClassNew]
        if   'PSSD1' == classNew:	printRoutine(pssd1GffFile, pssd1FaFile, geneId, gene)
        elif 'PSSD2' == classNew:	printRoutine(pssd2GffFile, pssd2FaFile, geneId, gene)
        elif 'DUP' == classNew:		printRoutine(dupGffFile, dupFaFile, geneId, gene)
        elif 'FRAG' == classNew:
            # fragments could be chopped-up pssd pgenes or singe non-processed exons or isolated real exons [zz]
            				printRoutine(fragGffFile, fragFaFile, geneId,gene)
        elif classNew.startswith('GENE'):
            # real gene candidates [zz]
            				printRoutine(realGffFile, realFaFile, geneId, gene)
        elif 'FP' == classNew: 		pass
        elif re.search(r'Numt', classNew, re.I): pass
        elif 'OR' == classNew: 		pass
        else:
            sys.stderr.write('this should not happen 2\n %s\n' % classNew)
            sys.exit(1)


    allGffFile.close()
    allFaFile.close()

    # bug in original closing non-existent files?
    pssd1GffFile.close()
    pssd1FaFile.close()
    pssd2GffFile.close()
    pssd2FaFile.close()

    fragGffFile.close()
    fragFaFile.close()
    dupGffFile.close()
    dupFaFile.close()
    realGffFile.close()
    realFaFile.close()

    tbf =  openOrFail('./pgenes/%s.parent.bounds.txt' % chr, 'w')
    for geneId in geneIds:
        print >>tbf, geneId, gToTrueGeneBounds[geneId]
    tbf.close()
