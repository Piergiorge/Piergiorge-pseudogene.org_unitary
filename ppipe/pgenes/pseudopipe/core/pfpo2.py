#!/usr/bin/env python

# informal parser for the output of tfasty, which we assume has the form

## header
## >>> comparison data
## [; name value]*
## [>> match data
##  [; name value]*
##  [[a-z]+\n]*
##  > seq data
##  [; name value]*
##  [[a-z]+\n]*
## ]*

#

import re, sys

nameRE = re.compile(r'; ([^:]+): ')

def readBindings(file):
    # process group of ; bindings
    b = {}
    for l in file:
        if ';' != l[0]: break
        try:
            mo = nameRE.match(l)
            name = mo.groups()[0]
            remainder = l[mo.span()[1]:-1]
        except:
            sys.stderr.write('Do not understand: ' + l)
            return None

        if b.has_key(name):
            if 'pg_name' != name and 'pg_ver' != name:
                sys.stderr.write('Warning, duplicated name: ' + l)

        b[name] = remainder

    return (b, l)

def readSequence(first, file):
    seq = [first]
    for l in file:
        if '>' == l[0]: break
        seq.append(l)

    return (seq, l)

def processComparison(first, file):
    # read header
    header = [first]
    for l in file:
        if l.startswith('The best scores are'): break
        header.append(l)
    else:
        sys.stderr.write('End of file reading header.')
        return None

    summaryData = []
    for l in file:
        if '\n' == l: break
        summary = l.split()
        # 5th field can be ( ddd) or (dddd). if the first case, split yields '(' and 'ddd)'. Fix this.
        if '(' == summary[5]:
            summary.pop(5)
            summary[5] = '( '+summary[5]
        summaryData.append(summary)
        
    l = file.next()
    if not l.startswith('>>>'):
        sys.stderr.write('Out of sync at line: ' + l)
        return None

    # this appears to work for both 3.1 and 3.4 output formats
    proChroPos = summaryData[0][0]

    # read comparison name/values
    (comparisonBindings, l) = readBindings(file)
    if not l:
        sys.stderr.write('End of file while reading comparison parameters.\n')
        return None

    matches = []
    while 1:
        if l.startswith('>>><<<'): break
        
        #read match name/values
        if not l.startswith('>>'):
            sys.stderr.write('Out of sync at line: ' + l)
            return None
        match = l[2:]
        (matchBindings, l) = readBindings(file)
        if not l:
            sys.stderr.write('End of file while reading match parameters.\n')
            return None

        # read query sequence
        if not l.startswith('>'):
            sys.stderr.write('Out of sync at line: ' + l)
            return None
        qSeqLabel = l[1:]
        (qSeqBindings, l) = readBindings(file)
        if not l:
            sys.stderr.write('End of file while reading query sequence parameters.\n')
            return None
        (qSeq, l) = readSequence(l, file)
        if not l:
            sys.stderr.write('End of file while reading query sequence.\n')
            return None

        # read matching sequence
        if not l.startswith('>'):
            sys.stderr.write('Out of sync at line: ' + l)
            return None
        mSeqLabel = l[1:]
        (mSeqBindings, l) = readBindings(file)
        if not l:
            sys.stderr.write('End of file while reading matching sequence parameters.\n')
            return None
        (mSeq, l) = readSequence(l, file)
        if not l:
            sys.stderr.write('End of file while reading matching sequence.\n')
            return None

        matches.append((match, matchBindings, qSeqLabel, qSeqBindings, qSeq, mSeqLabel, mSeqBindings, mSeq))
        
    # read trailer
    trailer = []
    for l in file:
        if l.startswith('Function used was'): break
        trailer.append(l)
    else:
        sys.stderr.write('End of file while reading trailer.\n')
        return None

    return (proChroPos, header, summaryData, matches, trailer)

def parseTFastyOutput(file):
    comparisons = {}
    for l in file:
        c = processComparison(l, file)
        if len(c[3]) == 0: continue # no matches reported.
        pro = c[0].split('===',1)[0] # uglyish hack to deal with zz naming conventions. last two bits are chr and loc, the rest is the name.
        if comparisons.has_key(pro):
            comparisons[pro][c[0]] = c
        else:
            comparisons[pro] = {c[0]: c}

    return comparisons
