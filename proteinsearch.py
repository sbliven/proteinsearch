#!/usr/bin/env python
"""
@author Spencer Bliven <sbliven@ucsd.edu>
"""

import sys
import os
import optparse
#import Bio
#from Bio.Seq import Seq
#from Bio import motifs
#from Bio.Alphabet.IUPAC import ExtendedIUPACProtein
try:
    import progressbar
except ImportError:
    sys.stderr.write("Warning: You don't seem to have progressbar installed. Install it or use -q\n")

import esmre

# From http://www.daimi.au.dk/~mailund/suffix_tree.html
#from suffix_tree import GeneralisedSuffixTree

def disambiguateRe(seq, addoriginal=True):
    """disambiguate(string) -> string

    returns a regular expression string replacing ambiguous residues (B,J,Z,X)
    with their matching character sets.

    If addoriginal is true (the default), also allow the ambiguous character to match.

    Example:
    disambiguate("ABERJF") => "A[BDN]ER[JIL]F"
    """
    regex = seq.upper()
    stdaas = "ACDEFGHIKLMNPQRSTVWY"
    if addoriginal:
        regex = regex.replace("B","[BDN]")
        regex = regex.replace("J","[JIL]")
        regex = regex.replace("Z","[ZEQ]")
        regex = regex.replace("X","[X%s]"%stdaas)
    else:
        regex = regex.replace("B","[DN]")
        regex = regex.replace("J","[IL]")
        regex = regex.replace("Z","[EQ]")
        regex = regex.replace("X","[%s]"%stdaas)
    return regex

# work around bug in esmre parsing of regular expressions
def disambiguate(seq, addoriginal=True,max_x=-1):
    """disambiguate(string) -> generator() of strings

    Expands ambiguous residues in a sequence to all possible matching strings.
    Expanded residues are:

    "J" -> ["I","L"]
    "B" -> ["D","N"]
    "Z" -> ["E","Q"]
    "X" -> A-Z \ BJOUXZ (20 standard AAs)

    These letters are each expanded independently, so words with multiple X can
    result in quite large iterators.

    If addall is true, also include the ambiguous character in the output (eg J -> JIL)

    max_x gives the maximum number of 'X' characters which can be
    disambiguated, to limit combinatorial explosion. Negative numbers indicate
    unlimited 'X's
    """
    if len(seq) < 1:
        yield ""
        return

    stdaas = "ACDEFGHIKLMNPQRSTVWY"
    nonstdaas = "OU"

    #not included in the output
    ignoredchars = " \t\n.-_+'\"\\/"

    first = seq[0].upper()
    other = seq[1:]

    if first in stdaas or first in nonstdaas:
        for pos in disambiguate(other,addoriginal,max_x):
            yield first + pos
    elif first == "J":
        for pos in disambiguate(other,addoriginal,max_x):
            yield "I" + pos
            yield "L" + pos
            if addoriginal:
                yield first + pos
    elif first == "B":
        for pos in disambiguate(other,addoriginal,max_x):
            yield "D" + pos
            yield "N" + pos
            if addoriginal:
                yield first + pos
    elif first == "Z":
        for pos in disambiguate(other,addoriginal,max_x):
            yield "E" + pos
            yield "Q" + pos
            if addoriginal:
                yield first + pos
    elif first == "X":
        if max_x < 0:
            # infinite x; normal case
            for pos in disambiguate(other,addoriginal,max_x):
                for aa in stdaas:
                    yield aa + pos
                if addoriginal:
                    yield first + pos
        elif max_x > 0:
            # finite x; decrement
            for pos in disambiguate(other,addoriginal,max_x-1):
                for aa in stdaas:
                    yield aa + pos
            if addoriginal:
                for pos in disambiguate(other,addoriginal,max_x):
                    yield first + pos
        else:
            # no x; only do original
            if addoriginal:
                for pos in disambiguate(other,addoriginal,max_x):
                    yield first + pos
            else:
                return
    elif first in ignoredchars:
        for pos in disambiguate(other,addoriginal,max_x):
            yield pos
    else:
        raise ValueError( ("Unknown character '%s'" %first) )


def trimWord(word):
    if word[-2:] == "'s":
        word = word[-2:]

    return word
if __name__ == "__main__":

    #print list(disambiguate(""))
    #print list(disambiguate("B"))
    #print list(disambiguate("J"))
    #print list(disambiguate("JF"))
    #print list(disambiguate("AJBF"))
    #print list(disambiguate("AGFJKB"))
    #print list(disambiguate("AXCXDB",True,0))
    #print list(disambiguate("AXCXDB",False,0))
    assert( len(list(disambiguate("AXCXDB",True,0)))    == 3 )
    assert( len(list(disambiguate("AXCXDB",False,0)))   == 0 )
    assert( len(list(disambiguate("AXDB",True,1)))      == 21*3 )
    assert( len(list(disambiguate("AXDB",False,1)))     == 20*2 )
    assert( len(list(disambiguate("AXCXDB",True,1)))    == 20*1*3+1*20*3+ 3)
    assert( len(list(disambiguate("AXCXDB",False,1)))   == 0 )
    assert( len(list(disambiguate("AXCXDB",True,2)))    == 21*21*3)
    assert( len(list(disambiguate("AXCXDB",False,2)))   == 20*20*2)

    parser = optparse.OptionParser( usage="usage: python %prog [options] vocab sequences" )
    parser.add_option("-v","--verbose", help="Long messages",
        dest="verbose",default=False, action="store_true")
    parser.add_option("-l","--min_len",help="Minimum word length (default 5)",
        dest="min_len",default=5, type="int")
    parser.add_option("-m","--matches",help="Output file for matches, or '-' for stdout",
        dest="matchfile",default=None)
    parser.add_option("-c","--counts", help="Output file for counts, or '-' for stdout",
        dest="countfile",default=None)
    parser.add_option("-q","--no-progress",help="Display progress bar",
        dest="progress",default=True,action="store_false")
    parser.add_option("--exact",help="Exact matches; don't disambiguate B,J,Z,and X",
        dest="exactmatch",default=False,action="store_true")
    parser.add_option("--max-x",help="Maximum number of 'x' in the word file to allow ambiguously",
        dest="max_x",default=-1,type="int")
    parser.add_option("--trim",help="Trim word suffixes like possives",
        dest="trim",default=False,action="store_true")
    parser.add_option("--sort-length",help="Sort countfile by word length",
        dest="sortlength",default=False,action="store_true")
    (options, args) = parser.parse_args()

    if len(args) != 2:
        parser.print_usage()
        parser.exit("Error: Expected 2 argument, but found %d"%len(args) )

    vocabfilename, sequencesfilename = args

    # read vocab
    vocab = esmre.Index()

    words = 0
    patterns = 0

    with open(vocabfilename,'r') as vocabfile:
        if options.progress:
            widgets=['Reading Vocab       ',progressBar.Percentage(),progressBar.Bar()]
            size = os.stat(vocabfilename).st_size or 1
            pbar = progressBar.ProgressBar(widgets=widgets,maxval=size).start()

        #vocabulary = [line.strip().upper() for line in vocabfile if len(line.strip())>=options.min_len]
        for line in vocabfile:
            if options.progress:
                pbar.update(vocabfile.tell())

            word = line.strip()
            if options.trim:
                word = trimWord(word)
            if len(word) >= options.min_len:
                words += 1
                if options.exactmatch:
                    ambiguous = [word.upper()]
                else:
                    ambiguous = disambiguate(word,max_x=options.max_x)
                for unambiguous in ambiguous:
                    patterns += 1
                    vocab.enter(unambiguous,word)

    if options.progress:
        pbar.finish()
        print "Read %d patterns from %d valid words in vocabulary" % (patterns,words)

    #vocabtree = GeneralisedSuffixTree(vocabulary)
    #vocabmotif = motifs.create( [ Seq(word, ExtendedIUPACProtein) for word in vocabulary] )


    # Open output files
    matchfile = None
    if options.matchfile:
        if options.matchfile == "-":
            matchfile = sys.stdout
        else:
            matchfile = open(options.matchfile,'w')
            matchfile.write("Word\tPDB.Chain\n")

    countfile = None
    if options.countfile:
        if options.countfile == '-':
            countfile = sys.stdout
        else:
            countfile = open(options.countfile,'w')

    found_words = {}

    with open(sequencesfilename,'r') as sequencesfile:

        if options.progress:
            widgets=['Finding matches     ',progressBar.Percentage(),progressBar.Bar()]
            size = os.stat(sequencesfilename).st_size or 1
            pbar = progressBar.ProgressBar(widgets=widgets,maxval=size).start()

        headerline = sequencesfile.readline()
        line_num = 2
        for line in sequencesfile:
            if options.progress:
                pbar.update(sequencesfile.tell())

            line = line.strip()[1:-1]
            fields = line.split('","')

            if len(fields) < 6: #should be 10 total
                continue

            pdb = fields[0]
            chain = fields[1]
            sequence = fields[5].upper()
            peptide = "Polypeptide(L)" == fields[8]

            #search for words
            if peptide and len(sequence) > options.min_len:
                #O(n^2) algorithm is easier than suffix tree
                #    for word in vocabulary:
                #        if word in sequence:
                #            found_words[word] = found_words.get(word,0)+1

                result = vocab.query(sequence)

                for word in result:
                    if matchfile is not None:
                        matchfile.write("%s\t%s.%s\n" % (word, pdb,chain) )
                    if countfile is not None:
                        found_words[word] = found_words.get(word,0)+1


            line_num += 1
            #if line_num > 100:
            #break

    if options.progress:
        pbar.finish()

    if matchfile is not None and matchfile != sys.stdout:
        matchfile.close()

    # Output counts
    if countfile is not None:
        if options.progress:
            widgets=['Writing counts      ',progressBar.Percentage(),progressBar.Bar()]
            size = len(found_words)
            pbar = progressBar.ProgressBar(widgets=widgets,maxval=max(size,1)).start()
            i = 0

        word_counts = found_words.items()

        if options.sortlength:
            word_counts.sort(key=lambda x:(len(x[0]),x[1]),reverse=True)
        else:
            word_counts.sort(key=lambda x:(x[1],len(x[0])),reverse=True)

        for word_count in word_counts:
            if options.progress:
                pbar.update(i)
                i += 1

            countfile.write( "%s\t%d\n" % word_count )

        if countfile != sys.stdout:
            countfile.close()

        if options.progress:
            pbar.finish()


