#!/usr/bin/env python
"""
@author Spencer Bliven <sbliven@ucsd.edu>
"""

import sys
import os
import optparse
import re
from cStringIO import StringIO
from proteinsearch import disambiguateRe


def outputPML(word,structures,prefix=None,out=sys.stdout):
    """outputPML(string, [(pdb,chain)...], handle) -> None
    """
    if len(structures) < 1:
        raise ValueError("No structures found")

    out.write("delete *\n")

    if prefix is None:
        prefix = word

    out.write("# Align structures containing '%s'\n" % word)
    out.write("fetch %s, async=0\n" % " ".join(set(["%s%s"%pdbchain for pdbchain in structures])))

    # Define selections
    for pdb,chain in structures:
        #out.write("create {0}.{1}, {0} and chain {1}\n".format(pdb,chain) )
        out.write("findseq {0}, {1}{2}, {3}_{1}{2}, firstOnly=1\n".format(word,pdb,chain,prefix))
    selections = ["{0}_{1}{2}".format(prefix,pdb,chain) for pdb,chain in structures]
    out.write("select %s, %s\n" % (prefix, " ".join(selections) ) )

    # superimpose on first structure
    for pdb,chain in structures[1:]:
        out.write("pair_fit {0}_{1}{2} and n. ca, {0}_{3}{4} and n. ca\n".format(prefix,pdb,chain,structures[0][0],structures[0][1]))

    # View
    out.write("""
pretty_protein
color_obj
set cartoon_transparency, .5
show sticks
hide everything, elem h

set antialias, 2
bg_color white
set label_font_id, 10
set label_font_id, 10
#set label_size, -2
set label_size, 56
set label_position, (0,0,5) #labels in front of atoms
set cartoon_transparency, .5
#set label_color, lightorange
set label_color, black
set label_outline_color, black
""")

    for sele in selections:
        out.write("create {0}, {0}\n".format(sele) )
    out.write("label {0} and n. ca, one_letter[resn]\n".format(selections[0]) )

    out.write("orient %s\n"%prefix)

    for pdb,chain in structures:
        out.write("set stick_transparency, .75, {0}{1}\n".format(pdb,chain))
    
    #out.write("ray\n")



if __name__ == "__main__":
    parser = optparse.OptionParser( usage="usage: python %prog [options] matchfile word" )
    parser.add_option("-v","--verbose", help="Long messages",
        dest="verbose",default=False, action="store_true")
    parser.add_option("--exact",help="Exact matches; don't disambiguate B,J,Z,and X",
        dest="exactmatch",default=False,action="store_true")
    (options, args) = parser.parse_args()

    if len(args) != 2:
        parser.print_usage()
        parser.exit("Error: Expected 2 argument, but found %d"%len(args) )

    matchfilename, word = args

    if options.exactmatch:
        wordexpr = word
    else:
        wordexpr = disambiguateRe(word)

    wordre = re.compile("%s\s+(\w{4})\.(\w+)\s$" % word,re.IGNORECASE)

    structures = []
    with open(matchfilename,'r') as matchfile:

        for line in matchfile:
            # filter matching lines
            match = wordre.match(line)
            if match:
                pdb,chain = match.groups()
                structures.append((pdb,chain))

    outputPML(wordexpr,structures,word)

