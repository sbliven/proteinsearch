proteinsearch
=============

Scripts to search for specific peptide sequences in the PDB. For instance, what
if you want to find all english words which appear in the PDB?

This code was written for the analysis at [The Structure of
Spencer](http://acsweb.ucsd.edu/~sbliven/2014/05/the-structure-of-spencer/).


Finding matches: proteinsearch.py
---------------------------------

Two input files are required:

1. A query file, giving one word per line
2. The sequence files, in the form of a PDB tabular report.

For the query file I used the ispell english dictionary from the [wordlist
project](http://wordlist.sourceforge.net/) and concatenated all the
subdictionaries to allow both British and American spellings. I also used  the
[1000 most popular male & female baby
names](http://www.ssa.gov/cgi-bin/popularnames.cgi) from 1986 for variety.

The PDB lets you download a tabular report of any search. I used a search for
all structures (100,147 on May 14, 2014), and then generated a report including
the sequences. After a few minutes of waiting you get a quoted CSV file. The
script relies on the following columns:

    1. PDB ID
    2. Chain ID
    6. Sequence
    9. Macromolecular Type (filtered for the string 'Polypeptide(L)')

Yes, fasta files would make more sense. Submit a feature request.

**Dependencies**

The script uses [progressbar](https://pypi.python.org/pypi/progressbar/2.2) unless you use the -q option.

It also relies heavily on the [esmre](https://pypi.python.org/pypi/esmre/0.3.1)
library for fast string searching. It implements the nice
[Aho-Corasick](https://en.wikipedia.org/wiki/Aho%E2%80%93Corasick_string_matching_algorithm)
algorithm to search for all queries simultaneously.

Both dependencies should be easily installed using `pip`.


**Examples**

To search for everything in your query file, use the following command

```bash
proteinsearch.py --exact --min_len 0 \
    --sort-length -m matches.tsv -c counts.tsv \
    query.txt tabularResults.csv
```

The `--exact` and `-l 0` parameters tell it to use all input words literally,
without any filtering. The `-m` and `-c` options give output files. The match
file gives the query and PDB ID for all search results. The counts file lists
all successful queries and the number of places they matched (may be more than
one per PDB). The `--sort-length` option causes the counts file to output in
descending word length, as long words are generally more
surprising/interesting.

There are some more options to ease filtering of the query files, since
spellcheck libraries contain lots of junk. The `--trim` option removes
possessive suffixes, and the `--min_len` option ignores short queries (defaults
to at least 5 residues).

Finally, what about the 6 letters which don't stand for amino acids? Well, they
can be interpretted as ambiguous characters. Thus, the query 'BIT' would match
either 'DIT' or 'NIT' (or 'BIT', if the PDB author wasn't sure whether the
residue was aspartate or aspagine). This behavior is enabled by default.

```bash
proteinsearch.py --max-x=1 --trim \
    -m matches.tsv -c counts.tsv --sort-length \
    query.txt tabularResults.csv
```

Since 'X' can stand for any letter, the counts were being dominated by queries
with lots of 'X', such as 'xxxix' (hey, sometimes you want to spell check roman
numerals). So the `--max-x` option filters these out.


Viewing Results: matchtopml.py
==============================

To produce nice figures in pymol, `matchtopml.py` will create a PyMol script file (.pml) to align all the matches for a particular word. It needs the matches.tsv file from proteinsearch.py, as well as the query you're interested in (probably selected from counts.tsv).

```bash
matchtopml.py matches.tsv word > word.pml
```

You may need to customize the `outputPML` method to suit your preferences. In
addition to commands to fetch files, find the sequence, and align structures,
the script also does some extensive visual configuration particular to my pymol
setup.

The script relies on the [findseq](http://pymolwiki.org/index.php/Findseq)
(which you may already have from the [Pymol script
repository](https://github.com/Pymol-Scripts/Pymol-script-repo), as well as
some custom aliases I use (pretty_protein,
[color_obj](http://www.pymolwiki.org/index.php/Color_Objects), etc.).
