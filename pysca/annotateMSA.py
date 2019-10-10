#! /usr/bin/env python3
"""
The annotateMSA script provides utilities to automatically annotate sequence
headers (for a fasta file) with taxonomic information. Currently this can be
done in one of two ways:

    1) For PFAM alignments, annotations can be extracted from the file
       pfamseq.txt (please download from:
       ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/database_files/pfamseq.txt.gz)

    2) For Blast alignments, annotations can be added using the NCBI Entrez
       utilities provided by BioPython. In this case, an additional command
       line argument (--idType, see below) should specify whether sequences are
       identified by GI or accession nmbers. These numbers are then used to
       query NCBI for taxonomy information (note that this approach requires a
       network connection).

To extract GI or accession numbers, use the scripts alnParseGI.py or
alnParseAcc.py, respectively.

For both the Pfam and NCBI utilities, the process of sequence annotation *can
be slow* (on the order of hours, particularly for NCBI entrez with larger
alignments). However, the annotation process only needs to be run once per
alignment.

**Arguments**
    Input_MSA.fasta (some input sequence alignment)

**Keyword Arguments**
    -o, --output        Specify an output file, Output_MSA.an
    -a, --annot         Annotation method. Options are 'pfam' or 'ncbi'.
                        Default: 'pfam'
    -t, --idType        ID used to annotate FASTA input ('gi' or 'acc').
                        Default: 'acc'
    -l, --idList        This argument is necessary for the 'ncbi' method.
                        Specifies a file containing a list of GI numbers
                        corresponding to the sequence order in the alignment; a
                        number of "0" indicates that a GI number wasn't
                        assigned for a particular sequence.
    -g, --giList        Deprecated. Identical to '--idList' and kept to keep
                        the CLI consistent with older versions of pySCA.
    -p, --pfam_seq      Location of the pfamseq.txt file. Defaults to
                        path2pfamseq (specified at the top of scaTools.py)

**Examples**::

  ./annotateMSA.py ../data/PF00186_full.txt -o ../output/PF00186_full.an -a 'pfam'
  ./annotateMSA.py ../data/DHFR_PEPM3.fasta -o ../output/DHFR_PEPM3.an -a 'ncbi' -l ../data/DHFR_PEPM3.gi -t 'gi'

:By: Rama Ranganathan, Kim Reynolds
:On: 9.22.2014

Copyright (C) 2015 Olivier Rivoire, Rama Ranganathan, Kimberly Reynolds

This program is free software distributed under the BSD 3-clause license,
please see the file LICENSE for details.
"""

import time
import sys
import argparse
import os
from Bio import Entrez
import scaTools as sca

import settings

Entrez.email = settings.entrezemail  # PLEASE use your email! (see settings.py)

if __name__ == '__main__':
    # parse inputs
    parser = argparse.ArgumentParser()
    parser.add_argument("Input_MSA", help='input sequence alignment')
    parser.add_argument("-o", "--output", dest="output", default='sys.stdout',
                        help="Outputfile name. Default: sys.stdout")
    parser.add_argument("-a", "--annot", dest="annot", default='pfam',
                        help="Annotation method. Options are 'pfam' or 'ncbi'."
                             " Default: 'pfam'")
    parser.add_argument("-t", "--idType", dest="idType", default='acc',
                        help="Sequence identifier. Options are 'acc' or 'gi'."
                             " Default: 'acc'")
    parser.add_argument("-l", "--idList", dest="idList", default=None,
                        help="This argument is necessary for the 'ncbi' "
                             "method. Specifies a file containing a list of "
                             "GI or accession numbers corresponding to the "
                             "sequence order in the alignment; a number of 0 "
                             "indicates that one wasn't assigned for a "
                             "particular sequence.")
    parser.add_argument("-g", "--giList", dest="idList", default=None,
                        help="Command kept for compatibility with previous "
                             "versions. Use '-l' or '--idList' instead.")
    parser.add_argument("-p", "--pfam_seq", dest="pfamseq", default=None,
                        help="Location of the pfamseq.txt file. Defaults to "
                             "path2pfamseq (specified in settings.py)")
    parser.add_argument("-d", "--pfam_db", dest="pfamdb", default=None,
                        help="Location of the pfamseq.db file. Priority over "
                             "pfamseq.txt file. Defaults to path2pfamseqdb "
                             "(specified in settings.py)")
    parser.add_argument("-e", "--entrez_email", dest="email", default=None,
                        help="email address for querying Entrez web API")
    options = parser.parse_args()

    if (options.annot != 'pfam') & (options.annot != 'ncbi'):
        sys.exit("The option -a must be set to 'pfam' or 'ncbi' - other"
                 " keywords are not allowed.")

    if (options.annot == 'ncbi'):
        if (options.idList is None):
            sys.exit("To use NCBI entrez annotation, you must specify a file "
                     "containing a list of GI numbers (see the --idList "
                     "argument).")
        if not ((options.idType == "acc") or (options.idType == "gi")):
            sys.exit("To use NCBI Entrez annotation, specify a valid "
                     "sequence identifier type ('acc' or 'gi').")

    if options.annot == 'pfam':
        # Annotate a Pfam alignment
        if options.pfamdb is not None:  # default to db query over txt search
            sca.AnnotPfamDB(options.Input_MSA, options.output, options.pfamdb)
        elif options.pfamseq is not None:
            sca.AnnotPfam(options.Input_MSA, options.output, options.pfamseq)
        else:
            # If no database or text file supplied to annotateMSA, then default
            # to the files defined in settings.py.
            if os.path.exists(settings.path2pfamseqdb):
                sca.AnnotPfamDB(options.Input_MSA, options.output)
            elif os.path.exists(settings.path2pfamseq):
                sca.AnnotPfam(options.Input_MSA, options.output)
    elif options.annot == 'ncbi':
        # Annotate using GI numbers/NCBI entrez
        if options.email is None:
            sca.AnnotNCBI(options.Input_MSA, options.output, options.idList,
                          options.idType)
        else:
            sca.AnnotNCBI(options.Input_MSA, options.output, options.idList,
                          options.idType, options.email)
