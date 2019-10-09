#! /usr/bin/env python3
"""
A script to convert GI numbers in the header of a FASTA file to accession
numbers.

**Arguments**
    Input_MSA.fasta (the alignment to be processed)

**Keyword Arguments**
    --output             output file name, default: FilteredAln.fa
    --delim              delimiter separating header fields, default: "_"
    --email              email to associate with Entrez web API queries

:By: Kim Reynolds
:On: 6.5.2015

Copyright (C) 2015 Olivier Rivoire, Rama Ranganathan, Kimberly Reynolds

This program is free software distributed under the BSD 3-clause
license, please see the file LICENSE for details.
"""

import argparse
import time
import sys
from Bio import Entrez
import scaTools as sca

import settings

if __name__ == '__main__':
    # Parse inputs
    parser = argparse.ArgumentParser()
    parser.add_argument("alignment", help='Input Sequence Alignment')
    parser.add_argument("-o", "--output", dest="outputfile",
                        default='output.acc',
                        help="specify an outputfile name")
    parser.add_argument("-d", "--delim", dest="delim", default="_",
                        help="specify the field delimiter in the header")
    parser.add_argument("-e", "--entrez_email", dest="email", default=None,
                        help="email address for querying Entrez web API")

    options = parser.parse_args()

    if options.email is None:
        Entrez.email = settings.entrezemail
    else:
        Entrez.email = options.email

    headers, seqs = sca.readAlg(options.alignment)
    gis = [h.split(options.delim)[1] for h in headers]

    # Check that the GI numbers are valid.
    good_idx = []
    for i, gi in enumerate(gis):
        if not gi.isdigit():
            print("Bad GI '%s' at line %s. Omitting." % (gi, i))
        else:
            good_idx.append(i)

    # Filter out headers and sequeneces with invalid GI numbers.
    headers = [headers[idx] for idx in good_idx]
    seqs = [seqs[idx] for idx in good_idx]
    gis = [h.split(options.delim)[1] for h in headers]

    gi_blocksize = 200  # more GIs need to be submitted as a POST request
    gi_blocks = [gis[x:x + gi_blocksize]
                 for x in range(0, len(gis), gi_blocksize)]

    # Query the Entrez web API with GI numbers and store the retured accession
    # numbers in an array.
    acc_ids = []
    for gi_block in gi_blocks:
        handle = Entrez.efetch(db="protein", rettype="acc", id=gi_block)
        res = handle.readlines()
        handle.close()
        if len(res) == len(gi_block):
            acc_ids.extend([acc_id.strip() for acc_id in res])
        else:
            sys.exit("ERROR: Entrez returned different number of Accession IDs.")

    # Replace GI field with accession numbers in the headers and write the
    # updated alignment to disk.
    f = open(options.outputfile, 'w')
    for i, header in enumerate(headers):
        fields = header.split(options.delim)
        fields[0] = 'res'
        fields[1] = acc_ids[i]
        f.write('>%s\n' % (options.delim).join(fields))
        f.write('%s\n' % seqs[i])
    print("Done. Output written to %s." % options.outputfile)
