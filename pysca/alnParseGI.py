#! /usr/bin/env python3
"""
A script to parse GI numbers from the headers of an alignment with typical
Blast formatting.

**Arguments**
    Input_MSA.fasta (the alignment to be processed)

**Keyword Arguments**
    --output             output file name, default: FilteredAln.fa
    --delim              delimiter for fields in the header for each sequence,
                         default: '_'

:By: Kim Reynolds
:On: 6.5.2015

Copyright (C) 2015 Olivier Rivoire, Rama Ranganathan, Kimberly Reynolds

This program is free software distributed under the BSD 3-clause
license, please see the file LICENSE for details.
"""

import argparse
import sys
import scaTools as sca

if __name__ == '__main__':
    # Parse inputs
    parser = argparse.ArgumentParser()
    parser.add_argument("alignment", help='input sequence alignment')
    parser.add_argument("-o", "--output", dest="outputfile",
                        default='input.gi',
                        help="specify an outputfile name")
    parser.add_argument("-d", "--delim", dest="delim", default="_",
                        help="specify the field delimiter in the header")
    options = parser.parse_args()

    # Read in the MSA.
    headers, seqs = sca.readAlg(options.alignment)

    # Get index of GI number in the header fields.
    try:
        gi_idx = (headers[0].split(options.delim)).index('gi') + 1
    except BaseException as e:
        print("ERROR: %s" % e)
        sys.exit("GI field not found in %s." % options.alignment)

    gis = [h.split(options.delim)[gi_idx] for h in headers]

    good_gis = []
    for i, gi in enumerate(gis):
        if not gi.isnumeric():
            print("Warning: non-numeric GI '%s' at line %s." % (gi, i))
            good_gis.append(0)
        else:
            good_gis.append(gi)

    f = open(options.outputfile, 'w')
    for gi in good_gis:
        f.write('%s\n' % gi)
    f.close()
