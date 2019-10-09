#! /usr/bin/env python3
"""
A script to parse GI numbers from the headers of an alignment with typical
Blast formatting.

**Arguments**
    Input_MSA.fasta (the alignment to be processed)

**Keyword Arguments**
    --output             output file name, default: FilteredAln.fa

:By: Kim Reynolds
:On: 6.5.2015

Copyright (C) 2015 Olivier Rivoire, Rama Ranganathan, Kimberly Reynolds

This program is free software distributed under the BSD 3-clause
license, please see the file LICENSE for details.
"""

import argparse
import scaTools as sca

if __name__ == '__main__':
    # Parse inputs
    parser = argparse.ArgumentParser()
    parser.add_argument("alignment", help='Input Sequence Alignment')
    parser.add_argument("-o", "--output", dest="outputfile", default='GI_Num',
                        help="specify an outputfile name")
    parser.add_argument("-d", "--delim", dest="delim", default="_",
                        help="specify the field delimiter in the header")
    options = parser.parse_args()

    headers, seqs = sca.readAlg(options.alignment)

    # Get index of GI number in the header fields.
    try:
        gi_idx = (headers[0].split(options.delim)).index('gi') + 1
    except BaseException as e:
        print("ERROR: %s" % e)
        sys.exit("GI field not found in %s." % options.alignment)

    gis = [h.split(options.delim)[gi_idx] for h in headers]
    bad_gis = [gi for gi in gis if not gi.isnumeric()]
    for bad_gi in bad_gis:
        print("Omitting '%s': non-numeric GI." % bad_gi)
    gis = [gi for gi in gis if gi not in bad_gis]

    f = open(options.outputfile, 'w')
    for gi in gis:
        f.write('%s\n' % gi)
    f.close()
