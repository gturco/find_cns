#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
gff_to_bed some code from <https://github.com/tanghaibao/quota-alignment/blob/master/scripts/gff_to_bed.py>
Author: Haibao Tang
"""

"""
%prog gff_file

convert the gff format into the bed format 
"""
#TODO change names to regex and change cds for splicing with how coge does it

import os.path as op
import sys

try:
    from BCBio.GFF import GFFParser
except:
    print >>sys.stderr, "BCBio.GFF module not found, try download and install from`" \
            "<http://github.com/chapmanb/bcbb/tree/master/gff/BCBio/>"
    sys.exit(1)


def gff_to_bed(gff_file):

    bed_fh=sys.stdout
    parser = GFFParser()
    seqids = parser.parse(gff_file, None)

    for seqid in seqids:
        for feat in seqid.features:
            subf = feat.sub_features
            
            if feat.type in ("chromosome", "protein"): continue
            strand = "+" if ==else
            is_cds = any(f.type=="mRNA" or f.type=="CDS" for f in subf) and\
                    feat.type=="gene"
            if is_cds:
                cds = [(int(f.location.start),int(f.location.end)) for f in
                        subf[0].sub_features if f.type=="CDS"]
                cds.sort()
                cds_col1 = ','.join([str(end-start) for start,end in cds])
                cds_col2 = ','.join([str(start - feat.location.start) for start,end in cds])
                
                print >>bed_fh, "\t".join(str(x) for x in (seqid.id, feat.location.start, \
                        feat.location.end, feat.id,(feat.location.end -
                            feat.location.start),feat.strand,'.',-1,-1,len(cds),cds_col1,cds_col2))


if __name__ == "__main__":
    import optparse

    parser = optparse.OptionParser(__doc__)
    parser.add_option("--noncoding", dest="cds", action="store_false", 
            default=True, help="extract coding features?")
    (options, args) = parser.parse_args()

    if len(args) != 1:
        sys.exit(parser.print_help())

    gff_file = args[0]

    gff_to_bed(gff_file)

