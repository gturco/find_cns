import sys
import os
import os.path as op
import numpy as np
import commands
from shapely.geometry import Point, Polygon, LineString, MultiLineString
from flatfeature import Bed

from processing import Pool
pool = None


def get_pair(pair_file, fmt, qbed, sbed, seen={}):
    """ read a line and make sure it's unique handles
    dag, cluster, and pair formats."""
    skipped = open('data/skipped.txt', 'w')
    fh = open(pair_file)
    for line in open(pair_file):
        if line[0] == "#": continue
        line = line.strip().split("\t")
        if fmt == 'dag':
            assert len(line) > 5, line
            pair = line[1], line[5]
        elif fmt in ('cluster', 'qa', 'raw'):
            assert len(line) == 5, line
            pair = line[1], line[3]
        elif fmt == 'pair':
            if len(line) == 1:
                line = line.split(",")
            assert len(line) >= 2, "dont know how to handle %s" % line
            pair = line[0], line[1]

        if fmt in ('qa', 'raw'):
            pair = int(pair[0]), int(pair[1])
        pair = tuple(pair)
        if pair in seen:
            continue
        seen[pair] = True
        try:
            if isinstance(pair[0], (int, long)):
                yield qbed[pair[0]], sbed[pair[1]]
            else:
                yield qbed.d[pair[0]], sbed.d[pair[1]]
        except KeyError, IndexError:
            print >>skipped, "%s\t%s" % pair
            print >>sys.stderr, "skipped %s %s" % pair
            continue



def main(qbed, sbed, pairs_file, pair_fmt, ncpu=8):
    """main runner for finding cnss"""
    pool = Pool(options.ncpu)
    outfile = sys.stdout
    
    pairs = [True]
    _get_pair_gen = get_pair(pairs_file, pair_fmt, qbed, sbed)
    # need this for parallization stuff.
    def get_pair_gen():
        try: return _get_pair_gen.next()
        except StopIteration: return None

    while any(pairs):
        pairs = [get_pair_gen() for i in range(ncpu)]
        print >> outfile, pairs
        # writw pairs to filw or mysql here
        # this helps in parallelizing.
    return None

if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser("usage: %prog [options] ")
    parser.add_option("-n", dest="ncpu", help="parallelize to this many cores", type='int', default=8)
    parser.add_option("-q", dest="qfasta", help="path to genomic query fasta")
    parser.add_option("--qbed", dest="qbed", help="query bed file")
    parser.add_option("-s", dest="sfasta", help="path to genomic subject fasta")
    parser.add_option("--sbed", dest="sbed", help="subject bed file")
    parser.add_option("-p", dest="pairs", help="the pairs file. output from dagchainer")
    choices = ("dag", "cluster", "pair", 'qa', 'raw')
    parser.add_option("--pair_fmt", dest="pair_fmt", default='raw',
                      help="format of the pairs, one of: %s" % str(choices),
                      choices=choices)
    (options, _) = parser.parse_args()


    if not (options.qfasta and options.sfasta and options.sbed and options.qbed):
        sys.exit(parser.print_help())

    qbed = Bed(options.qbed, options.qfasta); qbed.fill_dict()
    sbed = Bed(options.sbed, options.sfasta); sbed.fill_dict()

    main(qbed, sbed, options.pairs, options.pair_fmt, options.ncpu)
