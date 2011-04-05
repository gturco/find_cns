import collections
from flatfeature import Bed
import sys
import os
import os.path as op
sys.path.insert(0, os.path.dirname(__file__))
from find_cns import get_pair
import pickle 

def cns_id(cns_dict):
    c = cns_dict

    return "|".join(map(str,
           (c['qaccn'], c['qchr'], c['qstart'], c['qstop'], c['qstrand'],
            c['saccn'], c['schr'], c['sstart'], c['sstop'], c['sstrand'])))

def cns_link(cns_dict, qorg, sorg, base=" http://synteny.cnr.berkeley.edu/CoGe/GEvo.pl?prog=blastn&autogo=1&"):
    d = cns_dict.copy()
    d['qorg'] = qorg
    d['sorg'] = sorg
    
    url = "'dsid1=%(qorg)s&chr1=%(qchr)s&x1=%(qstart)s&dr1up=2000&dr1down=2000&\
dsid1=%(sorg)s;chr2=%(schr)s;x2=%(sstart)s;dr2up=2000;dr2down=2000;num_seqs=2;hsp_overlap_limit=0;hsp_size_limit=0'"
    return base + (url % d) 

def get_cns_dict(cnsfile):
    cnss = collections.defaultdict(list)
    for line in open(cnsfile):
        if line[0] == "#":
            continue
        line = line.rstrip().split(",")
        qchr, qaccn, schr, saccn = line[:4]

        cnslocs = map(int, line[4:])
        if len(cnslocs) % 4: raise

        for i in range(0, len(cnslocs), 4):
            key = (qchr, schr, tuple(cnslocs[i:i + 4]))
            cnss[key].append((qaccn, saccn))
    return cnss
 
class CNS(object):
    __slots__ = ("qchr", "schr", "qstart", "qstop", "sstart", "sstop")
    def __init__(self, cnsinfo):
        self.qchr, self.schr, (self.qstart, self.qstop, self.sstart, self.sstop) = cnsinfo

def cns_fmt_dict(cns, qfeat, sfeat):
    # qaccn, qchr, qstart, qstop, qstrand, saccn, schr, sstart, sstop, sstrand, link
    d = dict(qaccn=qfeat['accn'], qchr=qfeat['seqid'],
        qstart=cns.qstart, qstop=cns.qstop, qstrand=qfeat['strand'],
        saccn=sfeat['accn'], schr=sfeat['seqid'],
        sstart=cns.sstart, sstop=cns.sstop,
        sstrand=sfeat['strand'])
    return d


def get_start_stops(feat, cns):
    """
    return a range dependent on the position of the cns
    relative to the feature
    """
    if cns[0] > cns[1]: cns = cns[1], cns[0]
    if feat['start'] < cns[0] and feat['end'] > cns[1]:
        # intronicns cnsns:
        return cns[0], cns[1]
    featm = (feat['start'] + feat['end']) / 2.
    cnsm = (cns[0] + cns[1]) /2.
    if featm < cnsm:
        return min(feat['end'], cns[0]), max(feat['end'], cns[0])
    return sorted([cns[1], feat['start']])
    
    

def get_nearby_features(feat, cns_start_stop, bed):
    p0, p1 = get_start_stops(feat, cns_start_stop)
    inters = bed.get_features_in_region(feat['seqid'], p0, p1)
    return [f for f in inters if f["accn"] != feat["accn"]]


def assign(cnsdict, qbed, sbed, qpair_map, spair_map):

    for cnsinfo, accns in cnsdict.iteritems():
        cns = CNS(cnsinfo)

        def dist_sort(a, b):
            # sort features by distance to a given cns.
            qda = get_start_stops(a[0], (cns.qstart, cns.qstop))
            qdb = get_start_stops(b[0], (cns.qstart, cns.qstop))
            qda = abs(qda[1] - qda[0])
            qdb = abs(qdb[1] - qdb[0])

            sda = get_start_stops(a[1], (cns.sstart, cns.sstop))
            sdb = get_start_stops(b[1], (cns.sstart, cns.sstop))
            sda = abs(sda[1] - sda[0])
            sdb = abs(sdb[1] - sdb[0])
            return cmp(sda + qda, sdb + qdb)


        feats = [] 
        for qaccn, saccn in accns:
            try:
                feats.append((qbed.d[qaccn], sbed.d[saccn]))
            except KeyError:
                print >>sys.stderr, "skipped non top-level features:", qaccn, saccn
                raise
                #continue
        feats.sort(cmp=dist_sort)

        for qfeat, sfeat in feats:
            qint = get_nearby_features(qfeat, (cns.qstart, cns.qstop), qbed)
            sint = get_nearby_features(sfeat, (cns.sstart, cns.sstop), sbed)

            if len(qint) + len(sint) > 3: continue

            nqretained = sum(1 for q in qint if q['accn'] in qpair_map)
            nsretained = sum(1 for s in sint if s['accn'] in spair_map)
            
            # look for intervening retained
            # a-b-CNS
            # a-z-CNS-b
            # so it was assinged to a, but b intervenes.
            if nqretained + nsretained > 0: continue

            yield cns, qfeat, sfeat
            # the first 1 had to be the closest...
            break
        
        
def main(cnsfile, qbed_file, sbed_file, pairsfile, pairs_fmt, qorg, sorg):
    qcns_file = qbed_file.replace(".bed", "_cns.gff")
    assert qcns_file != qbed_file
    qcns_gff = open(qcns_file, 'w')
    print >>qcns_gff, "##gff-version 3"
    if sbed_file != qbed_file:
        scns_file = sbed_file.replace(".bed", "_cns.gff")
        assert scns_file != sbed_file
        scns_gff = open(scns_file, 'w')
        print >>scns_gff, "##gff-version 3"
    else:
        scns_gff = qcns_gff
    
    qbed = Bed(qbed_file); qbed.fill_dict()
    sbed = Bed(sbed_file); sbed.fill_dict()


    cnsdict = get_cns_dict(cnsfile)
    qpair_map, spair_map = make_pair_maps(pairsfile, pairs_fmt, qbed, sbed)
    out = sys.stdout

    fmt = "%(cns_id)s,%(qaccn)s,%(qchr)s,%(qstart)i,%(qstop)i,%(qstrand)s," + \
                       "%(saccn)s,%(schr)s,%(sstart)i,%(sstop)i,%(sstrand)s,%(link)s"

    print >>out, "#" + fmt.replace("%(","").replace(")s","").replace(")i","")
    for cns, qfeat, sfeat in assign(cnsdict, qbed, sbed, qpair_map, spair_map):
        d = cns_fmt_dict(cns, qfeat, sfeat)
        d['cns_id'] = cns_id(d)
        if d['sstop'] < d['sstart']:
            d['sstop'], d['sstart'] = d['sstart'], d['sstop']
        d['link'] = cns_link(d, qorg, sorg)
        print >>out, fmt % d
        write_gff(d, qcns_gff, scns_gff)

def write_gff(d, qcns_gff, scns_gff):

    qfmt = \
        "%(qchr)s\t.\tcns\t%(qstart)i\t%(qstop)i\t.\t%(qstrand)s\t.\tID=q__%(cns_id)s;qaccn=%(qaccn)s;saccn=%(saccn)s"
    sfmt = \
        "%(schr)s\t.\tcns\t%(sstart)i\t%(sstop)i\t.\t%(sstrand)s\t.\tID=s__%(cns_id)s;qaccn=%(qaccn)s;saccn=%(saccn)s"

    print >>qcns_gff, qfmt %d
    print >>scns_gff, sfmt %d

def make_pair_maps(pair_file, fmt, qbed, sbed):
    """
    make dicts of q => s and s => q
    """
    qmap = collections.defaultdict(list) # key is query, value is a list of subject hits
    smap = collections.defaultdict(list)
    print >>sys.stderr, "pair file:", pair_file
    for pair in get_pair(pair_file, fmt, qbed, sbed):
        if pair is None: break
        (qname, sname) = pair
        qmap[qname].append(sname)
        smap[sname].append(qname)
    return qmap, smap



DEBUG=False

if __name__ == "__main__":

    import optparse
    parser = optparse.OptionParser()
    parser.add_option("--qbed", dest="qbed", help="bed file of the query")
    parser.add_option("--sbed", dest="sbed", help="bed file of the subject")
    parser.add_option("--cns", dest="cns", help="path to the cns file created by find_cns.py")
    parser.add_option("--pairs", dest="pairs", help="path pairs file")
    choices = ("dag", "cluster", "pair", "qa")
    parser.add_option("--pair_fmt", dest="pair_fmt", default='dag',
                      help="format of the pairs, one of: %s" % str(choices),
                      choices=choices)
    parser.add_option("--qorg", dest="qorg", type="int", help="dsid number in coge for the query (same as export to bed input)")
    parser.add_option("--sorg", dest="sorg", type="int", help="dsid number in coge for the subject (same as export to bed input)")

    (options, _) = parser.parse_args()

    if not (options.qbed and options.sbed and options.cns and options.pairs):
        sys.exit(parser.print_help())

    res = main(options.cns, options.qbed, options.sbed,  options.pairs, options.pair_fmt, options.qorg, options.sorg )


