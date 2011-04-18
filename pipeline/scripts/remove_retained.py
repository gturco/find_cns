from flatfeature import Bed
import collections
import numpy
from numpy import array
import sys
import os
import os.path as op
sys.path.insert(0, os.path.dirname(__file__))

def get_pair(pair_file, fmt, qbed, sbed, seen={}):
    """ read a line and make sure it's unique handles
    dag, cluster, and pair formats."""
    #skipped = open('data/skipped.txt', 'w')
    fh = open(pair_file)
    if fmt == 'pck':
        pck_file = pickle.load(fh)
        for row in pck_file:
            region = row
            accn = row['sfeat']
            qfeat = qbed.accn(accn)
            pair = qfeat, region
            yield pair
    else:        
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
                    yield pair
            except KeyError, IndexError:
                print >>skipped, "%s\t%s" % pair
                continue

           
def get_cns_dict(cnsfile):
    cnss = collections.defaultdict(list)
    for line in open(cnsfile):
        if line[0] == "#":continue
        line = line.rstrip().split(",")
        qaccn, qchr, saccn, saccn_l, saccn_r, schr = line[0:6]
        cnslocs = map(int, line[6:])
        if len(cnslocs) % 4: raise

        for i in range(0, len(cnslocs), 4):
            key = (qchr, schr, tuple(cnslocs[i:i + 4]))
            cnss[key].append((qaccn, saccn, saccn_l, saccn_r))
    return cnss
    

def cns_fmt_dict(cns, qfeat, saccn, saccn_l, saccn_r):
    # saccn, schr, qaccn, qaccnL, qaccnR, qstart, qstop,saccn, schr, sstart, sstop,link
    d = dict(accn=saccn, saccnL=saccn_l,saccnR=saccn_r,schr=cns.schr,
        sstart=cns.sstart, sstop=cns.sstop,
        qaccn=qfeat['accn'], qchr=qfeat['seqid'],
        qstart=cns.qstart, qstop=cns.qstop)
    return d
    
 
class CNS(object):
    __slots__ = ("qchr", "schr", "qstart", "qstop", "sstart", "sstop")
    def __init__(self, cnsinfo):
        self.qchr, self.schr, (self.qstart, self.qstop, self.sstart, self.sstop) = cnsinfo  
        
def make_pair_maps(pair_file, fmt, qbed, sbed):
    """make dicts of q => s and s => q"""
    qmap_tuple = []
    for pair in get_pair(pair_file, fmt, qbed, sbed):
        if pair is None: break
        (qname, sname) = pair
        smap_tuple.append((sname,qname))
        smap_tuple.append((qname,sname))
    return qmap_tuple
        

 
def nearest_feat(feats, cns_start):
    "imputs:a list of feature start and stop postions and a cns start or stop postion\
     output: the feature start or stop postion that is closest to the cbs "
    dist = array([[abs(x - cns_start),abs(y - cns_start)] for x,y in feats])
    dist_min = numpy.unravel_index(dist.argmin(), dist.shape)
    nearest_feat =  feats[dist_min[0]]
    nearest_feat_start_stop = nearest_feat[dist_min[1]]
    return dist, dist_min
    
    
def same_chr_feat(feat_list, qbed, cns):
    f_list = []
    for f in feat_list:
        feat_bed = qbed.accn(f)
        if feat_bed['seqid'] == cns.qchr:
            f_list.append(f)
    return f_list[0]
    

def assign(cnsdict, qbed, pairsfile, sbed, qpair_map):
    "returns cns only if they are the closet of the spair_map"
    
    for cnsinfo, accns in cnsdict.iteritems():
        cns = CNS(cnsinfo)
        for qaccn, saccn, saccn_l, saccn_r in accns:
            left_retained = [l for (k,l) in qpair_map if saccn_r[:-1] in k]
            left_retained_one = same_chr_feat(left_retained, qbed, cns)
            right_retained = [l for (k,l) in qpair_map if saccn_l[1:] in k]
            right_retained_one = same_chr_feat(right_retained, qbed, cns)
            left_feat,right_feat = qbed.accn(left_retained_one), qbed.accn(right_retained_one)
            qfeat1 = qbed.accn(qaccn)
            h_feats = [left_feat['accn'], right_feat['accn'], qfeat1['accn']]
            homeolog_feats = [(left_feat['start'],left_feat['end']), (right_feat['start'],right_feat['end']), (qfeat1['start'], qfeat['end'])]
            dist_array, dist_min = nearest_feat(homeolog_feats, cns.qstart)
            if dist_min[0] < 2: continue
            yield cns, saccn, saccn_l, saccn_r, qfeat

        
            
def main(cnsfile, qbed_file, pairsfile):
    qbed = Bed(qbed_file); qbed.fill_dict()
    sbed = qbed
    cnsdict = get_cns_dict(cnsfile)
    qpair_map = make_pair_maps(pairsfile, 'pair', qbed, sbed) # map of all pairs that q and sbed in file
    out = sys.stdout
    
    fmt = "%(accn)s,%(saccnL)s,%(saccnR)s,%(schr)s,%(qaccn)s,%(qchr)s," + \
                     "%(sstart)i,%(sstop)i,%(qstart)i,%(qstop)i" 
                     
    print >>out, "#" + fmt.replace("%(","").replace(")s","").replace(")i","")
    for cns, saccn, saccn_l, saccn_r, qfeat in assign(cnsdict, qbed, pairsfile, sbed, qpair_map): 
        d = cns_fmt_dict(cns, qfeat, saccn, saccn_l, saccn_r)
        print >>out, fmt % d

        
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
    (options, _) = parser.parse_args()

    res = main(options.cns, options.qbed,  options.pairs)

            
#main('/Users/gturco/rice_v6_cns_res/04_08_10/test_mine/testassign.txt', '/Users/gturco/rice_v6_cns_res/04_08_10/test_org/rice_v6.bed' , '/Users/gturco/rice_v6_cns_res/04_08_10/test_org/rice_v6_rice_v6.pairs.txt')