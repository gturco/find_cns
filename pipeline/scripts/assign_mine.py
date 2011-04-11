from flatfeature import Bed
import collections
import numpy
from numpy import array
import sys
import os
import os.path as op
sys.path.insert(0, os.path.dirname(__file__))
from find_cns import get_pair


def assign_url(qcns, qseqid, scns, sseqid, sfeat, pairsfile, sbed, qbed,
               base = "http://synteny.cnr.berkeley.edu/CoGe/GEvo.pl?prog=blastn&autogo=1&"):
    "lines up coge based on the cns postion"
    params = {'qcns' : scns, 'qseqid' : qseqid , 'scns' : qcns , 'sseqid' : sseqid}
    params['sfeat'] = get_homeolog(sfeat,pairsfile, sbed, qbed)
    inside = 'dsid1=43388&dsgid1=9109&chr1=%(qseqid)s&x1=%(scns)s&dr1up=15000&dr1down=15000&dsid2=43388&dsgid2=9109&chr2=%(sseqid)s&x2=%(qcns)s&dr2up=15000;dr2down=15000&\
accn3=%(sfeat)s;dsid3=34580;dsgid3=34580;dr3up=15000;dr3down=15000;num_seqs=3;hsp_overlap_limit=0;hsp_size_limit=0' %params
    url = base + inside
    return url

def get_homeolog(sfeat,pairsfile, sbed, qbed):
    for region, sregion in get_pair(pairsfile, 'pck', sbed, qbed):
        if region['sfeat'] == sfeat[3]:
            return region['ORG2_qfeat']
           
def get_cns_dict(cnsfile):
    cnss = collections.defaultdict(list)
    for line in open(cnsfile):
        if line[0] == "#":continue
        line = line.rstrip().split(",")
        accn, qaccn_l, qaccn_r, qchr, saccn, schr = line[0:6]
        cnslocs = map(int, line[6:])
        if len(cnslocs) % 4: raise
        
        for i in range(0, len(cnslocs), 4):
            key = (qchr, schr, tuple(cnslocs[i:i + 4]))
            cnss[key].append((accn, qaccn_l, qaccn_r ,saccn))
    return cnss
    
    
def cns_fmt_dict(cns, sfeat, accn, qaccn_l, qaccn_r):
    # accn, qaccnL, qaccnR, qchr, qstart, qstop,saccn, schr, sstart, sstop,link
    d = dict(accn=accn, qaccnL=qaccn_l,qaccnR=qaccn_r,qchr=cns.qchr,
        qstart=cns.qstart, qstop=cns.qstop,
        saccn=sfeat['accn'], schr=sfeat['seqid'],
        sstart=cns.sstart, sstop=cns.sstop)
    return d
    
 
class CNS(object):
    __slots__ = ("qchr", "schr", "qstart", "qstop", "sstart", "sstop")
    def __init__(self, cnsinfo):
        self.qchr, self.schr, (self.qstart, self.qstop, self.sstart, self.sstop) = cnsinfo  
        
def make_pair_maps(cnsfile):
    """make dicts of q => s and s => q"""
    smap = []
    for line in open(cnsfile):
        if line[0] == "#":continue
        line = line.rstrip().split(",")
        sname = line[4]
        smap.append(sname)
    return smap 
        
def get_nearby_features(feat, bed, p0, p1):
    "grabs all features inbtween two postions (chr#, start, stop)"
    if p0 < p1:
        inters = bed.get_features_in_region(feat['seqid'], p0, p1)
    else:
        inters = bed.get_features_in_region(feat['seqid'], p1, p0)
    return [f for f in inters if f["accn"] != feat["accn"]]   
 
def nearest_feat(feats, cns_start):
    "imputs:a list of feature start and stop postions and a cns start or stop postion\
     output: the feature start or stop postion that is closest to the cbs "
    dist = array([[abs(x - cns_start),abs(y - cns_start)] for x,y in feats])
    dist_min = numpy.unravel_index(dist.argmin(), dist.shape)
    nearest_feat =  feats[dist_min[0]]
    nearest_feat_start_stop = nearest_feat[dist_min[1]]
    return nearest_feat_start_stop
    

def assign(cnsdict, sbed, spair_map):
    "returns cns only if they are the closet of the spair_map"
    
    for cnsinfo, accns in cnsdict.iteritems():
        cns = CNS(cnsinfo)
        for accn, qaccn_l, qaccn_r, saccn in accns:
            sfeat = sbed.d[saccn]
            feats = [locs for locs in sfeat['locs']]
            p1 = nearest_feat(feats, cns.qstart)
            genes_inbetween = get_nearby_features(sfeat,sbed, p1, cns.sstart)
            nsretained = sum(1 for gene in genes_inbetween if gene['accn'] in spair_map)
            if nsretained > 0 : continue
            yield cns, accn, qaccn_l, qaccn_r, sfeat

        
            
def main(cnsfile, sbed_file, pairsfile, qorg, sorg):
    sbed = Bed(sbed_file); sbed.fill_dict()
    qbed = sbed
    cnsdict = get_cns_dict(cnsfile)
    spair_map = make_pair_maps(cnsfile)
    out = sys.stdout
    
    fmt = "%(accn)s,%(qaccnL)s,%(qaccnR)s,%(qchr)s,%(qstart)i,%(qstop)i," + \
                     "%(saccn)s,%(schr)s,%(sstart)i,%(sstop)i,%(link)s" 
                     
    print >>out, "#" + fmt.replace("%(","").replace(")s","").replace(")i","")
    for cns, accn, qaccn_l, qaccn_r, sfeat in assign(cnsdict, sbed, spair_map): 
        d = cns_fmt_dict(cns, sfeat, accn, qaccn_l, qaccn_r)
        d['link'] = assign_url(cns.qstart, cns.qchr, cns.sstart, cns.schr,sfeat, pairsfile, sbed, qbed)
        print cns.schr, cns.sstart
        print >>out, fmt % d
        
        
# DEBUG=False
# 
# if __name__ == "__main__":
# 
#     import optparse
#     parser = optparse.OptionParser()
#     parser.add_option("--qbed", dest="qbed", help="bed file of the query")
#     parser.add_option("--sbed", dest="sbed", help="bed file of the subject")
#     parser.add_option("--cns", dest="cns", help="path to the cns file created by find_cns.py")
#     parser.add_option("--pairs", dest="pairs", help="path pairs file")
#     choices = ("dag", "cluster", "pair", "qa")
#     parser.add_option("--pair_fmt", dest="pair_fmt", default='dag',
#                       help="format of the pairs, one of: %s" % str(choices),
#                       choices=choices)
#     parser.add_option("--qorg", dest="qorg", type="int", help="dsid number in coge for the query (same as export to bed input)")
#     parser.add_option("--sorg", dest="sorg", type="int", help="dsid number in coge for the subject (same as export to bed input)")
# 
#     (options, _) = parser.parse_args()
# 
#     res = main(options.cns, options.sbed,  options.pairs, options.qorg, options.sorg )

            
main('/Users/gturco/rice_v6_cns_res/04_08_10/test_mine/testassign.txt', '/Users/gturco/rice_v6_cns_res/04_08_10/test_org/rice_v6.bed' , '/Users/gturco/rice_v6_cns_res/04_08_10/test_mine/rice_v6_rice_v6.pairs.pck', 'qorg', 'sorg')