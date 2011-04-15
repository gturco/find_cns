from flatfeature import Bed
import collections
import numpy
from numpy import array
import sys
import os
import os.path as op
sys.path.insert(0, os.path.dirname(__file__))
from find_cns import get_pair

def assign_url(scns, sseqid, qcns, qseqid, qfeat, pairsfile, sbed, qbed, sorg, qorg, padding,
               base = "http://synteny.cnr.berkeley.edu/CoGe/GEvo.pl?prog=blastn&autogo=1&"):
    "lines up coge based on the cns postion"
    sfeat = get_homeolog(qfeat,pairsfile, sbed, qbed)
    inside = "dsgid1={5}&chr1={0}&x1={1}&dr1up={7}&dr1down={7}&dsgid2={6}&chr2={2}&x2={3}&dr2up={7};dr2down={7}&"\
             "accn3={4};dr3up=15000;dr3down=15000;num_seqs=3;hsp_overlap_limit=0;hsp_size_limit=0".format(sseqid, scns, qseqid, qcns, sfeat, sorg, qorg, padding)
    url = base + inside
    return url

def get_homeolog(qfeat,pck, sbed, qbed):
    for region, sregion in get_pair(pck, 'pck', sbed, qbed):
        if region['sfeat'] == qfeat[3]:
            return region['ORG2_qfeat']
           
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
    d = dict(saccn=saccn, saccnL=saccn_l,saccnR=saccn_r,schr=cns.schr,
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
    for pair in get_pair(pair_file,fmt, qbed, sbed):
        if pair is None: break
        (sname, qname) = pair
        qmap_tuple.append((qname,sname))
    return qmap_tuple
#         
# def get_nearby_features(feat, bed, p0, p1):
#     "grabs all features inbtween two postions (chr#, start, stop)"
#     if p0 < p1:
#         inters = bed.get_features_in_region(feat['seqid'], p0, p1)
#     else:
#         inters = bed.get_features_in_region(feat['seqid'], p1, p0)
#     return [f for f in inters if f["accn"] != feat["accn"]]   
 


def nearest_feat(feats, cns_start, cns_stop):
    "imputs:a list of feature start and stop postions and a cns start or stop postion\
     output: the feature start or stop postion that is closest to the cbs "
    dist = array([[abs(x - cns_start),abs(y - cns_start)] for x,y in feats])
    params = [[(x , cns_start),(y , cns_start),(x , cns_stop),(y , cns_stop)] for x,y in feats]
    dist_min = numpy.unravel_index(dist.argmin(), dist.shape)
    #params_row = params[dist_min[0]]
    #p0, p1 =  feats[dist_min[1]]
    return dist_min[0]
    
def same_chr_feat(feat_list, qbed, cns):
    f_list = []
    for f in feat_list:
        feat_bed = qbed.accn(f)
        if feat_bed['seqid'] == cns.qchr:
            f_list.append(f)    
    return f_list[0]
    

def assign(cnsdict, qbed, qpair_map):
    "finds the nearest qfeat for each duplicat cns qstart qstop, sstart, sstop pos"
    
    for cnsinfo, accns in cnsdict.iteritems():
        cns = CNS(cnsinfo)
        qfeats = []
        for qaccn, saccn, saccn_l, saccn_r in accns:
            left_retained = [k for (k,l) in qpair_map if saccn_r[:-1] in k]
            left_retained_one = same_chr_feat(left_retained, qbed, cns)
            right_retained = [k for (k,l) in qpair_map if saccn_l[1:] in k]
            right_retained_one = same_chr_feat(right_retained, qbed, cns)
            left_feat,right_feat = qbed.accn(left_retained_one), qbed.accn(right_retained_one)
            qfeat1 = qbed.accn(qaccn)
            h_feats = [left_feat['accn'], right_feat['accn'], qfeat1['accn']]
            homeolog_feats = [(left_feat['start'],left_feat['end']), (right_feat['start'],right_feat['end']), (qfeat1['start'], qfeat['end'])]
            dist_min = nearest_feat(homeolog_feats, cns.qstart, cns.qstop)
            if dist_min < 2: continue
            try:
                qfeats.append(qbed.d[qaccn],saccn,saccn_l,saccn_r)
            except KeyError:
                print >>sys.stderr, "skipped non top-level features:", qaccn , saccn
                raise
        qfeat_start_stops = [(qfeat['start'], qfeat['end']) for qfeat,saccn,saccn_l,saccn_r in qfeats]
        pos = nearest_feat(qfeat_start_stops, cns.qstart, cns.qstop)
        qfeat,saccn,saccn_l,saccn_r = qfeats[pos]
    # genes_inbetween = get_nearby_features(qfeat,qbed, p0, p1) # these are the gene inbtween the nearest feat
    #         nsretained = sum(1 for gene in genes_inbetween if gene['accn'] in qpair_map)
    #         if nsretained > 0 : continue # if a gene inbetween the two is in the qpairmap then remove that cns from that gene
        yield cns, saccn, saccn_l, saccn_r, qfeat

        
            
def main(cnsfile, qbed_file, sbed_file, pairsfile, pck, qorg, sorg, padding):
    qbed = Bed(qbed_file); qbed.fill_dict()
    sbed = Bed(sbed_file); sbed.fill_dict()
    cnsdict = get_cns_dict(cnsfile)
    qpair_map = make_pair_maps(pairsfile, 'pair', qbed, sbed)
    out = sys.stdout
    
    fmt = "%(saccn)s,%(saccnL)s,%(saccnR)s,%(schr)s,%(sstart)i,%(sstop)i," + \
                     "%(qaccn)s,%(qchr)s,%(qstart)i,%(qstop)i,%(link)s" 
                     
    print >>out, "#" + fmt.replace("%(","").replace(")s","").replace(")i","")
    for cns, saccn, saccn_l, saccn_r, qfeat in assign(cnsdict, qbed, qpair_map): 
        d = cns_fmt_dict(cns, qfeat, saccn, saccn_l, saccn_r)
        d['link'] = assign_url(cns.sstart, cns.schr, cns.qstart, cns.qchr,qfeat, pck, sbed, qbed, sorg, qorg, padding)
        print >>out, fmt % d
        
        
DEBUG=False

if __name__ == "__main__":

    import optparse
    parser = optparse.OptionParser()
    parser.add_option("--qbed", dest="qbed", help="bed file of the query")
    parser.add_option("--sbed", dest="sbed", help="bed file of the subject")
    parser.add_option("--cns", dest="cns", help="path to the cns file created by find_cns.py")
    parser.add_option("--pairs", dest="pairs", help="path pairs file")
    choices = ("dag", "cluster", "pair", "qa", "pck")
    parser.add_option("--pair_fmt", dest="pair_fmt", default='dag',
                      help="format of the pairs, one of: %s" % str(choices),
                      choices=choices)
    parser.add_option("--qorg", dest="qorg", type="int", help="dsid number in coge for the query (same as export to bed input)")
    parser.add_option("--sorg", dest="sorg", type="int", help="dsid number in coge for the subject (same as export to bed input)")
    parser.add_option("--pad", dest="pad", type='int', default=12000,
                      help="how far from the end of each gene to set the padding for the link of the cnss")
    parser.add_option("--pck", dest="pck", help="pairs_file for pck/cns/region")
    (options, _) = parser.parse_args()

    res = main(options.cns, options.qbed, options.sbed, options.pairs, options.pck, options.qorg, options.sorg, options.pad)

            
#main('/Users/gturco/rice_v6_cns_res/04_08_10/test_mine/testassign.txt', '/Users/gturco/rice_v6_cns_res/04_08_10/test_org/rice_v6.bed', '/Users/gturco/rice_v6_cns_res/04_08_10/test_org/rice_v6.bed' , '/Users/gturco/rice_v6_cns_res/04_08_10/test_mine/rice_v6_rice_v6.pairs.pck', '9109', '9109', 1000)