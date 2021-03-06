from flatfeature import Bed
import collections
import numpy
from numpy import array
import sys
sys.path.append("scripts/")
import os
import os.path as op
sys.path.insert(0, os.path.dirname(__file__))
from find_cns import get_pair
from assign_org import get_cns_dict, cns_fmt_dict
from assign_region import nearest_feat

# the assign used for maize, assigns only based on qfeat

def assign_url(scns, sseqid, qcns, qseqid, sorg, qorg, padding,
               base = "http://coge.iplantcollaborative.org/CoGe/GEvo.pl?prog=blastn&autogo=1&show_cns=1&"):
    "lines up coge based on the cns postion"
    inside = "dsgid1={4}&chr1={0}&x1={1}&dr1up={6}&dr1down={6}&dsgid2={5}&chr2={2}&x2={3}&dr2up={6};dr2down={6}&"\
             "num_seqs=2;hsp_overlap_limit=0;hsp_size_limit=0".format(sseqid, scns, qseqid, qcns, sorg, qorg, padding)
    url = base + inside
    return url

class CNS(object):
    __slots__ = ("qchr", "schr", "qstart", "qstop", "sstart", "sstop")
    def __init__(self, cnsinfo):
        self.qchr, self.schr, (self.qstart, self.qstop, self.sstart, self.sstop) = cnsinfo  


def assign(cnsdict, qbed, sbed):
    "finds the nearest qfeat for each duplicat cns qstart qstop, sstart, sstop pos"
    
    for cnsinfo, accns in cnsdict.iteritems():
        cns = CNS(cnsinfo) 
        qfeats = []
        for qaccn, saccn in accns:
            try:
                qfeats.append((qbed.d[qaccn],sbed.d[saccn]))
            except KeyError:
                print >>sys.stderr, "skipped non top-level features:", qaccn , saccn
                raise
        qfeat_start_stops = [(qfeat['start'], qfeat['end']) for qfeat,sfeat in qfeats]
        pos = nearest_feat(qfeat_start_stops, cns.qstart, cns.qstop)
        qfeat,sfeat = qfeats[pos]
        yield cns, qfeat, sfeat

        
            
def main(cnsfile, qbed_file, sbed_file, qorg, sorg, padding):
    qbed = Bed(qbed_file); qbed.fill_dict()
    sbed = Bed(sbed_file); sbed.fill_dict()
    cnsdict = get_cns_dict(cnsfile)
    out = sys.stdout
    
    fmt = "%(qaccn)s,%(qchr)s,%(qstart)i,%(qstop)i,%(qstrand)s," + \
                       "%(saccn)s,%(schr)s,%(sstart)i,%(sstop)i,%(sstrand)s,%(link)s"
                     
    print >>out, "#" + fmt.replace("%(","").replace(")s","").replace(")i","")
    for cns, qfeat, sfeat in assign(cnsdict, qbed, sbed): 
        d = cns_fmt_dict(cns, qfeat, sfeat)
	if d['sstop'] < d['sstart']:
            d['sstop'], d['sstart'] = d['sstart'], d['sstop']        
	d['link'] = assign_url(cns.sstart, cns.schr, cns.qstart, cns.qchr, sorg, qorg, padding)
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
    (options, _) = parser.parse_args()

    res = main(options.cns, options.qbed, options.sbed, options.qorg, options.sorg, options.pad)

            
#main('/Users/gturco/code/freeling_lab/find_cns_gturco/pipeline/tests/resources/cns_2.csv', '/Users/gturco/rice_maize/rice_v6.bed', '/Users/gturco/rice_maize/maize_v2.bed', '9109', '9109', 1000)
