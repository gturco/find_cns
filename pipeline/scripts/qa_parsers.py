import sys
sys.path.append("../../cns_pipeline/bin/quota-alignment/scripts/")
from cleanup import DupLine
try:
    from bed_utils import Bed
except ImportError:
    raise 'bed_utils.py pacakge not found edit the sys.path.append of qa_parsers.py  to location of your quota-alignment'

class RawBed(object):
    """takes line from habos raw file and converts to brents bed line"""

    def __init__(self,bed):
        self.bed = Bed(bed)
        self.order = self.bed.get_order()

    def raw_to_bed(self,raw_pos):
        """returns the bed file for the raw_pos"""
        bed_info = self.bed[raw_pos]
        d = {}
        d['start'] = bed_info.start
        d['end'] = bed_info.end
        d['seqid'] = bed_info.seqid
        d['accn'] = bed_info.accn
        args = bed_info.stuff
        d['strand'] = args[1]
        #d['locs'] = loc_conv(args[-2],args[-1])
        return d

    def accn_to_raw(self,accn):
        "returns the raw line inputs for accn"
        pos =self.order[accn][0]
        seqid = self.order[accn][1].seqid
        return pos,seqid





class ParsePairs(object):
    """takes line from pairs returns the qaccn and saccn"""
    ##TODO add dag and cluster

    def __init__(self,line):
        self.line = line.strip().split('\t')

    def dag(self):
        assert len(self.line) > 5, self.line
        pairs = self.line[1], self.line[5]
        return pairs

    def cluster(self):
        assert len(self.line) == 5, self.line
        pair = self.line[1], self.line[3]
        return pair

    def qa(self):
        pairs = self.raw()
        return pairs

    def raw(self,qnolocaldups,snolocaldups):
        pairs = self.cluster()
        pairs = map(int,pairs)
        qbed = RawBed(qnolocaldups)##### slow because this parth should be
        ###moved to the init statemnt before the for loop..
        sbed = RawBed(snolocaldups)
        return qbed.raw_to_bed(pairs[0])['accn'],sbed.raw_to_bed(pairs[1])['accn']

    def pair(self):
        if len(self.line) == 1:
            line = self.line[0].split(",")
        assert len(self.line) >= 2, "dont know how to handle %s" % self.line
        pair = self.line[0], self.line[1]
        return pair


def pairs_to_qa(pair_file,fmt,qbed_file,sbed_file,raw_file):
    """takes new localdups file and new pairs file to create qa file"""
    new_qa = open(raw_file,'wb')
    qbed= RawBed(qbed_file)
    sbed = RawBed(sbed_file)
    print >>sys.stderr, "write tandem-filtered bed file {0}".format(raw_file)
    for line in open(pair_file):
            if line[0] == "#" : continue
            pair_line = ParsePairs(line)
            qfeat,sfeat= getattr(pair_line, fmt)()
            qpos,qseqid = qbed.accn_to_raw(qfeat)
            spos,sseqid = sbed.accn_to_raw(sfeat)
            new_line = "{0}\t{1}\t{2}\t{3}\t50\n".format(qseqid,qpos,sseqid,spos)
            new_qa.write(new_line)
    new_qa.close()

def write_nolocaldups(bed_path,localdups_file,out_name):
    bed = Bed(bed_path)
    children = []
    for line in open(localdups_file):
        dups= DupLine(line)
        children += dups.children
    print >>sys.stderr, "write tandem-filtered bed file {0}".format(out_name)
    fh = open(out_name, "w")
    for i, row in enumerate(bed):
        if row['accn'] in children: continue
        print >>fh, row
    fh.close()

#write_nolocaldups("../data/dick_m_tair_10/dick_m.bed","../data/dick_m_tair_10/dick_m.localdups","../data/dick_m_tair_10/dick_m.nolocaldups.bed")
#pairs_to_qa("../data/sorghum_n_setaria_n/sorghum_n_setaria_n.pairs.txt","pair","../data/sorghum_n_setaria_n/sorghum_n.nolocaldups.bed.local","../data/sorghum_n_setaria_n/setaria_n.nolocaldups.bed.local","../data/sorghum_n_setaria_n/sorghum_n_setaria_n.filtered.raw")
