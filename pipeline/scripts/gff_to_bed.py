import re
from collections import defaultdict,deque
import sys
sys.path.append("coann/")
from random_noncoding_seq import recursive_merge_both

def parse_accn(regex,accn_col):
    """ if given regex expression will look for otherwise uses first name"""
    m = re.search(regex,accn_col)
    names = accn_col.split(";")
    pid_all = names[0].split("=")[1]
    pid = re.split(".transposable_element|.transposon|.gene|.pseudogene|.mRNA|.CDS",pid_all)[0]
    try:
        accn = m.group(0)
    except AttributeError:
        accn = pid
        #l = re.search("Os\d\dg\d{5}",m)
    return accn,pid


def merge_feats(feats):
    merge_same_feats = set(feats)
    feats = list(merge_same_feats)
    feats.sort()
    if len(feats) > 1:
        merge_overlapping_feats = recursive_merge_both(deque(feats))
        merge_overlapping_feats.sort()
        feats = merge_overlapping_feats

    return feats



def parse_gff(regex,gff):
    parents = {}
    feats = defaultdict(list)
    for line in open(gff):
        if line[0] == "#": continue
        args = line.strip().split("\t")
        seqid = args[0]
        start = int(args[3])
        end = int(args[4])
        strand = args[6]
        accn,pid = parse_accn(regex,args[8])
        ftype = args[2]
        if ftype == "CDS":
            feats[pid].append((start,end))
        elif ftype in ["transposable_element","gene","transposon","pseudogene"]:
            values = [seqid,start,end,strand,accn]
            parents[pid] = values
    return parents, feats

def gff_to_bed(parents,feats,bed):
    bed_fh = open(bed,"wb")
    for pid in parents.keys():
        seqid,start,end,strand,accn = parents[pid]
        largest_feats = merge_feats(feats[pid])
        feat_starts = [str(s-start) for (s,e) in largest_feats]
        feat_lens = [str((e-s)+1) for (s,e) in largest_feats]
        bed_line = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t.\t.\t.\t{6}\t{7}\t{8}\n".format(seqid,start-1,end,accn,(end-start),
                strand,len(largest_feats),",".join(feat_lens),",".join(feat_starts))
        bed_fh.write(bed_line)

def main(regex,gff,bed):
    parents,feats = parse_gff(regex,gff)
    gff_to_bed(parents,feats,bed)

if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser("usage: %prog [options] ")
    parser.add_option("--re", dest="regex", default="This_will_not_be_found",  help="regex expression for gene name")
    parser.add_option("--gff", dest="gff", help="gff file from CoGe")
    parser.add_option("--bed", dest="bed", help="name of bed file output")
    (options, _) = parser.parse_args()

    main(options.regex,options.gff,options.bed)

