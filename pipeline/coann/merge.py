from flatfeature import Bed
from collections import defaultdict
import sys
import heapq
from random_noncoding_seq import recursive_merge_both
from collections import deque

def parse_missed_genes(missed_genes_path):
    """parses co-anno output: matches.txt tab sep file
    output: 1:1 ratio missing_gene,qaccn_hit"""
    handle = open(missed_genes_path)
    fh = handle.read()
    missed_genes = []
    for line in fh.split("\n")[:-1]:
        missed_gene,qaccns = line.split("\t")
        for qaccn in qaccns.split(','):
            missed_genes.append ((missed_gene,qaccn))
    return missed_genes

def no_intervening_genes(feat,b_feat,bed):
    """retunrs true is there are no intervening genes between feat and b_feat
    NOTE feat < b_feat... sort before hand"""
    if feat[0] == b_feat[0] and feat[4] == b_feat[4]:
        if feat[2] >= b_feat[1]: return True
        ### want to skip non merged feats for now
        feats = bed.get_features_in_region(feat[0],feat[2]+1, b_feat[1])
        strands = [f["strand"] for f in feats]
        if feat[4] not in strands: return True
        elif len(feats) > 0: return False
        else: return True
    else: return False

def update_locs(old_hit,new_hit):
    """changes start and stops to samllest locs is smaller then utr"""
    new_gene = old_hit
    new_gene["locs"] = old_hit["locs"] + new_hit["locs"]
    new_gene["locs"].sort()
    locs_start = new_gene["locs"][0][0]
    locs_end = new_gene["locs"][-1][1]
    new_gene["start"] = heapq.nsmallest(1,[locs_start,new_gene["start"]])[0]
    new_gene["end"] = heapq.nlargest(1,[locs_end,new_gene["end"]])[0]
    return new_gene


###########TODO#################
def near_by_gene():
    pass
    ### will fix in brents code

def merge_feats(feat):
    merge_same_feats = set(feat["locs"])
    feat["locs"] = list(merge_same_feats)
    feat["locs"].sort()
    if len(feat["locs"]) > 1:
        merge_overlapping_feats =recursive_merge_both(deque(feat["locs"]))
        merge_overlapping_feats.sort()
        feat["locs"] = merge_overlapping_feats
        locs_start = feat["locs"][0][0]
        locs_end = feat["locs"][-1][1]
        feat["start"] = heapq.nsmallest(1,[locs_start,feat["start"]])[0]
        feat["end"] = heapq.nlargest(1,[locs_end,feat["end"]])[0]
    return feat


def write_new_bed(gene_list, old_bed, missed_genes,out_file):
    merge_fh = open(out_file,"wb")
    hit_list = [hit for hit,qaccn in missed_genes]
    for i,gene in enumerate(old_bed):
        if gene["accn"] in hit_list: continue
        new_line = Bed.row_string(gene)
        merge_fh.write("{0}\n".format(new_line))
    for i,new_gene in enumerate(gene_list):
        ### merge overlapping here
        updated_feat = gene_list[new_gene]
        if len(updated_feat["locs"]) > 1:
            updated_feat = merge_feats(updated_feat)
        new_line = Bed.row_string(updated_feat)
        merge_fh.write("{0}\n".format(new_line))
        
###sort -n -k 1 -k 2 rice_v6.all2.bed  > rice_v6.all3.bed
################################

def group_genes_in_bed(missed_genes,old_bed,new_bed):
    """ if found in bed append to gene and give gene qaccn:[(chr,start,stop),(chr,start,stop)] otherwise give
    the regulart (chr,start,stop)"""
    missed_genes_grouped = defaultdict(list)
    missed_genes_dict = {}
    for hit_accn, qaccn in missed_genes:
        try:
            ### if with in gene of the old be merge with old bed
            old_hit = old_bed.accn(hit_accn)
            new_hit = new_bed.accn(hit_accn)
            new_gene = update_locs(old_hit,new_hit)
            ##### remove from old bed this removes probblems in merge_hits
            ### add hits to new_bed
            hit = new_gene
            hit_info = (hit["seqid"],hit["start"],hit["end"],hit["accn"],hit["strand"])
            missed_genes_grouped[qaccn].append(hit_info)
            missed_genes_dict[hit['accn']] = hit
        except KeyError:
            try:
                hit = new_bed.accn(hit_accn)
                hit_info = (hit["seqid"],hit["start"],hit["end"],hit["accn"],hit["strand"])
                missed_genes_grouped[qaccn].append(hit_info)
                missed_genes_dict[hit['accn']] = hit
            except KeyError: continue
        #new_new_bed[hit['accn']] = hit.row_to_dict()
    return missed_genes_grouped, missed_genes_dict

def merge_hits(hits,old_bed,missed_genes_dict):
    """sort hits if the hits are on the same chr at a given distance with no
    intervening genes they are joined add new hits to bed"""
    missed_genes_grouped_dict = {}
    hits.sort(key=lambda h: (h[0],h[1]))
    missed_genes_grouped_dict[hits[0][3]] = missed_genes_dict[hits[0][3]]
    for i,hit in enumerate(hits[:-1]):
        b_hit = hits[i+1]
        intervening = no_intervening_genes(hit,b_hit,old_bed)
        if hit[0] == b_hit[0] and (b_hit[1] - hit[2]) <= 7500 and intervening:
            try:
                ### check to see if we already added it to the old bed
                missed_genes_grouped_dict[b_hit[3]] = update_locs(missed_genes_grouped_dict[hit[3]],missed_genes_dict[b_hit[3]])
                del missed_genes_grouped_dict[hit[3]]
            except KeyError:
                #if not in bed add to bed
                missed_genes_grouped_dict[b_hit[3]] = update_locs(missed_genes_dict[hit[3]],missed_genes_dict[b_hit[3]])
                del missed_genes_grouped_dict[hit[3]]
        else:
            ### add only hit_b to bed
            missed_genes_grouped_dict[b_hit[3]] = missed_genes_dict[b_hit[3]]
    return missed_genes_grouped_dict

##################################################
def main(missed_genes_path,old_bed,new_bed,out_file):
    missed_genes = parse_missed_genes(missed_genes_path)
    missed_genes_grouped,missed_genes_dict = group_genes_in_bed(missed_genes,old_bed,new_bed)
    #print "grouped genes"
    ### make sure this ^^^^ works ^^^^^ 
    new_genes_final = {}
    for qaccn in missed_genes_grouped:
        hits = missed_genes_grouped[qaccn]
        hit_set = set(hits)
        hits = list(hit_set)
        #non_overlapping = merge_overlapping(hits)
        grouped_hits = merge_hits(hits,old_bed,missed_genes_dict)
        new_genes_final.update(grouped_hits)
    write_new_bed(new_genes_final,old_bed,missed_genes,out_file)


if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser("usage: %prog [options] ")
    parser.add_option("--missed_bed", dest="new_bed", help="missed ORGA from ORGB bed file from coanno ")
    parser.add_option("--missed_matches", dest="missed_genes", help="missed ORGA from ORGB matches.txt file from coanno")
    parser.add_option("--old_bed", dest="old_bed", help="orginal bed file for ORG")
    parser.add_option("--out", dest="out_fh", help = "out_file: where the new merged bed should go")
    (options, _) = parser.parse_args()

    new_bed = Bed(options.new_bed)
    old_bed = Bed(options.old_bed)

    main(options.missed_genes,old_bed,new_bed,options.out_fh)
    
#merge_same_hits(Bed('data/athaliana_lyrata2/missed_lyrata_from_athaliana.bed'),'data/athaliana_lyrata2/missed_lyrata_from_athaliana.matches.txt',Bed('data/athaliana_lyrata2/lyrata.bed'))
#merge(Bed('data/athaliana_lyrata2/lyrata.bed'),Bed('data/athaliana_lyrata2/missed_from_lyrata.bed'),'data/athaliana_lyrata2/lyrata.all.bed')
