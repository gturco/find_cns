import random
from flatfeature import Bed
from pyfasta import Fasta
from collections import deque
from itertools import chain
"""input gene space of intrest and features of intrest in space
output list of random features in the gene space same size as other feature not
overlapping other ones"""

# add cns to bed file
# mask genome noncodding seq...

def gene_size_dict(tab_sep_file):
    """select accn, group_concat(qdiff) from
    rice_b_setaria64_gturco_2011_4_11app_real group by accn"""
    handle = open(tab_sep_file)
    fh = handle.read()
    gene_size_dict = {}
    for line in fh.split("\n")[:-1]:
        gene, cnss = line.split("\t")
        gene_size_dict[gene] = cnss.split(",")
    return gene_size_dict

def find(x,y):
    """returns ture is a set is not disjoint"""
    return x[1] >= y[0] or x[1] >= y[1]

def union(x,y):
    """joins two points if they are not disjoint"""
    if find(x,y):
        start = min(x[0],y[0])
        stop = max(x[1],y[1])
        return deque([(start,stop)])
    else:
        return deque([x,y])

def merged_feats(feats):
    prv_count = 1
    output = deque([])
    #print "feats2{0}".format(feats)
    while prv_count != len(feats) and len(feats) > 0:
        prv_count = len(feats)
        if len(feats) > 1:
            x = feats.popleft()
            y = feats.popleft()
            merged = union(x,y)
            output += merged
        else:
            if feats[0] not in output:
                output += feats
                feats.popleft()
       # print "output: {0}".format(output)
    return output

def recursive_merge(feats):
    new_feats = feats
    prv_count = 0
    while prv_count != len(new_feats) and len(new_feats) > 1:
        prv_count = len(new_feats)
        #print prv_count
        new_feats = merged_feats(new_feats)
    #print new_feats
    return new_feats

def recursive_merge_both(feats):
    evens = recursive_merge(feats)
    first_even = deque([evens[0]])
    evens.popleft()
    odds = recursive_merge(evens)
    if len(odds) > 1:
        first_odd = deque([odds[0]]) 
        odds.popleft()
        odds_two = recursive_merge(odds)
        list_both = list(first_even) + list(odds_two) + list(first_odd)
    else:
        list_both = list(first_even) + list(odds)
    list_both.sort()
    return list_both

def non_coding(feats,cns_size,search_start,search_end):
    """find the inverse of genes ... the non_coding seq.. inverse the list"""
    non_cd_regions = []
    cns_size = int(cns_size)
    off_by_one = [range(stop+1,(feats[i+1][0])-cns_size) for
                    i,(start,stop) in enumerate(feats[:-1]) if feats[i+1][0]-cns_size > stop]
    region_start = range(search_start, int(feats[0][0])-int(cns_size))
    off_by_one.append(region_start)
    if search_end - cns_size > feats[-1][1]:
        for n in range(feats[-1][1]+1,search_end):
            if n + cns_size <= search_end:
                off_by_one.append([n])
    return off_by_one

def random_noncoding(cns_gene_dict,bed):
    random_cnss = []
    handle = open("/Users/gturco/gene_cns_dict.txt","wb")
    for gene in cns_gene_dict.keys():
        gene = bed.accn(gene)
        search_start = gene["start"] - 12000
        search_end = gene["end"] + 12000

        feats = bed.get_features_in_region(gene["seqid"],search_start,search_end)
        #print search_start,search_end
        #print feats
        locs =[feat["locs"] for feat in feats]
        flatten_locs = chain.from_iterable(locs)
        flatten_locs = list(flatten_locs)
        start_stops =[(int(start),int(end)) for (start,end) in flatten_locs]
        start_stops_in_range = []
        for s_start,s_stop in start_stops:
            if s_start <= search_end and s_stop >= search_start:
                start_stops_in_range.append((s_start,s_stop))
        #print search_start,search_end
        #print start_stops_in_range
        start_stops_in_range.sort()
        #print start_stops
        handle.write("{0}\n".format(gene["accn"]))
        if len(start_stops_in_range) > 1:
            start_stops_m = recursive_merge_both(deque(start_stops_in_range))
        else:
            start_stops_m = start_stops_in_range
        merged_feats = list(start_stops_m)
        merged_feats.sort()
        #print merged_feats
        for cns_size in cns_gene_dict[gene["accn"]]:
            ### cns start sites taking cns size into accoun so no overlap with
            off_by_one = non_coding(merged_feats,cns_size,search_start,search_end)
            non_cd_regions = []

            for range_list in off_by_one:
                if len(range_list) > 0:
                    for number in range_list:
                        non_cd_regions.append(number)
            #print non_cd_regions
            if len(non_cd_regions) > 0:
                random_start_sites = random.sample(non_cd_regions,1)
                #print random_start_sites
                #w= "start:{0},end:{1},gene:{2},merged_feats:{3}\n".format(random_start_sites,cns_size,gene,merged_feats)
                #handle.write(w)
                #random_cns = [(gene["seqid"],int(start),start+int(cns_size))
                #        for start in random_start_sites]
            
                #[random_cnss.append(cns) for cns in random_cns]
                cns = (gene["accn"],gene["seqid"],random_start_sites[0],random_start_sites[0]+int(cns_size))
                random_cnss.append(cns)
    return random_cnss

def get_seq(random_cns_list, fasta_file, new_fasta_name):
    new_fasta = open("{0}".format(new_fasta_name),"wb")
    n = 0
    for random_cns in random_cns_list:
        n = n +1
        accn,seqid,start,end = random_cns
        f = Fasta(fasta_file)
        seq = f[seqid][start:end]
        if "X" in seq:
            print accn,seqid,start,end
        if len(seq) < 15 and len(seq) > 0:
            print "OH NO!!!!!!"
        w = ">cns{0}\n".format(n)
        seq_w = "{0}\n".format(seq)
        new_fasta.write(w)
        new_fasta.write(seq_w)

dict_size = gene_size_dict('/Users/gturco/Desktop/rice_sorg_size.tsv')
#dict_size = gene_size_dict("test_file")
x = random_noncoding(dict_size,Bed('/Users/gturco/data/paper3/rice_b_sorghum_v1.nolocaldups.with_new_cns_mask.bed'))
print len(x)
#####print x
get_seq(x,"/Users/gturco/data/paper3/rice_b.fasta","/Users/gturco/test.fasta")
##
##### seq for cns
#handle = open("/Users/gturco/data/paper3/rice_b_sorghum_v1.cns.assigned_real.csv")
#fh = handle.read()
#cns_list = []
#for line in fh.split("\n")[:-1]:
#    if line[0] == "#": continue
#    cns_id,accn,seqid,start,end,strand = line.split(",")[:6]
#    cns_list.append((seqid,int(start),int(end)))
#
#len(cns_list)
#get_seq(cns_list,"/Users/gturco/data/paper3/rice_b.fasta","/Users/gturco/test_cns.fasta")
##

