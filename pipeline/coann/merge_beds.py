#### atmit one... overlap preexcesting genes...
from merge import update_locs, merge_feats
from flatfeature import Bed

def write_bed(gene,merge_fh):
     new_line = Bed.row_string(gene)
     merge_fh.write("{0}\n".format(new_line))

def overlapping(start,end,strand,seqid,accns):
    overlaps = []
    for accn in accns:
        #print accn
        if accn['seqid'] == seqid and accn['start'] <= end and accn['end'] >= start and accn['strand'] == strand:
            overlaps.append(accn)
    return overlaps

def new_genes(old_bed,new_bed):
    """finds all the genes that were changed in the new bed file"""
    new_genes = []
    for new_accn in new_bed:
        try:
            old_accn = old_bed.accn(new_accn['accn'])
            if list(old_accn)[:7] == list(new_accn)[:7]:continue
            ### watch out for . -1 changes
            else: new_genes.append(new_accn)
        except KeyError:
            new_genes.append(new_accn)
    return new_genes

def write_old_bed(new_genes_1,new_genes_2,old_bed,merge_fh):
    new_genes = new_genes_1 + new_genes_2
    new_accns = [new_gene['accn'] for new_gene in new_genes]
    for gene in old_bed:
        if gene in new_accns: continue
        write_bed(gene,merge_fh)

def main(old_bed,new_bed1,new_bed2,merge_path):
    merge_fh = open(merge_path,'w')
    new_genes_1 = new_genes(old_bed,new_bed1)
    new_genes_2 = new_genes(old_bed,new_bed2)
    print len(new_genes_1), len(new_genes_2)
    all_overlap = []
    for new_gene in new_genes_1:
        #print new_gene['accn']
        ### does it overlapp with any of the other new genes....
        overlapping_genes = overlapping(new_gene['start'],new_gene['end'],new_gene['strand'],new_gene['seqid'],new_genes_2)
        if overlapping_genes == 0:
            write_bed(gene,merge_fh)
            continue
        ### append all overlaping accns
        all_overlap.extend(overlapping_genes)
        for overlapping_gene in overlapping_genes:
            new_gene = update_locs(new_gene,overlapping_gene)
        merged_gene = merge_feats(new_gene)
        write_bed(new_gene,merge_fh)
        #### if it does merge the numbers
    for new_gene2 in new_genes_2:
        if new_gene2 not in all_overlap: write_bed(new_gene2,merge_fh)
    write_old_bed(new_genes_1,new_genes_2,old_bed,merge_fh)

main(Bed("data/rice_v6_setaria64/rice_v6.bed"),Bed("data/rice_v6_setaria64/rice_v6.all2.bed"),Bed("data/rice_v6_sorghum_v1/rice_v6.all2.bed"),"test")
