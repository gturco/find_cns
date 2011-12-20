#### atmit one... overlap preexcesting genes...
from merge import update_locs, merge_feats
from flatfeature import Bed

def write_bed(gene,merge_fh):
     new_line = Bed.row_string(gene)
     merge_fh.write("{0}\n".format(new_line))

def overlapping(start,end,strand,seqid,accns):
    overlaps = []
    for accn in accns:
        if accn['seqid'] == seqid & accn['start'] <= end & accn['end'] >= start & accn['strand'] == strand:
            overlaps.append(accn)
    return overlaps

def new_genes(old_bed,new_bed,merge_fh):
    """finds all the genes that were changed in the new bed file"""
    new_genes = []
    for new_accn in new_bed:
        try:
            old_accn = old_bed.accn(new_accn['accn'])
            if list(old_bed)[:7] == list(new_accn)[:7]:continue
            ### watch out for . -1 changes
            else: new_genes.append(new_accn)
        except KeyError:
            new_genes.append(new_accn)
    return new_genes

def main(old_bed,new_bed1,new_bed2,merge_fh):
    new_genes_1 = new_genes(old_bed,new_bed1)
    new_genes_2 = new_genes(old_bed,new_bed2)
    print len(new_genes_1), len(new_genes_2)
    all_overlap = []
    for new_gene in new_genes_1:
        ### does it overlapp with any of the other new genes....
        overlapping_genes = overlapping(new_gene['start'],new_gene['end'],new_gene['strand'],new_gene['seqid'],new_genes_2)
        if overlapping_genes == 0: write_bed(gene,merge_fh)
        ### append all overlaping accns
        all_overlap.extend(overlapping_genes)
        for overlapping_gene in overlapping_genes:
            new_gene = update_locs(new_gene,overlaping_gene)
        merged_gene = merge_feats(new_gene)
        write_bed(new_gene,merge_fh)
        #### if it does merge the numbers
    for new_gene2 in new_genes_2:
        if new_gene2 not in all_overlap: write_bed(new_genes_2,merge_fh)
main()
