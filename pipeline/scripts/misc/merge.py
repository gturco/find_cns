from flatfeature import Bed
from intersection import Intersecter, Feature
import sys

def near_gene(haccn, gene_bed):
    """if within 1000 bp of an exsiting gene"""
    hit_start = min(haccn['start'], haccn['end'])
    hit_stop = max(haccn['start'], haccn['end'])
    bpmin = hit_start - 100
    bpmax = hit_stop + 100
    seqid = haccn['seqid']
    intervening_genes = True
    feats = gene_bed.get_features_in_region(str(seqid), bpmin, bpmax)
    accns = [i['accn'] for i in feats]
    if len(accns)==0 or haccn['accn'] in accns and len(accns)==1:
        intervening_genes = False
    return intervening_genes

def get_intervening_genes(start,end,seqid, org_bed,saccn):
    inbetween = org_bed.get_features_in_region(seqid, start, end)
    accns = [i['accn']for i in inbetween ]
    intervening_genes = True
    if len(accns)==0 or saccn in accns and len(accns)==1:
        intervening_genes = False
    if abs(end - start) < 7500:
        intervening_genes = True
    return intervening_genes

def merge_same_hits(missed, fh_match, org_bed):
    """ groups genes that hit more then once """
    d = {}
    handle = open(fh_match)
    matches = handle.read()
    org_bed_path = org_bed.path
    path = org_bed_path.split('/')
    dirc = '/'.join(path[:-1])
    org = path[-1]
    fh = open('{0}/missed_from_{1}'.format(dirc,org), "wb")
    for match in matches.split('\n')[:-1]:
        qaccn,saccn = match.split('\t')
        #create dictionary
        try:
            seqid = missed.accn(qaccn)['seqid']
            haccn = missed.accn(qaccn)
        except KeyError: continue
        #if near_gene(haccn,org_bed)==True: continue
        if (seqid,saccn) not in d.keys():
            #append whole dict to keys
            d[(seqid,saccn)]= missed.accn(qaccn)
        else:
            #else add locs to exsting one
            gene_start = min(d[(seqid,saccn)]['locs'])[0]
            gene_end = max(d[(seqid,saccn)]['locs'])[1]
            missed_end = missed.accn(qaccn)['locs'][0][1]
            missed_start = missed.accn(qaccn)['locs'][0][0]
            if missed_end < gene_start:
                # if no intervening genes and they are close together...
                intervening_genes = get_intervening_genes(missed_end,gene_start,seqid, org_bed, d[(seqid,saccn)]['accn'])
                if intervening_genes is False:
                    d[(seqid,saccn)]['locs'] =  d[(seqid,saccn)]['locs'] + missed.accn(qaccn)['locs']
                    d[(seqid,saccn)]['start'] = missed_start
                    if 'Sb' in qaccn:
		    	        d[seqid,saccn]['accn'] = qaccn
                else:
                    d[(seqid,qaccn)] = missed.accn(qaccn)
            elif gene_end < missed_start:
                intervening_genes = get_intervening_genes(gene_end,missed_start,seqid, org_bed,d[(seqid,saccn)]["accn"])
                if intervening_genes is False:
                    d[(seqid,saccn)]['locs'] =  d[(seqid,saccn)]['locs'] + missed.accn(qaccn)['locs']
                    d[(seqid,saccn)]['end'] = missed_end
                    if 'Sb' in qaccn:
                        d[seqid,saccn]['accn'] = qaccn
                else:
                    d[(seqid,qaccn)]= missed.accn(qaccn)
            else:
                d[(seqid,saccn)]['locs'] =  d[(seqid,saccn)]['locs'] + missed.accn(qaccn)['locs']
        
    for key in d.keys():
        new_row = d[key]['locs'].sort()
        row = d[key]
        print >>fh, Bed.row_string(row)



def merge(org_bed, missed, merge_file):
    """creates blast.all file and updates everything"""
    merge_fh = open(merge_file, "w")
    #cds_missed = missed[missed['ftype'] == 'CDS']
    #count = org_bed.shape[0] + missed[missed['ftype'] !='CDS'].shape[0]
    new_rows = []
    seen_accns = {}
    # CDS added to existing gene.
    for row_missed in missed:
        if row_missed['accn'] in seen_accns: continue
        try:
            org_bed_row = org_bed.accn(row_missed['accn'])
             # it's a CDS
        except KeyError:
            #its a new gene
            new_rows.append(row_missed)
            seen_accns[row_missed['accn']] = True
            continue
        locs_interval = Intersecter()
        [locs_interval.add_interval(Feature(start,stop)) for start,stop in org_bed_row['locs']]
        for missed_start,missed_end in row_missed['locs']:
            if len(locs_interval.find(missed_start,missed_end)) > 0:
#                print >>sys.stderr, org_bed_row['accn']
                locs_intersects = [(l.start,l.stop) for l in locs_interval.find(missed_start,missed_end)]
                [org_bed_row['locs'].remove(locs_intersect) for locs_intersect in locs_intersects]
                locs_intersects = set(locs_intersects)
		locs_intersects.add((missed_start,missed_end))
                locs_start = min([start for start,end in locs_intersects])
                locs_end = max([end for start,end in locs_intersects])
                org_bed_row['locs'] = org_bed_row['locs'] + [(locs_start,locs_end)]
                row_missed['locs'].remove((missed_start,missed_end))

        org_bed_row['locs'] = org_bed_row['locs'] + row_missed['locs']
        #print >>sys.stderr, "{0},{1}".format(row_missed['accn'], locs)
        org_bed_row['locs'].sort()
        org_bed_row['start'] = min(min([start for start,end in org_bed_row['locs']]), org_bed_row['start'])
        org_bed_row['end'] = max(max([end for start,end in org_bed_row['locs']]), org_bed_row['end'])
        new_rows.append(org_bed_row)
        seen_accns[org_bed_row['accn']] =True

    for org_bed_rw in org_bed:
        if org_bed_rw['accn'] not in seen_accns:
            new_rows.append(org_bed_rw)
            seen_accns[org_bed_rw['accn']] =True

    def row_cmp(a,b):
        return cmp(a['seqid'], b['seqid']) or cmp(a['start'], b['start'])


    new_rows.sort(cmp=row_cmp)
    #print >>merge_fh, "\t".join(Bed.names)
    for i, row in enumerate(new_rows):
        print >>merge_fh, Bed.row_string(row)

def main(missed,fh_match,org_bed):
    """first megers all hits to the same gene... then updates the entire bed
    file output: all_ORG.bed """
    
    merge_same_hits(missed, fh_match, org_bed)
    org_bed_path = org_bed.path
    path = org_bed_path.split('/')
    dirc = '/'.join(path[:-1])
    org = path[-1]
    missed2 = '{0}/missed_from_{1}'.format(dirc,org)
    merge_fh = "{0}/all_{1}".format(dirc,org)
    print missed2
    merge(org_bed, Bed(missed2),merge_fh)

if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser("usage: %prog [options] ")
    parser.add_option("--missed", dest="missed", help="missed ORGA from ORGB bed file from coanno ")
    parser.add_option("--match", dest="fh_match", help="missed ORGA from ORGB matches.txt file from coanno")
    parser.add_option("--org", dest="org_bed", help="orginal bed file for ORG")
    (options, _) = parser.parse_args()

    missed_bed = Bed(options.missed)
    org_bed = Bed(options.org_bed)

    main(missed_bed,options.fh_match,org_bed)
    
#merge_same_hits(Bed('data/athaliana_lyrata2/missed_lyrata_from_athaliana.bed'),'data/athaliana_lyrata2/missed_lyrata_from_athaliana.matches.txt',Bed('data/athaliana_lyrata2/lyrata.bed'))
#merge(Bed('data/athaliana_lyrata2/lyrata.bed'),Bed('data/athaliana_lyrata2/missed_from_lyrata.bed'),'data/athaliana_lyrata2/lyrata.all.bed')
