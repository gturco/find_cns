import sys
import os.path as op
from flatfeature import Bed
sys.path.append("../")
from common import parse_raw_cns
from common import read_cns_to_rna, read_cns_to_protein_exons
sys.path.insert(0, "/home/gturco/src/quota-alignment/scripts")
from bed_utils import RawLine
import collections

from scipy.spatial import cKDTree

def cns_to_str(cns):
    key = "%(qseqid)s,%(qaccn)s,%(sseqid)s,%(saccn)s,%(qstart)i,%(qend)i,%(sstart)i,%(send)i,%(eval)s"
    return key % cns
     
def fill_tree(new_genes_all):
    trees = {}
    for k, new_genes in new_genes_all.items():
        for p, info in new_genes:
            assert k == (p['qseqid'], p['sseqid'])
            if not k in trees:
                trees[k] = []
            pt = (p['qstart'] + p['qend'])/2, (p['sstart'] + p['send'])/2
            trees[k].append(pt)
    for k in trees:
        trees[k] = cKDTree(trees[k])
    return trees


def get_new(cns, trees, key, prs, dist, seen={}):
    """
    find any new genes that are close (dist) to the given
    cns for the given qchr, schr key.
    """
    tree = trees.get(key)
    if tree is None: raise StopIteration
    
    pt = (cns['qstart'] + cns['qend'])/2, (cns['sstart'] + cns['send'])/2
    dists, idxs = tree.query(pt, p=2, distance_upper_bound=dist, k=64)
    idxs = idxs[idxs != tree.n]
    if not len(idxs): raise StopIteration
    pr = prs[key]
    # get the new gene from the list.
    for idx in idxs:
        gnew, info = pr[idx]
        yield gnew, info


def main(qbed, sbed, cnsfile, dist, orthology_path):
    """
    here, we remove cnss that have been called proteins/rnas from 
    the cns list, and add them to the bed files.
    AND have to do the preliminary assignment of cnss that remain to the new-genes
    that _were_ cnss. the proper assignment is then handled in assign.py
    """
    ortho_trees = read_orthos_to_trees(orthology_path, qbed, sbed)

    name, ext = op.splitext(cnsfile)
    real_cns_fh = open("%s.real%s" % (name, ext), "w")
    print >>sys.stderr, "writing to:", real_cns_fh.name
    outdir = op.dirname(cnsfile)
    print >>real_cns_fh, "#qseqid,qaccn,sseqid,saccn,qstart,qend,sstart,send,eval"

    crna = read_cns_to_rna(outdir)
    cpro = read_cns_to_protein_exons(outdir)

    cns_items = list(parse_raw_cns(cnsfile))
    proteins = collections.defaultdict(list)
    rnas = collections.defaultdict(list)
    real_cns_items = []
    for cns_id, cns in cns_items:
        key = (cns['qseqid'], cns['sseqid'])
        if cns_id in cpro:
            proteins[key].append((cns, cpro[cns_id]))
        elif cns_id in crna:
            rnas[key].append((cns, crna[cns_id]))
        else:
            real_cns_items.append((cns_id, cns))

    p_trees = fill_tree(proteins)
    r_trees = fill_tree(rnas)

    def assign_new_names(prs, protein_or_rna):
        n = {}
        for seqid_pair, li in prs.iteritems():
            if not seqid_pair in n: n[seqid_pair] = []
            for gnew, info in li[:]:
                new_qname = "%(qseqid)s_%(qstart)i_%(qend)i_cns" % gnew
                new_sname = "%(sseqid)s_%(sstart)i_%(send)i_cns" % gnew
                # and give them both an id so we know they were a pair.
                new_qname += "_%s" % (protein_or_rna)
                new_sname += "_%s" % (protein_or_rna)
                try:
                    qstrand = qbed.d[gnew['qaccn']]['strand']
                    sstrand = sbed.d[gnew['saccn']]['strand']
                except:
                    print >>sys.stderr, gnew
                    raise
                gnew['qaccn'] = new_qname
                gnew['saccn'] = new_sname
                gnew['qstrand'] = qstrand
                gnew['sstrand'] = sstrand
                n[seqid_pair].append((gnew, info))
        return n
    nproteins = assign_new_names(proteins, "protein")
    nrnas = assign_new_names(rnas, "rna")


    cns_seen = {}
    # go through the remaining cnss, print and assign them to the new
    # genes (previously cnss) in within dist.
    for cns_id, cns in real_cns_items:
        print >>real_cns_fh, cns_to_str(cns)
        key = (cns['qseqid'], cns['sseqid'])
        
        for pnew, info in get_new(cns, p_trees, key, nproteins, dist + 1000):
            cns['qaccn'] = pnew['qaccn']
            cns['saccn'] = pnew['saccn']
            cns_str = cns_to_str(cns)
	    if cns_str in cns_seen: continue
            cns_seen[cns_str] = 1
            print >>real_cns_fh, cns_str

        for rnew, info in get_new(cns, r_trees, key, nrnas, dist + 1000):
            cns['qaccn'] = rnew['qaccn']
            cns['saccn'] = rnew['saccn']
   	    cns_str = cns_to_str(cns)
            if cns_str in cns_seen: continue
            cns_seen[cns_str] = 1
            print >>real_cns_fh, cns_str

    qbed_list, qnew_pairs = merge_bed(qbed, nproteins, nrnas, ortho_trees, 'q')
    # dont need to do the orthos 2x so send in empty dict.
    sbed_list, snew_pairs_unused = merge_bed(sbed, nproteins, nrnas, {}, 's')

    # if it's the same org, we add the new cnss again to the same we send in both lists.
    # print_bed handles the repeats.
    if qbed.path == sbed.path:
        qbed_new = sbed_new = print_bed(qbed_list + sbed_list, qbed.path)
    else:
        qbed_new = print_bed(qbed_list, qbed.path)
        sbed_new = print_bed(sbed_list, sbed.path)

    return qbed_new, sbed_new, qnew_pairs


def print_bed(flist, old_path):
    ipath, ext = op.splitext(old_path)
    path = "%s.with_new%s" % (ipath, ext)

    print >>sys.stderr,  "writing to: %s.with_new%s" % (ipath, ext)
    fh = open(path, 'wb')
    seen = {}

    for item in flist:
        # convert the locs to a tuple.
        #print >>sys.stderr, item
        item = list(item)
        item[6] = tuple(item[6])
        item = tuple(item)
        if item in seen: continue
        seen[item] = 1
        locs = item[6] # tuple(sorted([item[1], item[2]]))

        row = dict(accn=item[3], start=item[1], end=item[2], seqid=item[0],
                   locs=locs, score='.', strand=item[5], rgb='.', thickstart='.', thickend=".")
        print >>fh, Bed.row_string(row)
    fh.close()
    return Bed(path)

def write_new_pairs(old_paralogy_path, old_orthology_path, qbed, qbed_new, sbed, sbed_new, new_pairs):
    """
    need to create a new .paralogy/.orthology file that includes the added
    pairs that were once cns's.
    this is a bit complicated because we need to map the old indexes in paralogy to the gene names,
    then back to the new index.
    """
    qbed.fill_dict(); sbed.fill_dict()
    qbed_new.fill_dict(); sbed_new.fill_dict()

    qi_new_lookup = dict((accn, i) for i, accn in enumerate(qbed_new['accn']))
    si_new_lookup = dict((accn, i) for i, accn in enumerate(sbed_new['accn']))

    ppath, pext = op.splitext(old_paralogy_path)
    opath, oext = op.splitext(old_orthology_path)
    ppath = "%s.with_new%s" % (ppath, pext)
    opath = "%s.with_new%s" % (opath, oext)

    oseen, pseen = {}, {}

    for line in open(old_paralogy_path):
        if line[0] == "#": continue
        r = RawLine(line)
        q, s = qbed[r.pos_a], sbed[r.pos_b]
        r.pos_a = qi_new_lookup[q['accn']]
        r.pos_b = si_new_lookup[s['accn']]
        pseen[(r.pos_a, r.pos_b)] = str(r)
        #print >>pfh, line,

    for line in open(old_orthology_path):
        if line[0] == "#": continue
        r = RawLine(line)
        q, s = qbed[r.pos_a], sbed[r.pos_b]
        r.pos_a = qi_new_lookup[q['accn']]
        r.pos_b = si_new_lookup[s['accn']]
        assert q['accn'] == qbed_new[r.pos_a]['accn']
        assert s['accn'] == sbed_new[r.pos_b]['accn']
        oseen[(r.pos_a, r.pos_b)] = str(r)
        #print >>ofh, line,

    for p in new_pairs:
        assert p['qstart'] < p['qend'] and p['sstart'] < p['send'], (p, )
        qi = qi_new_lookup[p['qaccn']]
        si = si_new_lookup[p['saccn']]
        if (qi, si) in pseen: continue
        line_str = "%s\t%i\t%s\t%i\t%.3g" % \
                 (p['qseqid'], qi, p['sseqid'], si, 40)
        pseen[(qi, si)] = line_str
        #print >>pfh, line_str
        if p['is_ortho']:
            if (qi, si) in oseen: continue
            oseen[(qi, si)] = line_str
            #print >>ofh, line_str

    print >>sys.stderr, "writing orthology pairs to %s" % opath
    print >>sys.stderr, "writing paralogy pairs to %s" % ppath
    pfh = open(ppath, 'wb')
    ofh = open(opath, 'wb')
    for qi, si in sorted(oseen):
        print >>ofh, oseen[(qi, si)]
    for qi, si in sorted(pseen):
        print >>pfh, pseen[(qi, si)]
        
        
def merge_bed(bed, proteins, rnas, ortho_trees, q_or_s):
    """
    merge the existing bed file with the new genes in proteins and rnas
    also write a file outdir/protein_rna.anno to use when making genelist
    """
    new_anno = open("%s/%s_protein_rna.anno" % (op.dirname(bed.path), q_or_s), "w")
    MAX_DIST = 20000
    bedlist = [tuple(ff) for ff in bed]
    seen = {}
    new_pairs = []
    for pr in (proteins, rnas):
        for seqid_pair, li in pr.iteritems():
            tree = ortho_trees.get(seqid_pair, None)
            for pnew, info in li:
                pnew = pnew.copy()
                # sort start, stop.
                if pnew['sstart'] > pnew['send']:
                    pnew['send'], pnew['sstart'] = pnew['sstart'], pnew['send']
                anew = (pnew[q_or_s + 'seqid'],
                        pnew[q_or_s + 'start'],
                        pnew[q_or_s + 'end'],
                        pnew[q_or_s + 'accn'],
                        '.', # score
                        pnew[q_or_s + 'strand'],
                        ((pnew[q_or_s + 'start'], pnew[q_or_s + 'end']),),
                        'gene_cns'
                       )
                if anew in seen: continue
                bedlist.append(anew) 

                is_ortho = False
                if tree is not None:
                    pt = (pnew['qstart'] + pnew['qend'])/2, (pnew['sstart'] + pnew['send'])/2
                    dists, idxs = tree.query(pt, k=2, distance_upper_bound=MAX_DIST)
                    idxs = idxs[idxs != tree.n]
                    is_ortho = bool(len(idxs))
                pnew['is_ortho'] = is_ortho
                pnew['pair'] = pnew['qaccn' if q_or_s == 's' else 'saccn']
                # accn, 

                sline = "%s,%s" % (pnew[q_or_s + 'accn'], info.replace(",", ";"))
                if sline in seen: continue
                new_pairs.append(pnew)
                seen[sline] = True
                print >>new_anno, sline
    # sort by seqid, then start.
    print >>sys.stderr, "wrote protein, rna info to %s" % new_anno.name
    return sorted(bedlist, key=lambda a:(a[0], a[1])), new_pairs
       
    
def read_orthos_to_trees(forthos, qbed, sbed):
    trees = {}
    for line in open(forthos):
        if line[0] == "#": continue
        raw = RawLine(line)
        qpos = qbed[raw.pos_a]
        spos = sbed[raw.pos_b]
        key = (raw.seqid_a, raw.seqid_b)
        if not key in trees: trees[key] = []
        qpos = (qpos['start'] + qpos['end']) / 2
        spos = (spos['start'] + spos['end']) / 2
        trees[key].append((int(qpos), int(spos)))
    for k in trees:
        trees[k] = cKDTree(trees[k])
    return trees

if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser()
    parser.add_option("--qbed", dest="qbed", help="query bed file")
    parser.add_option("--sbed", dest="sbed", help="subject bed file")
    parser.add_option("--cns", dest="cns", help="path to raw cns")
    parser.add_option("--dist", dest="dist", type='int', help="max dist from gene to cns", default=12000)
    parser.add_option("--paralogy", dest="paralogy", help="path to paralogy file")
    parser.add_option("--orthology", dest="orthology", help="path to orthology file")

    options, args = parser.parse_args()    

    if not (options.sbed and options.qbed and options.cns, options.orthology):
        sys.exit(parser.print_help())

    qbed = Bed(options.qbed); qbed.fill_dict()
    sbed = Bed(options.sbed); sbed.fill_dict()
    
    qbed_new, sbed_new, new_pairs = main(qbed, sbed, options.cns, options.dist, options.orthology)
    write_new_pairs(options.paralogy, options.orthology, qbed, qbed_new, sbed, sbed_new, new_pairs) 
