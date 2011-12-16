import simplejson
import logging
import bblast
import collections
import operator
import os
import os.path as op
import sys
import string
from biostuff import BlastLine
from pyfasta import Fasta
import mask_features as mf
from flatfeature import Flat, Bed



logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger('coanno')


MATCH_JOIN = "#"

def run_mask(cfg):
    for g in ("genome_a", "genome_b"):
        mf.main(cfg[g]["flat"], cfg[g]["fasta"])
        log.debug("MASKED: %s" % cfg[g]['name'])

def update_cfg(cfg):
    """
    set some defaults.

    >>> cfg = {'default': {}}
    >>> update_cfg(cfg)
    >>> cfg['default']['min_overlap']
    30
    """

    if not "overlap" in cfg["default"]:
        cfg["default"]["min_overlap"] = MIN_OVERLAP
    if not "coverage" in cfg["default"]:
        cfg["default"]["min_pct_coverage"] = MIN_PCT_COVERAGE


def main(cfg_path):
    cfg = simplejson.load(open(cfg_path))
    update_cfg(cfg)

    run_mask(cfg)

    run_blasts(cfg)

    dispatch(cfg)
    if cfg["default"].get("reciprocal"):
        dispatch(cfg, flip=True)

MIN_OVERLAP = 30
MIN_PCT_COVERAGE = 0.40

def merge_overlapping(new_genes, min_overlap, min_pct_coverage, match_file):
    """
    require at least `overlap` basepairs worth of overlap
    >>> list(merge_overlapping({}, 30, 0.4))
    []

    """
    match_fh = open(match_file, "w")

    merged = 0
    # since we have to back track sometimes
    # keep track of stuff on other strand.
    #dbg = open("genesets.txt", "w")

    for achr, genes in sorted(new_genes.iteritems()):
        genes = sorted([g for _name, g in genes.items()], key=operator.itemgetter('start'))
        i, j = 0, 1
        used_genes = set([])
        seen = {}
        while i < len(genes):

            # it was probably seen on the other strand.
            if i in used_genes: i += 1; j = i + 1; continue
            g = genes[i].copy()
            matches = g.pop('match')
            new_stop = g['end']

            last_gene = None # keep track of gene on other strand.

            #  PARAMETER: require 30bp overlap.
            while j < len(genes) and g['end'] > genes[j]['start'] + min_overlap:
                # can only merge genes with genes, and cds with cds.
                if j in used_genes or genes[j]['type'] != g['type']: 
                    j += 1
                    continue

                minstart = g['start']
                maxstop = max(new_stop, genes[j]['end'])
                maxstart = genes[j]['start']
                minstop = min(new_stop, genes[j]['end'])
                coverage = (minstop - maxstart) / (1.0 * maxstop - minstart)
                if genes[j]['strand'] != g['strand'] and coverage < min_pct_coverage:
                    last_gene = j
                    j += 1
                    continue
                new_stop = maxstop
                used_genes.update([j])
                matches += "," + genes[j].pop('match')
                j += 1
                merged += 1
          
            g['end'] = new_stop
            key = g['start'], g['end']
            if key in seen: continue
            seen[key] = None
            fix_name(g)
            assert len(g['accn']) <= 80, (len(g['accn']), g)

            print >>match_fh, "%s\t%s" % (g['accn'], matches)


            yield g
            used_genes.update([i])
            i = last_gene or j
            j = i + 1
    print >>sys.stderr, "merged %i genes" % merged

def fix_name(g):
    """the name is like: org_chr_start_stop but we may have changed start, stop

    >>> f = dict(name='name', start=2, end=3, strand='+', chr='chr3', 
    ...           type='CDS', attrs={'match':'amatch', 'ID': 'at2g26540_4_5'})
    >>> fix_name(f)
    >>> f
    {'start': 2, 'chr': 'chr3', 'name': 'at2g26540_2_3', 'type': 'CDS', 'end': 3, 'strand': '+', }
    
    `"""
    if g['type'] == 'CDS': return g['accn']
    new_name = g['accn']
    new_name = new_name.rstrip(string.digits).rstrip("_").rstrip(string.digits) + "%i_%i"
    new_name %= (g['start'], g['end'])
    g['accn'] = new_name

def dispatch(cfg, flip=False):
    """
    >>> cfg = simplejson.loads(open("tests/data/athaliana_athaliana.json").read())
    >>> update_cfg(cfg)
    >>> dispatch(cfg)

    """
    if flip is False:
        akey, bkey = "genome_a", "genome_b"
    else:
        akey, bkey = "genome_b", "genome_a"

    afasta = cfg[akey]["fasta"]
    bfasta = cfg[bkey]["fasta"]
    odir = cfg["default"]["out_dir"]
    assert os.path.exists(odir), "Need to create directory: %s" % odir
    min_len = cfg["default"]["min_len"]
    min_pct_coverage = cfg["default"]["min_pct_coverage"]

    ext = op.splitext(cfg[akey]['flat'])[1].lstrip(".")
    out_flat = "missed_%s_from_%s.%s" % (cfg[bkey]["name"], cfg[akey]["name"], ext)
    out_flat = os.path.join(cfg["default"]["out_dir"], out_flat)
    out_fh = open(out_flat, "w")
    

    a, b = fastas_for_features_vs_masked_genomic(afasta, bfasta)
    a_bnon_blast = bblast.get_blast_file(a, b, odir)
    a_b_blast = bblast.get_blast_file(afasta.replace(".fa", ".features.fa"),
                                      bfasta.replace(".fa", ".features.fa"), odir)
    new_genes = collections.defaultdict(dict)

    Klass = Bed if ext == "bed" else Flat

    aflat = Klass(cfg[akey]["flat"], cfg[akey]["fasta"])
    bflat = Klass(cfg[bkey]["flat"], cfg[bkey]["fasta"])

    for new_gene in find_missed(cfg[bkey]["name"],
                            aflat, bflat,
                            a_bnon_blast, a_b_blast, min_pct_coverage, flip=flip, min_len=min_len):
        #if new_gene['name'] in seen: continue
        n, achr = new_gene['accn'], new_gene['seqid']
        # found it 2x, add the match to the list.
        if n in new_genes[achr]:
            continue
        else:
            new_genes[achr][n] = new_gene
    # TODO: merge overlapping. only need to check chr, start, stops. so
    # sortable without index.
    match_file = "%s/missed_%s_from_%s.matches.txt" % (cfg['default']['out_dir'], 
                                           cfg[bkey]['name'], 
                                           cfg[akey]['name'])
    merged_genes = merge_overlapping(new_genes, cfg['default']['min_overlap'],
                                     cfg['default']['min_pct_coverage'], match_file)

    merged_genes = exclude_genes_in_high_repeat_areas(merged_genes, bfasta)
    if Klass == Flat:
        print >>out_fh, "\t".join(Klass.names)
    for i, new_gene in enumerate(merged_genes):
        new_gene['locs'] = [(new_gene['start'], new_gene['end'])]
        new_gene['score'] = new_gene['rgb'] = new_gene['thickend'] = new_gene['thickstart'] = "."
        print >>out_fh, Klass.row_string(new_gene)
    out_fh.close()
    log.debug("created %i new features in %s. with matches written to %s." \
                      % (i, out_flat, match_file))

    merge_file = "%s.all.%s" % (os.path.splitext(cfg[bkey]['flat'])[0], ext)
    log.debug("writing merged .%s file with new features to %s" % (ext, merge_file))
    merge(bflat, Klass(out_flat), merge_file, Klass)


def merge(main, missed, merge_file, Klass):
    merge_fh = open(merge_file, "w")
    #cds_missed = missed[missed['ftype'] == 'CDS']
    #count = main.shape[0] + missed[missed['ftype'] != 'CDS'].shape[0]
    new_rows = []
    seen_accns = {}
    # CDS added to existing gene.
    for row_missed in missed:
        if row_missed['accn'] in seen_accns: continue
        try:
            main_row = main.accn(row_missed['accn'])
            # it's a CDS
        except KeyError:
            # it's a new gene
            new_rows.append(row_missed)
            #seen_accns[row_missed['accn']] = True
            continue
        main_row['locs'] = main_row['locs'] + row_missed['locs']
        main_row['locs'].sort()
        new_rows.append(main_row)
        seen_accns[main_row['accn']] = True

    for row in (r for r in main if not r['accn'] in seen_accns):
        new_rows.append(row)

    def row_cmp(a, b):
        return cmp(a['seqid'], b['seqid']) or cmp(a['start'], b['start'])

    new_rows.sort(cmp=row_cmp)
    if Klass == Flat:
        print >>merge_fh, "\t".join(Klass.names)
    for i, row in enumerate(new_rows):
        if Klass == Flat:
            row['id'] = i + 1
        print >>merge_fh, Klass.row_string(row)


def exclude_genes_in_high_repeat_areas(merged_genes, bfasta):
    #print "FASTA:", afasta

    f = Fasta(bfasta)
    skipped = 0
    for gene in merged_genes:
        # get the total sequence length.
        seq = str(f[gene['seqid']][gene['start']:gene['end']])
        tot = len(seq)
        # and the lenght of real sequence.
        seq = seq.upper().replace('N', '').replace('X','')
        # if it's not > 80% good sequence, just skip it.
        if float(len(seq))/tot < .85: skipped+=1; continue
        yield gene
    log.info("removed %i (otherwise) new genes in masked areas" % skipped)


def dist(a, b):
    """
    >>> class O():
    ...    def __init__(self, start, stop):
    ...        self.sstart, self.sstop = start, stop

    >>> a, b, c = O(12, 15), O(10, 11), O(16, 17)
    >>> dist(a, b)
    1
    >>> dist(a, c)
    1
    >>> dist(b, c)
    5
    """

    d1 = abs(a.sstart - b.sstart)
    d2 = abs(a.sstop - b.sstart)
    d3 = abs(a.sstop - b.sstop)
    d4 = abs(a.sstart - b.sstop)
    return min(d1, d2, d3, d4)

def partition(slist, max_dist=1000):
    """
    return a list of lists where the hits are sorted by location and things 
    farther than max_dist apart are in separate lists.

    >>> class O():
    ...    def __init__(self, start, stop):
    ...        self.sstart, self.sstop = start, stop
    ...    def __repr__(self): return "O(%i, %i)" % (self.sstart, self.sstop)

    >>> partition([1])
    [[1]]

    >>> slist = [O(2, 5), O(1, 10), O(6, 7), O(160, 170), O(171, 180)]
    >>> partition(slist, max_dist=100)
    [[O(1, 10), O(2, 5), O(6, 7)], [O(160, 170), O(171, 180)]]

    >>> slist = [O(2, 5), O(1, 10), O(6, 7), O(26, 27), O(160, 170), O(171, 180)]
    >>> partition(slist, max_dist=10)
    [[O(1, 10), O(2, 5), O(6, 7)], [O(26, 27)], [O(160, 170), O(171, 180)]]

    """
    if not isinstance(slist, list):
        slist = list(slist)
    if len(slist) == 1: return [slist]
    slist.sort(key=operator.attrgetter('sstart'))
    lists = []
    current_list = [slist[0]]
    for a in slist[1:]:
        for b in current_list:
            if dist(a, b) < max_dist:
                current_list.append(a)
                break
        else:
            lists.append(current_list[:])
            current_list = [a]
    lists.append(current_list)
    return lists

def find_missed(sorg, qflat, sflat, q_snon_blast, q_s_blast,
                min_pct_coverage, min_len=30, flip=False):
    """ e.g.:
        >>> find_missed("papaya", "grape.flat", "papaya.flat",
        ...         "grape.features_vs_papaya.genomic.masked.blast",
        ...         "grape.features_vs_papaya.features.blast", 0.30) # doctest: +ELLIPSIS
        <generator object find_missed at ...>

    to find papaya mised exons. and:
        >>> find_missed("grape", "papaya.flat", "grape.flat",
        ...         "papaya.features_vs_grape.genomic.masked.blast",
        ...         "grape.features_vs_papaya.features.blast", 0.30, flip=True) # doctest: +ELLIPSIS
        <generator object find_missed at ...>

    to find grape missed exons.

    min coverage means if it finds a bunch of little spread out hits, it
    discards them unless they cover at least (for example) 0.4 == 40% of the
    total area. so hits like [105, 120], [190, 205] only cover 
       15 + 15 / (205 - 105) = 30%
    """

    # grouped_by_q has all the subject genomic hits mapped to the query feature.
    # and the subject start, stop are chromosomal positions.
    name = sorg + "_%(seqid)s_%(start)i_%(end)i"
    grouped_by_q = grouper(q_snon_blast)
    #print len(grouped_by_q)


    for qname, sdict in sorted(grouped_by_q.iteritems()):
        qfeat = qflat.accn(qname)
        qstrand = qfeat['strand'] == '-' and -1 or 1

        for schr, big_slist in sorted(sdict.iteritems()):
            slists = partition(big_slist)

            cover = 0.0
            for slist in slists:
                sstrand = slist[0].sstart < slist[0].sstop and 1 or -1
                if sstrand == 1:
                    cover += sum(x.sstop - x.sstart for x in slist)
                    sstart = min([x.sstart for x in slist])
                    sstop  = max([x.sstop for x in slist])
                else:
                    sstart = min([x.sstop for x in slist])
                    sstop  = max([x.sstart for x in slist])
                    cover += sum(x.sstart - x.sstop for x in slist)
                if abs(sstop - sstart) < min_len: continue

                # the hsps have to not be toooooo sparse
                if cover / (sstop - sstart) < min_pct_coverage: continue
                ##sstrand *= qstrand
                ##TODO create test file to confrim blastall strands
                sname = name % dict(seqid=schr, start=sstart, end=sstop)

                feat = dict(accn=sname , start=sstart , end=sstop
                           , seqid=schr , type="gene"
                           , strand= sstrand== 1 and "+" or "-", match=qfeat['accn'])

                # check if it's inside an existing gene. in which case, 
                # call it a cds and give the the same name as the parent.
                try:
                    # have to check here since we added extra for the case where we're inside an intron of an existing gene.
                    parents = [s for s in
                            sflat.get_features_in_region(feat['seqid'], feat['start'] -1000, feat['end'] + 1000)]
                except:
                    # this seqid doesnt have any features.
                    assert sflat[sflat['seqid'] == feat['seqid']] == []
                    yield feat
                    continue

                if len(parents) == 0:
                    yield feat
                else:
                    for parent in parents:
                        # have to check here since we added extra for the case where we're inside an intron of an existing gene.
                        parents['locs'].sort()
                        if parent['strand'] != feat['strand']: continue
                        if parent['locs'][0][0] < feat['start']: continue
                        if parent['locs'][-1][1] > feat['end']: continue
                        # when doing self-self dont want to add an annotation based
                        # on the same gene.
                        if parent['accn'] == qfeat['accn']:
                            continue

                        # here we see if it hit an already annotated CDS on another gene.
                        # http://toxic/CoGe/GEvo.pl?prog=blastn;spike_len=0;accn1=Bradi1g00480;fid1=35400621;dsid1=40124;dsgid1=1607;chr1=Bd1;dr1up=10000;dr1down=10000;ref1=1;accn2=Sb01g000240;fid2=19610158;dsid2=34580;dsgid2=93;chr2=1;dr2up=10000;dr2down=10000;ref2=1;num_seqs=2;hsp_overlap_limit=0;hsp_size_limit=0
                        do_break  = False
                        for start, end in parent['locs']:
                            l = end - start
                            if start - 0.1 * l <= feat['start'] and end + 0.1 * l >= feat['end']:
                                do_break = True
                        if do_break: break


                        feat['accn'] = parent["accn"]
                        feat['type'] = 'CDS'
                        # match
                        feat['strand'] = parent['strand']
                        #del feat['attrs']['match']
                        break
                    yield feat
                

def grouper(blast_file):
    """\
    group all subjects to a single query. so for
        grape.features_vs_papaya.genomic.masked.blast
    group all the papaya hits to the grape query"""
    g = collections.defaultdict(dict)
    for sline in open(blast_file):
        b = BlastLine(sline)
        # this removes low-copy transposons (length > 200, percent_id > 98)
        if b.pctid > 98.0 and b.hitlen > 200: continue
        if not b.subject in g[b.query]: g[b.query][b.subject] = []
        g[b.query][b.subject].append(b)
    return g

def fastas_for_features_vs_masked_genomic(afasta, bfasta):
    """
    >>> fastas_for_features_vs_masked_genomic("a.fasta", "b.fasta")
    ('a.features.fasta', 'b.genomic.masked.fasta')
    """
    a = afasta.replace(".fa", ".features.fa")
    b = bfasta.replace(".fa", ".genomic.masked.fa")
    return a, b

def run_blasts(config, test=False):
    """
    >>> cfg = simplejson.loads(open("tests/data/athaliana_athaliana.json").read())
    >>> run_blasts(cfg, test=True)

    """
    blast_cfg = config["blast"]
    blast_cfg.update({"m": 8, "p": "blastn"})

    out_dir = config["default"]["out_dir"]
    blast_cfg["o"] = out_dir
    blastall = config["default"]["blast_path"]
    # first blast a features to b genomic with features masked.
    agenes, bfasta =\
                fastas_for_features_vs_masked_genomic(config["genome_a"]["fasta"],
                                                      config["genome_b"]["fasta"])
    blast_cfg["i"] = agenes
    blast_cfg["d"] = bfasta

    if not test:
        bblast.blast(blast_cfg, blastall=blastall)

    bgenes, afasta =\
                fastas_for_features_vs_masked_genomic(config["genome_b"]["fasta"],
                                                      config["genome_a"]["fasta"])

    # then vice-versa
    if config["default"].get("reciprocal"):
        blast_cfg["i"] = bgenes
        blast_cfg["d"] = afasta

        if not test:
            bblast.blast(blast_cfg, blastall=blastall)
    

    # then genes to genes
    #blast_cfg["i"] = agenes
    #blast_cfg["d"] = bgenes
    #bblast.blast(blast_cfg, blastall=blastall)
    


if __name__ == "__main__":
    main(sys.argv[1])
