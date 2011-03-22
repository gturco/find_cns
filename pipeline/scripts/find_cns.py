# warrning make sure you delete the past rice_rice_split file before re-running if u uploaded a new file!!!!!!
import sys
import os
import os.path as op
import numpy as np
import commands
from shapely.geometry import Point, Polygon, LineString, MultiLineString
from flatfeature import Bed
import pickle 

from processing import Pool
pool = None


EXPON = 0.90


def retained_cnss(qfeat, sfeat, fbed, sfastas, cnss, mask='T'):
    """makes a dict of the seq_3 start and end and fasta along with the cns start and end for bl2seq
    returns a list of the hight scoring cns in seq3"""
    accn = qfeat['ORG2_qfeat']
    feat = fbed.accn(accn)
    feat_start = feat['start'] - 15000
    feat_stop = feat['end'] + 15000
    feat_fastas= get_masked_fastas(fbed)
    feat_fasta = feat_fastas[feat['seqid']]
    sfasta = sfastas[sfeat['seqid']]
    
    bl2seq = "/usr/bin/bl2seq " \
           "-p blastn -D 1 -E 2 -q -2 -r 1 -G 5 -W 7 -F %s " % mask + \
           " -Y 812045000 -d 26195 -e 2.11 -i %(feat_fasta)s -j %(sfasta)s \
              -I %(feat_start)d,%(feat_stop)d -J %(sstart)d,%(sstop)d | grep -v '#' \
            | grep -v 'WARNING' | grep -v 'ERROR' "
    
    for cns in cnss:
       print cns
       cns_start = cns[2]
       cns_stop  = cns[3]
       cmd = bl2seq % dict(feat_fasta=feat_fasta, sfasta=sfasta, feat_start=feat_start,
                           sstart=cns_start, feat_stop=feat_stop, sstop=cns_stop)
       print cmd
       retained_cns = (commands.getoutput(cmd))
       for line in retained_cns.split("\n"):
           if "WARNING:" in line: continue
           if "ERROR" in line: continue
           line = line.split("\t")
           seq3_cns = map(int, line[6:8])
           print >> sys.stdout, seq3_cns
           if len(seq3_cns) == 0: continue
           url = url_params(cns, qfeat['seqid'], sfeat['seqid'], feat['seqid'], seq3_cns)
           fcnss = sys.stdout
           print >> fcnss, "# qaccn,[qleft_gene,qright_gene],qseqid,saccn,sseqid,cns,url"#"#qseqid,qaccn,sseqid,saccn,[qstart,qend,sstart,send...]"
           print >> fcnss, "%s,[%s,%s],%s,%s,%s,%s,%s" %  (qfeat['accn'], qfeat['qleft_gene'], qfeat['qright_gene'], qfeat['seqid'], sfeat['accn'], sfeat['seqid'],cns, url)
           


def assign_url(params,
               base = "http://synteny.cnr.berkeley.edu/CoGe/GEvo.pl?prog=blastn&autogo=1&"):
    "lines up coge based on the cns postion"
    inside = 'dsid1=43388&dsgid1=9109&chr1=%(qseqid)s&x1=%(seq1)s&dr1up=1000&dr1down=1000&dsid2=43388&dsgid2=9109&chr2=%(sseqid)s&x2=%(seq2)s&dr2up=1000;dr2down=1000&\
dsid3=34580;dsgid3=34580;chr3=%(fseqid)s;x3=%(seq3)s;dr3up=1000;dr3down=1000;num_seqs=3;hsp_overlap_limit=0;hsp_size_limit=0' %params
    url = base + inside
    return url


def url_params(cns, qseqid, sseqid, fseqid, seq_3):
    params={}
    params['qseqid'] = qseqid
    params['sseqid'] = sseqid
    params['fseqid'] = fseqid
    params['seq1']= cns[0]
    params['seq2'] = cns[2]
    params['seq3'] = seq_3[0]
    url = assign_url(params)
    return url

def get_feats_in_space(locs, ichr, bpmin, bpmax, bed):
    """ locs == [start, stop]
    bpmin is the lower extent of the window, bpmax ...
    """
    assert bpmin < bpmax, (locs, ichr, bpmin, bpmax)
    feats = bed.get_features_in_region(str(ichr), bpmin, bpmax)
    feats = [f for f in feats if not (f['start'] == locs[0] and f['end'] == locs[1])]
    if len(feats) != 0:
        assert feats[0]['seqid'] == str(ichr)
    return [(f['start'], f['end'], f['accn']) for f in feats]

def parse_blast(blast_str, orient, qfeat, sfeat, qbed, sbed):
    blast = []
    slope = orient

    qgene = [qfeat['start'], qfeat['end']]
    sgene = [sfeat['start'], sfeat['end']]
    # qcds = qfeat['locs']
    # scds = sfeat['locs']


    sgene = sgene[::slope]
    # center = sum(qgene)/2., sum(sgene)/2.
    # 
    # EXP = EXPON
    # if abs(abs(qgene[1] - qgene[0]) - abs(sgene[1] - sgene[0])) > 3000:
    #     EXP = 0.94
    # 
    # 
    # #intercept = (sgene[0] + sgene[1])/2.  - slope * (qgene[0] + qgene[1])/2.
    # intercept = center[1] - slope * center[0]
    # rngx = qgene[1] - qgene[0]
    # rngy = abs(sgene[1] - sgene[0])
    # 
    # x = np.linspace(qgene[0] - pad, qgene[1] + pad, 50)
    # y = slope * x + intercept
    # 
    # 
    # xb = x + -slope * rngx/3. + -slope * np.abs(x - center[0])**EXP
    # yb = y + rngy/3. + np.abs(y - center[1])**EXP
    # 
    # xy = x + slope * rngx/3. + slope * np.abs(x - center[0])**EXP
    # yy = y - rngy/3. - np.abs(y - center[1])**EXP
    # 
    # if slope == 1:
    #     xall = np.hstack((xy[::-1], xb[::slope], xy[-1]))
    #     yall = np.hstack((yy[::-1],yb, yy[-1]))
    # if slope == -1:
    #     xall = np.hstack((xy, xb[::-1], xy[0]))
    #     yall = np.hstack((yy,yb[::-1], yy[0]))
    sstrat , sstop = grab_flanking_region(sfeat, qfeat)
    feats_nearby = {}
    feats_nearby['q'] = get_feats_in_space(qgene, qfeat['seqid'], qfeat['start'] ,qfeat['end'], qbed) # changed so that if looks for genes within region
    feats_nearby['s'] = get_feats_in_space(sgene, sfeat['seqid'], sstrat, sstop, sbed) #looks for genes in bowtie.....
    
    
    
    # genespace_poly = Polygon(zip(xall, yall))
    
    for sub in ('q', 's'):
        if len(feats_nearby[sub]) !=0:
            feats_nearby[sub] = MultiLineString([[(0, c0),(0, c1)] for c0, c1, fname in feats_nearby[sub]])
        else:
            feats_nearby[sub] = None
    
    cnss = set([])
    
    qgene_poly = LineString([(0.0, qgene[0]), (0.0, qgene[1])])
    sgene_poly = LineString([(0.0, sgene[0]), (0.0, sgene[1])])
    intronic_removed = 0

    for line in blast_str.split("\n"):
        if "WARNING:" in line: continue
        if "ERROR" in line: continue
        line = line.split("\t")
        locs = map(int, line[6:10])
        locs.extend(map(float, line[10:]))

        xx = locs[:2]
        yy = locs[2:4]

        # get rid of stuff on the wrong strand
        try:
            if slope == 1 and locs[2] > locs[3]: continue
            if slope == -1 and locs[2] < locs[3]: continue
        except:
            print >>sys.stderr, blast_str
            raise

        # to be saved. a hit must either be in an intron in both
        # genes, or in neither.

        ##########################################################
        # DEAL WITH INTRONIC cnss in the gene of interest.
        ##########################################################
        xls = LineString([(0, locs[0]), (0, locs[1])])
        yls = LineString([(0, locs[2]), (0, locs[3])])

        locs = tuple(locs) # make it hashable.
        # if not sgene_poly.intersects(yls):
        #     cnss.update((locs,))
        #     continue

        if  sgene_poly.intersects(yls):
            intronic_removed += 1
            continue

        ##########################################################

        ###############################################################
        # for all other genes, if it's in an intron, we dont keep it.
        ###############################################################
        intronic = False
        # get rid of stuff that overlaps another gene:
        for sub, (start, stop) in (('q', locs[:2]), ('s', locs[2:4])):
            feats = feats_nearby[sub]
            if feats is None: continue
            # the current hsp is overlapping another gene. we dont want that...
            if feats.contains(Point(0, start)) or feats.contains(Point(0, stop)):
                intronic = True
                break

        if intronic: continue #(if an intron dont updat cns... go back to results(forloop) and check another intron)

        ##########################################################

        # this is the bowtie.
        # if not genespace_poly.contains(LineString(zip(xx, yy))): continue #if it is in the bowtie update... otherwise get rid of it
        cnss.update((locs,))

    # cant cross with < 2 cnss.
    # get rid of the eval, bitscore stuff.
    if len(cnss) < 2: return [l[:4] for l in cnss]

    cnss = list(cnss)
    # need to flip to negative so the overlapping stuff still works.
    if orient == -1:
        cnss = list(cnss)
        for i, cns in enumerate(cnss):
            cns = list(cns)
            cns[2] *= - 1
            cns[3] *= - 1
            cnss[i] = tuple(cns)
        sgene[0] *= -1
        sgene[1] *= -1

    cnss = [l[:4] for l in remove_crossing_cnss(cnss, qgene, sgene)]
    if orient == -1:
        cnss = [(c[0], c[1], -c[2], -c[3]) for c in cnss]
    return cnss


def remove_overlapping_cnss(cnss):
    """for cases when there is nearly the same cns, but with 1
    basepair shfit up/down. that create many cnss stacked on top
    of each other. this reduces those down to one."""
    qcnss = [LineString([(0, cns[0]), (0, cns[1])]) for i, cns in enumerate(cnss)]
    scnss = [LineString([(0, cns[2]), (0, cns[3])]) for i, cns in enumerate(cnss)]

    remove = []
    for zcnss in (qcnss, scnss):
        for i, csi in enumerate(zcnss[:-1]):
            for _j, csj in enumerate(zcnss[i + 1:]):
                j = i + _j + 1 # cause enumerate starts at 0
                if csi.overlaps(csj):
                    if cnss[i][-2] < cnss[j][-2] or cnss[i][-1] > cnss[j][-2] or csi.y < csj.y:
                        remove.append(j)
                    else:
                        remove.append(i)
    remove = frozenset(remove)
    return [cns for i, cns in enumerate(cnss) if not i in remove]


def remove_crossing_cnss(cnss, qgene, sgene):
    diff = (sum(qgene)/2.) - (sum(sgene)/2.) # adjust subject so it's in same range as query
    cns_shapes = [LineString([((c[0] + c[1])/2., 0 ), ((c[2] + c[3])/2. + diff, 1000)]) for c in cnss]

    overlapping = len(cnss)
    cnss = remove_overlapping_cnss(cnss)
    overlapping -= len(cnss)
    cns_shapes = [LineString([((c[0] + c[1])/2., 0 ), ((c[2] + c[3])/2. + diff, 1000)]) for c in cnss]


    # and save a reference to the orginal cnss as that's the data we want.
    for cs, cns in zip(cns_shapes, cnss):
        cs.cns = cns
        # hold the number of times an hsp crosses any other.
        cs.cross_list = set([])
        # mark for removal.
        cs.do_remove = False


    for csi in cns_shapes:
        for csj in cns_shapes:
            if csi == csj: continue
            if csi.crosses(csj):
                csi.cross_list.update([csj])
                csj.cross_list.update([csi])

    ######################################################################
    # first remove anything that cross more than 5 other cnss.
    ######################################################################
    nremoved = 0
    any_removed = True
    while any_removed:
        # need this outer loop to refresh the sorting.
        cns_shapes = sorted(cns_shapes, reverse=True, cmp=lambda a, b: cmp(len(a.cross_list), len(b.cross_list)))[:]
        any_removed = False
        for i, cs in enumerate(cns_shapes):
            if len(cs.cross_list) > 3:
                # remove this from all other lists as it's a bad guy.
                for crossed in cs.cross_list:
                    crossed.cross_list.difference_update(set([cs]))
                cs.do_remove = True
                any_removed = True
                nremoved += 1
                del cns_shapes[i]
                break

    ######################################################################
    # then remove crosses one-by-one, keeping the < evalue, > bitscore.
    ######################################################################
    for csi in cns_shapes:
        if csi.do_remove or len(csi.cross_list) == 0: continue
        for csj in cns_shapes:
            if csj.do_remove or len(csj.cross_list) == 0: continue
            if csi.do_remove or len(csi.cross_list) == 0: continue
            if csi.crosses(csj):
                # access the assocated cns.
                # evalue: less is better       bitscore: more is better
                if csi.cns[-2] < csj.cns[-2] or csi.cns[-1] > csj.cns[-1]:
                    csj.do_remove = 1
                    map(lambda crossed: crossed.cross_list.difference_update(set([csj])), csj.cross_list)

                else:
                    csi.do_remove = 1
                    map(lambda crossed: crossed.cross_list.difference_update(set([csi])), csi.cross_list)
                    break

    for c in cns_shapes:
        if not c.do_remove: continue
        nremoved += 1
    return [c.cns for c in cns_shapes if not c.do_remove]

def grab_flanking_region(sfeat , flanking_genes):
    "grabs the start and end postion of the nearest gene to the left \
    and right of the sfeat"
    left_padding = (sfeat['start'] - 12000)
    right_padding = (sfeat['end'] + 12000)
    left_gene_end  = flanking_genes['sstart']
    right_gene_start = flanking_genes['send']
    if left_gene_end < left_padding and right_gene_start > right_padding:
        return left_padding, right_padding
    else:
        return left_gene_end, right_gene_start

                    
def get_pair(regions , sbed):
    "grabs the pairs from the region file"
    # pairs = []
    file= open(regions, "r")
    region_dict = pickle.load(file)
    for row in region_dict:
        region = row
        accn = row['sfeat']
        sfeat = sbed.accn(accn)
        pair = region, sfeat
        yield pair 
    #     pairs.append(pair)
    # return pairs

def get_masked_fastas(bed):
    """
    create the masked fasta files per chromosome. needed to run bl2seq.
    """
    f = bed.fasta.fasta_name
    fname = op.splitext(op.basename(f))[0]
    d = op.dirname(f) + "/%s_split" % fname
    try: os.mkdir(d)
    except OSError: pass

    fastas = {}
    for seqid, seq in bed.mask_cds():
        f = d + "/%s.fasta" % seqid
        fastas[seqid] = f
        if op.exists(f): continue
        fh = open(f, "wb")
        print >>fh, seq
        fh.close()
    return fastas

def main(qbed, sbed, fbed, pairs_file, mask='F', ncpu=8):
    """main runner for finding cnss"""
    pool = Pool(options.ncpu)


    bl2seq = "/usr/bin/bl2seq " \
           "-p blastn -D 1 -E 2 -q -2 -r 1 -G 5 -W 7 -F %s " % mask + \
           " -Y 812045000 -d 26195 -e 2.11 -i %(qfasta)s -j %(sfasta)s \
              -I %(qstart)d,%(qstop)d -J %(sstart)d,%(sstop)d | grep -v '#' \
            | grep -v 'WARNING' | grep -v 'ERROR' "

#    fcnss = sys.stdout
#    print >> fcnss, "# qaccn,res,urls"#"#qseqid,qaccn,sseqid,saccn,[qstart,qend,sstart,send...]"

    qfastas = get_masked_fastas(qbed)
    sfastas = get_masked_fastas(sbed) if qbed.filename != sbed.filename else qfastas



    pairs = [True]
    _get_pair_gen = get_pair(pairs_file , sbed)
    # need this for parallization stuff.
    def get_pair_gen():
        try: return _get_pair_gen.next()
        except StopIteration: return None

    while any(pairs):
        pairs = [get_pair_gen() for i in range(ncpu)]

        # this helps in parallelizing.
        def get_cmd(pair):
            if pair is None: return None
            qfeat, sfeat = pair
            
            #if qfeat['accn'] != "Bradi4g01820": return None
            #print >>sys.stderr, qfeat, sfeat

            qfasta = qfastas[qfeat['seqid']]
            sfasta = sfastas[sfeat['seqid']]

            qstart, qstop = qfeat['start'], qfeat['end']
            sstart, sstop = grab_flanking_region(sfeat, qfeat) # qfeat here is the final table with sfeat info from qfeat dict

            # assert qstop - qstart > 2 * pad or qstart == 1, (qstop, qstart)
            # assert sstop - sstart > 2 * pad or sstart == 1, (sstop, sstart)

            cmd = bl2seq % dict(qfasta=qfasta, sfasta=sfasta, qstart=qstart,
                                sstart=sstart, qstop=qstop, sstop=sstop)
            return cmd, qfeat, sfeat

        cmds = [c for c in map(get_cmd, [l for l in pairs if l]) if c]
        results = (r for r in pool.map(commands.getoutput, [c[0] for c in cmds]))
        #results = (r for r in map(commands.getoutput, [c[0] for c in cmds]))

        for res, (cmd, qfeat, sfeat) in zip(results, cmds):
            if not res.strip(): continue
            print >>sys.stderr,  "%s %s" % (qfeat["accn"], sfeat['accn']),
            orient = qfeat['strand'] == sfeat['strand'] and 1 or -1
            
            cnss =  parse_blast(res, orient, qfeat, sfeat, qbed, sbed)
            print >>sys.stderr, "(%i)" % len(cnss)
            if len(cnss) == 0: continue
            retained_cnss(qfeat, sfeat, fbed, sfastas, cnss, mask)
            
#            qname, sname = qfeat['accn'], sfeat['accn']
            
#       urls = url_params(cnss, qfeat['seqid'], sfeat['seqid'], qfeat['ORG2_qfeat'])
            
#            print >> fcnss, "%s,[%s,%s],%s,%s,%s,%s,%s" % (qname, qfeat['qleft_gene'], qfeat['qright_gene'], qfeat['seqid'], sname, sfeat['seqid'],
#                             ",".join(map(lambda l: ",".join(map(str,l)), cnss)), ",".join(urls))


    return None

if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser("usage: %prog [options] ")
    parser.add_option("-F", dest="mask", help="blast mask simple sequence [default: F]", default="F")
    parser.add_option("-n", dest="ncpu", help="parallelize to this many cores", type='int', default=8)
    parser.add_option("-q", dest="qfasta", help="path to genomic query fasta")
    parser.add_option("--qbed", dest="qbed", help="query bed file")
    parser.add_option("-s", dest="sfasta", help="path to genomic subject fasta")
    parser.add_option("--sbed", dest="sbed", help="subject bed file")
    parser.add_option("--fbed", dest="fbed", help="retained feauture bed file")
    parser.add_option("-f", dest="ffasta", help="retained feauture fasta file")
    parser.add_option("-p", dest="pairs", help="the pairs file. output from dagchainer")
    (options, _) = parser.parse_args()


    if not (options.qfasta and options.sfasta and options.sbed and options.qbed):
        sys.exit(parser.print_help())

    qbed = Bed(options.qbed, options.qfasta); qbed.fill_dict()
    sbed = Bed(options.sbed, options.sfasta); sbed.fill_dict()
    fbed = Bed(options.fbed, options.ffasta); sbed.fill_dict()
    assert options.mask in 'FT'

    main(qbed, sbed, fbed, options.pairs, options.mask, options.ncpu)
