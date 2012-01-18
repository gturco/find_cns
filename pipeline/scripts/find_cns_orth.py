import sys
import os
import os.path as op
import numpy as np
import commands
from shapely.geometry import Point, Polygon, LineString, MultiLineString
from flatfeature import Bed
import logging
from processing import Pool
#LOG_FILENAME = '/Users/gturco/find_cns.log'
#logging.basicConfig(filename=LOG_FILENAME,level=logging.INFO)
pool = None


EXPON = 0.90

def write_cns_bed(cns_bed,cns_path):
    """creates cns bed file for qcns"""
    write_file = open(cns_bed,"wb")
    for cns in csv.DictReader(open(cns_path)):
       new_line = '{0}\t{1}\t{2}\t{3}\t.\t{4}\t-1\t-1\t-1\t1\t{5}\t0'.format(cns["qchr"],
        int(cns["qstart"]), int(cns["qstop"]),
        cns["cns_id"],cns["qstrand"], abs(int(cns["qstart"]) -
        int(cns["qstop"])))
        write_file.write(new_line)
    write_file.close()

def replace_pairs_cns(cns_bed,pair_file_path):
    """ appends all postible cns with sfeat to the end """
    pair_file = open(pair_file_path,"a")
    pair_list = []
    for line in open(pair_file_path):
        if line[0] == "#": continue
        line = line.strip().split("\t")
        pair_list.append(line)
    for qaccn,saccn in pair_list:
        qfeat = qbed.accn(qaccn)
        cnss= cns_bed.get_features_in_region(qfeat["seqid"],qfeat["seqid"], abs(qfeat['start']-12000),qfeat['end'] + 12000)
        #"sed -i '/{0}\t{1}/ d'".format(qfeat,sfeat)
        for cns in cnss:
            new_pair = "{0}\t{1}\n".format(cns["accn"],saccn)
            pair_file.write(new_pair)
    pair_file.close()

def parse_blast(blast_str, orient, qfeat, sfeat, qbed, sbed, qpad, spad):
    blast = []
    slope = orient

    qgene = [qfeat['start'],qfeat['end']]
    sgene = [sfeat['start'], sfeat['end']]
    sgene = sgene[::slope]
    
    cnss = set([])
    for line in blast_str.split("\n"):
        if "WARNING:" in line: continue
        if "ERROR" in line: continue
        line = line.split("\t")
        if float(line[-1]) < 29.5: continue #finds 15/15 match
        locs = map(int, line[6:10])
        locs.extend(map(float, line[10:]))

        # get rid of stuff on the wrong strand
        try:
            if slope == 1 and locs[2] > locs[3]: continue
            if slope == -1 and locs[2] < locs[3]: continue
        except:
            print >>sys.stderr, blast_str
            raise

        cnss.update((locs,))

    cnss = list(cnss)
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
                    logging.info("overlapping:{0}".format(csi))
                if csi.contains(csj):
                    if cnss[i][-2] < cnss[j][-2] or cnss[i][-1] > cnss[j][-2] or cnss[i][1] < cnss[j][3]:
                        remove.append(j)
                    else:
                        remove.append(i)
                    logging.info("intersecting:{0}".format(cnss[i]))
                    #print >> sys.stderr, csi
    remove = frozenset(remove)
    return [cns for i, cns in enumerate(cnss) if not i in remove]


def remove_crossing_cnss(cnss, qgene, sgene):
    diff = qgene[0] - sgene[0] # adjust subject so it's in same range as query
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


def get_pair(pair_file, fmt, qbed, sbed, seen={}):
    """ read a line and make sure it's unique handles
    dag, cluster, and pair formats."""
    skipped = open('data/skipped.txt', 'w')
    fh = open(pair_file)
    for line in open(pair_file):
        if line[0] == "#": continue
        line = line.strip().split("\t")
        if fmt == 'dag':
            assert len(line) > 5, line
            pair = line[1], line[5]
        elif fmt in ('cluster', 'qa', 'raw'):
            assert len(line) == 5, line
            pair = line[1], line[3]
        elif fmt == 'pair':
            if len(line) == 1:
                line = line[0].split(",")
            assert len(line) >= 2, "dont know how to handle %s" % line
            pair = line[0], line[1]

        if fmt in ('qa', 'raw'):
            pair = int(pair[0]), int(pair[1])
        pair = tuple(pair)
        if pair in seen:
            continue
        seen[pair] = True
        try:
            if isinstance(pair[0], (int, long)):
                yield qbed[pair[0]], sbed[pair[1]]
            else:
                yield qbed.d[pair[0]], sbed.d[pair[1]]
        except KeyError, IndexError:
            print >>skipped, "%s\t%s" % pair
            print >>sys.stderr, "skipped %s %s" % pair
            continue

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

def main(qbed, sbed,cns_bed, pairs_file, qpad, spad, pair_fmt, blast_path, mask='F', ncpu=8):
    """main runner for finding cnss"""
    pool = Pool(ncpu)


    bl2seq = "%s " % blast_path + \
           "-p blastn -D 1 -E 2 -q -2 -r 1 -G 5 -W 7 -F %s " % mask + \
           " -e %(e_value).2f -i %(qfasta)s -j %(sfasta)s \
              -I %(qstart)d,%(qstop)d -J %(sstart)d,%(sstop)d | grep -v '#' \
            | grep -v 'WARNING' | grep -v 'ERROR' "

    fcnss = sys.stdout
    print >> fcnss, "#qseqid,qaccn,sseqid,saccn,[qstart,qend,sstart,send,evalue...]"

    qfastas = get_masked_fastas(qbed)
    sfastas = get_masked_fastas(sbed) if qbed.filename != sbed.filename else qfastas

    pairs = [True]
    _get_pair_gen = get_pair(pairs_file, pair_fmt, cns_bed, sbed)
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

            qstart, qstop = max(qfeat['start'] - qpad, 1), qfeat['end'] + qpad
            sstart, sstop = max(sfeat['start'] - spad, 1), sfeat['end'] + spad

            assert qstop - qstart > 2 * qpad or qstart == 1, (qstop, qstart)
            assert sstop - sstart > 2 * spad or sstart == 1, (sstop, sstart)
            
            #m = qstop - qstart
            #n = sstop - sstart
            #e_value = m*n*(2**(-28.51974)) # bit score above 15/15 noise
            #assert e_value > 0

            cmd = bl2seq % dict(qfasta=qfasta, sfasta=sfasta, qstart=qstart,
                                sstart=sstart, qstop=qstop, sstop=sstop,
                                e_value=30)
            return cmd, qfeat, sfeat

        cmds = [c for c in map(get_cmd, [l for l in pairs if l]) if c]
        results = (r for r in pool.map(commands.getoutput, [c[0] for c in cmds]))
        #results = (r for r in map(commands.getoutput, [c[0] for c in cmds]))

        for res, (cmd, qfeat, sfeat) in zip(results, cmds):
            if not res.strip(): continue
            print >>sys.stderr,  "%s %s" % (qfeat["accn"], sfeat['accn']),
            orient = qfeat['strand'] == sfeat['strand'] and 1 or -1

            cnss = parse_blast(res, orient, qfeat, sfeat, cns_bed, sbed, qpad, spad)
            print >>sys.stderr, "(%i)" % len(cnss)
            if len(cnss) == 0: continue

            qname, sname = qfeat['accn'], sfeat['accn']
            print >> fcnss, "%s,%s,%s,%s,%s" % (qfeat['seqid'], qname, sfeat['seqid'], sname,
                             ",".join(map(lambda l: ",".join(map(str,l)),cnss)))

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
    parser.add_option("-p", dest="pairs", help="the pairs file. output from dagchainer")
    choices = ("dag", "cluster", "pair", 'qa', 'raw')
    parser.add_option("--pair_fmt", dest="pair_fmt", default='raw',
                      help="format of the pairs, one of: %s" % str(choices),
                      choices=choices)
    parser.add_option("--qpad", dest="qpad", type='int', default=12000, help="how far from the end of the query gene to look for cnss")
    parser.add_option("--spad", dest="spad", type='int', default=26000, help="how far from the end of the subject gene to look for cnss")
    parser.add_option("--blast_path", dest="blast_path", type='string', help="path to bl2seq")
    parser.add_option("--cns_bed",dest="cns_bed",help="cns bed path")
    parser.add_option("--cns",dest="cns", help="final cns.csv assigned file")
    (options, _) = parser.parse_args()


    if not (options.qfasta and options.sfasta and options.sbed and options.qbed):
        sys.exit(parser.print_help())

    write_cns_bed(options.cns_bed,options.cns)

    qbed = Bed(options.qbed, options.qfasta); qbed.fill_dict()
    sbed = Bed(options.sbed, options.sfasta); sbed.fill_dict()
    cns_bed = Bed(options.cns_bed)
    assert options.mask in 'FT'

    replace_pairs_cns(cns_bed,options.pair_file)
    main(qbed, sbed,cns_bed, options.pairs, options.qpad, options.spad, options.pair_fmt, options.blast_path, options.mask, options.ncpu)
