from itertools import product
import commands
from flatfeature import Bed
from find_cns import parse_blast
from processing import Pool
### need to download processing

def parse_dups(localdup_path):
    localdup_file = open(localdup_path)
    dup_dic ={}
    for line in localdup_file:
        dups = line.split("\t")
        parent = dups[0]
        child = dups[:]
        dup_dic[parent] = child
    return dup_dic
# create all possible combinations for that pair if len = 3

def get_dups(qfeat_dups,sfeat_dups,qbed,sbed):
    for qaccn,saccn in product(qfeat_dups,sfeat_dups):
        qfeat = qbed.accn(qaccn)
        sfeat = sbed.accn(saccn)
        yield qfeat,sfeat

# add bed files..
def get_pairs(pair_file,fmt,dup_dict):
    skipped = open('data/skipped_dups.txt', 'w')
    fh = open(pair_file)
    dups = []
    for line in fh:
        if line[0] == "#" : continue
        line = line.strip().split("\t")
        if fmt == 'dag':
            assert len(line) > 5, line
            pair = line[1], line[5]
        elif fmt in ('cluster', 'qa', 'raw'):
            assert len(line) == 5, line
            pair = line[1], line[3]
        elif fmt == 'pair':
            if len(line) == 1:
                line = line.split(",")
            assert len(line) >= 2, "dont know how to handle %s" % line
            pair = line[0], line[1]

        if fmt in ('qa', 'raw'):
            pair = int(pair[0]), int(pair[1])
        pair = tuple(pair)
        if pair in seen: continue
        seen[pair] = True
        if pair[0] in dup_dict.keys() or pair[1] in dup_dict.keys():
            dups.append((pair[0],pair[1]))
        else: continue

def get_all_dups(dup_dict,feat):
    try:
        dups = dup_dict[feat]
    except KeyError:
        dups = [feat]
    return dups

def main(qdups_path,sdups_path,pair_file,fmt,qbed,sbed,qpad,spad,blast_path,mask='F',ncpu=8):
    pool = Pool(ncpu)
    bl2seq = "%s " % blast_path + \
            "-p blastn -D 1 -E 2 -q -2 -r 1 -G 5 -W 7 -F %s " % mask + \
            " -e %(e_value).2f -i %(qfasta)s -j %(sfasta)s \
             -I %(qstart)d,%(qstop)d -J %(sstart)d,%(sstop)d | grep -v '#' \
             | grep -v 'WARNING' | grep -v 'ERROR' "

    qdup_dict = parse_dups(qdups_path)
    sdup_dict = parse_dups(sdups_path)
    dups = get_pairs(pair_file,fmt,dup_dict)
    for (qfeat,sfeat) in dup_pairs:
        qfeat_dups = get_all_dups(qdup_pairs,qfeat)
        sfeat_dups = get_all_dups(sdup_pairs,sfeat)

        pairs = [True]
        _get_dups_gen = get_dups(qfeat_dups,sfeat_dups,qbed,sbed)

        def get_dups_gen():
            try: return _get_dups_gen.next()
            except StopIteration: return None
        while any(pairs):
            cnss_dups = []
            pairs = [get_dups_gen() for i in range(ncpu)]

            def get_cmd(pair):
                if pair is None: return None
                qfeat, sfeat = pair
                qfasta = qfastas[qfeat['seqid']]
                sfasta = sfastas[sfeat['seqid']]

                qstart, qstop = max(qfeat['start'] - pad, 1), qfeat['end'] + pad
                sstart, sstop = max(sfeat['start'] - pad, 1), sfeat['end'] + pad

                assert qstop - qstart > 2 * pad or qstart == 1, (qstop, qstart)
                assert sstop - sstart > 2 * pad or sstart == 1, (sstop, sstart)

                cmd = bl2seq % dict(qfasta=qfasta, sfasta=sfasta, qstart=qstart,
                        sstart=sstart, qstop=qstop, sstop=sstop)
                return cmd, qfeat, sfeat
            cmds = [c for c in map(get_cmd, [l for l in pairs if l]) if c]
            results = (r for r in pool.map(commands.getoutput, [c[0] for c in cmds]))

            for res, (cmd, qfeat, sfeat) in zip(results, cmds):
                if not res.strip(): continue
                print >>sys.stderr,  "%s %s" % (qfeat["accn"], sfeat['accn'])
                orient = qfeat['strand'] == sfeat['strand'] and 1 or -1
                cnss = parse_blast(res, orient, qfeat, sfeat, qbed, sbed, pad)
                print >>sys.stderr, "(%i)" % len(cnss)
                cnss_dups["{0}_{1}".format(qfeat["accn"],sfeat["qaccn"])] = cnss
        #### get largest pair
        ### rewrite_file

if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser("usage: %prog [options] ")
    parser.add_option("-F", dest="mask", help="blast mask simple sequence [default: F]", default="F")
    parser.add_option("-n", dest="ncpu", help="parallelize to this many cores", type='int', default=8)
    parser.add_option("--qbed", dest="qbed", help="query bed file WITH localdups")
    parser.add_option("--sbed", dest="sbed", help="subject bed file WITH localdups")
    parser.add_option("-p", dest="pairs", help="the pairs file. output from dagchainer")
    choices = ("dag", "cluster", "pair", 'qa', 'raw')
    parser.add_option("--pair_fmt", dest="pair_fmt", default='raw',
            help="format of the pairs, one of: %s" % str(choices),
            choices=choices)
    parser.add_option("--qpad", dest="qpad", type='int', default=12000, help="how far from the end of the query gene to look for cnss")
    parser.add_option("--spad", dest="spad", type='int', default=12000, help="how far from the end of the subject gene to look for cnss")
    parser.add_option("--blast_path", dest="blast_path", type='string', help="path to bl2seq")
    parser.add_option("--qdups", dest="qdups", type='string', help="path to query localdup_file")
    parser.add_option("--sdups", dest="sdups", type='string', help="path to subject localdup_file")
    (options, _) = parser.parse_args()

    qbed = Bed(options.qbed)
    sbed = Bed(options.sbed)
    assert options.mask in 'FT'
    
    main(options.qdups,options.sdups,options.pairs,options.pair_fmt,qbed,sbed,options.qpad,options.spad,options.blast_path,options.mask,options.ncpu)
