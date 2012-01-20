from itertools import product
import commands
from bed_utils import Bed as Orderbed
from flatfeature import Bed
from find_cns import parse_blast, get_masked_fastas
from processing import Pool
import sys
#from tempfile import mkstemp
#from shutil import move
from os import remove,close
#####sed '/Os12g41090\tSb08g020690/ d' -i r

### need to download processing

def parse_dups(localdup_path):
    localdup_file = open(localdup_path)
    dup_dic ={}
    for line in localdup_file:
        dups = line.strip().split("\t")
        parent = dups[0]
        child = dups[:]
        dup_dic[parent] = child
    localdup_file.close()
    return dup_dic
# create all possible combinations for that pair if len = 3

def get_dups(qfeat_dups,sfeat_dups,qbed,sbed):
    for qaccn,saccn in product(qfeat_dups,sfeat_dups):
        qfeat = qbed.accn(qaccn)
        sfeat = sbed.accn(saccn)
        if len(list(product(qfeat_dups,sfeat_dups))) > 9: yield None
        else: yield qfeat,sfeat

# add bed files..
def get_pairs(pair_file,fmt,qdup_dict,sdup_dict):
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
        if pair[0] in qdup_dict.keys() or pair[1] in sdup_dict.keys():
            dups.append((pair[0],pair[1]))
        else: continue
    return dups

def get_all_dups(dup_dict,feat):
    try:
        dups = dup_dict[feat]
    except KeyError:
        dups = [feat]
    return dups

####need to remove line from:#######
####################################
### pairs file sed 's/gene1\tgene2/genenew1\tgenenew2/' -i file_name ###### search and replace!!!!!!!
### raw.filtered file sed 's/qchr\tqpos\tsch\tspos/qchr\t/qfeat' -i file_name search and replace!!!!
### cns.txt file '/qchr,qaccn,sch,saccn/ d' -i file_name order doesnt matter
### .nolocaldups '/accn/ d' -i file_name -sort sort -n -k 1 -k 2 rice_v6.all2.bed  > rice_v6.all3.bed
### check sort for them

def update_pairs(qfeat,sfeat,qparent,sparent,pair_file):
    """replaces the old pair with the new head loaldup pairs"""
    #### new pairs file search and replace
    search_replace_pairs = "sed 's/{0}\t{1}/{2}\t{3}/' -i {4} ".format(qparent,sparent,qfeat,sfeat,pair_file)
    commands.getstatusoutput(search_replace_pairs)

def update_nolocaldups(bed,qfeat,sfeat,qnolocaldups_path,snolocaldups_path):
    """ removes the old head localdups and appends the new local dups"""
    qnolocaldups = open(qnolocaldups_path,'a')
    snolocaldups = open(snolocaldups_path,'a')
    qline = "{0}\n".format(bed.row_string(qfeat))
    sline = "{0}\n".format(bed.row_string(sfeat))
    qnolocaldups.write(qline)
    snolocaldups.write(sline)
    qnolocaldups.close()
    snolocaldups.close()
  
def remove_cnss_line(qbed,sbed,qparent,sparent,cns_file):
    """ removes any cnss with the old parent dup """
    remove_cnss = "sed '/{0},{1},{2},{3}/ d' -i {4}".format(qbed.accn(qparent)['seqid'],qparent,sbed.accn(sparent)['seqid'],sparent,cns_file)
    commands.getstatusoutput(remove_cnss)

def localdup_file(qparent,sparent,qfile,sfile,neworder):
    """ replaces the orginal parent local dup with the new order"""
    
    qdups = [qdup for cns_number,q_start,s_start,qdup,sdup,cns in neworder]
    sdups = [sdup for cns_number,q_start,s_start,qdup,sdup,cns in neworder]
    qreplace = "\t".join(qdups)
    sreplace = "\t".join(sdups)
    qsearch_replace = "sed -i 's/^{0}.*/{1}/g' {2}".format(qparent,qreplace,qfile)
    ssearch_replace = "sed -i 's/^{0}.*/{1}/g'{2}".format(sparent,sreplace,sfile)
    commands.getstatusoutput(qsearch_replace)
    commands.getstatusoutput(ssearch_replace)

def pairs_to_qa(pair_file,qbed_file,sbed_file):
    """takes new localdups file and new pairs file to create qa file"""
    ###-sort sort -n -k 1 -k 2 rice_v6.all2.bed  > rice_v6.all3.bed 
    header = pairs_file.split(".")[0]
    new_qa = open("{0}.raw.filtered".format(header),"wb")
    qbed = Orderbed(qbed_file)
    sbed = Orderbed(sbed_file)
    qorder = qbed.get_order()
    sorder = sbed.get_order()
    fh = open(pair_file)
    dups = []
    for line in fh:
        if line[0] == "#" : continue
        line = line.strip().split("\t")
        qfeat,sfeat = line
        qpos = qorder[qfeat][0]
        qchr = qorder[qfeat][1].seqid
        spos = sorder[sfeat][0]
        schr = sorder[sfeat][1].seqid
        new_line = "{0}\t{1}\t{2}\t{3}\t50\n".format(qchr,qpos,schr,spos)
        new_qa.write(new_line)
    fh.close()
    new_qa.close()

def main(cns_file,qdups_path,sdups_path,pair_file,fmt,qbed,sbed,qpad,spad,blast_path,mask='F',ncpu=8):
    pool = Pool(ncpu)
    bl2seq = "%s " % blast_path + \
            "-p blastn -D 1 -E 2 -q -2 -r 1 -G 5 -W 7 -F %s " % mask + \
            " -e %(e_value).2f -i %(qfasta)s -j %(sfasta)s \
             -I %(qstart)d,%(qstop)d -J %(sstart)d,%(sstop)d | grep -v '#' \
             | grep -v 'WARNING' | grep -v 'ERROR' "
    
    qfastas = get_masked_fastas(qbed)
    sfastas = get_masked_fastas(sbed) if qbed.filename != sbed.filename else qfastas

    qdup_dict = parse_dups(qdups_path)
    sdup_dict = parse_dups(sdups_path)
    dups = get_pairs(pair_file,fmt,qdup_dict,sdup_dict)
    small_dups = []
    large_dups = []
    for qfeat,sfeat in dups:
        if len(get_all_dups(qdup_dict,qfeat)) < 4 and len(get_all_dups(sdup_dict,sfeat)) < 4:
            small_dups.append((qfeat,sfeat))
        else:
            large_dups.append((qfeat,sfeat))

    qnolocaldups_path =  qbed.path.split(".")[0] + ".nolocaldups.bed"
    snolocaldups_path = sbed.path.split(".")[0] + ".nolocaldups.bed"
    for (qparent,sparent) in small_dups:
        remove_cnss_line(qbed,sbed,qparent,sparent,cns_file)
        if qparent not in [qdup for qdup,sdup in large_dups]:
            remove_qparent = "sed '/{0}/ d' -i {1}".format(qparent,qnolocaldups_path)
            x = commands.getstatusoutput(remove_qparent)
        if sparent not in [sdup for qdup,sdup in large_dups]:
            remove_sparent = "sed '/{0}/ d' -i {1}".format(sparent,snolocaldups_path)
            y = commands.getstatusoutput(remove_sparent)

    fcnss = open(cns_file, 'a')
    
    for (qparent,sparent) in small_dups:
        ### remove qparent,sparent from file
        qfeat_dups = get_all_dups(qdup_dict,qparent)
        sfeat_dups = get_all_dups(sdup_dict,sparent)
        cnss_size = []
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

                qstart, qstop = max(qfeat['start'] - qpad, 1), qfeat['end'] + qpad
                sstart, sstop = max(sfeat['start'] - spad, 1), sfeat['end'] + spad

                assert qstop - qstart > 2 * qpad or qstart == 1, (qstop, qstart)
                assert sstop - sstart > 2 * spad or sstart == 1, (sstop, sstart)

                cmd = bl2seq % dict(qfasta=qfasta, sfasta=sfasta,
                        qstart=qstart,sstart=sstart, qstop=qstop, sstop=sstop,
                        e_value=30)
                
                return cmd, qfeat, sfeat
            cmds = [c for c in map(get_cmd, [l for l in pairs if l]) if c]
            results = (r for r in pool.map(commands.getoutput, [c[0] for c in cmds]))

            for res, (cmd, qfeat, sfeat) in zip(results, cmds):
                if not res.strip(): continue
                print >>sys.stderr,  "%s %s" % (qfeat["accn"], sfeat['accn'])
                orient = qfeat['strand'] == sfeat['strand'] and 1 or -1
                cnss = parse_blast(res, orient, qfeat, sfeat, qbed, sbed, qpad,spad)
                print >>sys.stderr, "(%i)" % len(cnss)
                cnss_fmt = ",".join(map(lambda l: ",".join(map(str,l)),cnss))
                cnss_size.append((len(cnss)*-1,qfeat["start"],sfeat["start"],qfeat["accn"],sfeat["accn"],cnss_fmt))
        if len(cnss_size) == 0 : continue
        cnss_size.sort()
        cns_number,qfeat_start, sfeat_start,qaccn,saccn,largest_cnss = cnss_size[0]
        qfeat = qbed.accn(qaccn)
        sfeat = sbed.accn(saccn)
        print >>sys.stderr, "FINALL{0},{1},{2}".format(qaccn,saccn,cns_number)
        if cns_number == 0 :
            update_nolocaldups(qbed,qbed.accn(qparent), sbed.accn(sparent), qnolocaldups_path, snolocaldups_path)
        else:
            localdup_file(qparent,sparent,qdups_path,sdups_path,cnss_size)
            fcnss.write("%s,%s,%s,%s,%s\n" % (qfeat['seqid'], qaccn,sfeat['seqid'], saccn,largest_cnss))
            update_pairs(qfeat["accn"],sfeat["accn"],qparent,sparent,pair_file)
            update_nolocaldups(qbed, qfeat, sfeat, qnolocaldups_path, snolocaldups_path)
    fcnss.close()
    sort_qdups = "sort -n -k 1 -k 2 {0} -o {0}".format(qnolocaldups_path)
    commands.getstatusoutput(sort_qdups)
    sort_sdups = "sort -n -k 1 -k 2 {0} -o {0}".format(snolocaldups_path)
    commands.getstatusoutput(sort_sdups)
    pairs_to_qa(pair_file, qnolocaldups_path, snolocaldups_path)

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
    parser.add_option("-q", dest="qfasta", help="path to genomic query fasta")
    parser.add_option("-s", dest="sfasta", help="path to genomic subjectfasta")
    parser.add_option("--qpad", dest="qpad", type='int', default=12000, help="how far from the end of the query gene to look for cnss")
    parser.add_option("--spad", dest="spad", type='int', default=12000, help="how far from the end of the subject gene to look for cnss")
    parser.add_option("--blast_path", dest="blast_path", type='string', help="path to bl2seq")
    parser.add_option("--qdups", dest="qdups", type='string', help="path to query localdup_file")
    parser.add_option("--sdups", dest="sdups", type='string', help="path to subject localdup_file")
    parser.add_option("--cns_file",dest="cns_file", type='string', help="path to cns file cns.txt")
    (options, _) = parser.parse_args()

    
    qbed = Bed(options.qbed, options.qfasta); qbed.fill_dict()
    sbed = Bed(options.sbed, options.sfasta); sbed.fill_dict()
    assert options.mask in 'FT'
    
    main(options.cns_file,options.qdups,options.sdups,options.pairs,options.pair_fmt,qbed,sbed,options.qpad,options.spad,options.blast_path,options.mask,options.ncpu)
