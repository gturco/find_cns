import sys
sys.path.append("scripts/post_processing/")
import commands
from processing import Pool
from shutil import copyfile
from os import path, getcwd
from itertools import product, chain
from collections import defaultdict

from flatfeature import Bed
from find_cns import parse_blast, get_masked_fastas,get_cmd
from qa_parsers import pairs_to_qa,ParsePairs,write_nolocaldups
from orderedset import OrderedSet
from cleanup import DupLine

def condense_rdups(parents,q):
    """condensed the rdups of the same group meaning if
    a:(a,b),(a,b),(A,B),(A,b) condense to (a,b),(A,b) of best b"""
    repeats_cond = defaultdict(list)
    ### for each group get key with largest sum
    #repeat_group = defaultdict(list)
    for repeat in parents.keys():
        repeat_group = defaultdict(list) 
	for cnss in parents[repeat]:
            #print >> sys.stderr, parents[repeat]
            cns_number,qfeat_start, sfeat_start,qaccn,saccn,largest_cnss = cnss
            #### sort cns by qaccn key
            isoform = qaccn if q else saccn
	    repeat_group[isoform].append(cnss)
        for isoforms in repeat_group.keys():
            repeat_group[isoforms].sort()
            #print >> sys.stderr,"key {1} isoooooo {0}".format(repeat_group[isoforms],isoforms)
            isoparents = (repeat_group[isoforms][0], repeat)
	    repeats_cond[isoforms].append(isoparents)
    #print>> sys.stderr, repeats_cond
    return repeats_cond

def same_alpha(rdup,parents):
	qparent = []
	for parent in parents:
		qparent.append((rdup == parent[0]))
        if len(set(qparent)) == 1: return True
 	else: return  False

def best_repeats(rdups):
    """finds best dup when dups are repeated"""
    best_rep = defaultdict(list)
    for rdup in rdups.keys():
        #print rdup
        parents = rdups[rdup]
        #print parents
        if not same_alpha(rdup,parents) : continue 
        q=(rdup == parents.keys()[0][0])
        cond_reps =condense_rdups(parents,q)
        isoform_tuple = []
        for isoform in cond_reps.keys():
	    iso=  cond_reps[isoform]
            cns_total = [cns_number for (cns_number,qfeat_start,
            sfeat_start,qaccn,saccn,largest_cnss), parent in iso]
            isoform_tuple.append((sum(cns_total),iso))
        isoform_tuple.sort()
	#print isoform_tuple
        for i,p in isoform_tuple[0][1]:
        	best_rep[p] = (i)
    		#print "PPPPPPPPPP", p
    #print best_rep
    return best_rep


def get_large_dups(pairs,qdup_dict,sdup_dict):
    large_dups = []
    for qparent,sparent in pairs:
        qfeat_dups = get_all_dups(qdup_dict,qparent)
        sfeat_dups = get_all_dups(sdup_dict,sparent)
        if len(list(product(qfeat_dups,sfeat_dups))) > 9:
            large_dups.append((qparent,sparent))
    return large_dups

def skip_pair(qparent,sparent,rdups,rdups_both,ldups):
    """return ture if should be skiped because its a large dup seen more then
    once, too large of a combo, or both the query and subject are repeated more
    then once"""
    ldups_list = list(chain(*ldups))
    rdups_both = list(chain(*rdups_both))
    if qparent in rdups and qparent in ldups_list:return True
    elif sparent in rdups and sparent in ldups_list: return True
    elif qparent in rdups_both: return True
    elif sparent in rdups_both: return True
    elif (qparent,sparent) in ldups: return True
    else: return False


def parse_dups(localdup_path):
    localdup_file = open(localdup_path)
    dup_dic ={}
    for line in localdup_file:
        dup = DupLine(line)
        dup_dic[dup.parent] = dup.dups
    localdup_file.close()
    return dup_dic

def get_dups(qfeat_dups,sfeat_dups,qbed,sbed):
    """creates all poss combinations of qaccn,saccn if less then 9"""
    for qaccn,saccn in product(qfeat_dups,sfeat_dups):
        qfeat = qbed.accn(qaccn)
        sfeat = sbed.accn(saccn)
        if len(list(product(qfeat_dups,sfeat_dups))) > 9: yield None
        else: yield qfeat,sfeat

def get_pairs(pair_file,fmt,qdup_dict,sdup_dict):
    """parses pair file only returns pairs with atlest one localdup"""
    dups = []
    seen = []
    rdups =[]
    for line in open(pair_file):
        if line[0] == "#" : continue
        pair_line = ParsePairs(line)
        pair = getattr(pair_line,fmt)()
        pair = tuple(pair)
        if pair[0] in seen: rdups.append(pair[0])
        if pair[1] in seen: rdups.append(pair[1])
        if pair[0] in qdup_dict.keys() or pair[1] in sdup_dict.keys():
            dups.append((pair[0],pair[1]))
            if pair[0] in qdup_dict.keys(): seen.append(pair[0])
            if pair[1] in sdup_dict.keys(): seen.append(pair[1])
        else: continue
    return dups,rdups

def get_all_dups(dup_dict,feat):
    """ if local dup return all dup otherwise just return name of feat"""
    try:
        dups = dup_dict[feat]
    except KeyError:
        dups = [feat]
    return dups
############### update files ####################
def write_new_dups(npair_file,ncns_file,nqlocaldups,nslocaldups,cnss_size,qparent,sparent,qfeat,sfeat,qdups,sdups):
    """ reseach and replace cns file, localdups and pairs file"""
    cns_number,qfeat_start, sfeat_start,qaccn,saccn,largest_cnss = cnss_size[0]
    update_pairs(qfeat['accn'],sfeat['accn'],qparent,sparent,npair_file)
    if abs(cns_number) > 0:
        update_cnss_line(qfeat,sfeat,qparent,sparent,largest_cnss,ncns_file)
    qdup = list(qdups[qparent])  if qparent in qdups.keys() else [qparent]
    sdup = list(sdups[sparent])  if sparent in sdups.keys() else [sparent]
    #print "sssss", qdup,qparent,sparent
    write_localdup_file(qparent,sparent,nqlocaldups,nslocaldups,qaccn,saccn,qdup,sdup)

def make_copy_of_file(file1):
    """makes a copy of the file renaming it .local returning the new file
    path"""
    file2 = "{0}/{1}.local".format(path.dirname(file1),path.basename(file1))
    copyfile(file1,file2)
    print >>sys.stderr, "write tandem-filtered bed file {0}".format(file2)
    return file2

def update_pairs(qaccn,saccn,qparent,sparent,pair_file):
    """replaces the old pair with the new head loaldup pairs"""
    #### new pairs file search and replace
    search_replace_pairs = "sed 's/{0}\t{1}/{2}\t{3}/' -i {4}".format(qparent,sparent,qaccn,saccn,pair_file)
    out = commands.getstatusoutput(search_replace_pairs)
    assert out[0] == 0 

def update_cnss_line(qfeat,sfeat,qparent,sparent,largest_cnss,ncns_file):
    """ removes any cnss with the old parent dup and adds new one """
    rm_cns = '^{0},{1},{2},{3}.*'.format(qfeat['seqid'],qparent,sfeat['seqid'],sparent)
    removed = commands.getstatusoutput("sed '/{0}/ d' -i {1}".format(rm_cns,ncns_file))
    cns_line = '{0},{1},{2},{3},{4}'.format(qfeat['seqid'],qfeat['accn'],sfeat['seqid'],sfeat['accn'],largest_cnss)
    add_cns = commands.getstatusoutput("cat {0} | sed '$a {1}' -i {0}".format(ncns_file,cns_line))
    assert add_cns[0] == 0

def write_localdup_file(qparent,sparent,qfile,sfile,nqparent,nsparent,qdup,sdup):
    """ replaces the orginal parent local dup with the new order"""
    #print nqparent, qdup
    qdup.remove(nqparent)
    qnew = [nqparent] + qdup
    #print sdup,nsparent
    sdup.remove(nsparent) 
    snew = [nsparent] + sdup

    #### orderset: no repeats in order
    qreplace = "\t".join(qnew)
    sreplace = "\t".join(snew)
    qsearch_replace = "sed -i 's/^{0}.*/{1}/g' {2}".format(qparent,qreplace,qfile)
    
    ##### if starts somewhere else then appends it to the end of the line
    ssearch_replace = "sed -i 's/^{0}.*/{1}/g' {2}".format(sparent,sreplace,sfile)
    commands.getstatusoutput(qsearch_replace)
    assert commands.getstatusoutput(ssearch_replace)[0] == 0

###############################################################

def main(cns_file,qdups_path,sdups_path,pair_file,fmt,qbed,sbed,qpad,spad,blast_path,mask='F',ncpu=8):
    pool = Pool(ncpu)
    bl2seq = "%s " % blast_path + \
            "-p blastn -D 1 -E 2 -q -2 -r 1 -G 5 -W 7 -F %s " % mask + \
            " -e %(e_value).2f -i %(qfasta)s -j %(sfasta)s \
             -I %(qstart)d,%(qstop)d -J %(sstart)d,%(sstop)d | grep -v '#' \
             | grep -v 'WARNING' | grep -v 'ERROR' "

    qfastas = get_masked_fastas(qbed)
    sfastas = get_masked_fastas(sbed) if qbed.filename != sbed.filename else qfastas
    
    ################# file paths #####################
    qnolocaldups_path =  qbed.path.split(".")[0] + ".nolocaldups.bed"
    snolocaldups_path = sbed.path.split(".")[0] + ".nolocaldups.bed"
    qlocaldups_path = qbed.path.split(".")[0] + ".localdups"
    slocaldups_path = sbed.path.split(".")[0] + ".localdups"
    npair_file,nqlocaldups,nslocaldups, ncns_file = map(make_copy_of_file,[pair_file,qlocaldups_path,slocaldups_path,cns_file])
    ##########################################
    
    qdups = parse_dups(qdups_path)
    sdups = parse_dups(sdups_path)
    dups,rdups = get_pairs(pair_file,fmt,qdups,sdups)
    print len(dups), len(rdups)
    ldups = get_large_dups(dups,qdups,sdups)

    rdups_dic = defaultdict(dict)
    rdups_both = [(qparent,sparent) for qparent,sparent in dups if qparent in rdups and sparent in rdups]
    for (qparent,sparent) in dups:
        if skip_pair(qparent,sparent,rdups,rdups_both,ldups):continue
        cnss_size = []
        qfeat_dups = get_all_dups(qdups,qparent)
        sfeat_dups = get_all_dups(sdups,sparent)
        pairs = [True]
        _get_dups_gen = get_dups(qfeat_dups,sfeat_dups,qbed,sbed)

        def get_dups_gen():
            try: return _get_dups_gen.next()
            except StopIteration: return None
        while any(pairs):
            cnss_dups = []
            pairs = [get_dups_gen() for i in range(ncpu)]
            ###this is for parellization#########
            spad_map = [spad] * len(pairs)
            qpad_map = [qpad] * len(pairs)
            sfastas_map = [sfastas] * len(pairs)
            qfastas_map = [qfastas] * len(pairs)
            bl2seq_map =  [bl2seq] * len(pairs)
            ###################################
            cmds = [c for c in map(get_cmd, [l for l in pairs if l],
                bl2seq_map,qfastas_map,sfastas_map,qpad_map,spad_map) if c]
            results = (r for r in pool.map(commands.getoutput, [c[0] for c in cmds]))
            for res, (cmd, qfeat, sfeat) in zip(results, cmds):
                orient = qfeat['strand'] == sfeat['strand'] and 1 or -1 
                if not res.strip(): cnss = []
                else: cnss = parse_blast(res, orient, qfeat, sfeat, qbed, sbed, qpad,spad)
                print >>sys.stderr, "(%i)" % len(cnss)
                cnss_fmt = ",".join(map(lambda l: ",".join(map(str,l)),cnss))
                cnss_size.append((len(cnss)*-1,qfeat["start"],sfeat["start"],qfeat["accn"],sfeat["accn"],cnss_fmt))
            pairs = [pairs[-1]]
        ######################################################################
        if qparent in rdups:
            if (qparent,sparent) in rdups_dic[qparent].keys(): logging.info((qparent,sparent))
            rdups_dic[qparent].update({(qparent,sparent):cnss_size})
        elif sparent in rdups:
            if (qparent,sparent) in rdups_dic[sparent].keys(): logging.info((qparent,sparent))
            rdups_dic[sparent].update({(qparent,sparent):cnss_size})
        else:
            cnss_size.sort()
            cns_number,qfeat_start,sfeat_start,qaccn,saccn,largest_cnss = cnss_size[0]
            qfeat = qbed.accn(qaccn)
            sfeat = sbed.accn(saccn)
            print >>sys.stderr, "FINAL: {0},{1},{2}".format(qaccn,saccn,cns_number)
            write_new_dups(npair_file,ncns_file,nqlocaldups,nslocaldups,cnss_size,qparent,sparent,qfeat,sfeat,qdups,sdups)
    
    best_reps = best_repeats(rdups_dic)
    for dparents in best_reps.keys():
	#print dparents
        qparent,sparent = dparents
        #print parents,best_reps[parents]
        ### one or list? cnss[0]?
        cns_number,qfeat_start, sfeat_start,qaccn,saccn,largest_cnss = best_reps[dparents]
        write_new_dups(npair_file,ncns_file,nqlocaldups,nslocaldups,[best_reps[dparents]],qparent,sparent,qfeat,sfeat,qdups,sdups)

    write_nolocaldups(qbed.path,nqlocaldups,"{0}.nolocaldups.bed.local".format(qbed.path.split(".")[0]))
    write_nolocaldups(sbed.path,nslocaldups,"{0}.nolocaldups.bed.local".format(sbed.path.split(".")[0]))
    pairs_to_qa(npair_file,'pair',"{0}.nolocaldups.bed.local".format(qbed.path.split(".")[0]),"{0}.nolocaldups.bed.local".format(sbed.path.split(".")[0]),"{0}.raw.filtered.local".format(options.pairs.split(".")[0]))

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
    

    qnolocaldups_path =  qbed.path.split(".")[0] + ".nolocaldups.bed"
    snolocaldups_path =  sbed.path.split(".")[0] + ".nolocaldups.bed"
    #pairs_to_qa("{0}.local".format(options.pairs),'pair',"{0}.nolocaldups.bed.local".format(qbed.path.split(".")[0]),"{0}.nolocaldups.bed.local".format(sbed.path.split(".")[0]),"{0}.raw.filtered.local".format(options.pairs.split(".")[0]))

    import logging
    LOG_FILENAME = path.dirname(options.qfasta) + "dup_rdups.log"
    logging.basicConfig(filename=LOG_FILENAME,level=logging.INFO)



    main(options.cns_file,options.qdups,options.sdups,options.pairs,options.pair_fmt,qbed,sbed,options.qpad,options.spad,options.blast_path,options.mask,options.ncpu)
