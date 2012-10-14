import sys
import os
import os.path as op
import numpy as np
import commands
from find_cns import get_masked_fastas,get_cmd,get_pair,get_feats_nearby,get_genespace,remove_crossing_cnss
from shapely.geometry import Point, Polygon, LineString, MultiLineString
from flatfeature import Bed
from pyfasta import Fasta
import re

from processing import Pool
pool = None


def parse_blast(blast_str, orient, qfeat, sfeat, qbed, sbed, qpad, spad, unmasked_fasta):
    blast = []
    slope = orient

    qgene = [qfeat['start'], qfeat['end']]
    sgene = [sfeat['start'], sfeat['end']]

    sgene = sgene[::slope]
    center = sum(qgene)/2., sum(sgene)/2.



    intercept = center[1] - slope * center[0]
    x = np.linspace(qgene[0] - qpad, qgene[1] + qpad, 50)
    y = slope * x + intercept

    feats_nearby = get_feats_nearby(qgene,sgene,qfeat,sfeat,x,y,qbed,sbed)    
    qgene_space_poly,qgene_poly,sgene_space_poly,sgene_poly = get_genespace(qfeat,sfeat,qgene,sgene)
    intronic_removed = 0
    
    cnss = set([])
    for line in blast_str.split("\n"):
        if "WARNING:" in line: continue
        if "ERROR" in line: continue
        if line == '': continue
        line = line.split("\t")
        if float(line[-1]) < 29.5: continue #finds 15/15 match
       # if float(line[-1]) < 33.4: continue #finds 17/17 match
        locs = map(int, line[6:10])
        locs.extend(map(float, line[10:]))

        xx = locs[:2]
        yy = locs[2:4]

        
        #######################################################
        # MAIZE BOWTIE : JUST 5 PRIME 3 PRIME
        #######################################################
        qcenter = sum(qgene)/2
        scenter = sum(sgene)/2 * orient
        qcns_center = sum(xx)/2
        scns_center = sum(yy)/2 * orient
        if  scns_center > scenter and qcns_center < qcenter: continue
        if qcns_center > qcns_center and scns_center < scenter : continue
        
        
        # to be saved. a hit must either be in an intron in both
        # genes, or in neither.

        ##########################################################
        # DEAL WITH INTRONIC cnss in the gene of interest.
        ##########################################################
        xls = LineString([(0, locs[0]), (0, locs[1])])
        yls = LineString([(0, locs[2]), (0, locs[3])])

        locs = tuple(locs) # make it hashable.
        if qgene_poly.intersects(xls) and sgene_poly.intersects(yls):
            cnss.update((locs,))
            continue
        # has to be both or neither.
        if qgene_space_poly.intersects(xls) or sgene_space_poly.intersects(yls):
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
        if intronic: continue

        ##########################################################
        cnss.update((locs,))

    # cant cross with < 2 cnss.
    # get rid of the eval, bitscore stuff.
    if len(cnss) < 2: return [(c[0], c[1], c[2], c[3],c[-1]) for c in cnss]
    cnss = list(cnss)
    ####################################################################################
    #########split cns into groups based on inversion, seq marks in maize ##########
    #################################################################################
    def group_cns(cnss, group):
      """input list of cns and list of groups , this puts the cns in a dictionary fmt key = group
      values = cns that fall within range of group"""
      for cns in cnss:
        if cns[2] in range(group[0],group[1]): # group start and end pos
          key = group
          cns_groups.setdefault(key, []).append(cns)

    cns_groups = {}
    inversion_groups = find_inversions(unmasked_fasta, sfeat, spad)
    [group_cns(cnss, group) for group in inversion_groups] # creates dict where key = group value is appended cns
    # for each goup of cns values run the followiung
    cns_by_group = []
    for key in cns_groups.keys():
      # # first group, groups into smaller groups on strand
      values = cns_groups[key]

      opp_strand = []
      same_strand = []
      for cns in values:
          if slope == 1 and cns[2] > cns[3]:
              opp_strand.append(cns)
          elif slope == -1 and cns[2] < cns[3]:
              opp_strand.append(cns)
          else:
              same_strand.append(cns)
      # need to flip to negative so the overlapping stuff still works.
      if orient == -1:
          same_strand = map(change_orient, same_strand)
          opp_strand = map(change_orient, opp_strand)
          sgene[0] *= -1
          sgene[1] *= -1
      if abs(sgene[1]) in range(key[0], key[1]): # if the cns fall in same group as gene we know its same stand  as gene and dont need to run rest
        cnss_same_strand = [(c[0], c[1], c[2], c[3],c[-1]) for c in remove_crossing_cnss(same_strand, qgene, sgene)]
        map(cns_by_group.append, cnss_same_strand)
      else:
        cnss_same_strand = [(c[0], c[1], c[2], c[3],c[-1]) for c in remove_crossing_cnss(same_strand, qgene, sgene)]
        cnss_opp_strand = cns_opp_strand(opp_strand, qgene, sgene) # alternitive for cns on opp strand
        if len(cnss_same_strand) < len(cnss_opp_strand):
          map(cns_by_group.append, cnss_opp_strand)
        else: # what about if they are the same, use non reverse complment
          map(cns_by_group.append, cnss_same_strand)
    if orient == -1:
        cns_by_group = [(c[0], c[1], -c[2], -c[3],c[-1]) for c in cns_by_group]

    return cns_by_group
    
##########################################################################################################
def find_inversions(unmasked_fasta, sfeat, spad):
  """finds all the invered/cut of seq regions in maize and creates a start stop tuple for each region
  input: unmasked fasta, sfeat and padding used for searching for NNNNNNs """
  chrm = sfeat['seqid']
  chr_seq = unmasked_fasta[chrm]
  sstart, sstop = max(sfeat['start'] - spad, 1), sfeat['end'] + spad
  search_seq = chr_seq[sstart:sstop]
  inverstion_positions = [(sstart,sstart)]
  for m in re.finditer(r"N{50,1000}", search_seq):
#    print '%d-%d: %s' % (m.start(), m.end(), m.group(0))
    inverstion_positions.append((m.start()+sstart, m.end()+sstart))
  inverstion_positions.append((sstop, sstop))
  create_regions = [(inverstion_positions[i][1], inverstion_positions[i+1][0]) for i in range(0, (len(inverstion_positions)-1))]
  return create_regions

def change_orient(cns):
    return (cns[0], cns[1], cns[2] * -1, cns[3] * -1, cns[4], cns[5])

def cns_opp_strand(cnss, qgene, sgene):
    cnss = list(cnss)
    cnss = map(change_orient,cnss)
    sgene[0] *= -1
    sgene[1] *= -1

    cnss = [(c[0], c[1], c[2], c[3],c[-2]) for c in remove_crossing_cnss(cnss, qgene, sgene)]
    cnss_fixed = [(c[0], c[1], -c[2], -c[3],c[-1]) for c in cnss]
    return cnss_fixed


def main(qbed, sbed, pairs_file, qpad, spad, unmasked_fasta, pair_fmt,blast_path, mask='F', ncpu=8):
    """main runner for finding cnss"""
    pool = Pool(ncpu)
    
    bl2seq = "%s " % blast_path + \
            "-p blastn -D 1 -E 2 -q -2 -r 1 -G 5 -W 7 -F %s " % mask + \
            " -e %(e_value).2f -i %(qfasta)s -j %(sfasta)s \
            -I %(qstart)d,%(qstop)d -J %(sstart)d,%(sstop)d | grep -v '#' \
            | grep -v 'WARNING' | grep -v 'ERROR' "


    fcnss = sys.stdout
    print >> fcnss,
    "#qseqid,qaccn,sseqid,saccn,[qstart,qend,sstart,send,bitscore...]"

    qfastas = get_masked_fastas(qbed)
    sfastas = get_masked_fastas(sbed) if qbed.filename != sbed.filename else qfastas

    pairs = [True]
    _get_pair_gen = get_pair(pairs_file, pair_fmt, qbed, sbed)
    # need this for parallization stuff.
    
    def get_pair_gen():
        try: return _get_pair_gen.next()
        except StopIteration: return None

    while any(pairs):
        pairs = [get_pair_gen() for i in range(ncpu)]
        # this helps in parallelizing.
	spad_map = [spad] * len(pairs)
        qpad_map = [qpad] * len(pairs)
        sfastas_map = [sfastas] * len(pairs)
        qfastas_map = [qfastas] * len(pairs)
        bl2seq_map =  [bl2seq] * len(pairs)
	####################################       
 
	cmds = [c for c in map(get_cmd, [l for l in pairs if
                l],bl2seq_map,qfastas_map,sfastas_map,qpad_map,spad_map) if c]
	results = (r for r in pool.map(commands.getoutput, [c[0] for c in cmds]))

        for res, (cmd, qfeat, sfeat) in zip(results, cmds):
            if not res.strip(): continue
            print >>sys.stderr,  "%s %s" % (qfeat["accn"], sfeat['accn']),
            orient = qfeat['strand'] == sfeat['strand'] and 1 or -1
            cnss = parse_blast(res, orient, qfeat, sfeat, qbed, sbed, qpad, spad, unmasked_fasta)
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
    parser.add_option("--qpad", dest="qpad", type='int', default=12000,
                      help="how far from the end of the query gene to look for cnss")
    parser.add_option("--spad", dest="spad", type='int', default=26000,
                    help="how far from the end of the subject gene to look for cnss")
    parser.add_option("--UMfasta", dest="unmasked_fasta", help="path to unmasked fasta file file")
    parser.add_option("--blast_path", dest="blast_path", type='string', help="path to bl2seq")
    (options, _) = parser.parse_args()


    if not (options.qfasta and options.sfasta and options.sbed and options.qbed):
        sys.exit(parser.print_help())

    qbed = Bed(options.qbed, options.qfasta); qbed.fill_dict()
    sbed = Bed(options.sbed, options.sfasta); sbed.fill_dict()
    unmasked_fasta = Fasta(options.unmasked_fasta)
    assert options.mask in 'FT'

    main(qbed, sbed, options.pairs, options.qpad, options.spad, unmasked_fasta, options.pair_fmt, options.blast_path, options.mask, options.ncpu)
