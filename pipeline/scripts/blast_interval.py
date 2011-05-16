import re
#import rpy2.robjects as robjects
from numpy import array, median, mean
from flatfeature import Bed
import sys
import os
import os.path as op
import commands
from bx.intervals.intersection import Interval, IntervalTree


def find_median_interval(file_path):
  fh = open(file_path)
  interval_diff = fh.read() 
  interval_diff_list = interval_diff.split('\n')
  interval_ints = map(int, interval_diff_list)
  interval_pos = map(abs, interval_ints)
  a = array(interval_pos)
  print mean(a), median(a), max(interval_pos), min(interval_pos)
  #print robjects.r.assign('x_from_python', a)
  # # =print x
  #print robjects.r.median(a)
# find_median_interval('/Users/gturco/sorg_rice_maize_interval.txt')
# 
# 
# get pairs:
# qfeat = cns start, cns_stop + padding, chr, strand?
# sfeat = maize homelog + padding
# 
# mask_fasta_files
# 
# blast_cmd
# bed.cds_fasta()

def parse_blast(blast_str, cns_dict):
  """takes a blast file and cns_pair and sees if the query cns intersects with 
  any of the cns found"""
  scns_inteval = IntervalTree()
  for line in blast_str.split("\n"):
    if "WARNING" in line: continue
    if "ERROR" in line: continue
    line = line.split("\t")
    locus = map(int, line[6:10])
    locus.extend(map(float, line[10:]))
    
    s_start, s_end = locus[:2]
    s_start = min(int(s_start), int(s_end))
    s_end = max(int(s_start), int(s_end))
    scns_inteval.insert_interval(Interval(s_start, s_end))

  q_start = min(int(cns_dict['start']), int(cns_dict['end']))
  q_end = max(int(cns_dict['start']), int(cns_dict['end'])) 
  intersecting_cns = scns_inteval.find(q_start, q_end)
  return intersecting_cns
  
def get_fastas(bed, masked = True):
    """"
    if mask is Ture it masked all the cds and prints out a a masked fasta for
    each chr otherwise if flase it prints out a unmasked fasta for each chr
    """  
    f = bed.fasta.fasta_name
    fname = op.splitext(op.basename(f))[0]
    print fname
    d = op.dirname(f) + "/%s_split" % fname
    try: os.mkdir(d)
    except OSError: pass
    
    fastas = {}
    if masked == False:
      for seqid in bed.fasta.keys():
        f = d + "/%s.fasta" % seqid
        fastas[seqid] = f
        if op.exists(f): continue#os.remove(f)
        fh = open(f, "wb")
        print >>fh, bed.fasta[seqid]
        fh.close()
      return fastas
    elif masked == True:
      for seqid, seq in bed.mask_cds():
        f = d + "/%s.fasta" % seqid
        fastas[seqid] = f
        if op.exists(f): continue
        fh = open(f, "wb")
        print >>fh, seq
        fh.close()
      return fastas

def main(cns_bed, ortho_bed, cns_dict, qpad, spad, blast_path, mask):
  "imput a cns_dict and otholog_dict, cns_bed "
  
  bl2seq = "%s " % blast_path + \
          "-p blastn -D 1 -E 2 -q -2 -r 1 -G 5 -W 7 -F %s " % mask + \
          " -e %(e_value).2f -i %(qfasta)s -j %(sfasta)s \
             -I %(qstart)d,%(qstop)d -J %(sstart)d,%(sstop)d | grep -v '#' \
           | grep -v 'WARNING' | grep -v 'ERROR' "

  qfastas = get_fastas(cns_bed, True)
  sfastas = get_fastas(ortho_bed, False)

  def get_cmd(sfeat, qfeat, qpad, spad, mask):
        #if qfeat['accn'] != "Bradi4g01820": return None
        #print >>sys.stderr, qfeat, sfeat
        
        sfasta = sfastas[sfeat['seqid']]
        #qfasta = '/Users/gturco/trobble_shot/rice_v6_split/4.fasta'
        qfasta = qfastas[qfeat['seqid']]
        

        qstart, qstop = max(qfeat['start'] - qpad, 1), qfeat['end'] + qpad
        sstart, sstop = max(sfeat['start'] - spad, 1), sfeat['end'] + spad

        m = qstop - qstart
        n = sstop - sstart
      
        e_value = m*n*(2**(-28.51974)) # bit score above 15/15 noise
        assert e_value > 0

        cmd = bl2seq % dict(qfasta=qfasta, sfasta=sfasta, qstart=qstart,
                            sstart=sstart, qstop=qstop, sstop=sstop, e_value=e_value)
        return cmd
        
  sfeat = ortho_bed.accn('GRMZM2G094328')
  x = get_cmd(sfeat,cns_dict, qpad, spad, 'T')
  results = commands.getoutput(x)
  interval = parse_blast(results, cns_dict) 
  if len(interval) >= 1:
     print cns_dict['start'], cns_dict['end'], interval
  


if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser("usage: %prog [options] ")
    parser.add_option("-F", dest="mask", help="blast mask simple sequence [default: F]", default="F")
    parser.add_option("--ac", dest="qfasta", help="path to genomic query fasta")
    parser.add_option("--qbed", dest="qbed", help="query bed file")
    parser.add_option("--bd", dest="sfasta", help="path to genomic subject fasta")
    parser.add_option("--sbed", dest="sbed", help="subject bed file")
    parser.add_option("--qpad", dest="qpad", type='int', default=12000, help="how far from the end of the query gene to look for cnss")
    parser.add_option("--spad", dest="spad", type='int', default=26000, help="how far from the end of the subject gene to look for cnss")
    parser.add_option("--blast_path", dest="blast_path", type='string', help="path to bl2seq")
    (options, _) = parser.parse_args()
    
    # if not (options.qfasta and options.sfasta and options.sbed and options.qbed):
    #     print 'hi'
    #     sys.exit(parser.print_help())
        
    print options.qbed, options.qfasta, options.blast_path, options.mask
    cns_dict = {'start':31210231, 'end':31210254, 'seqid':'4'}
    qbed = Bed(options.qbed, options.qfasta); qbed.fill_dict()
    sbed = Bed(options.sbed, options.sfasta); sbed.fill_dict()
    assert options.mask in 'FT'

    main(qbed, sbed, cns_dict, options.qpad, options.spad, options.blast_path, options.mask)
