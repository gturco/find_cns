import sys
import numpy
from shapely.geometry import Point, Polygon, LineString, MultiLineString
from flatfeature import Bed
import csv
from collections import Counter
import datetime
from pandas import read_csv
from pyfasta import Fasta
from Bio.SeqUtils import Seq

def get_fasta(qfasta,sfasta,cns):
    "return the cns seqence for query and subject"
    qchr = qfasta[cns['qseqid']]
    schr = sfasta[cns['sseqid']]
    qstart,qend,sstart,send = map(int,[cns['qstart'],cns['qstop'],cns['sstart'],cns['sstop']])
    qseq = qchr[qstart,qend]
    sseq = schr[sstart,send]
    if cns['qstrand'] == '-':
        qseq = str(Seq(seq).reverse_complement())
    if cns['strand'] == '-':
        sseq = str(Seq(seq).reverse_complement())
    return qseq, sseq

def cns_link(qaccn,saccn, qdsid, sdsid, qpad,spad, base="http://synteny.cnr.berkeley.edu/CoGe/GEvo.pl?prog=blastn&autogo=1&show_cns=1&"):
  
    url = "dsgid1={0}&accn1={1}&dr1up={2}&dr1down={2}&\
dsgid2={3};accn2={4};dr2up={5};dr2down={5};num_seqs=2;hsp_overlap_limit=0;hsp_size_limit=0".format(qdsid,qaccn,qpad,sdsid,saccn,spad)
    return base + url

def group_locations(cns_dict):
    grouped_locations = []
    for cns in cns_dict:
        key = "{0}__{1}_{2}".format(cns["qaccn"],cns["saccn"],cns["type"])
        grouped_locations.append(key)
    grouped_locations_dic = Counter(grouped_locations)
    return grouped_locations_dic

def group_cns_number(cns_dict):
    grouped_cns = []
    for cns in cns_dict:
        key = "{0}__{1}".format(cns["qaccn"],cns["saccn"])
        grouped_cns.append(key)
    grouped_cns_number = Counter(grouped_cns)
    return grouped_cns_number

def cns_to_dic(cns,fmt):
    """ takes cns file in sql or csv format and creates a pck"""
    cns_file = list(csv.DictReader(open(cns)))
    return cns_file

def write_to_file_grouped(cns_grouped_number,grouped_locations_dic,cns_path,qdsid,sdsid):
    tdate = str(datetime.date.today())
    org1_org2 = cns_path.split("/")[1]
    write_file = open("data/{1}/{1}.cnslist-{0}.csv".format(tdate,org1_org2),"wb")
    header = "qaccn,saccn,number_of_cns,5_distal,5_proximal,5_UTR,intron,3_UTR,3_proximal,3_distal,url\n"
    write_file.write(header)
 
    for cns in cns_grouped_number:
        qaccn,saccn = cns.split('__')
        url = cns_link(qaccn,saccn, qdsid, sdsid,15000,15000)
	new_line = "{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10}\n".format(qaccn,saccn,cns_grouped_number[cns],
                grouped_locations_dic["{0}_5-distal".format(cns)],
        grouped_locations_dic["{0}_5-proximal".format(cns)],
        grouped_locations_dic["{0}_5-UTR".format(cns)],
        grouped_locations_dic["{0}_intron".format(cns)],
        grouped_locations_dic["{0}_3-UTR".format(cns)],
        grouped_locations_dic["{0}_3-proximal".format(cns)],
        grouped_locations_dic["{0}_3-distal".format(cns)],url)
        write_file.write(new_line)
    write_file.close()

def write_one_to_file(qfasta,sfasta,cns_dict,fmt,out_fh):
    """ imports into sql if fmt is pck otherwise writes file cns.location"""
    write_file = open(out_fh,"wb")
    if fmt == "csv":
        for cns in cns_dict:
            #string_values = map(str,cns.values())
            #all_values = "\t".join(string_values)
            qseq,sseq = get_fasta(qfasta,sfasta,cns)
            new_line = "{0},{1},{2}\n".format(cns["#cns_id"],cns["type"],qseq,sseq)
            write_file.write(new_line)
        write_file.close()
    pos_fasta = read_csv(out_fh,index_col=0)
    cnsf = read_csv(cns_file,index_col=0)
    cnslist = cns.join(pos_fasta)
    cnslist_fh ==
    cnslist.to_csv(cnslist_fh,sep=',')


def main(cns_path, fmt, query_bed_path, subject_bed_path):
  cns_dic = cns_to_dic(cns_path,fmt)
  query_bed = Bed(query_bed_path)
  subject_bed = Bed(subject_bed_path)
  utr_dict = {}
  for cns in cns_dic:
    cns['qstop'] = int(cns['qstop'])
    cns['qstart'] = int(cns['qstart'])
    cns['sstop'] = int(cns['sstop'])
    cns['sstart'] = int(cns['sstart'])
 
    qfeat = query_bed.accn(cns['qaccn'])
    sfeat = subject_bed.accn(cns['saccn']) 
    qgene_space_start = min(qfeat['locs'])[0]
    qgene_space_end = max(qfeat['locs'])[1]
    qgene_space_poly = LineString([(0.0, qgene_space_start), (0.0, qgene_space_end)])
    qgene_poly = LineString([(0.0, qfeat['start']), (0.0, qfeat['end'])])
    sgene_poly = LineString([(0.0, sfeat['start']), (0.0, sfeat['end'])])
    # if intron of one dont need to check other
    qcns = LineString([(0,cns['qstart']),(0,cns['qstop'])])
    scns = LineString([(0,cns['sstart']),(0,cns['sstop'])])
    cns_type(cns,qgene_space_poly, qgene_poly, sgene_poly, scns, qcns,qgene_space_start,qfeat)
    create_utr_list(utr_dict,qfeat, cns,"q")
    create_utr_list(utr_dict,sfeat, cns,"s")
  for cns in cns_dic:
    if cns['type'] == "5-prox_dist":
      qgene_start = min(utr_dict[cns['qaccn']])
      qgene_stop =  max(utr_dict[cns['qaccn']])
      # sstart = min(utr_dict[cns['saccn']])
      # sstop =  max(utr_dict[cns['saccn']])
      five_diff_pos = abs(qgene_start - cns["qstop"])
      five_diff_neg = abs(qgene_stop - cns["qstart"])
      if five_diff_pos <=1000 and cns["qstrand"] == "+" or five_diff_neg <=1000 and cns["qstrand"] == "-":
        cns["type"] = "5-proximal"
      elif five_diff_pos >1000 and cns["qstrand"] == "+" or five_diff_neg >1000 and cns["qstrand"] == "-":
        cns["type"] = "5-distal"
    elif cns['type'] == "3-prox_dist":
      qgene_start = min(utr_dict[cns['qaccn']])
      qgene_stop =  max(utr_dict[cns['qaccn']])
      three_diff_pos =  abs(cns["qstart"] - qgene_stop)
      three_diff_neg =  abs(cns["qstop"] - qgene_start)
      if three_diff_pos <=1000 and cns["qstrand"] == "+" or three_diff_neg <=1000 and cns["qstrand"] == "-":
        cns["type"] = "3-proximal"
      elif three_diff_pos > 1000 and cns["qstrand"] == "+" or three_diff_neg > 1000 and cns["qstrand"] == "-":
        cns["type"] = "3-distal"
  return cns_dic

def group_cns(cns_dic,cns_path,qid,sid):
    grouped_locations_dic = group_locations(cns_dic)
    cns_grouped_number = group_cns_number(cns_dic)
    write_to_file_grouped(cns_grouped_number,grouped_locations_dic,cns_path,qid,sid)

#write_to_file(cns,fmt)
    
def cns_type(cns, qgene_space_poly, qgene_poly, sgene_poly, scns, qcns,qgene_space_start,qfeat):
  if qgene_space_poly.intersects(qcns):
    cns["type"] = "intron"
  elif qgene_poly.intersects(qcns) or sgene_poly.intersects(scns):
    if qfeat['strand'] == "+" and cns['qstart'] < qgene_space_start or qfeat['strand'] == "-" and cns['qstart'] > qgene_space_start:
      cns['type'] = '5-UTR'
    elif qfeat['strand'] == "+" and cns['qstart'] > qgene_space_start or qfeat['strand'] == "-" and cns['qstart'] < qgene_space_start:
      cns['type'] = '3-UTR'
  elif qfeat['strand'] == "+" and cns['qstart'] > qgene_space_start or qfeat['strand'] == "-" and cns['qstart'] < qgene_space_start: 
    cns["type"] = "3-prox_dist"
  elif qfeat['strand'] == "+" and cns['qstart'] < qgene_space_start or qfeat['strand'] == "-" and cns['qstart'] > qgene_space_start: 
    cns["type"] = "5-prox_dist"
  return cns
        
def create_utr_list(utr_dict,feat, cns, letter):
  if feat['accn'] not in utr_dict.keys():
    utr_dict[feat["accn"]] = [feat['start'], feat["end"]]
  if cns['type'] == "UTR":
    utr_dict[feat['accn']].append((cns['{0}start'.format(letter)]))
    utr_dict[feat['accn']].append((cns['{0}end'.format(letter)]))
  
"postional 5' 3 'distance from utr - cns...."
"5' or 3' if < then gene_poly = 5 prime.... if grter then gene space poly then 3 prinme.."
  #   create intron list

    # qgene_poly = LineString([(0,qgene['start']),(0,cns['end'])])
    # sgene_poly = LineString([(0,sgene['start']),(0,cns['end'])])
    # 
if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser()
    parser.add_option("--qbed", dest="qbed", help="bed file of the query")
    parser.add_option("--sbed", dest="sbed", help="bed file of the subject")
    parser.add_option("--cns", dest="cns", help="path to the cns file created by find_cns.py")
    parser.add_option("--fmt", dest="fmt", help="fmt of file csv file or sql file")
    parser.add_option("--qdsgid",dest="qdsgid",help="dataset group id form coge for query org")
    parser.add_option("--sdsgid",dest="sdsgid",help="dataset group id from coge database for subject org")
    parser.add_option("--qfasta",dest="qfasta_path",help="fasta file for query cns ")
    parser.add_option("--sfasta",dest="sfasta_path",help="fasta file for subject cns")


    (options, _) = parser.parse_args()

    x= main(options.cns,options.fmt,options.qbed,options.sbed)
    qfasta = Fasta(options.qfasta_path)
    sfasta = Fasta(options.sfasta_path)
    write_one_to_file(qfasta, sfasta, x,"csv","{0}.location_indvi".format(options.cns))
    group_cns(x,options.cns,options.qdsgid,options.sdsgid)

  # x= main("/Users/gturco/code/freeling_lab/find_cns_gturco/pipeline/scripts/post_processing/find_cns_cns_test.pck","/Users/gturco/code/freeling_lab/find_cns_gturco/pipeline/scripts/post_processing/query_test.bed","/Users/gturco/code/freeling_lab/find_cns_gturco/pipeline/scripts/post_processing/subject_test.bed")
  ##### REMINDER CHANGE TABLE INSERTS INTO!! ##############################################
 # x= main("/Users/gturco/data/paper3/paper3_rice_b_sorghum_v1_gturco_2011_4_11app_real.pck","/Users/gturco/data/paper3/rice_b_sorg.nolocaldups.with_new.all.bed","/Users/gturco/data/paper3/sorghum_v1.nolocaldups.with_new.all.bed")
  #x=main("/Users/gturco/find_cns_thaliana_v10_thaliana_v10_gturco_2011_22_11app.pck","/Users/gturco/tair_10.nolocaldups.with_new.all.bed","/Users/gturco/tair_10.nolocaldups.with_new.all.bed")
  #load_in_table("rice_b_sorghum_v1_gturco_2011_4_11app_real_grouped")
  #x = main("../../data/rice_b_setaria64/rice_b_setaria64.cns.assigned_real.csv","csv","../../data/rice_b_setaria64/rice_b.nolocaldups.with_new.all.bed","../../data/rice_b_setaria64/setaria64.nolocaldups.with_new.all.bed")
  #group_cns(x,"../../data/rice_b_setaria64/rice_b_setaria64.cns.assigned_real.csv")
  #utr_present("/Users/gturco/data/find_cns_3_UTR.pck","/Users/gturco/data/rice_v6.nolocaldups.with_new.all.bed", 3)
  #utr_present("/Users/gturco/data/find_cns_5_UTR.pck","/Users/gturco/data/rice_v6.nolocaldups.with_new.all.bed", 5)
  #need to remove all cns rna ones.....
