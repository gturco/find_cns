import numpy
from shapely.geometry import Point, Polygon, LineString, MultiLineString
from flatfeature import Bed
import pickle
import MySQLdb

def main(cns_pck, query_bed_path, subject_bed_path):
  cns_handle = open(cns_pck)
  cns_pickle = pickle.load(cns_handle)
  query_bed = Bed(query_bed_path)
  subject_bed = Bed(subject_bed_path)
  utr_dict = {}
  for cns in cns_pickle:
    qfeat = query_bed.accn(cns['qaccn'])
    sfeat = subject_bed.accn(cns['saccn']) 
    qgene_space_start = min(qfeat['locs'])[0]
    qgene_space_end = max(qfeat['locs'])[1]
    qgene_space_poly = LineString([(0.0, qgene_space_start), (0.0, qgene_space_end)])
    qgene_poly = LineString([(0.0, qfeat['start']), (0.0, qfeat['end'])])
    sgene_poly = LineString([(0.0, sfeat['start']), (0.0, sfeat['end'])])
    # if intron of one dont need to check other
    qcns = LineString([(0,cns['qcns_start']),(0,cns['qcns_end'])])
    scns = LineString([(0,cns['scns_start']),(0,cns['scns_end'])])
    cns_type(cns,qgene_space_poly, qgene_poly, sgene_poly, scns, qcns,qgene_space_start,qfeat)
    create_utr_list(utr_dict,qfeat, cns,"q")
    create_utr_list(utr_dict,sfeat, cns,"s")
  for cns in cns_pickle:
    if cns['type'] == "5-prox_dist":
      qgene_start = min(utr_dict[cns['qaccn']])
      qgene_stop =  max(utr_dict[cns['qaccn']])
      # sstart = min(utr_dict[cns['saccn']])
      # sstop =  max(utr_dict[cns['saccn']])
      five_diff_pos = abs(qgene_start - cns["qcns_end"])
      five_diff_neg = abs(qgene_stop - cns["qcns_start"])
      if five_diff_pos <=1000 and cns["qstrand"] == "+" or five_diff_neg <=1000 and cns["qstrand"] == "-":
        cns["type"] = "5-proximal"
      elif five_diff_pos >1000 and cns["qstrand"] == "+" or five_diff_neg >1000 and cns["qstrand"] == "-":
        cns["type"] = "5-distal"
    elif cns['type'] == "3-prox_dist":
      qgene_start = min(utr_dict[cns['qaccn']])
      qgene_stop =  max(utr_dict[cns['qaccn']])
      three_diff_pos =  abs(cns["qcns_start"] - qgene_stop)
      three_diff_neg =  abs(cns["qcns_end"] - qgene_start)
      if three_diff_pos <=1000 and cns["qstrand"] == "+" or three_diff_neg <=1000 and cns["qstrand"] == "-":
        cns["type"] = "3-proximal"
      elif three_diff_pos > 1000 and cns["qstrand"] == "+" or three_diff_neg > 1000 and cns["qstrand"] == "-":
        cns["type"] = "3-distal"
    import_into_mysql(cns)

def import_into_mysql(cns):
  db = MySQLdb.connect(host="127.0.0.1", user="root", db = "find_cns")
  cursor = db.cursor()
  stmt = "INSERT INTO cns_postion_info (cns_id, qaccn, qseqid,qcns_start,qcns_end,qstrand, saccn, sseqid, scns_start, scns_end, sstrand, urls, postion) values('{0}','{1}',{2},{3},{4},'{5}','{6}',{7},{8},{9},'{10}','{11}','{12}')".format(cns['cns_id'], cns['qaccn'], cns["qseqid"],cns["qcns_start"],cns["qcns_end"],cns["qstrand"], cns["saccn"], cns["sseqid"], cns["scns_start"], cns["scns_end"], cns["sstrand"], cns["urls"], cns["type"])     
  print stmt
  cursor.execute(stmt)
    
def cns_type(cns, qgene_space_poly, qgene_poly, sgene_poly, scns, qcns,qgene_space_start,qfeat):
  if qgene_space_poly.intersects(qcns):
    cns["type"] = "intron"
  elif qgene_poly.intersects(qcns) or sgene_poly.intersects(scns):
    if qfeat['strand'] == "+" and cns['qcns_start'] < qgene_space_start or qfeat['strand'] == "-" and cns['qcns_start'] > qgene_space_start:
      cns['type'] = '5-UTR'
    elif qfeat['strand'] == "+" and cns['qcns_start'] > qgene_space_start or qfeat['strand'] == "-" and cns['qcns_start'] < qgene_space_start:
      cns['type'] = '3-UTR'
  elif qfeat['strand'] == "+" and cns['qcns_start'] > qgene_space_start or qfeat['strand'] == "-" and cns['qcns_start'] < qgene_space_start: 
    cns["type"] = "3-prox_dist"
  elif qfeat['strand'] == "+" and cns['qcns_start'] < qgene_space_start or qfeat['strand'] == "-" and cns['qcns_start'] > qgene_space_start: 
    cns["type"] = "5-prox_dist"

        
def create_utr_list(utr_dict,feat, cns, letter):
  if feat['accn'] not in utr_dict.keys():
    utr_dict[feat["accn"]] = [feat['start'], feat["end"]]
  if cns['type'] == "UTR":
    utr_dict[feat['accn']].append((cns['{0}cns_start'.format(letter)]))
    utr_dict[feat['accn']].append((cns['{0}cns_end'.format(letter)]))

  

  
"postional 5' 3 'distance from utr - cns...."
"5' or 3' if < then gene_poly = 5 prime.... if grter then gene space poly then 3 prinme.."
  #   create intron list

    # qgene_poly = LineString([(0,qgene['start']),(0,cns['end'])])
    # sgene_poly = LineString([(0,sgene['start']),(0,cns['end'])])
    # 
if __name__ == "__main__":
  pass
  # import optparse
  # parser = optparse.OptionParser()
  # parser.add_option("--qbed", dest="qbed", help="bed file of the query")
  # parser.add_option("--sbed", dest="sbed", help="bed file of the subject")
  # parser.add_option("--cns", dest="cns", help="path to the cns file created by find_cns.py")
  # 
  # (options, _) = parser.parse_args()
  # 

  # x= main("/Users/gturco/code/freeling_lab/find_cns_gturco/pipeline/scripts/post_processing/find_cns_cns_test.pck","/Users/gturco/code/freeling_lab/find_cns_gturco/pipeline/scripts/post_processing/query_test.bed","/Users/gturco/code/freeling_lab/find_cns_gturco/pipeline/scripts/post_processing/subject_test.bed")
  #x= main("/Users/gturco/data/find_cns_rice_v6_sorg_v1_gturco_2011_06_14app.pck","/Users/gturco/data/rice_v6.nolocaldups.with_new.all.bed","/Users/gturco/data/sorghum_v1.nolocaldups.with_new.all.bed")
  
  # 


#need to remove all cns rna ones.....
