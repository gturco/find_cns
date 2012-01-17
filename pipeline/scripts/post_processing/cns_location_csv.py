import sys
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
  db = MySQLdb.connect(host="127.0.0.1", user="root", db = "paper3")
  cursor = db.cursor()
  stmt = "INSERT INTO cns_postion_info (cns_id, qaccn,qseqid,qcns_start,qcns_end,qstrand, saccn, sseqid, scns_start, scns_end,sstrand, urls, postion) values('{0}','{1}',{2},{3},{4},'{5}','{6}','{7}',{8},{9},'{10}','{11}','{12}')".format(cns['cns_id'], cns['qaccn'], cns["qseqid"],cns["qcns_start"],cns["qcns_end"],cns["qstrand"], cns["saccn"], cns["sseqid"], cns["scns_start"], cns["scns_end"], cns["sstrand"], cns["urls"], cns["type"])     
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

def utr_present(cns_pck,query_bed_path, UTR):
  "checks to see if qaccn has utr region"
  db = MySQLdb.connect(host="127.0.0.1", user="root", db = "rice_gene_table")
  cursor = db.cursor()
  cns_handle = open(cns_pck)
  cns_pickle = pickle.load(cns_handle)
  query_bed = Bed(query_bed_path)
  for cns in cns_pickle:
    qfeat = query_bed.accn(cns['qaccn'])
    if qfeat['strand'] == "+":
      end = qfeat['end']
      start = qfeat["start"]
    else:
      end = qfeat['start']
      start = qfeat["end"]
    if UTR == 3:
      if end == min(qfeat['locs'])[0] or end == max(qfeat['locs'])[1]:
        stmt = "update MUSA_GENE_LIST_copy set MUSA_GENE_LIST_copy.3_UTR = 'ND' where MUSA_GENE_LIST_copy.Rice_MSU6_genes = '{0}'".format(cns['qaccn'])
        print stmt
        cursor.execute(stmt)
    elif UTR == 5:
      if start == min(qfeat['locs'])[0] or start == max(qfeat['locs'])[1]:
        stmt = "update MUSA_GENE_LIST_copy set MUSA_GENE_LIST_copy.5_UTR = 'ND' where MUSA_GENE_LIST_copy.Rice_MSU6_genes = '{0}'".format(cns['qaccn'])
        print stmt
        cursor.execute(stmt)

def load_in_table(table_name):
    "preprocessing create table cns_postion_info_grouped  as (select qaccn, \
    saccn, count(postion) as postion_number , postion from cns_postion_info \
    group by qaccn, saccn, postion) index qaccn, saccn"
    db = MySQLdb.connect(host="127.0.0.1", user="root", db = "paper3")
    cursor = db.cursor()
    colm_dict = {"3-distal":"3_distal",
            "3-proximal":"3_proximal",
            "3-UTR":"3_UTR",
            "5-distal":"5_distal",
            "5-proximal":"5_proximal",
            "5-UTR":"5_UTR",
            "intron": "intron"}
    keys =['5-distal','5-proximal','5-UTR','intron','3-UTR','3-proximal','3-distal']
    for key in keys:
        value = colm_dict[key]
        add_col = "ALTER TABLE {0} ADD COLUMN {1} int".format(table_name,value)
        cursor.execute(add_col)
        update = "UPDATE {0}, cns_postion_info_grouped SET {0}.{1} =cns_postion_info_grouped.postion_number WHERE {0}.qaccn =cns_postion_info_grouped.qaccn AND cns_postion_info_grouped.saccn = {0}.saccn and cns_postion_info_grouped.postion = '{2}'".format(table_name,value,key)
        #print >>sys.stderr,  update
        cursor.execute(update)
        update_null = "UPDATE {0} SET {0}.{1} = 0 where {0}.{1} is null".format(table_name, value)
        cursor.execute(update_null)
    cursor.close()
    db.commit()
    db.close()


  
"postional 5' 3 'distance from utr - cns...."
"5' or 3' if < then gene_poly = 5 prime.... if grter then gene space poly then 3 prinme.."
  #   create intron list

    # qgene_poly = LineString([(0,qgene['start']),(0,cns['end'])])
    # sgene_poly = LineString([(0,sgene['start']),(0,cns['end'])])
    # 
if __name__ == "__main__":
  #pass
  # import optparse
  # parser = optparse.OptionParser()
  # parser.add_option("--qbed", dest="qbed", help="bed file of the query")
  # parser.add_option("--sbed", dest="sbed", help="bed file of the subject")
  # parser.add_option("--cns", dest="cns", help="path to the cns file created by find_cns.py")
  # 
  # (options, _) = parser.parse_args()
  # 

  # x= main("/Users/gturco/code/freeling_lab/find_cns_gturco/pipeline/scripts/post_processing/find_cns_cns_test.pck","/Users/gturco/code/freeling_lab/find_cns_gturco/pipeline/scripts/post_processing/query_test.bed","/Users/gturco/code/freeling_lab/find_cns_gturco/pipeline/scripts/post_processing/subject_test.bed")
  ##### REMINDER CHANGE TABLE INSERTS INTO!! ##############################################
 # x= main("/Users/gturco/data/paper3/paper3_rice_b_sorghum_v1_gturco_2011_4_11app_real.pck","/Users/gturco/data/paper3/rice_b_sorg.nolocaldups.with_new.all.bed","/Users/gturco/data/paper3/sorghum_v1.nolocaldups.with_new.all.bed")
  #x=main("/Users/gturco/find_cns_thaliana_v10_thaliana_v10_gturco_2011_22_11app.pck","/Users/gturco/tair_10.nolocaldups.with_new.all.bed","/Users/gturco/tair_10.nolocaldups.with_new.all.bed")
  load_in_table("rice_b_sorghum_v1_gturco_2011_4_11app_real_grouped")
  #utr_present("/Users/gturco/data/find_cns_3_UTR.pck","/Users/gturco/data/rice_v6.nolocaldups.with_new.all.bed", 3)
  #utr_present("/Users/gturco/data/find_cns_5_UTR.pck","/Users/gturco/data/rice_v6.nolocaldups.with_new.all.bed", 5)
  #need to remove all cns rna ones.....
