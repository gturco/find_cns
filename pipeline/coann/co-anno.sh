#!/bin/sh
#on new synteny, run these. and put .bed and .fasta files in data/ directory.
#  perl export_to_bed.pl -fasta_name rice_v6.fasta -dsg 8163 -name_re "^Os\d\dg\d{5}$" > rice_v6.bed
# perl export_to_bed.pl -fasta_name sorghum_v1.4.fasta -dsg 95 > sorghum_v1.4.bed

# perl export_to_bed.pl -fasta_name brachy_v1.fasta -dsg 8120 > brachy_v1.bed
# use OrganismView in CoGe to look up other dsgs.


ORGA=$1
ORGB=$2
QUOTA=$3
BLAST_DIR=$4

#############################################
# dont edit below here
#############################################
DIR=data/${ORGA}_${ORGB}/
echo $DIR
python coann/create_json.py \
      --query ${ORGA} \
      --subject ${ORGB} \
      --out_dir ${DIR} \
      --blast_path ${BLAST_DIR}



coannotate.py $DIR/${ORGA}_${ORGB}.json
python coann/merge.py \
        --missed_bed ${DIR}/missed_${ORGA}_from_${ORGB}.bed \
        --missed_matches ${DIR}/missed_${ORGA}_from_${ORGB}.matches.txt \
        --old_bed ${DIR}/${ORGA}.bed \
        --out ${DIR}/${ORGA}.all.bed
python coann/merge.py \
        --missed_bed ${DIR}/missed_${ORGB}_from_${ORGA}.bed \
        --missed_matches ${DIR}/missed_${ORGB}_from_${ORGA}.matches.txt \
        --old_bed ${DIR}/${ORGB}.bed \
        --out ${DIR}/${ORGB}.all.bed
