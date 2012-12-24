#!/bin/sh
BLAST_DIR=~/blast-2.2.25/bin/
#on new synteny, run these. and put .bed and .fasta files in data/ directory.
#  perl export_to_bed.pl -fasta_name rice_v6.fasta -dsg 8163 -name_re "^Os\d\dg\d{5}$" > rice_v6.bed
# perl export_to_bed.pl -fasta_name sorghum_v1.4.fasta -dsg 95 > sorghum_v1.4.bed

# perl export_to_bed.pl -fasta_name brachy_v1.fasta -dsg 8120 > brachy_v1.bed
# use OrganismView in CoGe to look up other dsgs.


ORGA=rice_j
ORGB=setaria_n
QUOTA=1:1
NCPU=8
#############################################
# dont edit below here
#############################################
DIR=data/${ORGA}_${ORGB}/
echo finding cns...
python scripts/localdup.py \
	--qdups $DIR/${ORGA}.localdups --qbed $DIR/${ORGA}.bed \
	--sdups $DIR/${ORGB}.localdups --sbed $DIR/${ORGB}.bed \
  -q $DIR/${ORGA}.fasta -s $DIR/${ORGB}.fasta \
  -p $DIR/${ORGA}_${ORGB}.pairs.txt \
  --cns_file $DIR/${ORGA}_${ORGB}.cns.txt \
  -F T \
        -n 8 \
        --qpad 12000 \
        --spad 12000 \
        --blast_path ${BLAST_DIR}/bl2seq \
        --pair_fmt pair

