# on new synteny, run these. and put .bed and .fasta files in data/ directory.
#  perl export_to_bed.pl -fasta_name rice_v6.fasta -dsg 8163 -name_re "^Os\d\dg\d{5}$" > rice_v6.bed
# perl export_to_bed.pl -fasta_name sorghum_v1.4.fasta -dsg 95 > sorghum_v1.4.bed
# perl export_to_bed.pl -fasta_name brachy_v1.fasta -dsg 8120 > brachy_v1.bed
# use OrganismView in CoGe to look up other dsgs.


ORGA=rice_v6
ORGB=rice_v6
QUOTA=1:1
NCPU=8

#############################################
# dont edit below here
#############################################
DIR=data/${ORGA}_${ORGB}/

#sh quota.sh $DIR/${ORGA} $DIR/${ORGB} $QUOTA $NCPU
python scripts/find_cns.py \
        -q $DIR/${ORGA}.fasta --qbed $DIR/${ORGA}.bed \
        -s $DIR/${ORGB}.fasta --sbed $DIR/${ORGB}.bed \
        -p $DIR/${ORGA}_${ORGB}.pairs.txt \
        -F T \
        -n 8 \
        --qpad 12000 \
        --spad 12000 \
        --blast_path ~/src/blast-2.2.25/bin/bl2seq \
        --pair_fmt pair > $DIR/${ORGA}_${ORGB}.cns.txt

python scripts/assign.py \
      --qbed $DIR/${ORGA}.nolocaldups.bed \
      --sbed $DIR/${ORGB}.nolocaldups.bed \
      --cns $DIR/${ORGA}_${ORGB}.cns.txt \
      --pairs $DIR/${ORGA}_${ORGB}.pairs.txt \
      --qorg 9109 \
      --sorg 95 \
      --pad 15000 \
      --pair_fmt pair > $DIR/${ORGA}_${ORGB}.cns.assigned.csv


# load orga
#python scripts/load_simpledb.py \
#    --db data/db/bsr.db \
#    --prefix $DIR/${ORGA} \
#    --comparison ${O:w
#RGA}_${ORGB} \
#    --qors q \
#   --assigned-cns $DIR/${ORGA}_${ORGB}.cns.assigned.csv
#echo "loaded orga"
# load orgb
#python scripts/load_simpledb.py \
#    --db data/db/bsr.db \
#    --prefix $DIR/${ORGB} \
#    --comparison ${ORGA}_${ORGB} \
#    --qors s \
#    --assigned-cns $DIR/${ORGA}_${ORGB}.cns.assigned.csv

