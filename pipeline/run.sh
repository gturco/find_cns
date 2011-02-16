# on new synteny, run these. and put .bed and .fasta files in data/ directory.
#  perl export_to_bed.pl -fasta_name rice_v6.fasta -dsg 8163 -name_re "^Os\d\dg\d{5}$" > rice_v6.bed
# perl export_to_bed.pl -fasta_name sorghum_v1.4.fasta -dsg 95 > sorghum_v1.4.bed
# perl export_to_bed.pl -fasta_name brachy_v1.fasta -dsg 8120 > brachy_v1.bed
# use OrganismView in CoGe to look up other dsgs.


ORGA=thaliana_v9
ORGB=thaliana_v9
QUOTA=1:1
NCPU=8

#############################################
# dont edit below here
#############################################
DIR=data/${ORGA}_${ORGB}/
#
sh quota.sh $DIR/${ORGA} $DIR/${ORGB} $QUOTA $NCPU

python scripts/assign_pair.py \
        -q $DIR/${ORGA}.fasta --qbed $DIR/${ORGA}.bed \
        -s $DIR/${ORGB}.fasta --sbed $DIR/${ORGB}.bed \
        -p $DIR/${ORGA}_${ORGB}.pairs.txt \
        --pair_fmt pair > $DIR/${ORGA}_${ORGB}.bed_pairs.txt

