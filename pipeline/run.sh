#!/bin/sh
BLAST_DIR=../cns_pipline/blast/bin/
#on new synteny, run these. and put .bed and .fasta files in data/ directory.
#  perl export_to_bed.pl -fasta_name rice_v6.fasta -dsg 8163 -name_re "^Os\d\dg\d{5}$" > rice_v6.bed
# perl export_to_bed.pl -fasta_name sorghum_v1.4.fasta -dsg 95 > sorghum_v1.4.bed

# perl export_to_bed.pl -fasta_name brachy_v1.fasta -dsg 8120 > brachy_v1.bed
# use OrganismView in CoGe to look up other dsgs.


ORGA=sorghum_v1
ORGB=brachy
QUOTA=1:1
NCPU=8
#### look on coge for dsgid information
qDSGID=95
sDSGID=9109
#############################################
# dont edit below here
#############################################
DIR=data/${ORGA}_${ORGB}/

source cns_pipeline/bin/activate

echo checking for unannotated proteins......
sh coann/co-anno.sh ${ORGA} ${ORGB} $QUOTA $BLAST_DIR
echo finding syntenic regions...
sh quota.sh $DIR/${ORGA} $DIR/${ORGB} $QUOTA $NCPU
echo finding cns...
python scripts/find_cns.py \
	-q $DIR/${ORGA}.fasta --qbed $DIR/${ORGA}.all.bed \
	-s $DIR/${ORGB}.fasta --sbed $DIR/${ORGB}.all.bed \
        -p $DIR/${ORGA}_${ORGB}.pairs.txt \
        -F T \
        -n 8 \
        --qpad 12000 \
        --spad 12000 \
        --blast_path ${BLAST_DIR}/bl2seq \
        --pair_fmt pair > $DIR/${ORGA}_${ORGB}.cns.txt 

python scripts/localdup.py \
       -q $DIR/${ORGA}.fasta --qbed $DIR/${ORGA}.all.bed \
       -s $DIR/${ORGB}.fasta --sbed $DIR/${ORGB}.all.bed \
        -p $DIR/${ORGA}_${ORGB}.pairs.txt \
       --cns_file $DIR/${ORGA}_${ORGB}.cns.txt \
        -F T \
        -n 8 \
        --qpad 12000 \
        --spad 12000 \
        --blast_path ${BLAST_DIR}/bl2seq \
        --pair_fmt pair \
       --qdups $DIR/${ORGA}.all.localdups \
       --sdups $DIR/${ORGB}.all.localdups

python scripts/post_processing/cns_to_fasta.py \
                -c $DIR/${ORGA}_${ORGB}.cns.txt.local \
                --qfasta $DIR/${ORGA}.genomic.masked.fasta \
                --sfasta $DIR/${ORGB}.genomic.masked.fasta \
                --qorg ${ORGA} \
                --sorg ${ORGB} \
                --min_len=18 \
                > $DIR/${ORGA}_${ORGB}.cns.fasta
##echo removing cns that have hits in arabidopsis as rna or protein
##wget -O data/at_protein.fasta ftp://ftp.arabidopsis.org/home/tair/Sequences/blast_datasets/TAIR10_blastsets/TAIR10_pep_20101214
##wget -O data/os_protein.fasta ftp://ftp.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_6.1/all.dir/all.pep
#bblast.py -p blastx -d data/at_protein.fasta -i $DIR/${ORGA}_${ORGB}.cns.fasta -e 0.01 -m 8 -a ${NCPU} -o $DIR/at_protein.blast
#bblast.py -p blastx -d data/os_protein.fasta -i $DIR/${ORGA}_${ORGB}.cns.fasta -e 0.01 -m 8 -a ${NCPU} -o $DIR/os_protein.blast
#
#python scripts/post_processing/find_exons.py \
#                 --cns $DIR/${ORGA}_${ORGB}.cns.txt.local \
#                 -o $DIR \
#                 $DIR/at_protein.blast $DIR/os_protein.blast
#
###NEED TO EDIT find_rna.py SO IT looks for the correct cns fasta file
### have to modify below file to be valid(er) gff3 and remove the chromosome types from teh body.
##wget -O data/thaliana_v10.gff ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes_transposons.gff
##wget -O data/thaliana_v10.description ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR10_genome_release/TAIR10_functional_descriptions
##wget -O data/thaliana_v10.fasta ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas
##perl -pi -e "s/>(.).*/>\$1/" data/thaliana_v10.fasta
##perl -pi -e "s/.{3}//" data/thaliana_v10.gff
#
# #add ncpu to run.py!!!!
#python scripts/post_processing/find_rna.py \
#	 -g data/thaliana_v9.gff \
#         -f data/thaliana_v9.fasta \
#         -b $DIR/${ORGA}_${ORGB}_cns_vs_at_rnas.blast \
#     	 -q ${ORGA}  \
#      	 -s ${ORGB} \
#   	 -o $DIR \
#   	 -d data/thaliana_v9.description
##
#cp $DIR/${ORGA}_${ORGB}.raw.filtered.local  $DIR/${ORGA}_${ORGB}.raw2.filtered.local

#python scripts/post_processing/shuffle_protein_cns.py \
#    --qbed $DIR/${ORGA}.nolocaldups.bed.local \
#    --sbed $DIR/${ORGB}.nolocaldups.bed.local \
#    --cns  $DIR/${ORGA}_${ORGB}.cns.txt.local \
#    --paralogy  $DIR/${ORGA}_${ORGB}.raw.filtered.local \
#    --orthology $DIR/${ORGA}_${ORGB}.raw2.filtered.local \
#    --pairs $DIR/${ORGA}_${ORGB}.pairs.txt.local \ 
# creates: $DIR/${ORGA}_${ORGB}.quota.with_new.orthology
#


#python scripts/post_processing/assign.py \
#      --qbed $DIR/${ORGA}.nolocaldups.bed.with_new.local \
#      --sbed $DIR/${ORGB}.nolocaldups.bed.with_new.local \
#      --cns $DIR/${ORGA}_${ORGB}.cns.txt.real.local \
#      --pairs $DIR/${ORGA}_${ORGB}.pairs.txt.local \
#      --qdsgid 11821 \
#      --sdsgid 8120 \
#      --qpad 15000 \
#      --spad 15000 \
#      --pair_fmt pair > $DIR/${ORGA}_${ORGB}.cns.assigned.csv.local

#
python scripts/post_processing/assign.py \
      --qbed $DIR/${ORGA}.nolocaldups.bed.local \
      --sbed $DIR/${ORGB}.nolocaldups.bed.local \
      --cns $DIR/${ORGA}_${ORGB}.cns.txt.local \
      --pairs $DIR/${ORGA}_${ORGB}.pairs.txt.local \
      --qdsgid $qDSGID \
      --sdsgid $sDSGID\
      --qpad 15000 \
      --spad 15000 \
      --pair_fmt pair > $DIR/${ORGA}_${ORGB}.cns.assigned.csv.local

#
#python scripts/cns_to_fasta.py \
#                -c $DIR/${ORGA}_${ORGB}.cns.real.txt \
#                --qfasta $DIR/${ORGA}.genomic.masked.fasta \
#                --sfasta $DIR/${ORGB}.genomic.masked.fasta \
#                --qorg ${ORGA} \
#                --sorg ${ORGB} \
#                > $DIR/${ORGA}_${ORGB}.cns_real.fasteI
#
