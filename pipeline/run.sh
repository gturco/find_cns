# on new synteny, run these. and put .bed and .fasta files in data/ directory.
#  perl export_to_bed.pl -fasta_name rice_v6.fasta -dsg 8163 -name_re "^Os\d\dg\d{5}$" > rice_v6.bed
# perl export_to_bed.pl -fasta_name sorghum_v1.4.fasta -dsg 95 > sorghum_v1.4.bed
# perl export_to_bed.pl -fasta_name brachy_v1.fasta -dsg 8120 > brachy_v1.bed
# use OrganismView in CoGe to look up other dsgs.


ORGA=rice_v6
ORGB=sorghum_v1
QUOTA=1:1
NCPU=8

#############################################
# dont edit below here
#############################################
DIR=data/${ORGA}_${ORGB}/

#coannotate:
     # this creates .all.bed files.
#     python scripts/create_json.py \
#       --query $ORGA \
#       --subject $ORGB
#     coannotate.py $DIR/${ORGA}_${ORGB}.json
#     # this will overwrite the genomic.masked.fasta from the coannotate.py call above.
#     # so just delete them to avoid asking.
#     mkdir $DIR/old
#     mv $DIR/${ORGA}.bed $DIR/old/
#     mv  $DIR/${ORGB}.bed  $DIR/old/
#     mv  $DIR/${ORGA}.all.bed $DIR/${ORGA}.bed  
#     mv $DIR/${ORGB}.all.bed $DIR/${ORGB}.bed
#
#sh quota.sh $DIR/${ORGA} $DIR/${ORGB} $QUOTA $NCPU
#python scripts/find_cns.py \
#        -q $DIR/${ORGA}.fasta --qbed $DIR/${ORGA}.bed \
#        -s $DIR/${ORGB}.fasta --sbed $DIR/${ORGB}.bed \
#        -p $DIR/${ORGA}_${ORGB}.pairs.txt \
#        -F T \
#        -n 8 \
#        --qpad 12000 \
#        --spad 12000 \
#        --blast_path ~/src/blast-2.2.25/bin/bl2seq \
#        --pair_fmt pair > $DIR/${ORGA}_${ORGB}.cns.txt

#postprocessing......  
#cns_sequence:
#python scripts/cns_to_fasta.py \
#                -c $DIR/${ORGA}_${ORGB}.cns.txt \
#                --qfasta $DIR/${ORGA}.genomic.masked.fasta \
#                --sfasta $DIR/${ORGB}.genomic.masked.fasta \
#                --qorg ${ORGA} \
#                --sorg ${ORGB} \
#                --min_len=18 \
#                > $DIR/${ORGA}_${ORGB}.cns.fasta
#proteins_and_rna:
#### THIS is CDS/protein stuff.
#wget -O data/at_protein.fasta ftp://ftp.arabidopsis.org/home/tair/Sequences/blast_datasets/TAIR9_blastsets/TAIR9_pep_20090619
#wget -O data/os_protein.fasta ftp://ftp.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_6.1/all.dir/all.p
#bblast.py -p blastx -d data/at_protein.fasta -i $DIR/${ORGA}_${ORGB}.cns.fasta -e 0.01 -m 8 -a 6 -o $DIR/at_protein.blast
#bblast.py -p blastx -d data/os_protein.fasta -i $DIR/${ORGA}_${ORGB}.cns.fasta -e 0.01 -m 8 -a 6 -o $DIR/os_protein.blast

#python scripts/find_exons.py \
#                 -q ${ORGA}\
#                 -s ${ORGB}\
#                 -o $DIR \
#                 $DIR/at_protein.blast $DIR/os_protein.blast
#NEED TO EDIT find_rna.py SO IT looks for the correct cns fasta file
#python scripts/find_rna.py -g data/thaliana_v9.gff \
#         -f data/thaliana_v9.fasta \
#         -b $DIR/${ORGA}_${ORGB}_cns_vs_at_rnas.blast \
#     -q ${ORGA}  \
#       -s ${ORGB} \
#   -o $DIR \
#   -d data/thaliana_v9.description

#getrna:
 # have to modify below file to be valid(er) gff3 and remove the chromosome types from teh body.
 #wget -O data/thaliana_v9.gff ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR9_genome_release/TAIR9_gff3/TAIR9_GFF3_genes_transposons.gff
 #wget -O data/thaliana_v9.description ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR9_genome_release/TAIR9_functional_descriptions
 #wget -O data/sativa_v6.1.description ftp://ftp.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_6.1/all.dir
 #wget -O data/thaliana_v9.fasta ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR9_genome_release/TAIR9_chr_all.fas
 # perl -pi -e "s/>(.).*/>\$1/" data/thaliana_v9.fasta
 # THIS is NON-cds stuff.

#CREATE PARAPLOGY AND ORTHOLOGY WITH RAW FILTERED FILE
#shuffle_protein_cns:
#  # add protein, rna cnss to the gene flat file
#  # and remove them from the _cns.txt
#  # this will create *.with_new.bed, _cns.real.txt, # with_new.*
#  python scripts/shuffle_protein_cns.py \
#    --qbed $DIR/${ORGA}.nolocaldups.bed \
#    --sbed $DIR/${ORGB}.nolocaldups.bed \
#    --cns  $DIR/${ORGA}_${ORGB}.cns.txt \
#    --paralogy  $DIR/${ORGA}_${ORGB}.paralogy \
#    --orthology $DIR/${ORGA}_${ORGB}.orthology \
#    # creates: $DIR/${ORGA}_${ORGB}.quota.with_new.orthology
#
python scripts/assign.py \
      --qbed $DIR/${ORGA}.nolocaldups.with_new.bed \
      --sbed $DIR/${ORGB}.nolocaldups.with_new.bed \
      --cns $DIR/${ORGA}_${ORGB}.cns.real.txt \
      --pairs $DIR/${ORGA}_${ORGB}.pairs.txt \
      --qdsid 9109 \
      --sdsid 95 \
      --qpad 15000 \
      --spad 15000 \
      --pair_fmt pair > $DIR/${ORGA}_${ORGB}.cns.assigned.csv
#

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

