QUOTA_DIR=/home/gturco/src/quota-alignment/
## these are set from the command-line do not edit.
#qprefix=$1
#sprefix=$2 # if your files are rice_v6.fasta and rice_v6.bed, this will be rice_v6
quota=1:1
#cpus=$4
#dirc=$5
ORGA=rice_j
ORGB=sorghum_n
dirc=data/${ORGA}_${ORGB}
qprefix=${dirc}/${ORGA}
sprefix=${dirc}/${ORGB}
JOIN_PREFIX=${dirc}/${ORGA}_${ORGB}
QDSID=47668
SDSID=48263
#SDSID=47667



##########################################################
##########################################################
# dont edit below here unless you know what you're doing #
##########################################################
##########################################################
JOIN_PREFIX=`dirname ${qprefix}`/`basename ${qprefix}`_`basename ${sprefix}`

###### move orginal localdups info used for finding pairs to pairs_data
#mkdir ${dirc}/pairs_data/
#mv * ${dirc}/*.localdups ${dirc}/pairs_data/
#mv  ${dirc}/*raw* ${dirc}/pairs_data/
#mv  ${dirc}/*.nolocaldups.bed ${dirc}/quota_data/
#mv ${JOIN_PREFIX}.pairs.txt ${dirc}/quota_data/
#
###### re-run local dups with max of 3 dups
#python ${QUOTA_DIR}/scripts/blast_to_raw.py \
#            --qbed ${qprefix}.bed \
#            --sbed ${sprefix}.bed \
#            --cscore 0.5 \
#            --no_strip_names \
#            --tandem_Nmax 4 \
#            ${JOIN_PREFIX}.blast
##
##### self opt in quota_align
#python ${QUOTA_DIR}/quota_align.py \
#            --format raw --merge \
#            --Dm 20 --min_size 4 \
#            --quota $quota \
#            ${JOIN_PREFIX}.raw
#

########### run make genelist ########
python scripts/post_processing/make_genelist_gt.py \
   --qflat_all ${qprefix}.bed \
   --sflat_all ${sprefix}.bed \
   --qflat_new ${qprefix}.nolocaldups.with_new.local \
   --sflat_new ${sprefix}.nolocaldups.with_new.local \
        --qorg ${ORGA} \
        --sorg ${ORGB} \
        --qdsid ${QDSID} \
        --sdsid ${SDSID}  \
        --qdups ${qprefix}.localdups \
        --sdups ${sprefix}.localdups \
        --qlocal ${qprefix}.localdups.local \
        --slocal ${sprefix}.localdups.local \
        --datasheet ${JOIN_PREFIX}.cns.assigned_real.csv \
        --orthology ${JOIN_PREFIX}.raw.with_new.filtered \
         --paralogy ${JOIN_PREFIX}.raw.with_new.filtered \
####### if post processing not ran run raw.filtered.local

######### run cns_location ########
#python scripts/post_processing/cns_location_csv.py \
#  --qbed ${qprefix}.nolocaldups.with_new.bed \
#  --sbed ${sprefix}.nolocaldups.with_new.bed \
#  --cns ${JOIN_PREFIX}.cns.assigned_real.csv \
#  --fmt csv \
#python scripts/post_processing/cns_location_csv.py \
#  --qbed ${qprefix}.nolocaldups.bed \
#  --sbed ${sprefix}.nolocaldups.bed \
#  --cns ${JOIN_PREFIX}.cns.assigned_real.csv \
#  --fmt csv \
#
