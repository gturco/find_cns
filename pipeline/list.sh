QUOTA_DIR=../cns_pipeline/bin/quota-alignment/

## these are set from the command-line do not edit.
ORGA=$1
ORGB=$2 # if your files are rice_v6.fasta and rice_v6.bed, this will be rice_v6
quota=$3
QDSGID=$4
SDSGID=$5
dirc=data/${ORGA}_${ORGB}
qprefix=${dirc}/${ORGA}
sprefix=${dirc}/${ORGB}
JOIN_PREFIX=${dirc}/${ORGA}_${ORGB}
#QDSID=47668
#SDSID=48263



##########################################################
##########################################################
# dont edit below here unless you know what you're doing #
##########################################################
##########################################################

###### re-run local dups with max of 3 dups
python ${QUOTA_DIR}/scripts/blast_to_raw.py \
            --qbed ${qprefix}.all.bed \
            --sbed ${sprefix}.all.bed \
            --cscore 0.5 \
            --no_strip_names \
            --tandem_Nmax 4 \
            ${JOIN_PREFIX}.blast

##### self opt in quota_align
python ${QUOTA_DIR}/quota_align.py \
            --format raw --merge \
            --Dm 20 --min_size 4 \
            --quota $quota \
            ${JOIN_PREFIX}.raw


########### run make genelist ########
python scripts/post_processing/make_genelist_gt.py \
   --qflat_all ${qprefix}.all.bed \
   --sflat_all ${sprefix}.all.bed \
   --qflat_new ${qprefix}.all.nolocaldups.bed.with_new.local \
   --sflat_new ${sprefix}.all.nolocaldups.bed.with_new.local \
        --qorg ${ORGA} \
        --sorg ${ORGB} \
        --qdsgid ${QDSGID} \
        --sdsgid ${SDSGID}  \
        --qdups ${qprefix}.all.localdups \
        --sdups ${sprefix}.all.localdups \
        --qlocal ${qprefix}.all.localdups.local \
        --slocal ${sprefix}.all.localdups.local \
        --datasheet ${JOIN_PREFIX}.cns.assigned.csv.local \
        --orthology ${JOIN_PREFIX}.raw.with_new.filtered \
         --paralogy ${JOIN_PREFIX}.raw.with_new.filtered \
####### if post processing not ran run raw.filtered.local

######### run cns_location ########
python scripts/post_processing/cns_location_csv.py \
  --qbed ${qprefix}.all.nolocaldups.bed.with_new.local\
  --sbed ${sprefix}.all.nolocaldups.bed.with_new.local \
  --cns ${JOIN_PREFIX}.cns.assigned.csv.local \
  --fmt csv \
