QUOTA_DIR=/home/gturco/src/quota-alignment/


# these are set from the command-line do not edit.
qprefix=$1
sprefix=$2 # if your files are rice_v6.fasta and rice_v6.bed, this will be rice_v6
quota=$3 # something like 2:1
cpus=$4


##########################################################
##########################################################
# dont edit below here unless you know what you're doing #
##########################################################
##########################################################

JOIN_PREFIX=`dirname ${qprefix}`/`basename ${qprefix}`_`basename ${sprefix}`

python -c "from flatfeature import Bed; b = Bed('${qprefix}.bed', '${qprefix}.fasta'); b.cds_fasta(outfile='${qprefix}.features.fasta')";
python -c "from flatfeature import Bed; b = Bed('${sprefix}.bed', '${sprefix}.fasta'); b.cds_fasta(outfile='${sprefix}.features.fasta')";
#
#
LASTZ=lastz
python scripts/blastz.py -i ${qprefix}.features.fasta \
                 -d ${sprefix}.features.fasta \
                 -a $cpus \
                 --path $LASTZ \
                 -o ${JOIN_PREFIX}.blast


python ${QUOTA_DIR}/scripts/blast_to_raw.py \
            --qbed ${qprefix}.bed \
            --sbed ${sprefix}.bed \
            --cscore 0.5 \
            --no_strip_names \
            --tandem_Nmax 20 \
            ${JOIN_PREFIX}.blast

python ${QUOTA_DIR}/quota_align.py \
            --format raw --merge \
            --Dm 20 --min_size 4 \
            --quota $quota \
            ${JOIN_PREFIX}.raw

python ${QUOTA_DIR}/scripts/qa_plot.py --qbed ${qprefix}.nolocaldups.bed \
                                   --sbed ${sprefix}.nolocaldups.bed \
                                    ${JOIN_PREFIX}.raw.filtered

python ${QUOTA_DIR}/scripts/blast_plot.py --qbed ${qprefix}.bed \
                                   --sbed ${sprefix}.bed \
                                    ${JOIN_PREFIX}.blast

python ${QUOTA_DIR}/scripts/qa_to_pairs.py --qbed ${qprefix}.nolocaldups.bed \
                                       --sbed ${sprefix}.nolocaldups.bed \
                                       ${JOIN_PREFIX}.raw.filtered > ${JOIN_PREFIX}.pairs.txt
