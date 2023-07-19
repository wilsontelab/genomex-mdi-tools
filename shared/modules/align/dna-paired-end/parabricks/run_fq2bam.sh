# determine how to access Parabricks pbrun
# default to using an existing installation, including one provided by user call to $PARABRICKS_PREP_COMMAND
# which might run "module load XXX", create an alias compatible with the server, etc.
# alternatively, $PARABRICKS_PREP_COMMAND defaults to "module load singularity"
# with singularity used to download and run the parabricks container
echo "locating pbrun"
if [ "$PARABRICKS_PREP_COMMAND" != "" ]; then 
    $PARABRICKS_PREP_COMMAND
    checkPipe
fi
if command -v pbrun &> /dev/null; then 
    PB_RUN=pbrun
else 
    if [[ "$PARABRICKS_IMAGE" == "" || "$PARABRICKS_IMAGE" == "NA" ]]; then
        CONTAINER=nvcr.io/nvidia/clara/clara-parabricks:$PARABRICKS_VERSION_TAG
        CONTAINER_DIR=$MDI_DIR/containers/parabricks
        export PARABRICKS_IMAGE=$CONTAINER_DIR/clara-parabricks-$PARABRICKS_VERSION_TAG.sif
        if [ ! -e $PARABRICKS_IMAGE ]; then 
            echo "downloading container: $CONTAINER"
            mkdir -p $CONTAINER_DIR
            singularity pull $PARABRICKS_IMAGE docker://$CONTAINER
        fi
    fi
    SINGULARITY_EXEC="singularity exec --cleanenv --nv --bind $PB_TMP_DIR:/workdir --pwd /workdir --env TCMALLOC_MAX_TOTAL_THREAD_CACHE_BYTES=268435456 $PARABRICKS_IMAGE"
    PB_RUN="$SINGULARITY_EXEC pbrun"
fi

# # set the sort options; these require Parabricks v4.1.0 or greater
# --align-only has a bug and causes segmentation fault in v4.1.1
# if [ "$BAM_SORT" = "name" ]; then
#     END_AT="align-only"
# elif [ "$BAM_SORT" = "coordinate" ]; then
#     END_AT="no-markdups"
# else
#     echo "--bam-sort both is not compatible with --n-gpu > 0"
#     exit 1
# fi

# prepare the working directory for the parabricks container
echo "preparing the working directory, including genome files"
cd $PB_TMP_DIR
BWA_GENOME_FASTA_NAME=`basename $BWA_GENOME_FASTA`
cp ${BWA_GENOME_FASTA}* input # only reliable way to ensure container can access them genome
OUT_FILE_NAME=`basename $NAME_BAM_FILE`
LOG_FILE=$LOG_FILE_PREFIX.fq2bam.txt
LOG_FILE_NAME=`basename $LOG_FILE`

# run the alignment in two stages
#   soft-clip supplemental (--bwa-options must use "=")
#   suppress secondary
#   job runs in bound $TMP_WRK_DIR
echo "running Parabricks fq2bam"

# echo
# echo "================================================================================"
# echo "BEFORE PARABRICKS, IN CONTAINER: $PARABRICKS_IMAGE"
# echo "pwd"
# $SINGULARITY_EXEC pwd
# echo "ls -lh input"
# $SINGULARITY_EXEC ls -lh input
# echo "nvidia-smi"
# $SINGULARITY_EXEC nvidia-smi
# echo "python3 --version"
# $SINGULARITY_EXEC python3 --version
# echo "================================================================================"
# echo

LOW_MEM_FLAG=""
if [[ "$LOW_MEMORY" != "0" && "$LOW_MEMORY" != "" ]]; then
    LOW_MEM_FLAG="--low-memory"
fi
run_fq2bam(){
    echo "================================================================================"
    echo "aligning $1 reads"
    $PB_RUN fq2bam \
    --logfile output/$1.$LOG_FILE_NAME \
    --out-bam output/$1.$OUT_FILE_NAME \
    $2 \
    --ref input/$BWA_GENOME_FASTA_NAME \
    --bwa-options=-Y \
    --tmp-dir tmp \
    --filter-flag 256 \
    --gpuwrite --gpusort $LOW_MEM_FLAG \
    --no-markdups # see above, --align-only buggy; if/when fixed, expose END_AT and adjust final output below    
                  # however, --no-markdups also seems buggy as the step appears to run even with this flag set
    checkPipe
}
run_fq2bam unmerged "--in-fq input/unmerged_1.fastq.gz input/unmerged_2.fastq.gz --fix-mate" 
run_fq2bam   merged "--in-se-fq input/merged.fastq.gz" 
rm -rf input

# echo
# echo "================================================================================"
# echo "AFTER PARABRICKS, ON NODE"
# echo "ls -lh output"
# ls -lh output
# echo "================================================================================"
# echo
# # ================================================================================
# # AFTER PARABRICKS, ON NODE
# # ls -lh output
# # total 220M
# # -rw-rw-r-- 1 wilsonte wilsonte 4.5K Jul 13 16:59 merged.NA12878_dev.fq2bam.txt
# # -rw-rw-r-- 1 wilsonte wilsonte 3.8K Jul 13 16:59 merged.NA12878_dev.hg38.name_chrs.txt
# # -rw-rw-r-- 1 wilsonte wilsonte  14M Jul 13 16:59 merged.NA12878_dev.hg38.name.cram
# # -rw-rw-r-- 1 wilsonte wilsonte 2.9K Jul 13 16:59 merged.NA12878_dev.hg38.name.cram.crai
# # -rw-rw-r-- 1 wilsonte wilsonte 4.6K Jul 13 16:58 unmerged.NA12878_dev.fq2bam.txt
# # -rw-rw-r-- 1 wilsonte wilsonte 4.0K Jul 13 16:58 unmerged.NA12878_dev.hg38.name_chrs.txt
# # -rw-rw-r-- 1 wilsonte wilsonte 206M Jul 13 16:58 unmerged.NA12878_dev.hg38.name.cram
# # -rw-rw-r-- 1 wilsonte wilsonte  20K Jul 13 16:58 unmerged.NA12878_dev.hg38.name.cram.crai
# # ================================================================================

# export files from node to permanent drive
echo "combining and copying out final bam/cram file"
cat output/unmerged.$LOG_FILE_NAME output/merged.$LOG_FILE_NAME > $LOG_FILE
if [ "$BAM_SORT" = "name" ]; then
    # samtools cat --threads $N_CPU $CRAM_OUTPUT_OPTIONS output/unmerged.$OUT_FILE_NAME output/merged.$OUT_FILE_NAME | 
    # slurp -s 500M -o $NAME_BAM_FILE
    # checkPipe
    samtools merge   --threads $N_CPU --reference $BWA_GENOME_FASTA -u -o - output/unmerged.$OUT_FILE_NAME output/merged.$OUT_FILE_NAME |
    samtools collate --threads $N_CPU -O $CRAM_OUTPUT_OPTIONS - $TMP_DIR_WRK/collate.tmp | # cannot use fast mode, it filters out merged reads
    slurp -s 500M -o $NAME_BAM_FILE
    checkPipe
elif [ "$BAM_SORT" = "coordinate" ]; then
    samtools merge --threads $N_CPU $CRAM_OUTPUT_OPTIONS --write-index -o $COORDINATE_BAM_FILE output/unmerged.$OUT_FILE_NAME output/merged.$OUT_FILE_NAME
    checkPipe
else
    echo "--bam-sort both is not compatible with --n-gpu > 0"
    exit 1
fi

# cleanup
rm -fr $TMP_DIR_WRK/*
