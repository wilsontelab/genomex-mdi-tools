# action:
#     set environment variables to guide a subsequent read alignment
# expects:
#     source $MODULES_DIR/genome/set_genome_vars.sh
# optional:
#     $USE_CRAM  [default: creates .bam file]
# usage:
#     source $MODULES_DIR/align/set_alignment_vars.sh

# set the alignment log files
export FASTP_LOG_PREFIX=$LOG_FILE_PREFIX.fastp
export BWA_LOG_FILE=$LOG_FILE_PREFIX.bwa.log
export BWA_REALIGN_LOG_FILE=$LOG_FILE_PREFIX.bwa.realign.log
export MINIMAP_LOG_FILE=$LOG_FILE_PREFIX.minimap.log

# set the product bam/cram files
if [[ "$USE_CRAM" = "" || "$USE_CRAM" = "0" || "$USE_CRAM" = "null" ]]; then
    export NAME_BAM_FILE=$DATA_GENOME_PREFIX.name.bam
    export COORDINATE_BAM_FILE=$DATA_GENOME_PREFIX.coordinate.bam
    export COORDINATE_BAM_INDEX=$COORDINATE_BAM_FILE.bai
    export NAME_REALIGNED_BAM_FILE=$DATA_GENOME_PREFIX.name.realigned.bam
    export COORDINATE_REALIGNED_BAM_FILE=$DATA_GENOME_PREFIX.coordinate.realigned.bam
    export CRAM_OUTPUT_OPTIONS=""
else
    export NAME_BAM_FILE=$DATA_GENOME_PREFIX.name.cram
    export COORDINATE_BAM_FILE=$DATA_GENOME_PREFIX.coordinate.cram
    export COORDINATE_BAM_INDEX=$COORDINATE_BAM_FILE.crai
    export NAME_REALIGNED_BAM_FILE=$DATA_GENOME_PREFIX.name.realigned.cram
    export COORDINATE_REALIGNED_BAM_FILE=$DATA_GENOME_PREFIX.coordinate.realigned.cram
    export CRAM_OUTPUT_OPTIONS="--output-fmt CRAM --reference $GENOME_FASTA"
fi

# set the product PAF file, for minimap2 long reads
if [ "$PAF_FILE_TYPE" == "" ]; then
    export NAME_PAF_FILE=$DATA_GENOME_PREFIX.name.paf.gz
    export COORDINATE_PAF_FILE=$DATA_GENOME_PREFIX.coordinate.paf.gz
else 
    export NAME_PAF_FILE=$DATA_GENOME_PREFIX.$PAF_FILE_TYPE.name.paf.gz
    export COORDINATE_PAF_FILE=$DATA_GENOME_PREFIX.$PAF_FILE_TYPE.coordinate.paf.gz
fi
