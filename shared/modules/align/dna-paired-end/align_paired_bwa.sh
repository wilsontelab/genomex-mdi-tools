# action:
#     align paired dna-seq read files to genome with read merging using fastp and bwa
#     if $N_GPU > 0, use Parabricks fq2bam for GPU acceleration
#     if requested in CPU mode, use minimap2 as a drop-in replacement for bwa mem
# expects:
#     source $MODULES_DIR/genome/set_genome_vars.sh
#     source $MODULES_DIR/source/set_read_file_vars.sh
#     source $MODULES_DIR/align/set_alignment_vars.sh
#     $MIN_INSERT_SIZE
# optional:
#     $ADAPTER_SEQUENCE [default: merge-inherent trimming only]
#     $FORCE_ALIGNMENT  [default: don't overwrite NAME_BAM_FILE]
#     $USE_CRAM         [default: creates .bam file]
#     $UMI_FILE         [default: no UMIs]
#     $UMI_SKIP_BASES
#     $MIN_QUAL         [default: no quality filtering]
#     $N_TERMINAL_BASES 
#     $N_GPU            for GPU acceleration
#     $USE_MINIMAP2     use minimap2 instead of bwa mem
# input:
#     if FASTQ files are found (.fastq.gz) they are used
#     otherwise searches for SRA (.sra) files that are converted to FASTQ in a stream
# output:
#     $NAME_BAM_FILE
#     if "$USE_CRAM" != "", $NAME_BAM_FILE will be .cram format

#------------------------------------------------------------------
# set the product bam/cram file; abort silently if exists and not forced
#------------------------------------------------------------------
if [[ "$FORCE_ALIGNMENT" != "" && "$FORCE_ALIGNMENT" != "0" && "$FORCE_ALIGNMENT" != "false" && -e $NAME_BAM_FILE ]]; then
    echo "forcing overwrite of bam/cram file: $NAME_BAM_FILE"
    rm -f $NAME_BAM_FILE
fi
if [ -e $NAME_BAM_FILE ]; then
    echo "bam/cram file already exists"

#------------------------------------------------------------------
# check for input sequence read files
#------------------------------------------------------------------
elif [[ "$FASTQ_FILE1" != "" && "$FASTQ_FILE2" = "" ]]; then
    echo "missing fastq file(s); expected paired .fastq.gz files"
    exit 1
elif [[ "$FASTQ_FILE1" = "" && "$SRA_FILES" = "" ]]; then
    echo "missing sequence read file(s); expected two paired .fastq.gz files or a set of .sra files"
    exit 1  
else

#------------------------------------------------------------------
# initialize the alignment process
#------------------------------------------------------------------

# provide log feedback
echo "aligning paired end reads to genome $GENOME with pre-alignment merging"
echo "  read length: $READ_LEN"
if [ "$FASTQ_FILE1" = "" ]; then
    echo "  SRA files:" 
    echo "$SRA_FILES" | perl -ne 'print "    ".join("\n    ", split(" ", $_)),"\n"'
else
    echo "  read #1 fastq files:"
    echo "$FASTQ_FILE1" | perl -ne 'print "    ".join("\n    ", split(" ", $_)),"\n"'
    echo "  read #2 fastq files:"
    echo "$FASTQ_FILE2" | perl -ne 'print "    ".join("\n    ", split(" ", $_)),"\n"'
fi
echo "  genome: $BWA_GENOME_FASTA" 

# if requested, perform "extra" adapter trimming (in addition to merge-inherent trimming)
if [[ "$ADAPTER_SEQUENCE" != "NA" && "$ADAPTER_SEQUENCE" != "null" && "$ADAPTER_SEQUENCE" != "" ]]; then
    export ADAPTER_SEQUENCE="--adapter_sequence $ADAPTER_SEQUENCE"
else
    export ADAPTER_SEQUENCE=""
fi

# if requested, align paired-end reads as two single reads, i.e., without aligner pairing
if [ "$SUPPRESS_SMART_PAIRING" = "" ]; then
    export SMART_PAIRING="-p"
else
    export SMART_PAIRING=""
fi

#------------------------------------------------------------------
# set branching to either standard bwa, minimap2 or the gpu-accelerated parabricks fq2bam
#------------------------------------------------------------------
if [ "$N_GPU" != "0" ]; then
    echo "--n-gpu > 0, using Parabricks GPU acceleration"
    ALIGN_SCRIPT1="perl $SHARED_MODULE_DIR/parabricks/split_fastq_to_file.pl"
    ALIGN_SCRIPT2="bash $SHARED_MODULE_DIR/parabricks/run_fq2bam.sh"
    export PB_TMP_DIR=$TMP_DIR_WRK/parabricks
    mkdir -p $PB_TMP_DIR/input
    mkdir -p $PB_TMP_DIR/output
    mkdir -p $PB_TMP_DIR/tmp
    export UNMERGED_FILE_1="$PB_TMP_DIR/input/unmerged_1.fastq.gz"
    export UNMERGED_FILE_2="$PB_TMP_DIR/input/unmerged_2.fastq.gz"
    export MERGED_FILE="$PB_TMP_DIR/input/merged.fastq.gz"
    N_FASTP_THREADS=$N_CPU # although usually fastp is not the rate-limiting step
elif [[ "$USE_MINIMAP2" != "" && "$USE_MINIMAP2" != "0" ]]; then
    if [[ "$BANDWIDTH" == "" || "$BANDWIDTH" == "NA" ]]; then
        export BANDWIDTH_LOG="aligner default"
        export BANDWIDTH=""
    else
        export BANDWIDTH_LOG="$BANDWIDTH"
        export BANDWIDTH="-r $BANDWIDTH"
    fi
    echo "  bandwidth: $BANDWIDTH_LOG" 
    echo "--n-gpu==0, --use-minimap2 set, using CPU-based minimap2"
    ALIGN_SCRIPT1="bash $SHARED_MODULE_DIR/minimap2/run_minimap2.sh"
    ALIGN_SCRIPT2="bash $SHARED_MODULE_DIR/bwa/sort_and_index.sh" # yes, the same for minimap2 and bwa
    N_FASTP_THREADS=3 # the fastp default
else
    if [[ "$BANDWIDTH" == "" || "$BANDWIDTH" == "NA" ]]; then
        export BANDWIDTH_LOG="aligner default"
        export BANDWIDTH=""
    else
        export BANDWIDTH_LOG="$BANDWIDTH"
        export BANDWIDTH="-w $BANDWIDTH"
    fi
    echo "  bandwidth: $BANDWIDTH_LOG" 
    echo "--n-gpu==0, using CPU-based bwa mem"
    ALIGN_SCRIPT1="bash $SHARED_MODULE_DIR/bwa/run_bwa.sh"
    ALIGN_SCRIPT2="bash $SHARED_MODULE_DIR/bwa/sort_and_index.sh"
    N_FASTP_THREADS=3 # the fastp default
fi

#------------------------------------------------------------------
# process reads and align to genome
#------------------------------------------------------------------

# pull quality reads from various source types to a consistent interleaved format
perl $SHARED_MODULE_DIR/prepare_fastq.pl |

# use fastp for one-pass adapter trimming, read merging and quality filtering
# large numbers of threads do not improve the speed
# dup evaluation is memory intensive and not needed here
fastp \
--thread $N_FASTP_THREADS \
--stdin --interleaved_in --stdout \
--dont_eval_duplication \
--length_required $MIN_INSERT_SIZE $ADAPTER_SEQUENCE \
--merge --include_unmerged --correction \
--html $FASTP_LOG_PREFIX.html --json $FASTP_LOG_PREFIX.json \
--report_title \"$DATA_NAME\" 2>/dev/null | 

# tweak the way read pair merge status is reported in QNAME line
perl $SHARED_MODULE_DIR/adjust_merge_tags.pl |

# align to genome using BWA; soft-clip supplementary
$ALIGN_SCRIPT1
checkPipe

export -f checkPipe
$ALIGN_SCRIPT2

#------------------------------------------------------------------
# index reads for fast recovery of unaltered inputs, if requested
#------------------------------------------------------------------
if [ "$CREATE_FASTQ_INDEX" != "" ]; then
    echo "indexing input read pairs"
    tabix -s 1 -b 2 -e 2 -f $DATA_FILE_PREFIX".indexed_reads.bgz"
fi

echo "done"

fi
