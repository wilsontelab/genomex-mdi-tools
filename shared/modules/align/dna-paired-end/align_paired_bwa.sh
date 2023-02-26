# action:
#     align paired dna-seq read files to genome with read merging using fastp and bwa
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
    ADAPTER_SEQUENCE="--adapter_sequence $ADAPTER_SEQUENCE"
else
    ADAPTER_SEQUENCE=""
fi

# if requested, align paired-end reads as two single reads, i.e., without aligner pairing
if [ "$SUPPRESS_SMART_PAIRING" = "" ]; then
    SMART_PAIRING="-p"
else
    SMART_PAIRING=""
fi
    
#------------------------------------------------------------------
# process reads and align to genome
#------------------------------------------------------------------

# pull reads from various source types to a consistent interleaved format
perl $SHARED_MODULE_DIR/prepare_fastq.pl |

# use fastp for one-pass adapter trimming, read merging and quality filtering
# large numbers of threads do not improve the speed
# dup evaluation is memory intensive and not needed here
fastp \
--stdin --interleaved_in --stdout \
--dont_eval_duplication \
--length_required $MIN_INSERT_SIZE $ADAPTER_SEQUENCE \
--merge --include_unmerged --correction \
--html $FASTP_LOG_PREFIX.html --json $FASTP_LOG_PREFIX.json \
--report_title \"$DATA_NAME\" 2>/dev/null |

# tweak the way read pair merge status is reported in QNAME line
perl $SHARED_MODULE_DIR/adjust_merge_tags.pl |

# align to genome using BWA; soft-clip supplementary
bwa mem $SMART_PAIRING -Y -t $N_CPU $BWA_GENOME_FASTA - 2>$BWA_LOG_FILE |

# convert to bam/cram and add mate information while still name sorted
samtools fixmate -@ $N_CPU -m $CRAM_OUTPUT_OPTIONS - - |
slurp -s 100M -o $NAME_BAM_FILE
checkPipe

#------------------------------------------------------------------
# handle bam/cram file coordinate sorting, if requested
#------------------------------------------------------------------
if [[ "$BAM_SORT" = "coordinate" || "$BAM_SORT" = "both" ]]; then
    echo "sorting alignments by coordinate"
    slurp -s 500M $NAME_BAM_FILE |
    samtools sort $CRAM_OUTPUT_OPTIONS --threads $N_CPU -m $SORT_RAM_PER_CPU_INT -T $TMP_FILE_PREFIX.samtools.sort - |
    slurp -s 500M -o $COORDINATE_BAM_FILE
    checkPipe
    if [ "$BAM_SORT" = "coordinate" ]; then
        rm -rf $NAME_BAM_FILE
    fi
fi

echo "done"

fi
