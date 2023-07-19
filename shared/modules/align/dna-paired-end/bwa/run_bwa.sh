# handle read alignment by bwa
# called as needed by align_paired_bwa.sh

bwa mem $SMART_PAIRING -Y -t $N_CPU $BWA_GENOME_FASTA 2>$BWA_LOG_FILE - |

# convert to bam/cram and add mate information while still name sorted
samtools fixmate -@ $N_CPU -m $CRAM_OUTPUT_OPTIONS - - |
slurp -s 100M -o $NAME_BAM_FILE
