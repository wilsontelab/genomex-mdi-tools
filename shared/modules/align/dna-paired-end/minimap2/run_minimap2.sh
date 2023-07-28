# handle short-read alignment by minimap2 as a drop-in replacement to bwa mem
# minimap2 is 3x as fast as "about as accurate" as bwa mem for paired reads >100 bp
# minimap2 implicitly handles admixed single and paired reads in an interleaved stream
# called as needed by align_paired_bwa.sh

minimap2 -ax sr -t $N_CPU -Y --secondary=no $BWA_GENOME_FASTA - 2>$MINIMAP_LOG_FILE |

# convert to bam/cram and add mate information while still name sorted
samtools fixmate -@ $N_CPU -m $CRAM_OUTPUT_OPTIONS - - |
slurp -s 100M -o $NAME_BAM_FILE
