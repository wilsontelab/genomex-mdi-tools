# handle short-read alignment by minimap2 as a drop-in replacement to bwa mem
# minimap2 is advertised as 3x as fast and "about as accurate" as bwa mem for paired reads >100 bp
# minimap2 implicitly handles admixed single and paired reads in an interleaved stream
# called as needed by align_paired_bwa.sh

minimap2 -ax sr -t $N_CPU -Y $BANDWIDTH $BWA_GENOME_FASTA - 2>$MINIMAP_LOG_FILE |

# convert to bam/cram and add mate information while still name sorted
samtools fixmate -@ $N_CPU -m $CRAM_OUTPUT_OPTIONS - - |
slurp -s 100M -o $NAME_BAM_FILE

# -x STR
# Preset []. This option applies multiple options at the same time. It should be applied before other options because options applied later will overwrite the values set by -x.
# sr
# Short single-end reads without splicing (-k21 -w11 --sr --frag=yes -A2 -B8 -O12,32 -E2,1 -b0 -r100 -p.5 -N20 -f1000,5000 -n2 -m25 -s40 -g100 -2K50m --heap-sort=yes --secondary=no).
