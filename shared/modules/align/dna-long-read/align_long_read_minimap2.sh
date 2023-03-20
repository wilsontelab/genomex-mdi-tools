# action:
#     align long read files to genome using minimap2
# expects:
#     source $MODULES_DIR/genome/set_genome_vars.sh
#     source $MODULES_DIR/source/set_read_file_vars.sh
#     source $MODULES_DIR/align/set_alignment_vars.sh
# optional:
#     $FORCE_ALIGNMENT  [default: don't overwrite output file]
# input:
#     a set of FASTQ files
# output:
#     $NAME_PAF_FILE = PAF format, with CIGAR strings

#------------------------------------------------------------------
# set the product bam/cram file; abort silently if exists and not forced
#------------------------------------------------------------------
if [[ "$FORCE_ALIGNMENT" != "" && "$FORCE_ALIGNMENT" != "0" && "$FORCE_ALIGNMENT" != "false" && -e $NAME_PAF_FILE ]]; then
    echo "forcing overwrite of cram file: $NAME_PAF_FILE"
    rm -f $NAME_PAF_FILE
fi
if [ -e $NAME_PAF_FILE ]; then
    echo "alignment file already exists"

#------------------------------------------------------------------
# check for input sequence read files
#------------------------------------------------------------------
elif [ "$FASTQ_FILES" = "" ]; then
    echo "missing sequence read file(s); expected INPUT_DIR/*.fastq.gz"
    exit 1  
else

# ------------------------------------------------------------------
# create and save the appropriate minimap2 alignment index
# ------------------------------------------------------------------
if [ ! -f "$MINIMAP2_INDEX" ]; then
    echo "create minimap2 index"
    echo "  "$MINIMAP2_INDEX
    mkdir -p `dirname $MINIMAP2_INDEX`
    minimap2 -x $ALIGNMENT_MODE -t 3 -d $MINIMAP2_INDEX $GENOME_FASTA
    checkPipe
fi

#------------------------------------------------------------------
# provide log feedback
#------------------------------------------------------------------
echo "aligning long reads to genome $GENOME, minimap2 mode '$ALIGNMENT_MODE'"
echo "  input: $INPUT_DIR"
echo "  genome: $GENOME_FASTA" 
echo "  output: $NAME_PAF_FILE"

#------------------------------------------------------------------
# process reads and align to genome; soft clip supplemental
# output routinely includes both supplemental and secondary, and many unmapped
#   N     FLAG
#   31003 0    # aligned, top and bottom strands
#   30824 16
#   18990 4    # unaligned
#    7802 256  # secondary, top and bottom strands
#    7596 272
#    1915 2048 # supplemental, top and bottom strands
#    1870 2064
# --cs=short if it becomes necessary to get base-level information (but it is of dubious utility)
#------------------------------------------------------------------
ALN_CPU=$(( N_CPU - 1 )) # minimap2 uses 1 additional core for IO

perl $SHARED_MODULE_DIR/index_fastq.pl | # for sequence retrieval during SV extraction, to avoid carrying all big sequences in PAF
minimap2 -x $ALIGNMENT_MODE -t $ALN_CPU --secondary=no -c $MINIMAP2_INDEX - 2>$MINIMAP_LOG_FILE |
pigz -p $N_CPU -c | 
slurp -s 100M -o $NAME_PAF_FILE
checkPipe

echo "done"

fi

# Usage: minimap2 [options] <target.fa>|<target.idx> [query.fa] [...]
# Options:
#   Indexing:
#     -H           use homopolymer-compressed k-mer (preferrable for PacBio)
#     -k INT       k-mer size (no larger than 28) [15]
#     -w INT       minimizer window size [10]
#     -I NUM       split index for every ~NUM input bases [4G]
#     -d FILE      dump index to FILE []
#   Mapping:
#     -f FLOAT     filter out top FLOAT fraction of repetitive minimizers [0.0002]
#     -g NUM       stop chain enlongation if there are no minimizers in INT-bp [5000]
#     -G NUM       max intron length (effective with -xsplice; changing -r) [200k]
#     -F NUM       max fragment length (effective with -xsr or in the fragment mode) [800]
#     -r NUM[,NUM] chaining/alignment bandwidth and long-join bandwidth [500,20000]
#     -n INT       minimal number of minimizers on a chain [3]
#     -m INT       minimal chaining score (matching bases minus log gap penalty) [40]
#     -X           skip self and dual mappings (for the all-vs-all mode)
#     -p FLOAT     min secondary-to-primary score ratio [0.8]
#     -N INT       retain at most INT secondary alignments [5]
#   Alignment:
#     -A INT       matching score [2]
#     -B INT       mismatch penalty (larger value for lower divergence) [4]
#     -O INT[,INT] gap open penalty [4,24]
#     -E INT[,INT] gap extension penalty; a k-long gap costs min{O1+k*E1,O2+k*E2} [2,1]
#     -z INT[,INT] Z-drop score and inversion Z-drop score [400,200]
#     -s INT       minimal peak DP alignment score [80]
#     -u CHAR      how to find GT-AG. f:transcript strand, b:both strands, n:don't match GT-AG [n]
#   Input/Output:
#     -a           output in the SAM format (PAF by default)
#     -o FILE      output alignments to FILE [stdout]
#     -L           write CIGAR with >65535 ops at the CG tag
#     -R STR       SAM read group line in a format like '@RG\tID:foo\tSM:bar' []
#     -c           output CIGAR in PAF
#     --cs[=STR]   output the cs tag; STR is 'short' (if absent) or 'long' [none]
#     --MD         output the MD tag
#     --eqx        write =/X CIGAR operators
#     -Y           use soft clipping for supplementary alignments
#     -t INT       number of threads [3]
#     -K NUM       minibatch size for mapping [500M]
#     --version    show version number
#   Preset:
#     -x STR       preset (always applied before other options; see minimap2.1 for details) []
#                  - map-pb/map-ont - PacBio CLR/Nanopore vs reference mapping
#                  - map-hifi - PacBio HiFi reads vs reference mapping
#                  - ava-pb/ava-ont - PacBio/Nanopore read overlap
#                  - asm5/asm10/asm20 - asm-to-ref mapping, for ~0.1/1/5% sequence divergence
#                  - splice/splice:hq - long-read/Pacbio-CCS spliced alignment
#                  - sr - genomic short-read mapping
