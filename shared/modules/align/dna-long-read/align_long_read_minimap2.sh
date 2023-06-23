# action:
#     align long read files to genome using minimap2
# expects:
#     source $MODULES_DIR/genome/set_genome_vars.sh
#     source $MODULES_DIR/source/set_read_file_vars.sh
#     source $MODULES_DIR/align/set_alignment_vars.sh
# optional:
#     $FORCE_ALIGNMENT   [default: don't overwrite output file]
#     $MINIMAP2_ACCURACY ["high", "low"; default if missing: high, which places cigar/cg tag in PAF]
# input (in precedence order):
#     a single unaligned sam.gz file
#     a single unaligned bam file
#     a set of FASTQ files, or
# output:
#     $NAME_PAF_FILE = PAF format, +/- CIGAR strings depending on $MINIMAP2_ACCURACY

#------------------------------------------------------------------
# set the product PAF file; abort silently if exists and not forced
#------------------------------------------------------------------
if [[ "$FORCE_ALIGNMENT" != "" && "$FORCE_ALIGNMENT" != "0" && "$FORCE_ALIGNMENT" != "false" && -e $NAME_PAF_FILE ]]; then
    echo "forcing overwrite of PAF file: $NAME_PAF_FILE"
    rm -f $NAME_PAF_FILE
fi
if [ -e $NAME_PAF_FILE ]; then
    echo "alignment file already exists"

#------------------------------------------------------------------
# check for input sequence read files
#------------------------------------------------------------------
elif [[ "$USAM_FILES" == "" && "$UBAM_FILES" == "" && "$FASTQ_FILES" == "" ]]; then
    echo "missing sequence read file(s); expected INPUT_DIR/*.unaligned.sam.gz, INPUT_DIR/*.unaligned.bam or INPUT_DIR/*.fastq.gz"
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

# ------------------------------------------------------------------
# set the input files and parser
# ------------------------------------------------------------------
if [ "$USAM_FILES" != "" ]; then
    INPUT_TYPE="usam" # convert to fastq as required by minimap2 
    PARSE_READS1="zcat $USAM_FILES"
    PARSE_READS2="samtools fastq -"
elif [ "$UBAM_FILES" != "" ]; then
    INPUT_TYPE="ubam" # convert to fastq as required by minimap2 
    PARSE_READS1="samtools fastq $UBAM_FILES"
    PARSE_READS2="cat"
else 
    INPUT_TYPE="fastq"
    PARSE_READS1="perl $SHARED_MODULE_DIR/index_fastq.pl" # for sequence retrieval during SV extraction, to avoid carrying all big sequences in PAF
    PARSE_READS2="cat"
fi

# ------------------------------------------------------------------
# set the bandwidth
# ------------------------------------------------------------------
if [ "$BANDWIDTH" == "" ]; then
    BANDWIDTH_LOG="minimap2 default"
else
    BANDWIDTH_LOG="$BANDWIDTH"
    BANDWIDTH="-r $BANDWIDTH"
fi

# ------------------------------------------------------------------
# set the accuracy
# ------------------------------------------------------------------
if [ "$MINIMAP2_ACCURACY" == "low" ]; then
    CIGAR_FLAG=""
    MINIMAP2_ACCURACY="low (cigar not present in output PAF)"
else
    CIGAR_FLAG="-c"
    MINIMAP2_ACCURACY="high (cigar present in output PAF)"
fi

#------------------------------------------------------------------
# provide log feedback
#------------------------------------------------------------------
echo "aligning long reads to genome $GENOME with minimap2"
echo "  input type: $INPUT_TYPE"
echo "  input dir:  $INPUT_DIR"
echo "  genome:     $GENOME_FASTA" 
echo "  output:     $NAME_PAF_FILE"
echo "  mode:       $ALIGNMENT_MODE"
echo "  bandwith:   $BANDWIDTH_LOG"
echo "  accuracy:   $MINIMAP2_ACCURACY"

#------------------------------------------------------------------
# process reads and align to genome; soft clip supplemental
# --cs=short if it becomes necessary to get base-level information
# --no-long-join Disable the long gap patching heuristic. (which one?)
# --rmq=no|yes Use the minigraph chaining algorithm [no]
#------------------------------------------------------------------
ALN_CPU=$(( N_CPU - 1 )) # minimap2 uses 1 additional core for IO

$PARSE_READS1 |
$PARSE_READS2 |
minimap2 -x $ALIGNMENT_MODE -t $ALN_CPU --secondary=no $BANDWIDTH $CIGAR_FLAG $MINIMAP2_INDEX - 2>$MINIMAP_LOG_FILE |
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
