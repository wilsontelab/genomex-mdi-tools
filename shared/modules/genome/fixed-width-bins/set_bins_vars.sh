
# action:
#     set environment variables that define a set of fixed width bins
# expects:
#     that prepareBins has already been successfully run
#     $GENOME
#     $BIN_SIZE
#     $KMER_LENGTH
#     $N_ERRORS
# usage:
#     source $MODULES_DIR/genome/fixed-width-bins/set_bins_vars.sh

# set derivative genome variables
source $MODULES_DIR/genome/set_genome_vars.sh

# path prefixes
export FIXED_WIDTH_BINS_DIR=$MDI_GENOME_DIR/fixed_width_bins
export FIXED_WIDTH_BINS_PREFIX=$FIXED_WIDTH_BINS_DIR/$GENOME

# parsed genome files
export GAP_BAD_REGIONS_BED=$FIXED_WIDTH_BINS_PREFIX.gap+bad_regions.bed.gz
export GAP_BED_BILE=$FIXED_WIDTH_BINS_PREFIX.gap.bed.gz

# individual bin attribute files
export GC5_BED_BILE=$FIXED_WIDTH_BINS_PREFIX.gc5.size_$BIN_SIZE.bed.gz
export GENMAP_BED_FILE=$FIXED_WIDTH_BINS_PREFIX.genmap.size_$BIN_SIZE.k_$KMER_LENGTH.e_$N_ERRORS.bed.gz
export UMAP_BED_FILE=$FIXED_WIDTH_BINS_PREFIX.umap.size_$BIN_SIZE.k_$KMER_LENGTH.bed.gz

# composite bin attribute file (typically use this)
export BINS_SUFFIX=bins.size_$BIN_SIZE.k_$KMER_LENGTH.e_$N_ERRORS.bed.gz
export FIXED_BINS_BED_FILE=$FIXED_WIDTH_BINS_PREFIX.$BINS_SUFFIX

# set the single-sample output file corresponding to the input bins
export SAMPLE_BIN_COUNT_BED=$DATA_GENOME_PREFIX.count.$BINS_SUFFIX

# ensure that prepareBins had already created $FIXED_BINS_BED_FILE
if [ ! -e  $FIXED_BINS_BED_FILE ]; then
    echo "missing file: $FIXED_BINS_BED_FILE"
    exit 1
fi

# set the number of columns in the bins file (before sample additions)
export N_FIXED_BIN_COLUMNS=`zcat $FIXED_BINS_BED_FILE | head -n 1 | wc -w`
export FIXED_BIN_COLUMNS=`perl -e 'print join(",", 1..'$N_FIXED_BIN_COLUMNS')'`
export FIXED_BIN_SAMPLE_COLUMN=$(($N_FIXED_BIN_COLUMNS+1))

