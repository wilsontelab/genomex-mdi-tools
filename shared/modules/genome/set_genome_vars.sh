# action:
#     set environment variables for genome file access
# expects:
#     $GENOME
#     $GENOMES_DIR
#     sub-directories in $GENOMES_DIR previously created by genomex-mdi-tools/download, or manually        
# usage:
#     source $MODULES_DIR/genome/set_genome_vars.sh

# set the file prefix for genome-specific output files
export DATA_GENOME_PREFIX=$DATA_FILE_PREFIX.$GENOME
export PLOT_GENOME_PREFIX=$PLOT_PREFIX.$GENOME

# fasta files and aligner indices
export IGENOME_DIR=`echo $GENOMES_DIR/iGenomes/*/UCSC/$GENOME`
if [ "$ALIGNMENT_MODE" = "NA" ]; then
    export ALIGNMENT_MODE=$DEFAULT_ALIGNMENT_MODE
fi
if [ -e $IGENOME_DIR ]; then # use iGenomes download, if available
    export SPECIES=`echo $IGENOME_DIR | awk 'BEGIN{FS="/"}{print $(NF-2)}'`
    export GENOME_FASTA=$IGENOME_DIR/Sequence/WholeGenomeFasta/genome.fa
    export BWA_GENOME_FASTA=$IGENOME_DIR/Sequence/BWAIndex/genome.fa
    export CHROM_FASTA_DIR=$IGENOME_DIR/Sequence/Chromosomes
    export MINIMAP2_INDEX=$IGENOME_DIR/Sequence/minimap2/$GENOME.$ALIGNMENT_MODE.mmi # built by us, not iGenomes
else # otherwise, expect to find a custom genome installation in GENOMES_DIR/GENOME
    export CUSTOM_GENOME_DIR=$GENOMES_DIR/$GENOME
    export SPECIES=$GENOME
    export GENOME_FASTA=$CUSTOM_GENOME_DIR/$GENOME.fa
    export BWA_GENOME_FASTA=$GENOME_FASTA
    export CHROM_FASTA_DIR=$CUSTOM_GENOME_DIR/Chromosomes
    export MINIMAP2_INDEX=$CUSTOM_GENOME_DIR/minimap2/$GENOME.$ALIGNMENT_MODE.mmi
fi

# metadata directories
export GENOME_BINS_DIR=$GENOMES_DIR/bins/$GENOME
export GENOME_METADATA_DIR=$GENOMES_DIR/metadata/$GENOME
export GENOME_ENCODE_DIR=$GENOME_METADATA_DIR/ENCODE
export GENOME_UCSC_DIR=$GENOME_METADATA_DIR/UCSC
export GENOME_UMAP_DIR=$GENOME_METADATA_DIR/umap

# bad genome regions
export BAD_REGIONS_FILE=$GENOME_ENCODE_DIR/$GENOME-blacklist.v2.bed.gz

# genome size from canonical chromosomes, or all chromosomes if requested
if [[ "$USE_ALL_CHROMS" == "" || "$USE_ALL_CHROMS" == "0" ]]; then
    export GENOME_SIZE=`awk '$1!~/_/{s+=$2}END{print s}' $GENOME_FASTA.fai`
    export GENOME_CHROMS=`cat $GENOME_FASTA.fai | cut -f1 | grep -v _ | grep -v chrEBV` #  | sort | uniq
else
    export GENOME_SIZE=`awk '{s+=$2}END{print s}' $GENOME_FASTA.fai`
    export GENOME_CHROMS=`cat $GENOME_FASTA.fai | cut -f1` # | sort | uniq
fi
