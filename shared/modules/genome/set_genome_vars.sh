# action:
#     set environment variables for genome file access
# expects:
#     $GENOME
#     $GENOMES_DIR
#     sub-directories in $GENOMES_DIR previously created by genomex-mdi-tools/download
# usage:
#     source $MODULES_DIR/genome/set_genome_vars.sh

# set the file prefix for genome-specific output files
export DATA_GENOME_PREFIX=$DATA_FILE_PREFIX.$GENOME

# root paths to genome for different data sources
export IGENOME_DIR=`echo $GENOMES_DIR/iGenomes/*/UCSC/$GENOME`
export SPECIES=`echo $IGENOME_DIR | awk 'BEGIN{FS="/"}{print $(NF-2)}'`

# metadata directories
export GENOME_BINS_DIR=$GENOMES_DIR/bins/$GENOME
export GENOME_METADATA_DIR=$GENOMES_DIR/metadata/$GENOME
export GENOME_ENCODE_DIR=$GENOME_METADATA_DIR/ENCODE
export GENOME_UCSC_DIR=$GENOME_METADATA_DIR/UCSC
export GENOME_UMAP_DIR=$GENOME_METADATA_DIR/umap

# fasta files and aligner indices
export GENOME_FASTA=$IGENOME_DIR/Sequence/WholeGenomeFasta/genome.fa
export BWA_GENOME_FASTA=$IGENOME_DIR/Sequence/BWAIndex/genome.fa
export CHROM_FASTA_DIR=$IGENOME_DIR/Sequence/Chromosomes

# bad genome regions
export BAD_REGIONS_FILE=$GENOME_ENCODE_DIR/$GENOME-blacklist.v2.bed.gz

# genome size from canonical chromosomes
export GENOME_SIZE=`awk '$1!~/_/{s+=$2}END{print s}' $GENOME_FASTA.fai`

# chromosomes, including chr#, chrX, chrY, chrM
export GENOME_CHROMS=`cat $GENOME_FASTA.fai | cut -f1 | grep -v _ | grep -v chrEBV | sort | uniq`    
