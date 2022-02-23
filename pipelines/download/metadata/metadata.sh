#!/bin/bash

# set the directories
GENOME=$DATA_NAME
GENOME_METADATA_DIR=$OUTPUT_DIR/../../metadata/$GENOME
mkdir -p $GENOME_METADATA_DIR
checkPipe
cd $GENOME_METADATA_DIR
mkdir -p UCSC
mkdir -p ENCODE
mkdir -p umap

echo "attempting to download metadata files one at a time"

# download files
# do NOT check for existence, some may fail for some genomes
# preserve the source files names, but put into our structured paths
function wgetFile {
    FOLDER=$1   
    FILE=$2
    URL_BASE=$3
    FLAG=$4
    URL=$URL_BASE/$FILE
    echo
    echo $URL
    wget --no-verbose -O $FOLDER/$FILE $URL $FLAG
}

# genomes gaps
wgetFile UCSC gap.txt.gz https://hgdownload.cse.ucsc.edu/goldenPath/$GENOME/database 

# gc content
wgetFile UCSC $GENOME.gc5Base.wigVarStep.gz http://hgdownload.cse.ucsc.edu/goldenpath/$GENOME/bigZips

# bad genome regions
wgetFile ENCODE $GENOME-blacklist.v2.bed.gz https://github.com/Boyle-Lab/Blacklist/raw/master/lists

# umap mappability
KMER_LENGTHS="50 100" # 24 36 also available, larger files, rarely used
for KMER_LENGTH in $KMER_LENGTHS; do
    wgetFile umap k$KMER_LENGTH.umap.bedgraph.gz https://bismap.hoffmanlab.org/raw/$GENOME --no-check-certificate
done

echo
echo "done"
