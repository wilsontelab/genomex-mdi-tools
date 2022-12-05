#!/bin/bash

echo "parsing mappability in $BIN_SIZE bp genome bins"

# check input files
RUNS_FILE_PREFIX=$DATA_FILE_PREFIX.pemap.r_$READ_LENGTH.q_$MIN_MAPQ
if [ ! -e $RUNS_FILE_PREFIX.chr1.runs.gz ]; then
    echo "missing runs files"
    echo "please run pemap/align first"
    echo $RUNS_FILE_PREFIX.*.runs.gz
else

# set the composite output file with all canonical chromosomes present
BINS_FILE=$DATA_FILE_PREFIX.pemap.r_$READ_LENGTH.q_$MIN_MAPQ.s_$BIN_SIZE.bed
rm -rf $BINS_FILE # hard rm since using append in loop below

# set the working, canonical chromosomes
CANONICAL_CHROMS=`perl $ACTION_DIR/../align/get_chroms.pl`
for CHROM in $CANONICAL_CHROMS; do 
echo $CHROM
export CHROMOSOME=$CHROM

# collapse base-level runs to bins with average mappability over all bases
zcat $RUNS_FILE_PREFIX.$CHROM.runs.gz | 
perl $ACTION_DIR/average_bin_scores.pl >> $BINS_FILE
checkPipe

done

# compress
BINS_FILE_GZ=$BINS_FILE.gz
cat $BINS_FILE | gzip -c > $BINS_FILE_GZ
rm -rf $BINS_FILE

# log file snippet
echo
echo $BINS_FILE_GZ
zcat $BINS_FILE_GZ | head
echo "..."
zcat $BINS_FILE_GZ | tail

fi 
