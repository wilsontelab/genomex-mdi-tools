#!/bin/bash

# if needed, sort and index NAME_BAM_FILE to COORDINATE_BAM_FILE
if [ ! -f $COORDINATE_BAM_FILE ]; then
    echo "sorting alignments by coordinate"
    slurp -s 500M $NAME_BAM_FILE |
    samtools sort $CRAM_OUTPUT_OPTIONS --threads $N_CPU -m $SORT_RAM_PER_CPU_INT -T $TMP_FILE_PREFIX.samtools.sort - |
    slurp -s 500M -o $COORDINATE_BAM_FILE
    checkPipe
fi
if [ ! -f $COORDINATE_BAM_INDEX ]; then
    echo "indexing bam/cram file"
    samtools index $COORDINATE_BAM_FILE
    checkPipe
fi

# assemmble the pileup, one chromosome per thread
echo "assembling pileup"
perl $SHARED_MODULE_DIR/pileup.pl

# merges into a single pileup file for convenience
echo "merging pileup"
zcat $SHM_DIR_WRK/pileup.*.txt.gz | 
bgzip -c --threads $N_CPU |
slurp -s 500M -o $PILEUP_FILE
checkPipe

# note: this file has stripped the "chr" prefix to limit bgz file size
echo "indexing pileup"
tabix -s 1 -b 2 -e 2 --force $PILEUP_FILE
checkPipe
