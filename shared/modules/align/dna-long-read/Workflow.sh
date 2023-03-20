#!/bin/bash

# set derivative environment variables
export SHARED_SUITE=genomex-mdi-tools
export SHARED_MODULES_DIR=$SUITES_DIR/$SHARED_SUITE/shared/modules
export SHARED_MODULE_DIR=$SHARED_MODULES_DIR/align/dna-long-read
source $SHARED_MODULES_DIR/genome/set_genome_vars.sh

# set working directory to INPUT_DIR to avoid too-long argument list with multiple FASTQ
cd $INPUT_DIR
export FASTQ_FILES=`ls *.fastq.gz 2>/dev/null`

# always force CRAM output since ONT and CLR reads might have CIGAR strings too long for BAM
export USE_CRAM=TRUE
source $SHARED_MODULES_DIR/align/set_alignment_vars.sh

# set the sort parameters
source $SHARED_MODULES_DIR/utilities/shell/create_temp_dir.sh
export SORT_RAM_PER_CPU_INT=$(($RAM_PER_CPU_INT - 1000000000))

# align reads to genome and create temporary output files
runWorkflowStep 1 align $SHARED_MODULE_DIR/align_long_read_minimap2.sh

# clean up
rm -r $TMP_DIR_WRK

