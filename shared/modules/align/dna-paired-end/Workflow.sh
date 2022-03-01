#!/bin/bash

# set derivative environment variables
export SHARED_SUITE=genomex-mdi-tools
export SHARED_MODULES_DIR=$SUITES_DIR/$SHARED_SUITE/shared/modules
export SHARED_MODULE_DIR=$SHARED_MODULES_DIR/align/dna-paired-end
source $SHARED_MODULES_DIR/genome/set_genome_vars.sh
source $SHARED_MODULES_DIR/source/set_read_file_vars.sh
source $SHARED_MODULES_DIR/align/set_alignment_vars.sh

# set the sort parameters
source $SHARED_MODULES_DIR/utilities/shell/create_temp_dir.sh
export SORT_RAM_PER_CPU_INT=$(($RAM_PER_CPU_INT - 1000000000))

# align reads to genome and create temporary output files
runWorkflowStep 1 align $SHARED_MODULE_DIR/align_paired_bwa.sh

# clean up
rm -r $TMP_DIR_WRK
