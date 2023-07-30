#!/bin/bash

# set derivative environment variables
export SHARED_SUITE=genomex-mdi-tools
export SHARED_MODULES_DIR=$SUITES_DIR/$SHARED_SUITE/shared/modules
export SHARED_MODULE_DIR=$SHARED_MODULES_DIR/pileup
source $SHARED_MODULES_DIR/genome/set_genome_vars.sh
source $SHARED_MODULES_DIR/align/set_alignment_vars.sh

# create temp directories
source $SHARED_MODULES_DIR/utilities/shell/create_shm_dir.sh
source $SHARED_MODULES_DIR/utilities/shell/create_temp_dir.sh
export SORT_RAM_PER_CPU_INT=$(($RAM_PER_CPU_INT - 1000000000))
rm -f $SHM_DIR_WRK/*
rm -f $TMP_DIR_WRK/*

# name the output file
export PILEUP_FILE=$DATA_GENOME_PREFIX.pileup.txt.bgz

# run the pileup to create a file readable by browser track base_pileup
runWorkflowStep 1 pileup $SHARED_MODULE_DIR/pileup.sh

# clean up
rm -r $SHM_DIR_WRK
rm -r $TMP_DIR_WRK
