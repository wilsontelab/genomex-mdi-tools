#!/bin/bash

# set derivative environment variables
export SHARED_SUITE=genomex-mdi-tools
export SHARED_MODULES_DIR=$SUITES_DIR/$SHARED_SUITE/shared/modules
export SHARED_MODULE_DIR=$SHARED_MODULES_DIR/nanopore/basecall
source $SHARED_MODULES_DIR/genome/set_genome_vars.sh

# create temp directories
source $SHARED_MODULES_DIR/utilities/shell/create_shm_dir.sh
source $SHARED_MODULES_DIR/utilities/shell/create_temp_dir.sh
POD5_BUFFER_DIR=$SHM_DIR_WRK
if [ "$POD5_BUFFER" = "tmp" ]; then POD5_BUFFER_DIR=$TMP_DIR_WRK; fi

# parse MDI to Dorado options in preparation for basecalling
source $SHARED_MODULE_DIR/parse_dorado_options.sh

# convert ONT read files from POD5/FAST5 to BAM, i.e., call bases (and  maybe align reads)
runWorkflowStep 1 basecall $SHARED_MODULE_DIR/basecall.sh
