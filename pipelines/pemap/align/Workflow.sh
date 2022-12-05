#!/bin/bash

# initialize the genome
source $MODULES_DIR/genome/set_genome_vars.sh

# create a shared memory location for sequence retrieval
source $MODULES_DIR/utilities/shell/create_shm_dir.sh

# and a directory for sorting
source $MODULES_DIR/utilities/shell/create_temp_dir.sh

# discover cells and align them to genome one at a time, with parallelization
runWorkflowStep 1 align align.sh
