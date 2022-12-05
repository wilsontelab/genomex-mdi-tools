#!/bin/bash

# initialize the genome
source $MODULES_DIR/genome/set_genome_vars.sh

# discover cells and align them to genome one at a time, with parallelization
runWorkflowStep 1 bin bin.sh
