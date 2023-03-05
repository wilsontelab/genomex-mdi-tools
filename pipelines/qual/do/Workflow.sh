#!/bin/bash

source $MODULES_DIR/source/set_read_file_vars.sh

# perform the FASTQ file QUAL scan
runWorkflowStep 1 qual $ACTION_DIR/qual.sh
