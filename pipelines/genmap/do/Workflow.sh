#!/bin/bash

# validate genome files and create needed derivatives
checkWorkflowStep 1 makeMap Snakefile

# set the output directories
export GENOME=$DATA_NAME
export GENOME_DIR=$TASK_DIR
source $MODULES_DIR/genome/set_genome_vars.sh
export MAPS_DIR=$GENOME_DIR/maps
mkdir -p $MAPS_DIR
TARGET_FILE=maps/$GENOME.mappability.k_$KMER_LENGTH.e_$N_ERRORS.bedgraph

# run the nested file assembly
if [ "$SN_FORCEALL" != "" ]; then rm -rf $GENOME_DIR/.snakemake; fi
snakemake $SN_DRY_RUN $SN_FORCEALL \
--cores $N_CPU \
--snakefile $ACTION_DIR/Snakefile \
--directory $GENOME_DIR \
$TARGET_FILE
checkPipe

# report our results
echo
echo $GENOME_DIR
# tree $GENOME_DIR

finishWorkflowStep
