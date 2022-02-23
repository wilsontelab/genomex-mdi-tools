#!/bin/bash

# validate genome files and create needed derivatives
checkWorkflowStep 1 makeBinsFile Snakefile

# set the output directory to a sub-directory of the parent genome
export GENOME=$DATA_NAME
BINS_DIR=$TASK_DIR/fixed_width_bins
TARGET_FILE=$GENOME.bins.size_$BIN_SIZE.k_$KMER_LENGTH.e_$N_ERRORS.bed.gz

# run the nested file assembly
if [ "$SN_FORCEALL" != "" ]; then rm -rf $BINS_DIR/.snakemake; fi
snakemake $SN_DRY_RUN $SN_FORCEALL \
--cores $N_CPU \
--snakefile $ACTION_DIR/Snakefile \
--directory $BINS_DIR \
$TARGET_FILE
checkPipe

# report our results
TARGET_FILE=$BINS_DIR/$TARGET_FILE
echo
ls -l $BINS_DIR
echo
echo $TARGET_FILE
zcat $TARGET_FILE | head
echo "..."
zcat $TARGET_FILE | tail

finishWorkflowStep
