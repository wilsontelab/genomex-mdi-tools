#!/bin/bash

# process and set working paths
EXPANDED_INPUT_DIR=`echo ${INPUT_DIR}`
if [ "$DORADO_OUTPUT_TYPE" == "" ]; then
    FASTQ_FILE=${DATA_FILE_PREFIX}.fastq.gz
    USAM_FILE=${DATA_FILE_PREFIX}.unaligned.sam.gz
else 
    FASTQ_FILE=${DATA_FILE_PREFIX}.$DORADO_OUTPUT_TYPE.fastq.gz
    USAM_FILE=${DATA_FILE_PREFIX}.$DORADO_OUTPUT_TYPE.unaligned.sam.gz
fi

# begin log report
echo "calling bases"
echo "  Dorado version: "${DORADO_VERSION}
echo "  model:          "${ONT_MODEL}
echo "  input:          "${EXPANDED_INPUT_DIR}
echo "  pod5 buffer:    "${POD5_BUFFER_DIR}
if [ $DORADO_OUTPUT_FORMAT == "usam" ]; then
    echo "  output: "${USAM_FILE}
    rm -f $TMP_DIR_WRK/*.unaligned.sam.gz
    rm -f ${USAM_FILE}
else 
    echo "  output: "${FASTQ_FILE}
    rm -f $TMP_DIR_WRK/*.fastq.gz
    rm -f ${FASTQ_FILE}
fi

# set basecalling options
READ_IDS_FILE=""
if [ "$DORADO_READ_IDS" != "" ]; then # read ids file cannot be gzipped 
    if [ -f $DORADO_READ_IDS ]; then
        READ_IDS_FILE=`echo $DORADO_READ_IDS`
    else 
        READ_IDS_FILE=`echo ${DATA_FILE_PREFIX}.${GENOME}.*.qNames.txt`
    fi
    if [[ "$READ_IDS_FILE" == "" || ! -e $READ_IDS_FILE ]]; then
        echo "requested read ids file not found: $DORADO_READ_IDS"
        exit 1
    fi
    echo "  reads file: "${READ_IDS_FILE}   
    READ_IDS_FILE="--read-ids $READ_IDS_FILE"
fi
EMIT_MOVES=""
if [[ "$DORADO_EMIT_MOVES" == "true" ||  "$DORADO_EMIT_MOVES" == "TRUE" || "$DORADO_EMIT_MOVES" == "1" ]]; then
    EMIT_MOVES="--emit-moves"
fi

# initialize pod5 sources
cd ${EXPANDED_INPUT_DIR}
CHECK_COUNT=`ls -1 *.pod5 2>/dev/null | wc -l`
if [ "$CHECK_COUNT" == "0" ]; then
    POD5_FILES=(*.fast5) # support implicit fallback to fast5 instead of the preferred pod5
else
    POD5_FILES=(*.pod5)
fi
rm -rf $POD5_BUFFER_DIR/*

# functions for copy and calling a batch of pod5 files
do_batch_copy () {
    if [ "$COPY_IN_I" != "" ]; then
        BATCH_FILES=${POD5_FILES[@]:COPY_IN_I:POD5_BATCH_SIZE}
        COPY_DIR=$POD5_BUFFER_DIR/$COPY_IN_I
        mkdir $COPY_DIR  
        cp $BATCH_FILES $COPY_DIR
        checkPipe
    fi
    if [ "$COPY_OUT_I" != "" ]; then
        if [ $DORADO_OUTPUT_FORMAT == "usam" ]; then
            CALL_FILE=$TMP_DIR_WRK/$COPY_OUT_I.unaligned.sam.gz
            cat $CALL_FILE >> $USAM_FILE
        else 
            CALL_FILE=$TMP_DIR_WRK/$COPY_OUT_I.fastq.gz
            cat $CALL_FILE >> $FASTQ_FILE
        fi
        checkPipe
        rm $CALL_FILE        
    fi
}
run_dorado () {
    CALL_DIR=$POD5_BUFFER_DIR/$CALL_I
    if [ $DORADO_OUTPUT_FORMAT == "usam" ]; then  
        ${DORADO_EXECUTABLE} basecaller $EMIT_MOVES $READ_IDS_FILE $ONT_MODEL_DIR $CALL_DIR | 
        samtools view - | # strip header
        pigz -p $N_CPU -c > $TMP_DIR_WRK/$CALL_I.unaligned.sam.gz
    else 
        ${DORADO_EXECUTABLE} basecaller --emit-fastq $READ_IDS_FILE $ONT_MODEL_DIR $CALL_DIR |
        pigz -p $N_CPU -c > $TMP_DIR_WRK/$CALL_I.fastq.gz
    fi
    checkPipe
    rm -r $CALL_DIR    
}

# run basecaller one pod5 file batch at a time, working from /dev/shm to a local write buffer
# see: https://github.com/nanoporetech/dorado/issues/223
COPY_IN_I=0
echo "waiting for batch copy $COPY_IN_I"
do_batch_copy
CALL_I=$COPY_IN_I
for ((COPY_IN_I=POD5_BATCH_SIZE; COPY_IN_I < ${#POD5_FILES[@]}; COPY_IN_I+=POD5_BATCH_SIZE)); do 
    do_batch_copy &
    COPY_PID=$!
    run_dorado
    echo "waiting for batch copy $COPY_IN_I"
    wait $COPY_PID
    COPY_OUT_I=$CALL_I
    CALL_I=$COPY_IN_I
done
COPY_IN_I=""
do_batch_copy &
COPY_PID=$!
run_dorado
echo "finishing final file copy"
wait $COPY_PID
COPY_OUT_I=$CALL_I
do_batch_copy


# $ dorado --help
# Usage: dorado [options] subcommand

# Positional arguments:
# basecaller
# download
# duplex

# Optional arguments:
# -h --help               shows help message and exits
# -v --version            prints version information and exits
# -vv                     prints verbose version information and exits


# $ dorado basecaller --help
# Usage: dorado [-h] [--device VAR] [--read-ids VAR] [--max-reads VAR] [--min-qscore VAR] [--batchsize VAR] [--chunksize VAR] [--overlap VAR] [--num_runners VAR] [--modified-bases VAR...] [--modified-bases-models VAR] [--emit-fastq] [--emit-moves] model data

# Positional arguments:
#   model                         the basecaller model to run. 
#   data                          the data directory. 

# Optional arguments:
#   -h, --help                    shows help message and exits 
#   -v, --version                 prints version information and exits 
#   -v, --verbose          
#   -x, --device                  device string in format "cuda:0,...,N", "cuda:all", "metal" etc.. [default: "cuda:all"]
#   -l, --read-ids                A file with a newline-delimited list of reads to basecall. If not provided, all reads will be basecalled [default: ""]
#   -n, --max-reads               [default: 0]
#   --min-qscore                  [default: 0]
#   -b, --batchsize               if 0 an optimal batchsize will be selected [default: 0]
#   -c, --chunksize               [default: 10000]
#   -o, --overlap                 [default: 500]
#   -r, --num_runners             [default: 2]
#   --modified-bases              [nargs: 1 or more] 
#   --modified-bases-models       a comma separated list of modified base models [default: ""]
#   --emit-fastq           
#   --emit-moves

