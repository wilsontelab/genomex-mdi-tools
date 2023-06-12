#!/bin/bash

# process and set working paths
EXPANDED_INPUT_DIR=`echo ${INPUT_DIR}`
if [ "$DORADO_OUTPUT_TYPE" == "" ]; then
    FASTQ_FILE=${DATA_FILE_PREFIX}.fastq.gz
    UBAM_FILE=${DATA_FILE_PREFIX}.unaligned.bam
else 
    FASTQ_FILE=${DATA_FILE_PREFIX}.$DORADO_OUTPUT_TYPE.fastq.gz
    UBAM_FILE=${DATA_FILE_PREFIX}.$DORADO_OUTPUT_TYPE.unaligned.bam
fi

# begin log report
echo "calling bases"
echo "  Dorado version: "${DORADO_VERSION}
echo "  model:          "${ONT_MODEL}
echo "  input:          "${EXPANDED_INPUT_DIR}
echo "  pod5 buffer:    "${POD5_BUFFER_DIR}
if [ $DORADO_OUTPUT_FORMAT == "ubam" ]; then
    echo "  output: "${UBAM_FILE}
    rm -f $TMP_DIR_WRK/*.unaligned.bam
else 
    echo "  output: "${FASTQ_FILE}
    rm -f $TMP_DIR_WRK/*.fastq.gz
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

# run basecaller one pod5 file chunk at a time, working from /dev/shm to a local write buffer
# see: https://github.com/nanoporetech/dorado/issues/223
cd ${EXPANDED_INPUT_DIR}
POD5_FILES=(*.pod5)
rm -f $POD5_BUFFER_DIR/*.pod5
for ((i=0; i < ${#POD5_FILES[@]}; i+=POD5_BATCH_SIZE)); do 
    echo "pod5 chunk i = $i"
    POD5_CHUNK=${POD5_FILES[@]:i:POD5_BATCH_SIZE}
    cp $POD5_CHUNK $POD5_BUFFER_DIR
    if [ $DORADO_OUTPUT_FORMAT == "ubam" ]; then  
        ${DORADO_EXECUTABLE} basecaller $EMIT_MOVES $READ_IDS_FILE $ONT_MODEL_DIR $POD5_BUFFER_DIR | 
        samtools view -b - > $TMP_DIR_WRK/$i.unaligned.bam
        checkPipe
    else 
        ${DORADO_EXECUTABLE} basecaller --emit-fastq $READ_IDS_FILE $ONT_MODEL_DIR $POD5_BUFFER_DIR |
        pigz -p $N_CPU -c > $TMP_DIR_WRK/$i.fastq.gz
        checkPipe
    fi
    rm $POD5_BUFFER_DIR/*.pod5
done

# merge the output files from local disk to shared data drive
cd $TMP_DIR_WRK
echo "merging output"
if [ $DORADO_OUTPUT_FORMAT == "ubam" ]; then
    samtools cat *.unaligned.bam |
    slurp -s 100M -o ${UBAM_FILE}
    checkPipe
    rm *.unaligned.bam
else 
    zcat *.fastq.gz |
    pigz -p $N_CPU -c | 
    slurp -s 100M -o ${FASTQ_FILE}
    checkPipe
    rm *.fastq.gz
fi


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

