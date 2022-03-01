# action:
#     set environment variables that identify and characterize input read-pair files
# expects:
#     sequence read file(s) in one of these patterns: (used in this precedence order)
#         $INPUT_DIR/$INPUT_NAME/*.fastq.gz
#         $INPUT_DIR/$INPUT_NAME/*.sra
#         $INPUT_DIR/$DATA_NAME/*.fastq.gz
#         $INPUT_DIR/$DATA_NAME/*.sra
# usage:
#     source $MODULES_DIR/source/set_read_file_vars.sh

# set the sequence read input directory
if [[ "$INPUT_NAME" = "" || "$INPUT_NAME" = "null" ]]; then
    WORKING_INPUT_DIR=$INPUT_DIR/$DATA_NAME
else
    WORKING_INPUT_DIR=$INPUT_DIR/$INPUT_NAME
fi
if [ ! -d  $WORKING_INPUT_DIR ]; then
    echo "missing directory: $WORKING_INPUT_DIR"
    exit 1
fi

# set the sequence read input files
# set values derived from inputs
export FASTQ_FILES=`ls -d $WORKING_INPUT_DIR/*.fastq.gz 2>/dev/null`
if [ "$FASTQ_FILES" = "" ]; then
    export SRA_FILES=`ls -d $WORKING_INPUT_DIR/*.sra 2>/dev/null`
    if [ "$SRA_FILES" = "" ]; then
        echo "no .fastq.gz or .sra files found: $WORKING_INPUT_DIR"
        exit 1 
    fi
    export FASTQ_FILE1=""
    export FASTQ_FILE2=""
    export READ_LEN=`fastq-dump --stdout --split-files $SRA_FILES 2>/dev/null | head -n2 | tail -n1 | awk '{print length($0)}'`     
else
    export SRA_FILES=""
    export N_FASTQ_FILES=`echo $FASTQ_FILES | wc -w`
    if [ $N_FASTQ_FILES -gt 2 ]; then # with multiple files per read, must have the Illumina _Rx_ portion of the file name
        export FASTQ_FILE1=`ls -d $WORKING_INPUT_DIR/*.fastq.gz 2>/dev/null | grep _R1_ | tr "\n" " "`
        export FASTQ_FILE2=`ls -d $WORKING_INPUT_DIR/*.fastq.gz 2>/dev/null | grep _R2_ | tr "\n" " "`
    else # otherwise, if only two files, they can be named anything, really
        export FASTQ_FILE1=`echo $FASTQ_FILES | cut -d " " -f1`
        export FASTQ_FILE2=`echo $FASTQ_FILES | cut -d " " -f2`        
    fi
    export READ_LEN=`zcat $FASTQ_FILE1 | head -n2 | tail -n1 | awk '{print length($0)}'`     
fi
