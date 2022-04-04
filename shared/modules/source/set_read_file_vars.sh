# action:
#     set environment variables that identify and characterize input read-pair files
# expects:
#     sequence read file(s) in one of these patterns: (used in this precedence order)
#         $INPUT_DIR/$INPUT_NAME/*.fastq.gz  i.e., one named folder per sample
#         $INPUT_DIR/$INPUT_NAME_*.fastq.gz  i.e., one common folder with prefixed file names
#         $INPUT_DIR/$DATA_NAME/*.fastq.gz
#         $INPUT_DIR/$DATA_NAME_*.fastq.gz
#         $INPUT_DIR/$INPUT_NAME/*.sra  i.e., SRA files instead of FASTQ
#         $INPUT_DIR/$INPUT_NAME_*.sra
#         $INPUT_DIR/$DATA_NAME/*.sra
#         $INPUT_DIR/$DATA_NAME_*.sra
# usage:
#     source $MODULES_DIR/source/set_read_file_vars.sh

# set the sequence read input directory/prefix (supports wildcards)
if [[ "$INPUT_NAME" = "" || "$INPUT_NAME" = "null" ]]; then
    export INPUT_NAME=$DATA_NAME 
fi
WORKING_INPUT_DIR=$INPUT_DIR/$INPUT_NAME"/"
WORKING_INPUT_PREFIX=$INPUT_DIR/$INPUT_NAME"_"

# set the sequence read input files
export FASTQ_FILES=`ls -d $WORKING_INPUT_DIR*.fastq.gz 2>/dev/null`
if [ "$FASTQ_FILES" = "" ]; then
    export FASTQ_FILES=`ls -d $WORKING_INPUT_PREFIX*.fastq.gz 2>/dev/null`
fi
if [ "$FASTQ_FILES" = "" ]; then
    export SRA_FILES=`ls -d $WORKING_INPUT_DIR*.sra 2>/dev/null`
fi
if [[ "$FASTQ_FILES" = "" && "$SRA_FILES" = "" ]]; then
    export SRA_FILES=`ls -d $WORKING_INPUT_PREFIX*.sra 2>/dev/null`
fi

# set values derived from inputs
if [ "$FASTQ_FILES" != "" ]; then
    export SRA_FILES=""
    export N_FASTQ_FILES=`echo $FASTQ_FILES | wc -w`
    if [ $N_FASTQ_FILES -gt 2 ]; then # with multiple files per read, must have the Illumina _Rx_ portion of the file name
        export FASTQ_FILE1=`ls -d $FASTQ_FILES 2>/dev/null | grep _R1_ | tr "\n" " "`
        export FASTQ_FILE2=`ls -d $FASTQ_FILES 2>/dev/null | grep _R2_ | tr "\n" " "`
    else # otherwise, if only two files, they can be named anything, really
        export FASTQ_FILE1=`echo $FASTQ_FILES | cut -d " " -f1`
        export FASTQ_FILE2=`echo $FASTQ_FILES | cut -d " " -f2`        
    fi
    export READ_LEN=`zcat $FASTQ_FILE1 | head -n2 | tail -n1 | awk '{print length($0)}'`  
elif [ "$SRA_FILES" != "" ]; then
    export FASTQ_FILE1=""
    export FASTQ_FILE2=""
    export READ_LEN=`fastq-dump --stdout --split-files $SRA_FILES 2>/dev/null | head -n2 | tail -n1 | awk '{print length($0)}'`     
else
    echo -e "no .fastq.gz or .sra files found\n$INPUT_DIR + $INPUT_NAME"
    exit 1 
fi
