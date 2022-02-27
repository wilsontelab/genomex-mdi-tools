# action:
#     parse and check the requested input name-sorted bam file
# expects:
#     source $GENOMEX_MODULES_DIR/align/set_alignment_vars.sh
# optional:
#     $BAM_FILE, to override $NAME_BAM_FILE
# usage:
#     source $GENOMEX_MODULES_DIR/source/check_name_bam_file.sh

# allow users to provide their own name-sorted bam file
export IS_USER_BAM=0
if [ "$BAM_FILE" != "" ]; then
    export IS_USER_BAM=1
    export NAME_BAM_FILE=$BAM_FILE
fi

# make sure the input bam file exists
if [ ! -e $NAME_BAM_FILE ]; then
    echo -e "bam file not found:\n    $NAME_BAM_FILE"
    exit 1
fi
