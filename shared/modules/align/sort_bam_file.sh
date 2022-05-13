# action:
#     coordinate-sort one or more input name-sorted bam file(s), in series
# expects:
#     source $GENOMEX_MODULES_DIR/genome/set_genome_vars.sh
#     source $GENOMEX_MODULES_DIR/align/set_alignment_vars.sh
#     source $MODULES_DIR/utilities/shell/create_temp_dir.sh
# input:
#     $NAME_BAM_FILES, or, if unset, $NAME_BAM_FILE
# output:
#     coordinate-sorted bam/cram files with ".name." replaced with ".coordinate." in file name

# get and check the files to process
if [ "$NAME_BAM_FILES" = "" ]; then
    export NAME_BAM_FILES=$NAME_BAM_FILE
fi
if [ "$NAME_BAM_FILES" = "" ]; then
    echo -e "ERROR: expected variable: NAME_BAM_FILES or NAME_BAM_FILE"
    exit 1
fi

# if not already done, coordinate-sort and index each name-sorted bam file
SAMTOOLS_RAM_PER_CPU_INT=$(($RAM_PER_CPU_INT - 1000000000))
for NBF in $NAME_BAM_FILES; do
    CBF=`echo $NBF | sed 's/name\.realigned/coordinate\.realigned/g'`

    if [ ! -f $CBF ]; then
        echo "sorting: $NBF"
        samtools sort $CRAM_OUTPUT_OPTIONS \
        --threads $N_CPU \
        -m $SAMTOOLS_RAM_PER_CPU_INT \
        -T $TMP_FILE_PREFIX.samtools.sort \
        -o $CBF $NBF
        checkPipe

        samtools index -@ $N_CPU $CBF
        checkPipe
    else
        echo "exists: $CBF" 
    fi
done

echo "done"
