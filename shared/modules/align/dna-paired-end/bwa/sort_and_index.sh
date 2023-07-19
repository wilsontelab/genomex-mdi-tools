#------------------------------------------------------------------
# handle bam/cram file coordinate sorting, if requested
#------------------------------------------------------------------
if [[ "$BAM_SORT" = "coordinate" || "$BAM_SORT" = "both" ]]; then

    echo "sorting alignments by coordinate"
    slurp -s 500M $NAME_BAM_FILE |
    samtools sort $CRAM_OUTPUT_OPTIONS --threads $N_CPU -m $SORT_RAM_PER_CPU_INT -T $TMP_FILE_PREFIX.samtools.sort - |
    slurp -s 500M -o $COORDINATE_BAM_FILE
    checkPipe

    echo "indexing coordinate bam file"
    samtools index $COORDINATE_BAM_FILE
    checkPipe

    if [ "$BAM_SORT" = "coordinate" ]; then
        rm -rf $NAME_BAM_FILE
    fi
fi
