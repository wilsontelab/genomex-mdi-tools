# action:
#     pre-scan FASTQ files to summarize their read QUAL distributions
# expects:
#     source $MODULES_DIR/source/set_read_file_vars.sh
#     $N_READS
#     $N_TERMINAL_BASES
# input:
#     if FASTQ files are found (.fastq.gz) they are used
#     otherwise searches for SRA (.sra) files that are converted to FASTQ in a stream
# output:
#     $PLOT_PREFIX/qual-plots.png

export INTERIM_FILE=$DATA_FILE_PREFIX.interim.txt.gz

echo "calculating average QUALs per read"
perl $ACTION_DIR/pull-quals.pl |
gzip -c > $INTERIM_FILE
checkPipe

echo "making QUAL distribution plots"
Rscript $ACTION_DIR/plot-quals.R
checkPipe

rm -f $INTERIM_FILE

echo "done"
