
RUN_IDS=`echo $RUN_IDS | sed 's/,/ /g'`

# download files with essentially no limit on file size
echo "downloading one run file at a time"
echo $RUN_IDS | tr ' ' '\n' | xargs -I SRR prefetch --max-size 1000000000 --output-file $TASK_DIR/$EXPERIMENT_NAME/SRR.sra SRR
checkPipe

echo "done"
