#!/bin/bash

# set the directories
IGENOMES_DIR=$OUTPUT_DIR/../iGenomes
if [ ! -e $IGENOMES_DIR ]; then mkdir $IGENOMES_DIR; fi
cd $IGENOMES_DIR
checkPipe

# parse the target genome urls
URLS=`echo $URLS | sed 's/,/ /g'`

# download genome tar.gz files
echo "downloading one genome at a time"
for URL in $URLS; do
    echo
    echo $URL
    wget --no-verbose --no-check-certificate $URL
    checkPipe
    tar xzf *.tar.gz
    checkPipe
    unlink *.tar.gz
    checkPipe
done

echo
echo "done"
