#!/bin/bash

echo "running self-realignment of genome"

# check for prior output files and give feedback
RUNS_FILE_PREFIX=$DATA_FILE_PREFIX.pemap.r_$READ_LENGTH.q_$MIN_MAPQ
EXISTING_FILES=`ls -l $RUNS_FILE_PREFIX.*.runs.gz 2>/dev/null`
if [ "$EXISTING_FILES" != "" ]; then
    echo "one or more output file(s) already exist"
    echo "manually delete them to force reprocessing"
    echo $RUNS_FILE_PREFIX.*.runs.gz
fi

# load genome fasta into RAM drive
export GENOME_FASTA_SHM=$SHM_DIR_WRK/genome.fa
cp $GENOME_FASTA $GENOME_FASTA_SHM
cp $GENOME_FASTA.fai $GENOME_FASTA_SHM.fai

# step through the working, canonical chromosomes
CANONICAL_CHROMS=`perl $ACTION_DIR/get_chroms.pl`
for CHROM in $CANONICAL_CHROMS; do 
echo $CHROM
export CHROMOSOME=$CHROM

# set output file; never overwrite (since it takes a long time to generate)
RUNS_FILE=$RUNS_FILE_PREFIX.$CHROM.runs.gz
if [ ! -e $RUNS_FILE ]; then

# pull faux read pairs
perl $ACTION_DIR/make_read_pairs.pl |
samtools faidx --length $READ_LENGTH --region-file - $GENOME_FASTA_SHM |
perl $ACTION_DIR/parse_faidx.pl | 

# align them back to the reference geome
bwa mem -p -Y -t $N_CPU $BWA_GENOME_FASTA - 2>/dev/null |

# count fragment endpoints that were aligned with high confidence
perl $ACTION_DIR/parse_bwa.pl | 
gzip -c | 
slurp -s 100M -o $RUNS_FILE
checkPipe

fi

done
