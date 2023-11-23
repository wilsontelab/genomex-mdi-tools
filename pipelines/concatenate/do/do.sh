#!/bin/bash

# collect genome information
export GENOME1=`echo ${GENOMES} | awk 'BEGIN{FS=":"}{print $1}'`
export GENOME2=`echo ${GENOMES} | awk 'BEGIN{FS=":"}{print $2}'`
export ANNOTATION1=`echo ${ANNOTATIONS} | awk 'BEGIN{FS=":"}{print $1}'`
export ANNOTATION2=`echo ${ANNOTATIONS} | awk 'BEGIN{FS=":"}{print $2}'`
export GENOME1_CHROMS=`echo ${ORDERED_CHROMS} | sed 's/,/ /g' | awk 'BEGIN{FS=":"}{print $1}'`
export GENOME2_CHROMS=`echo ${ORDERED_CHROMS} | sed 's/,/ /g' | awk 'BEGIN{FS=":"}{print $2}'` 
if [ "$FORCE_CREATE" == "0" ]; then FORCE_CREATE=""; fi
if [ "$FORCE_CREATE" != "" ]; then echo -e "\n--force-create is set; all files will be created anew\n"; fi

# parse and create additional directories
IGENOMES_DIR=${GENOMES_DIR}/iGenomes
CUSTOM_DIR=${GENOMES_DIR}/custom
GENOME1_DIR=`echo ${IGENOMES_DIR}/*/UCSC/${GENOME1}`
GENOME2_DIR=`echo ${IGENOMES_DIR}/*/UCSC/${GENOME2}`
CHROMS_PATH=Sequence/Chromosomes
GENOME_DIR=${CUSTOM_DIR}/${GENOME}
mkdir -p ${GENOME_DIR}/Chromosomes
mkdir -p ${GENOME_DIR}/minimap2
METADATA_DIR=${GENOMES_DIR}/metadata
METADATA_DIR_OUT=${GENOMES_DIR}/metadata/${GENOME}
mkdir -p ${METADATA_DIR_OUT}/ENCODE
mkdir -p ${METADATA_DIR_OUT}/UCSC

# make the composite sequence files
GENOME_FASTA=${GENOME_DIR}/${GENOME}.fa
echo
if [[ ! -f ${GENOME_FASTA} || ${FORCE_CREATE} ]]; then
    echo "creating composite genome.fa"
    rm -f ${GENOME_FASTA}
    for CHROM in ${GENOME1_CHROMS}; do
        CHROM_FILE_IN=${GENOME1_DIR}/${CHROMS_PATH}/${CHROM}.fa
        CHROM_FILE_OUT=`echo $CHROM_FILE_IN | sed s/.fa/_${GENOME1}.fa/`
        CHROM_FILE_OUT=`basename $CHROM_FILE_OUT`
        echo "    $CHROM_FILE_OUT"
        cat $CHROM_FILE_IN | 
        awk '{if($0~/>/){ print $0"_'${GENOME1}'"} else { print $0 }}' >> ${GENOME_FASTA}
        checkPipe
        cp $CHROM_FILE_IN ${GENOME_DIR}/Chromosomes/$CHROM_FILE_OUT
        checkPipe
    done 
    for CHROM in ${GENOME2_CHROMS}; do
        CHROM_FILE_IN=${GENOME2_DIR}/${CHROMS_PATH}/${CHROM}.fa
        CHROM_FILE_OUT=`echo $CHROM_FILE_IN | sed s/.fa/_${GENOME2}.fa/`
        CHROM_FILE_OUT=`basename $CHROM_FILE_OUT`
        echo "    $CHROM_FILE_OUT"
        cat $CHROM_FILE_IN | 
        awk '{if($0~/>/){ print $0"_'${GENOME2}'"} else { print $0 }}' >> ${GENOME_FASTA}
        checkPipe
        cp $CHROM_FILE_IN ${GENOME_DIR}/Chromosomes/$CHROM_FILE_OUT
        checkPipe
    done 
else 
    echo "already exists: ${GENOME_FASTA}"
fi

# create genome index files
echo
if [[ ! -f ${GENOME_FASTA}.fai || ${FORCE_CREATE} ]]; then
    echo "indexing ${GENOME} faidx"
    samtools faidx ${GENOME_FASTA}
    checkPipe
else 
    echo "${GENOME} faidx index already exists"
fi
echo
if [[ ! -f ${GENOME_FASTA}.bwt || ${FORCE_CREATE} ]]; then
    echo "indexing ${GENOME} bwa"
    bwa index ${GENOME_FASTA}
    checkPipe
else 
    echo "${GENOME} bwa index already exists"
fi

# concatenate gc5Base
GC5BASE_FILE=${METADATA_DIR_OUT}/UCSC/${GENOME}.gc5Base.wigVarStep.gz
echo
if [[ ! -f ${GC5BASE_FILE} || ${FORCE_CREATE} ]]; then
    echo "concatenating gc5Base.wigVarStep.gz"
    cat <(
        zcat ${METADATA_DIR}/${GENOME1}/UCSC/${GENOME1}.gc5Base.wigVarStep.gz | # variableStep chrom=chr1 span=5
        perl -e '
            $printing = 0;
            %chroms = map { $_ => 1 } split(" ", $ENV{GENOME1_CHROMS});
            while (my $line = <STDIN>){
                if($line =~ m/chrom=(\S+)/){
                    $printing = $chroms{$1};
                    $line =~ s/(chrom=chr\S+)/$1_$ENV{GENOME1}/;
                }
                $printing and print $line;
            }
        '
    )  <(
        zcat ${METADATA_DIR}/${GENOME2}/UCSC/${GENOME2}.gc5Base.wigVarStep.gz | # variableStep chrom=chr1 span=5
        perl -e '
            $printing = 0;
            %chroms = map { $_ => 1 } split(" ", $ENV{GENOME2_CHROMS});
            while (my $line = <STDIN>){
                if($line =~ m/chrom=(\S+)/){
                    $printing = $chroms{$1};
                    $line =~ s/(chrom=chr\S+)/$1_$ENV{GENOME2}/;
                }
                $printing and print $line;
            }
        '
    ) | gzip -c > ${GC5BASE_FILE}
    checkPipe
else 
    echo "${GENOME} gc5Base file already exists"
fi

# concatenate gaps
GAP_FILE=${METADATA_DIR_OUT}/UCSC/gap.txt.gz
echo
if [[ ! -f ${GAP_FILE} || ${FORCE_CREATE} ]]; then
    echo "concatenating gap.txt.gz"
    cat <(
        zcat ${METADATA_DIR}/${GENOME1}/UCSC/gap.txt.gz  | 
        perl -e '
            %chroms = map { $_ => 1 } split(" ", $ENV{GENOME1_CHROMS});
            while(<STDIN>){
                @f = split("\t");
                $chroms{$f[1]} and print $_;
            }
        ' | 
        perl -ne '$_ =~ s/(chr\S+)/$1_$ENV{GENOME1}/; print $_;'
    ) <(
        zcat ${METADATA_DIR}/${GENOME2}/UCSC/gap.txt.gz  | 
        perl -e '
            %chroms = map { $_ => 1 } split(" ", $ENV{GENOME2_CHROMS});
            while(<STDIN>){
                @f = split("\t");
                $chroms{$f[1]} and print $_;
            }
        ' | 
        perl -ne '$_ =~ s/(chr\S+)/$1_$ENV{GENOME2}/; print $_;'
    ) | gzip -c > ${GAP_FILE}
    checkPipe
else 
    echo "${GENOME} gap file already exists"
fi

# concatenate blacklist
BLACKLIST_FILE=${METADATA_DIR_OUT}/ENCODE/${GENOME}-blacklist.v2.bed.gz
echo
if [[ ! -f ${BLACKLIST_FILE} || ${FORCE_CREATE} ]]; then
    echo "concatenating blacklist.v2.bed.gz"
    cat <(
        zcat ${METADATA_DIR}/${GENOME1}/ENCODE/${GENOME1}-blacklist.v2.bed.gz  | 
        perl -e '
            %chroms = map { $_ => 1 } split(" ", $ENV{GENOME1_CHROMS});
            while(<STDIN>){
                @f = split("\t");
                $chroms{$f[0]} and print $_;
            }
        ' | 
        perl -ne '$_ =~ s/(chr\S+)/$1_$ENV{GENOME1}/; print $_;'
    ) <(
        zcat ${METADATA_DIR}/${GENOME2}/ENCODE/${GENOME2}-blacklist.v2.bed.gz  | 
        perl -e '
            %chroms = map { $_ => 1 } split(" ", $ENV{GENOME2_CHROMS});
            while(<STDIN>){
                @f = split("\t");
                $chroms{$f[0]} and print $_;
            }
        ' | 
        perl -ne '$_ =~ s/(chr\S+)/$1_$ENV{GENOME2}/; print $_;'
    ) | gzip -c > ${BLACKLIST_FILE}
    checkPipe
else 
    echo "${GENOME} blacklist file already exists"
fi

# create metadata file (this is always created anew)
GENOME_METADATA_FILE=${GENOME_DIR}/${GENOME}.yml
cat ${ACTION_DIR}/metadata_template.yml | 
sed -e 's/__GENOME__/'$GENOME'/g' \
    -e 's/__GENOME1__/'$GENOME1'/g' \
    -e 's/__GENOME2__/'$GENOME2'/g' \
    -e 's/__ANNOTATION1__/'$ANNOTATION1'/g' \
    -e 's/__ANNOTATION2__/'$ANNOTATION2'/g' > ${GENOME_METADATA_FILE}

# report results summaries
echo
echo ${GENOME_DIR}
tree -h ${GENOME_DIR}
echo
echo ${METADATA_DIR_OUT}
tree -h ${METADATA_DIR_OUT}
echo
echo ${GENOME_METADATA_FILE}
cat ${GENOME_METADATA_FILE}
echo 
