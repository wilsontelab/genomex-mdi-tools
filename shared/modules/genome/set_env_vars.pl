#!/usr/bin/perl
use strict;
use warnings;

# required by pipeline.yml to set genome-related environment variables at launch time

$ENV{DATA_GENOME_PREFIX} = "$ENV{DATA_FILE_PREFIX}.$ENV{GENOME}";
$ENV{PLOT_GENOME_PREFIX} = "$ENV{PLOTS_DIR}/$ENV{DATA_NAME}.$ENV{GENOME}";
$ENV{BAM_FILE}           = "$ENV{DATA_GENOME_PREFIX}.bam";
$ENV{INS_SIZES_FILE}     = "$ENV{DATA_GENOME_PREFIX}.insertSizes.txt";
$ENV{COUNTS_FILE}        = "$ENV{DATA_GENOME_PREFIX}.counts.txt";

1;
