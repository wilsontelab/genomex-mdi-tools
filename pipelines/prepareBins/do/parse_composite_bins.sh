#!/bin/bash

# for unclear reasons, chrY is sometimes missing from mappability files from https://bismap.hoffmanlab.org/
# thus, must base our output bins on gc5base and add potentially null mappability values (i.e., -1) for chrY

INTERSECT="bedtools intersect -a stdin" # -sorted <<< sorting not that helpful, and see note above
GROUP_BY="bedtools groupby"

# version-sort bins with GC in bed score column
perl $MODULES_DIR/genome/sort_chroms.pl | # a small stream, sorts fast without special handling
        
# add exclusions to windows as number of excluded bases (merge of gaps and bad regions)
$INTERSECT -wao -b <(zcat $GENOME.gap+bad_regions.bed.gz) |
$GROUP_BY -g 1,2,3,4,5,6 -c 10 -o sum |

# then also add gaps and bad region counts as non-merged base counts for clarity and separate handling
$INTERSECT -wao -b <(zcat $GENOME.gap.bed.gz | cut -f 1-3) |
$GROUP_BY -g 1,2,3,4,5,6,7 -c 11 -o sum |
$INTERSECT -wao -b <(zcat ../../../metadata/$GENOME/ENCODE/$GENOME-blacklist.v2.bed.gz | cut -f 1-3) |
$GROUP_BY -g 1,2,3,4,5,6,7,8 -c 12 -o sum |

# add umap mappability score
# no longer supported; we now just add a value of 1 to every bin
awk 'BEGIN{OFS="\t"}{print $0, 1}' |

# add genmap mappability score
# no need to group, mappability and GC windows are both desired bin spans at this point
$INTERSECT -loj -b <(zcat $GENOME.genmap.size_$BIN_SIZE.k_$KMER_LENGTH.e_$N_ERRORS.bed.gz) |
cut -f 1-10,15 | 

# add bin numbers and compress
awk 'BEGIN{OFS="\t"}{$4=NR; print $0}' |
gzip -c
