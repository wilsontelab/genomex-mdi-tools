#!/bin/bash

# aggregate mappability per $BIN_SIZE bp bins
# approach is to take a weighted average of non-gapped bases in a bin

# for umap, add back missing runs with 0 values (except for known gaps)
perl $MODULES_DIR/genome/add_missing_runs.bedgraph.pl |
bedtools intersect -v -a stdin -b <(zcat $GENOME.gap.bed.gz) | # remove expanded runs in gaps

# calculate value per bin
perl $ACTION_DIR/average_bin_scores.bedgraph.pl |
awk 'BEGIN{OFS="\t"}{
    $4 = NR;
    $5 = $5 >= 0 ? int($5 * 10000 + 0.5) / 10000 : -1;
    print $0;
}' |
gzip -c


#--------------------------------------------------------------------------------
# umap
#--------------------------------------------------------------------------------
# has missing spans for real, non-gap bases whose values are implied zeros (need to add back)
# values represent the fraction of k-mers crossing position i that are uniquely mappable
# values range from 0 to 1 per base
#--------------------------------------------------------------------------------
#track type=bedGraph name='umap M100' description='umap multi-read mappability of hg38 with k100'
#chr1 10008 10009  0.0
#chr1 10009 10010  0.01
#chr1 10010 10011  0.02
#chr1 10011 10012  0.03
#chr1 10012 10013  0.04
#chr1 10013 10014  0.05
#chr1 10014 10015  0.06
#chr1 10015 10016  0.07
#chr1 10016 10017  0.08
#chr1 10017 10018  0.09
#chr1 10018 10019  0.1
#chr1 10019 10020  0.11
#chr1 10020 10021  0.12
#chr1 10021 10022  0.13
#chr1 10022 10023  0.14
#chr1 10023 10024  0.15
#chr1 10024 10025  0.16
#chr1 10025 10026  0.17
#chr1 10026 10027  0.18
#chr1 10027 10028  0.19
#chr1 10028 10029  0.2
#chr1 10029 10030  0.21
#chr1 10030 10031  0.22
#chr1 10031 10032  0.23

#--------------------------------------------------------------------------------
# genmap
#--------------------------------------------------------------------------------
# contains values for all bases except those where the kmer crosses a genome gap
# values represent mappability of the kmer whose leftmost base is at position i
# values range from ~0 to 1 per base
#--------------------------------------------------------------------------------
#chr1 9999 10000 0.00238663
#chr1 10000 10001  0.00210084
#chr1 10001 10003  0.00209205
#chr1 10003 10004  0.00210084
#chr1 10004 10005  0.00209644
#chr1 10005 10007  0.00210084
#chr1 10007 10009  0.00209205
#chr1 10009 10010  0.00238663
#chr1 10010 10011  0.00242718
#chr1 10011 10037  1
#chr1 10037 10043  0.5
#chr1 10043 10050  0.333333
#chr1 10050 10062  1
#chr1 10062 10063  0.5
#chr1 10063 10438  1
#chr1 10438 10522  0.5
#chr1 10522 10611  1
#chr1 10611 10614  0.5
#chr1 10614 10617  1
#chr1 10617 10625  0.5
#chr1 10625 10627  0.04
#chr1 10627 10630  0.0384615
#chr1 10630 10631  0.037037
#chr1 10631 10632  0.0357143
#chr1 10632 10636  0.037037
