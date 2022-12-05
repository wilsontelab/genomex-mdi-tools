#!/bin/bash

# aggregate mappability per $BIN_SIZE bp bins
# approach is to take a weighted average of non-gapped bases in a bin

# add back missing runs with 0 values (except for known gaps)
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
#chr1 10011 10037  1            bases at these kmer positions map uniquely
#chr1 10037 10043  0.5          bases at these kmer positions have two genomic matches
#chr1 10043 10050  0.333333     bases at these kmer positions have three genomic matches, etc.
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
