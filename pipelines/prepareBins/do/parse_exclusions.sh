#!/bin/bash

# parse bad/missing genome regions to BED format

cut -f 1-3 | # 3-column bed, i.e., spans only
sort -k1,1 -k2,2n |
bedtools merge -i stdin | # combine overlapping gaps and bad regions into single spans
gzip -c
