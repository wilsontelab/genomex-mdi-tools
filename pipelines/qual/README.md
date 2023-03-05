---
title: qual
parent: Stage 1 Pipelines
has_children: true
nav_order: 100
---

## Summarize FASTQ file quality scores

'qual' provides a simple utility that pre-scans FASTQ files
to provide information on the distribution of QUAL scores.
It will help you intelligently set quality
filter options for genome alignment pipelines and modules. 

### Required inputs

The pipeline requires either one or two FASTQ files as input
provided via option `input-dir`. An equal
number of reads are taken from the head of each file.

## Quality assessment regions and outputs

The first value distribution is for the 
average Phred-scaled QUAL score over the entire read.

Additionally, two average QUAL values are calculated for
the first and last `--n-terminal-bases` of each read.

Finally, a ratio is calculated of the last over the first
`--n-terminal-bases` of each read, which reveals how much
reads decayed in quality as sequencing proceeded. The best
possible reads would have a high average quality and 
a high last/first ratio.

The output is a single image file with a composite of distribution plots.
