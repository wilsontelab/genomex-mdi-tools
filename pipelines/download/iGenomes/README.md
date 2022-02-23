---
title: iGenomes
parent: download
grand_parent: Stage 1 Pipelines 
has_children: false
nav_order: 3
---

## Illumina iGenomes

Illumina provides a convenient and consistent single source for 
downloading genome sequences.

- <https://support.illumina.com/sequencing/sequencing_software/igenome.html>

The 'mdi download iGenomes' action is simple. Just go to the web site above, 
right-click and copy the URL for the genome tar.gz file(s) you wish to download, 
and provide it as option '--urls'.

> Please note: most Wilson lab pipelines use the UCSC genome, so it is usually
 the one you want to download.

To download multiple genome file sets, enter a comma-delimited set of tar.gz
file urls as option '--urls'; they will be downloaded one at a time.

Option '--data-name' should be set to the
name of the genome, e.g., 'hg38', if you are downloading a single genome,
or a name that describes a set of multiple genomes, e.g., 'human'.

Option 'output-dir' must be set to 'GENOMES_DIR/download',
where GENOMES_DIR is the value you will subsequently
provide to other pipelines. If you would like to use the default genomes 
directory, it can be entered as:

```yml
optionFamilies:
    output:
        output-dir: $MDI_DIR/resources/genomes/$PIPELINE_NAME
        data-name:  bosTau8
    iGenomes:
        urls: http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Bos_taurus/UCSC/bosTau8/Bos_taurus_UCSC_bosTau8.tar.gz
```
