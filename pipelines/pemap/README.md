---
published: false
---

## pemap: Create a paired end mappability map of a genome

'pemap' creates a mappability map of a genome that is
designed to closely mimic the true mappability for paired-end
reads as seen by an aligner, including all of the aligner's
specific handling of willingness to pair, MAPQ, etc.

The approach is to create a set of pseudo-randomized
genomic paired end fragments of a specific read length
such that every genome base pair is used exactly once as the 
start and end of a fragment. Those fragments are passed to
BWA MEM and aligned. Alignments are subjected to filtering
for proper pairs and a given MAPQ and the number of starts 
and ends falling at each base is counted and divided by 2.

Result are reported as a bedgraph of runs of values that
are always 0, 0.5 or 1, depending on whether the start 
and/or end positions at that based proved mappable.

### Specifying Genome Name and Directory

For genmap, option 'data-name' must be set to the
UCSC format name of the reference genome, e.g., 'hg38'.

Option 'output-dir' must be set to '\<genomes-dir\>/pemap',
where \<genomes-dir\> is the value you will subsequently
provide to other pipelines and thus should already have
other sub-directories such as iGenomes, etc. If you would 
like to use the default genomes directory, it can be entered as:

```
optionFamilies:
    output:
        output-dir: $MDI_DIR/resources/genomes/$PIPELINE_NAME 
        data-name:  hg38
```

### Required Genome Inputs

Required input files provide information about the genome
in use and must include (remember that 'data-name' is
synonymous with 'genome', per above):

```
<output-dir>/../iGenomes/<species>/UCSC/<genome>/Sequence/WholeGenomeFasta/genome.fa
```

File 'genome.fa' is from the iGenomes reference repository.
You can obtain it here:

<https://support.illumina.com/sequencing/sequencing_software/igenome.html>

or, even better, use pipeline 'download' to install it, i.e., 
'mdi download iGenomes ...'.
