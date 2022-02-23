---
published: false
---

## genmap: Create a GenMap Mappability Profile

'genmap' runs the genmap utility to create a mappability
profile for a given kmer length and error rate. See:

<https://academic.oup.com/bioinformatics/article/36/12/3687/5815974>

The resulting maps can be used for prepareBins or other
pipelines that incorporate read mappability information.

### Specifying Genome Name and Directory

For genmap, option 'data-name' must be set to the
UCSC format name of the reference genome, e.g., 'hg38'.

Option 'output-dir' must be set to '\<genomes-dir\>/genmap',
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
