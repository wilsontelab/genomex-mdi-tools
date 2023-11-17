---
title: concatenate
parent: Stage 1 Pipelines
has_children: false
nav_order: 2
---

## Making composite genomes

The 'concatenate' pipeline creates a new, custom, composite
genome from two pre-existing UCSC genomes. The resulting
genome files are appropriate for use in experiments
where DNA from two species are known to be mixed together.

### Prerequisites

You must have already run pipeline actions 
`download iGenomes` and `download metadata`
as appropriate for the two source UCSC genomes.

In addition to the genome FASTA files themselves,
each of the following metadata files must exist for each genome. 
- metadata/\<genome\>/UCSC/\<genome\>.gc5Base.wigVarStep.gz
- metadata/\<genome\>/UCSC/gap.txt.gz
- metadata/\<genome\>/ENCODE/\<genome\>-blacklist.v2.bed.gz

If a file is missing for a genome, manually create it as an empty file if necessary, e.g.:
```sh
cat /dev/null | gzip -c > metadata/\<genome\>/ENCODE/\<genome\>-blacklist.v2.bed.gz
```

### Settings options

Set option `--genome` to the name of the composite genome to create.
The name "\<genome1\>_\<genome2\>" is recommended for clarity.

Option `--genomes` should be set as two colon-delimited UCSC genome names, e.g., `hg38:sacCer3`.

Option `--ordered-chroms` should be set as two comma-delimited ordered chromosomes lists
separated by a colon, e.g., `chr1,chr2:chrI,chrII`, in the same order as `--genomes`.

### Output

The new genome is created as a custom genome in folder
`<GENOMES_DIR>/custom/<GENOME>`. 

All chromosome names in all output files will be suffixed
with `_<genome>`, e.g., input chromosome `chr1` from genome `hg38`
will become `chr1_hg38`.

A metadata file with information about the genome will
be created as 
`<GENOMES_DIR>/custom/<GENOME>/<GENOME>.yml`.
You may wish to make manual edits to this metadata 
as appropriate for your genome and/or code that will use it.
