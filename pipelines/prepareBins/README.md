---
published: false
---

## prepareBins: Pre-Assemble Fixed-Width Genome Bins

'prepareBins' uses a series of genome input files to
pre-assemble a set of fixed-width (i.e., equal size) bins
(a.k.a. windows) along a reference genome. The content of
each been is annotated with values for gap content, bad
genome regions per ENCODE, mappability, and GC content.

The resulting output files can be used as inputs to
other pipelines that count reads or other features within
such bins.

'prepareBins' is independent of any sample-specific data.

### Specifying Genome Name and Directory

For prepareBins, option 'data-name' must be set to the
UCSC format name of the reference genome, e.g., 'hg38'.

Option 'output-dir' must be set to '\<genomes-dir\>/bins',
where \<genomes-dir\> is the value you will subsequently
provide to other pipelines and thus should already have
other sub-directories such as iGenomes, etc. If you would 
like to use the default genomes directory, it can be entered as:

```
optionFamilies:
    output:
        output-dir: $MDI_DIR/resources/genomes/bins 
        data-name:  hg38
```

### Required Genome Inputs

Required input files provide information about the genome
in use and must include (remember that 'data-name' is
synonymous with 'genome', per above):

```
<output-dir>/../metadata/<genome>/UCSC/gap.txt.gz
<output-dir>/../metadata/<genome>/UCSC/<genome>.gc5Base.wigVarStep.gz
<output-dir>/../metadata/<genome>/ENCODE/<genome>-blacklist.v2.bed.gz    
<output-dir>/../genmap/<genome>/maps/<genome>.mappability.k_<kmer-length>.e_<n-errors>.bedgraph
```

The genmap mappability file can be obtained by running the 'genmap' 
pipeline in this pipeline suite.

The various other metadata files can be downloaded or created manually,
or, better yet, using the 'metadata' action of the 'download' pipeline.

### Output file format

The output is a gzipped extended BED format file named:

```
${OUTPUT_DIR}/${DATA_NAME}/fixed_width_bins/${DATA_NAME}.bins.size_${BIN_SIZE}.k_${KMER_LENGTH}.e_${N_ERRORS}.bed.gz
```

with columns:

```bash
chrom     # chromosom
start     # 0-referenced, per BED specifications
end       # 1-referenced, per BED specifications
binN      # from 1 to N
gc        # fraction GC bases
strand    # unused, always "."
excluded  # number of excluded bases
gap       # number of gap bases
badRegion # number of bases in bad genome regions
umap      # unused, always 1
genmap    # fraction mappability, as per genmap
```

### Note regarding umap mappability files

umap mappability files are no longer supported in this version of
genomex-mdi-tools/prepareBins and beyond. Please use genmap mappability instead.

For file format consistency, the umap column is still present in the
bins file, but it always has a value of 1. Thus, it will not break
any downstream code that uses it, but output won't be adjusted
for mappability as you would expect.
