---
title: Metadata
parent: download
grand_parent: Stage 1 Pipelines 
has_children: false
nav_order: 4
---

## Genome metadata files

Various other genome metadata files are useful for genomic analysis 
pipelines. Action 'metadata' retrieves many of these files in a single
step in a path suitable for other pipelines such as 'prepareBins'.

For metadata, option '--data-name' must be set to the
UCSC format name of the reference genome, e.g., 'hg38'.

Option '--output-dir' must be set to 'GENOMES_DIR/download/metadata',
where GENOMES_DIR is the value you will subsequently
provide to other pipelines. If you would like to use the default genomes 
directory, it can be entered as:

```yml
optionFamilies:
    output:
        output-dir: $MDI_DIR/resources/genomes/$PIPELINE_NAME/$PIPELINE_ACTION
        data-name:  bosTau8
```

### Metadata file sources

File 'gap.txt.gz' is a UCSC file that lists the known gap
regions in a genome assembly, used to correct coverage expectations. 

- <https://hgdownload.cse.ucsc.edu/goldenPath/GENOME/database/>

File 'GENOME.gc5Base.wigVarStep.gz' is a UCSC file that
tabulates the GC content of the reference genome, used
during coverage normalization. 

- <http://hgdownload.cse.ucsc.edu/goldenpath/GENOME/bigZips/>

File 'GENOME-blacklist.v2.bed.gz' is an ENCODE file that 
lists known low quality regions of the genome that should be 
masked and ignored. 

- <https://pubmed.ncbi.nlm.nih.gov/31249361/>
- <https://github.com/Boyle-Lab/Blacklist/>

Please note that as of this version of genomex-mdi-tools/download/metadata
and beyond, umap mappability files are no longer supported. See further
notes in genomex-mdi-tools/prepareBins.
