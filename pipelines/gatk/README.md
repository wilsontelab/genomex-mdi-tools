---
published: false
---

## gatk: Simple wrapper around the Genome Analysis Toolkit

`gatk` is a simple wrapper around the Broad Institute's
Genome Analysis Toolkit. It uses the GATK4, which is
"100% open source and released under the Apache 2.0 license".

GATK4 brings together well-established tools from the GATK and 
Picard codebases under a streamlined framework, to enable selected 
tools to be run in a massively parallel way on local clusters or 
in the cloud using Apache Spark. 

For more information, see:

- https://github.com/broadinstitute/gatk 
- https://gatk.broadinstitute.org/hc/en-us

`genomex-mdi-tools/gatk` isn't really a pipeline, per se, 
it is a way to help install the toolkit using conda and 
run one-off commands. A more complex pipeline can be built 
around gatk (and other) tools by invoking the `genomex-mdi-tools/gatk`
conda module within your own named pipeline.

## Usage

### Open a command shell with a live GATK installation

To use GATK, you might first simply wish to 
create its conda environment:

```
mdi gatk conda --create
```

open a shell into that environment:

```
mdi -d gatk shell
```

and run any command at will:

```bash
# from the environment shell command prompt
$ gatk

 Usage template for all tools (uses --spark-runner LOCAL when used with a Spark tool)
    gatk AnyTool toolArgs

<snip>
```

### Run a single GATK command using the MDI wrapper

Alternatively, you can make use of MDI local environments and the like
by using the gatk pipeline to run a single GATK4 command:

```
mdi gatk do --command MyCommand --arguments "--xx 11 -zz 22"
```
