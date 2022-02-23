---
title: SRA
parent: download
grand_parent: Stage 1 Pipelines 
has_children: false
nav_order: 2
---

## Sequence Read Archive (SRA)

SRA files are publicly available and require no other information
to obtain than the target identifier. The name can get confusing,
but what you want to identify and provide to 'mdi download sra'
is a Run identifier(s), which starts with 'SRR'.

An experiment is a library type applied to a specific biosample.
An experiment might have multiple sequencing runs, and you will
typically want to retrieve them all, because they all represent
data from the same sequenced library.

As an example, exome sequencing applied to NA12878 is an Experiment
that might have been run twice on an Illumina sequencer, yielding
two Run files. You want to provide the Run ids, not the Experiment.

You will be asked to provide you own name for the experiment (you
can use the SRA/SRX name, or change it).

The pipeline will download every requested SRR file,
one file at a time, into $TASK_DIR/$EXPERIMENT_NAME.

### sra vs. fastq format files

The 'mdi download sra' command use 'prefetch' from the SRA Toolkit
to download the data  as .sra files.  Note: the pipeline does
NOT convert those files to fastq!  That is a waste of time and disk
space because you can simply convert to fastq in a stream when
you go to align the data. It takes much less disk space to store the
files in .sra format, and does not slow alignment down because the
alignment step (e.g., bwa) is slower than fastq conversion.

There is information in Illumina read names that is lost in the .sra
files you get from SRA; this information was permanently discarded by SRA
and cannot be recovered by converting to FASTQ or by any other means.
