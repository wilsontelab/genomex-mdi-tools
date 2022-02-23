---
title: BaseSpace
parent: download
grand_parent: Stage 1 Pipelines 
has_children: false
nav_order: 1
---

## Illumina BaseSpace

To download data from BaseSpace, follow the following general
steps. 

1) Obtain an Illumina BaseSpace account.

2) Access the developer dashboard (log in as needed):
https://developer.basespace.illumina.com/dashboard

3) Click 'Create a New Application' (if you don't already have one)

> Select "Web/Desktop - Other" as the application type.
Other values are up to you (e.g., the name; an example might be
"My File Download").

4) In the developer dashboard, select your new app.

5) Navigate to the "Credentials" tab and copy your Access Token -
you need to provide this value to
'mdi download basespace --access-token'.

6) In data.yml, use value 'NA' for 'project-id' or 'sample-names' to obtain lists
of projects, and then samples.  Then use those values for 'project-id'
and 'sample-names', where the latter is a comma-delimited list.

The pipeline will download every file for each requested sample, one
file at a time, into $TASK_DIR/$SAMPLE_NAME. Typically, the files you want and will get will be fastq.gz files.

For additional help, see:

- <https://developer.basespace.illumina.com/docs/content/documentation/faq/developer-faq>
- <https://developer.basespace.illumina.com/docs/content/documentation/rest-api/api-reference>
- <https://bioconductor.org/packages/release/bioc/vignettes/BaseSpaceR/inst/doc/BaseSpaceR.pdf>
- <https://bioconductor.org/packages/release/bioc/manuals/BaseSpaceR/man/BaseSpaceR.pdf>
