---
title: download
parent: Stage 1 Pipelines
has_children: true
nav_order: 1
---

## Resource downloader

'download' provides utilities that make it easy to obtain
sequence data files from common external sources such as
Illumina BaseSpace Sequence Hub and public data repositories.

In general, the action name is the name of the external
data source, e.g., 'mdi download basespace'.

Please see subsections / subfolders to learn how to configure
a job for each action type.

### Avoid IO bottlenecks, connection limits

In general, you should write data.yml to download one
file at a time; it is usually a bad idea to have multiple
processes downloading multiple files at once and doesn't usually
speed things up anyway when bandwidth is limiting. Also,
some sources might limit the number of concurrent connections
from any one client.

### A note of caution about drive space

This tool makes it easy to download
a LOT of data, make sure you have drive space to store
what you are downloading! Some commands have features to
check the total file size before downloading.
