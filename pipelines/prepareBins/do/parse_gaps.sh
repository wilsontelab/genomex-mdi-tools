#!/bin/bash

# parse genome gaps to BED format

cut -f 2,3,4,8 |
awk '{print $0"\t0\t."}' | 
gzip -c
