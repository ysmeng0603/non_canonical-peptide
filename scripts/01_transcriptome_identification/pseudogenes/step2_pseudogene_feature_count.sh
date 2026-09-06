#!/bin/bash
  featureCounts \
    -T 8 \
    -p \
    --countReadPairs \
    -s 0 \
    -a annotation/gencode.v49.pseudogene.clean.gtf \
    -o pseudogene_featureCounts.txt \
    -t exon \
    -g gene_id \
    *.bam