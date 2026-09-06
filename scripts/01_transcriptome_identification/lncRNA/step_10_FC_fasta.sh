#!/bin/bash
set -euo pipefail

echo "============================================================"
echo "Step 10: Extract FASTA sequences of differentially expressed lncRNAs"
echo "Start time: $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"

GTF=./lncRNA/StringTie/result/step_5_len_exon_filter/final_lncRNA_with_exons.gtf
GENOME=./lncRNA/StringTie/result/step_5_len_exon_filter/hg38.fa
OUT=./lncRNA/StringTie/result/step_10_FC_fasta
mkdir -p $OUT
# 
gffread \
  -w $OUT/lncRNA_transcripts.fa \
  -g $GENOME \
  $GTF
grep -c "^>" $OUT/lncRNA_transcripts.fa

echo "============================================================"
echo "Step 10 completed! End time: $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"