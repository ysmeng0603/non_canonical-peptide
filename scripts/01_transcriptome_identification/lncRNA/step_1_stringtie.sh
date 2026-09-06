#!/bin/bash
info_file="./lncRNA/StringTie/GSE89408_bam_path.txt"
#
ref_gtf="./reference/Homo_sapiens.GRCh38.113.gtf"
outdir="./lncRNA/StringTie/result/step_1_assembly"
# 
echo "start:$(date '+%Y-%m-%d %H:%M:%S')"
while IFS=$'\t' read -r bam_file sampleID; do
  if [[ -f "$bam_file" ]]; then
    echo "Processing $bam_file"
    stringtie "$bam_file" \
        -G "$ref_gtf" \
        -o "$outdir/${sampleID}.gtf" \
        -p 50 # 
    # 
  else
    echo "File not found: $bam_file"
  fi
done < "$info_file"

echo "finish:$(date '+%Y-%m-%d %H:%M:%S')

