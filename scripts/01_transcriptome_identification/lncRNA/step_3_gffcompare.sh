#!/bin/bash
#
echo "start:$(date '+%Y-%m-%d %H:%M:%S')"

# 
merged_gtf="./lncRNA/StringTie/result/step_2_merge/merged.gtf"
ref_gtf="./reference/Homo_sapiens.GRCh38.113.gtf"

#
outdir="./lncRNA/StringTie/result/step_3_gffcompare"
mkdir -p "$outdir"

# 运行 gffcompare
gffcompare \
    -r "$ref_gtf" \
    -o "$outdir/compare" \
    "$merged_gtf"

echo "=== gffcompare finish ==="
echo "ourdir:$outdir"
echo "finish:$(date '+%Y-%m-%d %H:%M:%S')"