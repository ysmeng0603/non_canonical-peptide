#!/bin/bash
set -euo pipefail   #

echo "=============================================="
echo "run jcast "
echo "begin: $(date '+%Y-%m-%d %H:%M:%S')"
echo "=============================================="

##############################################################################################
# rmats dir
rmats_dir="./circpeptide4/AS_rMATS_jcast/rMATS/result/"
#
gtf="./reference/gencode.v44.primary_assembly.annotation.gtf"
fasta="./reference/GRCh38.primary_assembly.genome.fa"
# 
out_dir="./circpeptide4/AS_rMATS_jcast/jcast/result/"

jcast \
  "$rmats_dir" \
  "$gtf" \
  "$fasta" \
  -o "$out_dir" \
  -r 1 \
  -q 0 0.05 \
  -m \
  -c
# qvalue


echo "end:$(date '+%Y-%m-%d %H:%M:%S')

