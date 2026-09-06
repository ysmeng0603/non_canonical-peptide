#!/bin/bash
set -euo pipefail  

echo "=============================================="
echo "run rMATS-turbo v4.3.0 "
echo "time:$(date '+%Y-%m-%d %H:%M:%S')"
echo "=============================================="
## ===== paths =====
RMATS_HOME="./software/rMATS.4.0.1/rMATS-turbo-Linux-UCS4"
GTF="./reference/gencode.v44.primary_assembly.annotation.gtf"

hc_bam_path="./AS_rMATS_jcast/rMATS/need_file/HC_bamPath.txt"
disease_bam_path="./AS_rMATS_jcast/rMATS/need_file/disease_bamPath.txt"

out_dir_1="./AS_rMATS_jcast/rMATS/result"
mkdir -p "$out_dir_1"

ulimit -n 6666

cd "$RMATS_HOME"

python rmats.py \
  --b1 "$hc_bam_path"\
  --b2 "$disease_bam_path" \
  --gtf "$GTF" \
  -t paired \
  --readLength 101 \
  --nthread 20 \
  --tstat 4 \
  --cstat 0.0001 \
  --libType fr-unstranded \
  --od "$out_dir_1"

echo "=============================================="
echo "finish: $(date '+%Y-%m-%d %H:%M:%S')"
echo "=============================================="
