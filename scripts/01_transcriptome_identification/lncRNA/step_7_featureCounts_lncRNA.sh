#!/bin/bash
set -euo pipefail

echo "============================================================"
echo "Step 7: Quantify known lncRNAs using featureCounts"
echo "Start time: $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"

# 1) 
lnc_known_gtf="./reference/gencode.v49.long_noncoding_RNAs.gtf"
# 2) BAM list
bam_list="./lncRNA/StringTie/GSE89408_bam_path.txt"
# 3) 
outdir="./lncRNA_known/StringTie/result/step_7_featureCounts"
mkdir -p "$outdir"

counts_out="$outdir/lncRNA_transcript_counts.txt" # 
logfile="$outdir/featureCounts.lncRNA.log"

echo "[1/3] Checking input files..."
[[ -s "$lnc_known_gtf" ]] || { echo "ERROR: GTF file does not exist or is empty"; exit 1; }
[[ -s "$bam_list" ]] || { echo "ERROR: BAM list does not exist or is empty"; exit 1; }
echo "[2/3] Reading BAM file list..."
mapfile -t BAM_FILES < <(cut -f1 "$bam_list")

echo "  BAM count: ${#BAM_FILES[@]}"

echo "[3/3] Running featureCounts (paired-end)..."

featureCounts \
  -T 60 \
  -p -B -C \
  -t transcript \
  -g transcript_id \
  -a "$lnc_known_gtf" \
  -o "$counts_out" \
  "${BAM_FILES[@]}" \
  > "$logfile" 2>&1


echo "→ featureCounts complete"
echo "→ Result file: $counts_out"
echo "→ Log file: $logfile"

echo "============================================================"
echo "Step 7 complete! End time: $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"
