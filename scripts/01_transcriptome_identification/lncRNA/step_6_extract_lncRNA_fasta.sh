#!/bin/bash
set -euo pipefail

echo "============================================================"
echo "Step 6: Extract transcript sequences (FASTA) from the final lncRNA GTF"
echo "Start time: $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"
# ref genome fasta
genome_fa="./reference/GRCh38.primary_assembly.genome.fa"

final_gtf="./lncRNA/StringTie/result/step_5_len_exon_filter/final_lncRNA_with_exons.gtf"
outdir="./lncRNA/StringTie/result/step_6_fasta"
mkdir -p "$outdir"

lnc_fa="$outdir/final_lncRNA_transcripts.fa"

echo "[1/2] Extracting transcript sequences using gffread..."
gffread "$final_gtf" -g "$genome_fa" -w "$lnc_fa"

echo "→ lncRNA FASTA sequences saved to: $lnc_fa"
echo "  Sequence count: $(grep -c '^>' "$lnc_fa")"

echo "============================================================"
echo "Step 6 completed!"
echo "End time: $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"
