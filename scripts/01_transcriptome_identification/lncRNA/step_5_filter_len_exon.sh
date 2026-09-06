#!/bin/bash
set -euo pipefail

echo "============================================================"
echo "Step 5: Filter standard lncRNAs by length >=200 nt and exon count >=2"
echo "Start time: $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"
# 
gtf_with_exons="./circpeptide4/lncRNA/StringTie/result/step_4_filter/lncRNA_candidates_with_exons.gtf"
# 
outdir="./circpeptide4/lncRNA/StringTie/result/step_5_len_exon_filter"
mkdir -p "$outdir"

stats_file="$outdir/lncRNA_len_exon_stats.txt"
pass_id="$outdir/final_lncRNA_ids.txt"
final_gtf="$outdir/final_lncRNA_with_exons.gtf"

echo "[1/3] Calculating the exon count and total length for each transcript..."
awk '
    $3 == "exon" {
        # 解析 transcript_id
        match($0, /transcript_id "([^"]+)"/, arr)
        tid = arr[1]
        if (tid != "") {
            len = $5 - $4 + 1  
            exon_cnt[tid]++
            exon_len[tid] += len
        }
    }
    END {
        for (t in exon_cnt) {
            print t "\t" exon_cnt[t] "\t" exon_len[t]
        }
    }
' "$gtf_with_exons" > "$stats_file"

echo "→ Statistics completed, number of lines: $(wc -l < "$stats_file")"
echo "   Format: transcript_id  exon_count  total_length"

echo "[2/3] Filtering transcripts by exon count >=2 and total length >=200..."
# Transcripts with exon count > 2 and total length > 200
awk '$2 >= 2 && $3 >= 200 {print $1}' "$stats_file" \
    | sort -u > "$pass_id"
# sort -u = remove duplicates and sort
echo "→ Number of lncRNA transcripts passing the filter: $(wc -l < "$pass_id")"

echo "[3/3] Extracting these transcripts and their exons from the GTF (generating the final GTF)..."

awk '
    NR==FNR { keep[$1]=1; next }
    {
        match($0, /transcript_id "([^"]+)"/, arr)
        tid = arr[1]
        if (tid in keep) print $0
    }
' "$pass_id" "$gtf_with_exons" > "$final_gtf"

echo "→ Number of lines in the final GTF: $(wc -l < "$final_gtf")"

echo "============================================================"
echo "Step 5 completed! Output directory: $outdir"
echo "   ├── lncRNA_len_exon_stats.txt     # Exon count and total length for each transcript"
echo "   ├── final_lncRNA_ids.txt          # lncRNA transcript_ids passing the filter"
echo "   └── final_lncRNA_with_exons.gtf   # Final standard lncRNAs (transcript + exon lines)"
echo "============================================================"
echo "End time: $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"
