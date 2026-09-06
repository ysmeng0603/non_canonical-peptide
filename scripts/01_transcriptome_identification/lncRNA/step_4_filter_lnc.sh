#!/bin/bash
set -euo pipefail  

echo "============================================================"
echo "# Step 4: Extract candidate lncRNA transcripts and all exons from the gffcompare results"
echo "start：$(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"
# ================ path ================
annotated="./lncRNA/StringTie/result/step_3_gffcompare/compare.annotated.gtf"
outdir="./lncRNA/StringTie/result/step_4_filter"
mkdir -p "$outdir"

# # Check the annotated file
if [[ ! -s "$annotated" ]]; then
    echo "ERROR: compare.annotated.gtf does not exist or is empty: $annotated"
    exit 1
fi
# ================ Step 4.1 ================
echo "[1/3] Extracting transcripts with class_code u, i, x, o, and j..."
awk '$3 == "transcript" && /class_code "[uixoj]"/' "$annotated" \
    > "$outdir/lncRNA_candidates.gtf"

cand="$outdir/lncRNA_candidates.gtf"

if [[ ! -s "$cand" ]]; then
    echo "ERROR: No candidate lncRNA transcripts were found! Please check the gffcompare results."
    exit 1
fi

echo "→ Candidate transcript count: $(wc -l < "$cand")"

# ================ Step 4.2 ================
echo "[2/3] Extracting all exons for these transcripts..."

awk '
    ARGIND==1 {
        if ($3 == "transcript") {
            match($0, /transcript_id "([^"]+)"/, a)
            if (a[1] != "") ids[a[1]]=1
        }
        next
    }
    {
        match($0, /transcript_id "([^"]+)"/, a)
        tid=a[1]
        if (tid in ids) print $0
    }
' "$cand" "$annotated" > "$outdir/lncRNA_candidates_with_exons.gtf"

with_exons="$outdir/lncRNA_candidates_with_exons.gtf"

if [[ ! -s "$with_exons" ]]; then
    echo "ERROR: Failed to extract exons. Please check whether transcript_id was parsed correctly."
    exit 1
fi

echo "→ Total number of transcript + exon lines: $(wc -l < "$with_exons")"
# ================ Step 4.3 ================
echo "[3/3] Generating the transcript_id list..."

grep 'transcript_id' "$cand" \
    | sed -E 's/.*transcript_id "([^"]+)".*/\1/' \
    | sort -u > "$outdir/lncRNA_transcript_ids.txt"

echo "→ transcript_id count: $(wc -l < "$outdir/lncRNA_transcript_ids.txt")"

echo "============================================================"
echo "All steps completed! End time: $(date '+%Y-%m-%d %H:%M:%S')"
echo "Output directory: $outdir"
echo "   ├── lncRNA_candidates.gtf"
echo "   ├── lncRNA_candidates_with_exons.gtf"
echo "   └── lncRNA_transcript_ids.txt"
echo "============================================================"
echo "Finished: $(date '+%Y-%m-%d %H:%M:%S')"