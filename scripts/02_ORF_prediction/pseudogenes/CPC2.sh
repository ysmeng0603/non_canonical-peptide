#!/bin/bash

# ============================================================
# Run CPC2 for pseudogene transcript sequences
# ============================================================

CPC2="./software/CPC2-beta/bin/CPC2.py"

input_dir="./circpeptide4/pseudogenes/result/step4_seq/DE_pseudogene_fasta/by_group"
output_dir="./circpeptide4/pseudogenes/result/step5_coding/result/CPC2"

groups=(
    "H_eRA"
    "H_RA"
    "RA_eRA"
)

for group in "${groups[@]}"; do

    python2.7 "$CPC2" \
        -i "$input_dir/${group}.pseudogene.transcripts.fasta" \
        -o "$output_dir/${group}_output.txt"

    echo "Processed: $group"

done