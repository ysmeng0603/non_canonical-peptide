#!/bin/bash

# ============================================================
# Run CPC2 for lncRNA sequences
# ============================================================

CPC2="./software/CPC2-beta/bin/CPC2.py"
input_dir="./circpeptide4/lncRNA_known/StringTie/result/step_10_FC_fasta"
output_dir="./circpeptide4/lncRNA_known/coding/CPAT_CPC2/result/CPC2"

groups=(
    "H_eRA"
    "H_RA"
    "RA_eRA"
)

for group in "${groups[@]}"; do

    python2.7 "$CPC2" \
        -i "$input_dir/${group}.fasta" \
        -o "$output_dir/${group}_output.txt"

    echo "Processed: $group"

done