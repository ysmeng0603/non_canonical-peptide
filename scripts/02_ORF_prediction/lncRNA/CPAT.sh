#!/bin/bash

# ============================================================
# Run CPAT for lncRNA sequences
# ============================================================

model="./software/CPAT/Human_logitModel.RData"
hexamer="./software/CPAT/Human_Hexamer.tsv"

input_dir="./circpeptide4/lncRNA_known/StringTie/result/step_10_FC_fasta"
output_base="./circpeptide4/lncRNA_known/coding/CPAT_CPC2/result/CPAT"

# Group name, input FASTA name, and output subdirectory
groups=(
    "H_eRA HC_eRA"
    "H_RA HC_RA"
    "RA_eRA RA_eRA"
)

for item in "${groups[@]}"; do

    read -r group output_group <<< "$item"

    cpat \
        -g "$input_dir/${group}.fasta" \
        -d "$model" \
        -x "$hexamer" \
        -o "$output_base/${output_group}/CPAT_${group}_output.txt"

    echo "Processed: $group"

done