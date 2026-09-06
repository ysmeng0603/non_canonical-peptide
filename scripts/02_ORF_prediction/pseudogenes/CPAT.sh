#!/bin/bash

# ============================================================
# Run CPAT for pseudogene transcript sequences
# ============================================================

model="./software/CPAT/Human_logitModel.RData"
hexamer="./software/CPAT/Human_Hexamer.tsv"

input_dir="./circpeptide4/pseudogenes/result/step4_seq/DE_pseudogene_fasta/by_group"
output_base="./circpeptide4/pseudogenes/result/step5_coding/result/CPAT"

groups=(
    "H_eRA HC_eRA"
    "H_RA HC_RA"
    "RA_eRA RA_eRA"
)

for item in "${groups[@]}"; do

    read -r group output_group <<< "$item"

    cpat \
        -g "$input_dir/${group}.pseudogene.transcripts.fasta" \
        -d "$model" \
        -x "$hexamer" \
        -o "$output_base/${output_group}/CPAT_${group}_output.txt"

    echo "Processed: $group"

done