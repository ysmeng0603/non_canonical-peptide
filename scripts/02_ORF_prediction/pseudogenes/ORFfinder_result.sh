#!/bin/bash

# ============================================================
# Run ORFfinder for lncRNA sequences
# ============================================================

ORFfinder="./software/ORFfinder"

input_dir="./circpeptide4/lncRNA_known/coding/CPAT_CPC2/result/filter/fasta"
rna_output_dir="./circpeptide4/lncRNA_known/coding/ORFfinder/result/fasta/result"
protein_output_dir="./circpeptide4/lncRNA_known/coding/ORFfinder/result/fasta/p_result"

groups=(
    "HC_eRA"
    "HC_RA"
    "RA_eRA"
)

for group in "${groups[@]}"; do

    input_fasta="$input_dir/${group}_ENST.fasta"

    # Nucleotide ORF sequences
    "$ORFfinder" \
        -in "$input_fasta" \
        -s 0 \
        -ml 30 \
        -out "$rna_output_dir/${group}_result.fasta" \
        -outfmt 1

    # Protein ORF sequences
    "$ORFfinder" \
        -in "$input_fasta" \
        -s 0 \
        -ml 30 \
        -out "$protein_output_dir/${group}_result.fasta" \
        -outfmt 0

    echo "Processed: $group"

done