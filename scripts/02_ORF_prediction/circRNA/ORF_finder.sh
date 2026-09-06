#!/bin/bash

# ============================================================
# Run ORFfinder for multiple groups
# ============================================================

ORFfinder="./software/ORF_Finder/ORFfinder"


# Function for running ORFfinder
run_orffinder() {

    input_dir="$1"
    output_dir="$2"
    p_output_dir="$3"

    # Process each FASTA file
    for file in "$input_dir"/*.fasta; do

        # Extract filename
        filename=$(basename -- "$file")

        # Set output file paths
        output_file="$output_dir/result_$filename"
        p_output_file="$p_output_dir/p_result_$filename"

        # Run ORFfinder: nucleotide output
        "$ORFfinder" \
            -in "$file" \
            -s 0 \
            -ml 30 \
            -out "$output_file" \
            -outfmt 1

        # Run ORFfinder: protein output
        "$ORFfinder" \
            -in "$file" \
            -s 0 \
            -ml 30 \
            -out "$p_output_file" \
            -outfmt 0

        echo "Processed: $filename"

    done
}


# 1. H_eRA_down
run_orffinder \
"./circRNA_peptide3/part2_BSJ/seq_repeat_4/H_eRA_down/data_seq4/repeat4" \
"./circRNA_peptide3/part2_BSJ/seq_repeat_4/H_eRA_down/result_H_eRA_down" \
"./circRNA_peptide3/part2_BSJ/seq_repeat_4/H_eRA_down/p_result_H_eRA_down"


# 2. H_eRA_up
run_orffinder \
"./circRNA_peptide3/part2_BSJ/seq_repeat_4/H_eRA_up/data_seq4/repeat4" \
"./circRNA_peptide3/part2_BSJ/seq_repeat_4/H_eRA_up/result_H_eRA_up" \
"./circRNA_peptide3/part2_BSJ/seq_repeat_4/H_eRA_up/p_result_H_eRA_up"


# 3. H_RA_down
run_orffinder \
"./circRNA_peptide3/part2_BSJ/seq_repeat_4/H_RA_down/data_seq4/repeat4" \
"./circRNA_peptide3/part2_BSJ/seq_repeat_4/H_RA_down/result_H_RA_down" \
"./circRNA_peptide3/part2_BSJ/seq_repeat_4/H_RA_down/p_result_H_RA_down"


# 4. H_RA_up
run_orffinder \
"./circRNA_peptide3/part2_BSJ/seq_repeat_4/H_RA_up/data_seq4/repeat4" \
"./circRNA_peptide3/part2_BSJ/seq_repeat_4/H_RA_up/result_H_RA_up" \
"./circRNA_peptide3/part2_BSJ/seq_repeat_4/H_RA_up/p_result_H_RA_up"


# 5. RA_eRA_down
run_orffinder \
"./circRNA_peptide3/part2_BSJ/seq_repeat_4/RA_eRA_down/data_seq4/repeat4" \
"./circRNA_peptide3/part2_BSJ/seq_repeat_4/RA_eRA_down/result_RA_eRA_down" \
"./circRNA_peptide3/part2_BSJ/seq_repeat_4/RA_eRA_down/p_result_RA_eRA_down"