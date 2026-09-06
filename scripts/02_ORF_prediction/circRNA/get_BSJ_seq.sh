#!/bin/bash

# ============================================================
# Extract BSJ sequences
# Take the last 10 nt + first 10 nt of each circRNA sequence
# ============================================================

# Function for extracting BSJ sequences
extract_bsj() {

    input_folder="$1"
    output_folder="$2"

    # Process each FASTA file
    for file in "$input_folder"/*.fasta; do

        # Extract filename without path
        filename=$(basename -- "$file")

        # Extract filename without extension
        filename_no_ext="${filename%.*}"

        # Set output file path
        output_file="$output_folder/BSJseq_$filename"

        # Extract the first 10 nucleotides
        head -n 2 "$file" | tail -n 1 | cut -c 1-10 > "first_10_$filename.txt"

        # Extract the last 10 nucleotides
        tail -n 1 "$file" | rev | cut -c 1-10 | rev > "last_10_$filename.txt"

        # Concatenate the last 10 nt and first 10 nt
        paste -d "" \
            "last_10_$filename.txt" \
            "first_10_$filename.txt" \
            > "$output_file"

        # Remove intermediate files
        rm "first_10_$filename.txt" "last_10_$filename.txt"

        echo "Processed: $filename"

    done
}


# 1. H_eRA_down
extract_bsj \
"./circRNA_peptide3/part2_BSJ/seq_repeat_4/H_eRA_down/data_seq4/select_raw" \
"./circRNA_peptide3/part2_BSJ/seq_repeat_4/H_eRA_down/BSJ_files"


# 2. H_eRA_up
extract_bsj \
"./circRNA_peptide3/part2_BSJ/seq_repeat_4/H_eRA_up/data_seq4/select_raw" \
"./circRNA_peptide3/part2_BSJ/seq_repeat_4/H_eRA_up/BSJ_files"


# 3. H_RA_down
extract_bsj \
"./circRNA_peptide3/part2_BSJ/seq_repeat_4/H_RA_down/data_seq4/select_raw" \
"./circRNA_peptide3/part2_BSJ/seq_repeat_4/H_RA_down/BSJ_files"


# 4. H_RA_up
extract_bsj \
"./circRNA_peptide3/part2_BSJ/seq_repeat_4/H_RA_up/data_seq4/select_raw" \
"./circRNA_peptide3/part2_BSJ/seq_repeat_4/H_RA_up/BSJ_files"


# 5. RA_eRA_down
extract_bsj \
"./circRNA_peptide3/part2_BSJ/seq_repeat_4/RA_eRA_down/data_seq4/select_raw" \
"./circRNA_peptide3/part2_BSJ/seq_repeat_4/RA_eRA_down/BSJ_files"