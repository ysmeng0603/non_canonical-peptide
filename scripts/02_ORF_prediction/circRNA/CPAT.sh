#!/bin/bash

# ============================================================
# Run CPAT for multiple groups
# ============================================================

# CPAT model files
model="./software/CPAT/Human_logitModel.RData"
hexamer="./software/CPAT/Human_Hexamer.tsv"


# Function for running CPAT
run_cpat() {

    input_file="$1"
    input_dir="$2"
    output_dir="$3"

    while IFS= read -r id; do

        cpat \
            -g "$input_dir/$id.fasta" \
            -d "$model" \
            -x "$hexamer" \
            -o "$output_dir/${id}_CPAT.txt"

        echo "Processed: $id"

    done < "$input_file"
}


# 1--- H_eRA_down
run_cpat \
"./circRNA_peptide3/part1_coding_prob/z_bash_need/ID/H_eRA_down_ID.txt" \
"./circRNA_peptide3/get_sequence/sequence_file/H_eRA_down/fasta" \
"./circRNA_peptide3/part1_coding_prob/CPAT/H_eRA_down"


# 2--- H_eRA_up
run_cpat \
"./circRNA_peptide3/part1_coding_prob/z_bash_need/ID/H_eRA_up_ID.txt" \
"./circRNA_peptide3/get_sequence/sequence_file/H_eRA_up/fasta" \
"./circRNA_peptide3/part1_coding_prob/CPAT/H_eRA_up"


# 3--- H_RA_down
run_cpat \
"./circRNA_peptide3/part1_coding_prob/z_bash_need/ID/H_RA_down_ID.txt" \
"./circRNA_peptide3/get_sequence/sequence_file/H_RA_down/fasta" \
"./circRNA_peptide3/part1_coding_prob/CPAT/H_RA_down"


# 4--- H_RA_up
run_cpat \
"./circRNA_peptide3/part1_coding_prob/z_bash_need/ID/H_RA_up_ID.txt" \
"./circRNA_peptide3/get_sequence/sequence_file/H_RA_up/fasta" \
"./circRNA_peptide3/part1_coding_prob/CPAT/H_RA_up"


# 5--- RA_eRA
run_cpat \
"./circRNA_peptide3/part1_coding_prob/z_bash_need/ID/RA_eRA_ID.txt" \
"./circRNA_peptide3/get_sequence/sequence_file/RA_eRA_down/fasta" \
"./circRNA_peptide3/part1_coding_prob/CPAT/RA_eRA_down"