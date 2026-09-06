#!/bin/bash
# Activate the corresponding environment before running
cd ./circRNA_peptide3/part1_coding_prob/IRES_finder/merge_fasta

python ./software/IRESfinder/IRESfinder.py -f ./H_eRA_down_merge.fasta -o ../Result/H_eRA_down_IRES.txt
python ./software/IRESfinder/IRESfinder.py -f ./H_eRA_up_merge.fasta -o ../Result/H_eRA_up_IRES.txt
python ./software/IRESfinder/IRESfinder.py -f ./H_RA_down_merge.fasta -o ../Result/H_RA_down_IRES.txt
python ./software/IRESfinder/IRESfinder.py -f ./H_RA_up_merge.fasta -o ../Result/H_RA_up_IRES.txt
python ./software/IRESfinder/IRESfinder.py -f ./RA_eRA_down_merge.fasta -o ../Result/RA_eRA_down_IRES.txt