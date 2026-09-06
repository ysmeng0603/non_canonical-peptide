#!/bin/bash
cd ./software/TABR-BERT
input_file="./TCR/TABR_BERT/need/input_healthy.csv"
output_dir="./TCR/TABR_BERT/result/healthy/healthy_output.csv"

python ./predict_tcr_pmhc_binding.py \
 --input "$input_file"  --output "$output_dir" \
  --healthy_tcr ./software/TABR-BERT/data/small_healthy_tcr.csv \
  --pseudo_sequence_dict ./software/TABR-BERT/data/mhcflurry.allele_sequences_homo.csv \
  --tcr_pmhc_model ./software/TABR-BERT/model/tcr_pmhc_model.pt \
  --tcr_model ./software/TABR-BERT/model/tcr_model.pt \
  --pmhc_model ./software/TABR-BERT/model/pmhc_model.pt \
  --GPUs 0

# github
# https://github.com/Freshwind-Bioinformatics/TABR-BERT

# Usage: predict_tcr_pmhc_binding.py [options]
#       --input STRING: The data to be predicted (*.csv) 
#                       Required columns: ["allele", "peptide", "cdr3"]
#       --healthy_tcr STRING: TCRs from healthy people for generating negative cases (*.csv)
#                             Required columns: "cdr3" 
#       --pseudo_sequence_dict STRING: allele name to pseudo sequence (*.csv)
#                                      Required columns: ["allele" "sequence"]   
#       --tcr_pmhc_model STRING: TCR-pMHC prediction model dir (*.pt)
#       --tcr_model STRING: TCR embedding model dir (*.pt)
#       --pmhc_model STRING: pMHC embedding model dir (*.pt)                           
#       --output STRING: output file dir (*.csv)

# Optional:
#       --batchsize INT: mini batchsize (default: 256)
#       --embedding_batchsize INT: mini batchsize of generation embedding (default: 256)
#       --pmhc_d_model INT: dimention of pmhc embedding (default: 256)
#       --tcr_d_model INT: dimention of pmhc embedding (default: 256)
#       --GPUs INT: num of GPUs used in this task [if you have GPU recommend 1, if not, recommend 0] (default: 0)