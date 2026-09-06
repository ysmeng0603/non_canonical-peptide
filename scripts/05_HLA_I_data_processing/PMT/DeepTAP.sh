#!/bin/bash

# ============================================================
# DeepTAP prediction
# ============================================================

# Input file
input_file="./MS_unspecific/circRNA/result/PMT/TAP/for_TAP.csv"

# Output directory
out_dir="./MS_unspecific/circRNA/result/PMT/TAP/result"

# DeepTAP installation directory
cd ./software/DeepTAP

# Classification model prediction
python deeptap.py -t cla -f "$input_file" -o "$out_dir"

# Regression model prediction
python deeptap.py -t reg -f "$input_file" -o "$out_dir"