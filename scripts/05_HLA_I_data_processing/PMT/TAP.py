import pandas as pd

# ============================================================
# Prepare peptide input for DeepTAP
# ============================================================

# Read NetChop results
for_DeepTAP = pd.read_csv(
    "./MS_unspecific/circRNA/result/PMT/netchop/result/filter_result/netChop_result_1.txt",
    sep="\t",
    header=0
)

# Extract unique peptide sequences
for_DeepTAP = (
    for_DeepTAP["Peptide"]
    .drop_duplicates()
    .to_frame(name="peptide")
)

# Export DeepTAP input
for_DeepTAP.to_csv(
    "./circRNA_peptide3/A_circ_2026/supple_2/Radar_plot/need_file/TAP/need_file/for_TAP.csv",
    index=False
)