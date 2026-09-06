#!/bin/bash

# ============================================================
# Run netMHCpan-4.1 predictions for MHC class I
# ============================================================

parameter="MHC1_file_paths.txt"
num_tasks=20

output_dir="./circRNA_peptide3/part4_netMHCpan/MHC_1/netMHCpan_1_result"

netmhcpan="./softs/netMHCpan-4.1/netMHCpan"
rdir="./softs/netMHCpan-4.1/Linux_x86_64"

synlist="$rdir/data/synlist.bin"
pseudo="$rdir/data/MHC_pseudo.dat"


# Count tasks excluding the header
num_lines=$(( $(wc -l < "$parameter") - 1 ))

# Calculate the number of batches
num_batches=$(( (num_lines + num_tasks - 1) / num_tasks ))


for ((i = 0; i < num_batches; i++)); do

    echo "[Batch $((i + 1))/$num_batches]"

    start_line=$((i * num_tasks + 2))
    end_line=$((i * num_tasks + num_tasks + 1))

    temp_file="temp_${parameter}"

    # Extract tasks for the current batch
    sed -n "${start_line},${end_line}p" "$parameter" > "$temp_file"

    while IFS=$'\t' read -r task_index task_file_path task_file_name task_HLA_1; do
        {
            output_file="${output_dir}/${task_file_name}.txt"

            echo "[START] ${task_index} | ${task_file_name}"

            # Skip completed tasks
            if [[ -f "$output_file" ]]; then

                echo "[SKIP] ${task_file_name} already completed"

            else

                "$netmhcpan" \
                    -rdir "$rdir" \
                    -syn "$synlist" \
                    -hlapseudo "$pseudo" \
                    -a "$task_HLA_1" \
                    -BA \
                    -f "$task_file_path" \
                    -xls \
                    -xlsfile "$output_file"

            fi

            echo "[DONE] ${task_index} | ${task_file_name}"

        } &

    done < "$temp_file"

    # Wait for the current batch to finish
    wait

    rm -f "$temp_file"

done

echo "All netMHCpan-4.1 tasks completed."