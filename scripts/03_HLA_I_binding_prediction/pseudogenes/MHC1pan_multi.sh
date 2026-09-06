#!/bin/bash

set -euo pipefail

# ============================================================
# Run netMHCpan-4.1 predictions for MHC class I
# Usage: ./MHC1pan_multi.sh MHC1_file_paths_only_up.txt 35
# ============================================================

parameter="$1"
num_tasks="$2"

output_dir="./pseudogenes/result/step8_netMHCpan/result/HLA_1/netMHCpan_1_result"

netmhcpan="./software/netMHCpan-4.1/Linux_x86_64/bin/netMHCpan"
rdir="./software/netMHCpan-4.1/Linux_x86_64"

synlist="$rdir/data/synlist.bin"
pseudo="$rdir/data/MHC_pseudo.dat"

tmp_dir="./pseudogenes/result/step8_netMHCpan/result/HLA_1/tmp"

mkdir -p "$tmp_dir" "$output_dir"


# ============================================================
# Run a single netMHCpan task
# ============================================================

run_netmhcpan() {

    IFS=$'\t' read -r \
        task_index \
        task_file_path \
        task_file_name \
        task_HLA_1 <<< "$1"

    local output_file="${output_dir}/${task_file_name}.xls"
    local running_flag="${output_dir}/${task_file_name}_running"

    echo "[START] ${task_index} | ${task_file_name}"

    # Skip completed tasks
    if [[ -f "$output_file" ]]; then
        echo "[SKIP] ${task_file_name} already completed"
        return 0
    fi

    # Prevent duplicate execution
    if ! (set -o noclobber; > "$running_flag") 2>/dev/null; then
        echo "[SKIP] ${task_file_name} already running"
        return 0
    fi

    # Remove the running flag when the task exits
    trap 'rm -f "$running_flag"' RETURN

    TMPDIR="$tmp_dir" "$netmhcpan" \
        -rdir "$rdir" \
        -syn "$synlist" \
        -hlapseudo "$pseudo" \
        -a "$task_HLA_1" \
        -BA \
        -f "$task_file_path" \
        -xls \
        -xlsfile "$output_file"

    echo "[DONE] ${task_index} | ${task_file_name}"
}


export -f run_netmhcpan

export output_dir
export netmhcpan
export rdir
export synlist
export pseudo
export tmp_dir


# ============================================================
# Run tasks in parallel
# ============================================================

echo "============================================================"
echo "Running netMHCpan-4.1"
echo "Task file: $parameter"
echo "Maximum parallel tasks: $num_tasks"
echo "Output directory: $output_dir"
echo "============================================================"

tail -n +2 "$parameter" |
    xargs -P "$num_tasks" -I {} bash -c 'run_netmhcpan "$1"' _ "{}"

echo "============================================================"
echo "All netMHCpan-4.1 tasks completed"
echo "============================================================"