#!/bin/bash

# ============================================================
# Run MaxQuant for three alternative-splicing comparisons
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

MAXQUANT_DLL="./software/MaxQuant_v2.6.1.0/bin/MaxQuantCmd.dll"

COMBINED_SRC="./MaxQuant_raw_data/PXD037581/combined"

BASE_OUT_DIR="./MS_unspecific/AS/result/contain_jcast/PXD037581"

mkdir -p "$BASE_OUT_DIR"


# ============================================================
# Define MaxQuant tasks
# ============================================================

declare -a TASKS=(
    "HC_eRA|contain_PXD037581_HC_eRA.xml|HC_eRA_result"
    "HC_RA|contain_PXD037581_HC_RA.xml|HC_RA_result"
    "RA_eRA|contain_PXD037581_RA_eRA.xml|RA_eRA_result"
)


# ============================================================
# Run MaxQuant
# ============================================================

success_count=0
total_count=${#TASKS[@]}

for task in "${TASKS[@]}"; do

    IFS="|" read -r GROUP_NAME XML_FILE OUT_SUBDIR <<< "$task"

    XML_PATH="${SCRIPT_DIR}/${XML_FILE}"
    OUT_DIR="${BASE_OUT_DIR}/${OUT_SUBDIR}"

    echo "============================================================"
    echo "Running: $GROUP_NAME"
    echo "XML: $XML_PATH"
    echo "Output: $OUT_DIR"
    echo "============================================================"

    # Skip the task if the XML file is missing
    if [[ ! -f "$XML_PATH" ]]; then
        echo "[ERROR] XML file not found: $XML_PATH"
        continue
    fi

    mkdir -p "$OUT_DIR"

    # Run MaxQuant
    if dotnet "$MAXQUANT_DLL" "$XML_PATH"; then
        echo "[DONE] MaxQuant completed: $GROUP_NAME"
    else
        echo "[ERROR] MaxQuant failed: $GROUP_NAME" >&2
        continue
    fi

    # Copy the MaxQuant combined output
    if cp -r "$COMBINED_SRC" "$OUT_DIR"/; then
        echo "[DONE] Combined output copied"
    else
        echo "[WARNING] Failed to copy combined output" >&2
    fi

    ((success_count++))

done


# ============================================================
# Summary
# ============================================================

echo "============================================================"
echo "MaxQuant batch analysis completed"
echo "Total tasks: $total_count"
echo "Successful tasks: $success_count"
echo "Failed tasks: $((total_count - success_count))"
echo "Output directory: $BASE_OUT_DIR"
echo "============================================================"