#!/usr/bin/env bash
set -euo pipefail

if (( $# > 2 )); then
    printf 'Usage: bash %s [approved_ppmi_data_directory] [local_output_directory]\n' "$0" >&2
    exit 2
fi

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
data_dir=${1:-01_Data/PPMI}
out_dir=${2:-03_Result/PPMI}

exec python3 "$script_dir/build_scan_diagnosis.py" \
    --scan-list "$data_dir/PPMI_selected_scans.csv" \
    --research-dx "$data_dir/Primary_Research_Diagnosis_23Nov2026.csv" \
    --clinical-dx "$data_dir/Clinical_Diagnosis_23Nov2026.csv" \
    --out-dir "$out_dir"
