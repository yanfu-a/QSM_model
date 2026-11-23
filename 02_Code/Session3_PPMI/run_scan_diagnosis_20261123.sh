#!/usr/bin/env bash
set -euo pipefail

if (( $# == 3 || $# > 4 )); then
    printf 'Usage: bash %s [data_directory] [output_directory] [research_csv clinical_csv]\n' "$0" >&2
    exit 2
fi

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "$script_dir/../.." && pwd)
data_dir=${1:-$repo_root/01_Data/PPMI}
out_dir=${2:-${PPMI_SCAN_DX_OUTPUT_DIR:-}}
if [[ -z $out_dir ]]; then
    printf 'Specify an approved output directory outside the Git checkout (argument 2 or PPMI_SCAN_DX_OUTPUT_DIR).\n' >&2
    exit 2
fi
research_csv=${3:-$data_dir/Primary_Research_Diagnosis_23Nov2026.csv}
clinical_csv=${4:-$data_dir/Clinical_Diagnosis_23Nov2026.csv}

exec python3 "$script_dir/build_scan_diagnosis.py" \
    --scan-list "$data_dir/PPMI_selected_scans.csv" \
    --research-dx "$research_csv" \
    --clinical-dx "$clinical_csv" \
    --out-dir "$out_dir"
