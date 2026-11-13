
#!/bin/bash
set -euo pipefail

# ==========================================================
# 1. External arguments
# ==========================================================

if [[ $# -lt 2 ]]; then
    echo "Usage:"
    echo "bash $0 START_INDEX END_INDEX [JOBS]"
    echo
    echo "Example:"
    echo "bash $0 1 100 4"
    exit 1
fi

START_INDEX="$1"
END_INDEX="$2"
JOBS="${3:-4}"

# Validate arguments
for value in "$START_INDEX" "$END_INDEX" "$JOBS"; do
    if ! [[ "$value" =~ ^[0-9]+$ ]]; then
        echo "ERROR: Arguments must be positive integers."
        exit 1
    fi
done

if [[ "$START_INDEX" -lt 1 ]] ||
   [[ "$END_INDEX" -lt "$START_INDEX" ]] ||
   [[ "$JOBS" -lt 1 ]]; then
    echo "ERROR: Invalid start, end or jobs."
    exit 1
fi

# ==========================================================
# 2. Environment
# ==========================================================

module load fsl

export FSLOUTPUTTYPE=NIFTI_GZ
export OMP_NUM_THREADS=1

# ==========================================================
# 3. Path configuration
# ==========================================================

INPUT_ROOT="/cwStorage/nodecw_group/FY_data/PPMI/T1_nifti"

OUTPUT_ROOT="/cwStorage/nodecw_group/FY_data/PPMI/T1_MNI2mm"

MNI_REF="$FSLDIR/data/standard/MNI152_T1_2mm.nii.gz"

MNI_BRAIN="$FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz"

FNIRT_CONFIG="$FSLDIR/etc/flirtsch/T1_2_MNI152_2mm.cnf"

# BET parameter
BET_F="${BET_F:-0.25}"

# Optional sequence filter
SERIES_FILTER="${SERIES_FILTER:-}"

# ==========================================================
# 4. Environment checks
# ==========================================================

for cmd in fslhd bet flirt fnirt applywarp; do
    command -v "$cmd" >/dev/null 2>&1 || {
        echo "ERROR: Cannot find $cmd"
        exit 1
    }
done

[[ -d "$INPUT_ROOT" ]] || {
    echo "ERROR: Input directory does not exist"
    exit 1
}

[[ -f "$MNI_REF" ]] || {
    echo "ERROR: MNI reference does not exist"
    exit 1
}

[[ -f "$MNI_BRAIN" ]] || {
    echo "ERROR: MNI brain reference does not exist"
    exit 1
}

[[ -f "$FNIRT_CONFIG" ]] || {
    echo "ERROR: FNIRT configuration does not exist"
    exit 1
}

mkdir -p "$OUTPUT_ROOT"

RUN_ID=$(date +%Y%m%d_%H%M%S)

STATUS_DIR="$OUTPUT_ROOT/_status_${RUN_ID}"

mkdir -p "$STATUS_DIR"

export INPUT_ROOT OUTPUT_ROOT STATUS_DIR
export MNI_REF MNI_BRAIN FNIRT_CONFIG BET_F

# ==========================================================
# 5. Build sorted NIfTI file list
# Compatible with older Bash versions
# ==========================================================

echo "Building NIfTI file list..."

FILES=()

if [[ -n "$SERIES_FILTER" ]]; then

    while IFS= read -r -d '' file; do
        FILES+=("$file")
    done < <(
        find "$INPUT_ROOT" \
            -type f \
            -name '*.nii.gz' \
            -path "*${SERIES_FILTER}*" \
            -print0 |
        sort -z
    )

else

    while IFS= read -r -d '' file; do
        FILES+=("$file")
    done < <(
        find "$INPUT_ROOT" \
            -type f \
            -name '*.nii.gz' \
            -print0 |
        sort -z
    )

fi

TOTAL_FILES=${#FILES[@]}

if [[ "$TOTAL_FILES" -eq 0 ]]; then
    echo "ERROR: No NIfTI files found."
    exit 1
fi

if [[ "$START_INDEX" -gt "$TOTAL_FILES" ]]; then
    echo "ERROR: START_INDEX exceeds total number of files."
    exit 1
fi

if [[ "$END_INDEX" -gt "$TOTAL_FILES" ]]; then
    echo "ERROR: END_INDEX exceeds total number of files."
    exit 1
fi

SELECTED_COUNT=$((END_INDEX - START_INDEX + 1))

SELECTED_FILES=(
    "${FILES[@]:$((START_INDEX - 1)):$SELECTED_COUNT}"
)

echo "Total NIfTI files: $TOTAL_FILES"
echo "Start index: $START_INDEX"
echo "End index: $END_INDEX"
echo "Selected files: ${#SELECTED_FILES[@]}"

# ==========================================================
# 6. Registration function
# ==========================================================

register_one() {

    local input_file="$1"

    local relative_path="${input_file#"$INPUT_ROOT"/}"

    local relative_dir
    relative_dir=$(dirname "$relative_path")

    local filename
    filename=$(basename "$input_file")

    local base="${filename%.nii.gz}"

    local output_dir="$OUTPUT_ROOT/$relative_dir"

    mkdir -p "$output_dir/transforms"

    local key
    key=$(printf '%s' "$relative_path" |
          md5sum | awk '{print $1}')

    # ------------------------------------------------------
    # Output files
    # ------------------------------------------------------

    local MNI_OUT="$output_dir/${base}_MNI2mm.nii.gz"

    local LINEAR_OUT="$output_dir/${base}_linear_brain_MNI2mm.nii.gz"

    local NATIVE_MASK_OUT="$output_dir/${base}_brain_mask_native.nii.gz"

    local MNI_MASK_OUT="$output_dir/${base}_brain_mask_MNI2mm.nii.gz"

    local AFFINE_OUT="$output_dir/transforms/${base}_to_MNI2mm_affine.mat"

    local WARP_OUT="$output_dir/transforms/${base}_to_MNI2mm_warpcoef.nii.gz"

    local LOG_FILE="$output_dir/${base}_registration.log"

    # ------------------------------------------------------
    # Skip completed files
    # ------------------------------------------------------

    if [[ -s "$MNI_OUT" ]] &&
       [[ -s "$LINEAR_OUT" ]] &&
       [[ -s "$NATIVE_MASK_OUT" ]] &&
       [[ -s "$MNI_MASK_OUT" ]] &&
       [[ -s "$AFFINE_OUT" ]] &&
       [[ -s "$WARP_OUT" ]]; then

        touch "$STATUS_DIR/${key}.skipped"

        echo "[SKIP] $relative_path"

        return 0
    fi

    echo "[START] $relative_path"

    local START_TIME
    START_TIME=$(date +%s)

    # ------------------------------------------------------
    # Temporary directory
    # ------------------------------------------------------

    local WORK_DIR

    WORK_DIR=$(mktemp -d "$output_dir/.work_${base}.XXXXXX")

    local BRAIN_IMG="$WORK_DIR/T1_brain.nii.gz"

    local BRAIN_MASK="$WORK_DIR/T1_brain_mask.nii.gz"

    local AFFINE_TMP="$WORK_DIR/T1_to_MNI_affine.mat"

    local LINEAR_TMP="$WORK_DIR/T1_linear.nii.gz"

    local WARP_TMP="$WORK_DIR/T1_to_MNI_warpcoef.nii.gz"

    local MNI_TMP="$WORK_DIR/T1_MNI2mm.nii.gz"

    local MNI_MASK_TMP="$WORK_DIR/T1_brain_mask_MNI2mm.nii.gz"

    echo "Input: $input_file" > "$LOG_FILE"

    # ------------------------------------------------------
    # Step 1. Header QC
    # ------------------------------------------------------

    local HEADER

    if ! HEADER=$(fslhd "$input_file" 2>>"$LOG_FILE"); then

        touch "$STATUS_DIR/${key}.failed"

        echo "[FAILED HEADER] $relative_path"

        return 0
    fi

    local QFORM_CODE
    local SFORM_CODE
    local DIM0
    local DIM4

    QFORM_CODE=$(awk '$1=="qform_code"{print $2}' <<< "$HEADER")

    SFORM_CODE=$(awk '$1=="sform_code"{print $2}' <<< "$HEADER")

    DIM0=$(awk '$1=="dim0"{print $2}' <<< "$HEADER")

    DIM4=$(awk '$1=="dim4"{print $2}' <<< "$HEADER")

    if [[ "${QFORM_CODE:-0}" -eq 0 ]] &&
       [[ "${SFORM_CODE:-0}" -eq 0 ]]; then

        echo "Invalid qform/sform" >> "$LOG_FILE"

        touch "$STATUS_DIR/${key}.failed"

        echo "[FAILED HEADER] $relative_path"

        return 0
    fi

    if [[ "${DIM0:-0}" -ne 3 ]] ||
       [[ "${DIM4:-0}" -ne 1 ]]; then

        echo "Only 3D images are accepted" >> "$LOG_FILE"

        touch "$STATUS_DIR/${key}.failed"

        echo "[FAILED DIMENSION] $relative_path"

        return 0
    fi

    # ------------------------------------------------------
    # Step 2. Brain extraction
    # Keep original T1 geometry unchanged.
    # ------------------------------------------------------

    if ! bet \
        "$input_file" \
        "$BRAIN_IMG" \
        -R \
        -f "$BET_F" \
        -m >> "$LOG_FILE" 2>&1; then

        touch "$STATUS_DIR/${key}.failed"

        echo "[FAILED BET] $relative_path"

        echo "Work directory: $WORK_DIR"

        return 0
    fi

    if [[ ! -s "$BRAIN_IMG" ]] ||
       [[ ! -s "$BRAIN_MASK" ]]; then

        touch "$STATUS_DIR/${key}.failed"

        echo "[FAILED BET OUTPUT] $relative_path"

        return 0
    fi

    # ------------------------------------------------------
    # Step 3. FLIRT linear registration
    # ------------------------------------------------------

    if ! flirt \
        -in "$BRAIN_IMG" \
        -ref "$MNI_BRAIN" \
        -omat "$AFFINE_TMP" \
        -out "$LINEAR_TMP" \
        -dof 12 \
        -cost corratio >> "$LOG_FILE" 2>&1; then

        touch "$STATUS_DIR/${key}.failed"

        echo "[FAILED FLIRT] $relative_path"

        echo "Work directory: $WORK_DIR"

        return 0
    fi

    # ------------------------------------------------------
    # Step 4. FNIRT nonlinear registration
    #
    # --in: original native T1
    # --aff: FLIRT affine matrix
    # --cout: nonlinear warp coefficients
    # --iout: final MNI T1
    # ------------------------------------------------------

    if ! fnirt \
        --in="$input_file" \
        --aff="$AFFINE_TMP" \
        --ref="$MNI_REF" \
        --config="$FNIRT_CONFIG" \
        --cout="$WARP_TMP" \
        --iout="$MNI_TMP" >> "$LOG_FILE" 2>&1; then

        touch "$STATUS_DIR/${key}.failed"

        echo "[FAILED FNIRT] $relative_path"

        echo "Work directory: $WORK_DIR"

        return 0
    fi

    # ------------------------------------------------------
    # Step 5. Apply nonlinear transformation to brain mask
    # ------------------------------------------------------

    if ! applywarp \
        --in="$BRAIN_MASK" \
        --ref="$MNI_REF" \
        --warp="$WARP_TMP" \
        --out="$MNI_MASK_TMP" \
        --interp=nn >> "$LOG_FILE" 2>&1; then

        touch "$STATUS_DIR/${key}.failed"

        echo "[FAILED APPLYWARP] $relative_path"

        echo "Work directory: $WORK_DIR"

        return 0
    fi

    # ------------------------------------------------------
    # Step 6. Output QC
    # ------------------------------------------------------

    local REF_DIMS
    local OUT_DIMS

    REF_DIMS=$(fslhd "$MNI_REF" |
        awk '$1=="dim1"{a=$2} $1=="dim2"{b=$2} $1=="dim3"{c=$2} END{print a,b,c}')

    OUT_DIMS=$(fslhd "$MNI_TMP" |
        awk '$1=="dim1"{a=$2} $1=="dim2"{b=$2} $1=="dim3"{c=$2} END{print a,b,c}')

    if [[ "$REF_DIMS" != "$OUT_DIMS" ]]; then

        echo "MNI dimensions mismatch" >> "$LOG_FILE"

        touch "$STATUS_DIR/${key}.failed"

        echo "[FAILED DIMENSION QC] $relative_path"

        return 0
    fi

    # ------------------------------------------------------
    # Step 7. Save outputs
    # ------------------------------------------------------

    cp "$BRAIN_MASK" "$NATIVE_MASK_OUT"

    mv "$MNI_TMP" "$MNI_OUT"

    mv "$LINEAR_TMP" "$LINEAR_OUT"

    mv "$MNI_MASK_TMP" "$MNI_MASK_OUT"

    mv "$AFFINE_TMP" "$AFFINE_OUT"

    mv "$WARP_TMP" "$WARP_OUT"

    rm -rf "$WORK_DIR"

    local END_TIME
    END_TIME=$(date +%s)

    local ELAPSED=$((END_TIME - START_TIME))

    echo "Elapsed time: ${ELAPSED} seconds" >> "$LOG_FILE"

    touch "$STATUS_DIR/${key}.success"

    echo "[SUCCESS] $relative_path (${ELAPSED}s)"

}

export -f register_one

# ==========================================================
# 7. Parallel execution
# ==========================================================

printf '%s\0' "${SELECTED_FILES[@]}" |
xargs -0 -r -n 1 -P "$JOBS" \
    bash -c 'register_one "$1"' _

# ==========================================================
# 8. Summary
# ==========================================================

SUCCESS=$(find "$STATUS_DIR" -type f -name '*.success' | wc -l)

FAILED=$(find "$STATUS_DIR" -type f -name '*.failed' | wc -l)

SKIPPED=$(find "$STATUS_DIR" -type f -name '*.skipped' | wc -l)

TOTAL=$((SUCCESS + FAILED + SKIPPED))

echo
echo "=========================================="
echo "Registration completed"
echo "=========================================="

echo "Available files: $TOTAL_FILES"
echo "Selected range: $START_INDEX - $END_INDEX"
echo "Selected files: $SELECTED_COUNT"
echo "Success: $SUCCESS"
echo "Failed: $FAILED"
echo "Skipped: $SKIPPED"

if [[ "$TOTAL" -gt 0 ]]; then
    awk -v s="$SUCCESS" -v t="$TOTAL" \
        'BEGIN {printf "Success rate: %.2f%%\n", s/t*100}'
fi

echo "Output: $OUTPUT_ROOT"
echo "Status: $STATUS_DIR"

echo "=========================================="
