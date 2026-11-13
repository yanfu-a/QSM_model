
#!/bin/bash
set -euo pipefail

# ==========================================================
# 1. Environment
# ==========================================================

module load fsl

export FSLOUTPUTTYPE=NIFTI_GZ

# 避免多进程与程序内部线程过度竞争
export OMP_NUM_THREADS=1

# ==========================================================
# 2. Path configuration
# ==========================================================

# 已经通过 dcm2niix 转换的 NIfTI
INPUT_ROOT="/cwStorage/nodecw_group/FY_data/PPMI/T1_nifti"

# MNI 2mm 输出目录
OUTPUT_ROOT="/cwStorage/nodecw_group/FY_data/PPMI/T1_MNI2mm"

# MNI152 2mm templates
MNI_REF="$FSLDIR/data/standard/MNI152_T1_2mm.nii.gz"

MNI_BRAIN="$FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz"

# FNIRT configuration
FNIRT_CONFIG="$FSLDIR/etc/flirtsch/T1_2_MNI152_2mm.cnf"

# ==========================================================
# 3. Processing parameters
# ==========================================================

# 并行任务数
JOBS="${JOBS:-4}"

# BET fractional intensity threshold
BET_F="${BET_F:-0.25}"

# 是否进行 robustfov
# 1 = 裁剪颈部
# 0 = 保留原始 FOV
USE_ROBUSTFOV="${USE_ROBUSTFOV:-1}"

# 是否进行 fslreorient2std
# 默认关闭，避免掩盖错误的方向信息
REORIENT="${REORIENT:-0}"

# 可选：指定序列名称关键词
# 空字符串表示处理所有 NIfTI
SERIES_FILTER="${SERIES_FILTER:-}"

# ==========================================================
# 4. Environment checks
# ==========================================================

for cmd in fslhd fslreorient2std robustfov bet \
           flirt fnirt applywarp; do

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
    echo "ERROR: MNI brain template does not exist"
    exit 1
}

[[ -f "$FNIRT_CONFIG" ]] || {
    echo "ERROR: FNIRT configuration does not exist"
    exit 1
}

[[ "$JOBS" =~ ^[1-9][0-9]*$ ]] || {
    echo "ERROR: JOBS must be a positive integer"
    exit 1
}

mkdir -p "$OUTPUT_ROOT"

RUN_ID=$(date +%Y%m%d_%H%M%S)

STATUS_DIR="$OUTPUT_ROOT/_status_$RUN_ID"

mkdir -p "$STATUS_DIR"

export INPUT_ROOT OUTPUT_ROOT STATUS_DIR
export MNI_REF MNI_BRAIN FNIRT_CONFIG
export BET_F USE_ROBUSTFOV REORIENT

# ==========================================================
# 5. Registration function
# ==========================================================

register_one() {

    local input_file="$1"

    # ------------------------------------------------------
    # A. Input information
    # ------------------------------------------------------

    local relative_path="${input_file#"$INPUT_ROOT"/}"

    local relative_dir
    relative_dir=$(dirname "$relative_path")

    local filename
    filename=$(basename "$input_file")

    local base="${filename%.nii.gz}"

    local output_dir="$OUTPUT_ROOT/$relative_dir"

    mkdir -p "$output_dir"
    mkdir -p "$output_dir/transforms"

    local key
    key=$(printf '%s' "$relative_path" |
          md5sum | awk '{print $1}')

    # ------------------------------------------------------
    # B. Output files
    # ------------------------------------------------------

    local MNI_OUT="$output_dir/${base}_MNI2mm.nii.gz"

    local BRAIN_OUT="$output_dir/${base}_brain_MNI2mm.nii.gz"

    local MASK_OUT="$output_dir/${base}_brain_mask_MNI2mm.nii.gz"

    local LINEAR_OUT="$output_dir/${base}_linear_brain_MNI2mm.nii.gz"

    local AFFINE_OUT="$output_dir/transforms/${base}_to_MNI2mm_affine.mat"

    local WARP_OUT="$output_dir/transforms/${base}_to_MNI2mm_warpcoef.nii.gz"

    local LOG_FILE="$output_dir/${base}_registration.log"

    # ------------------------------------------------------
    # C. Skip completed subjects
    # ------------------------------------------------------

    if [[ -s "$MNI_OUT" ]] &&
       [[ -s "$BRAIN_OUT" ]] &&
       [[ -s "$MASK_OUT" ]] &&
       [[ -s "$LINEAR_OUT" ]] &&
       [[ -s "$AFFINE_OUT" ]] &&
       [[ -s "$WARP_OUT" ]]; then

        touch "$STATUS_DIR/${key}.skipped"

        echo "[SKIP] $relative_path"

        return 0
    fi

    echo
    echo "[START] $relative_path"

    # ------------------------------------------------------
    # D. Temporary working directory
    # ------------------------------------------------------

    local WORK_DIR

    WORK_DIR=$(mktemp -d "$output_dir/.work_${base}.XXXXXX")

    local REORIENT_IMG="$WORK_DIR/T1_reorient.nii.gz"

    local CROPPED_IMG="$WORK_DIR/T1_cropped.nii.gz"

    local BRAIN_IMG="$WORK_DIR/T1_brain.nii.gz"

    local BRAIN_MASK="$WORK_DIR/T1_brain_mask.nii.gz"

    local AFFINE_TMP="$WORK_DIR/T1_to_MNI_affine.mat"

    local LINEAR_TMP="$WORK_DIR/T1_linear.nii.gz"

    local WARP_TMP="$WORK_DIR/T1_to_MNI_warpcoef.nii.gz"

    local MNI_TMP="$WORK_DIR/T1_MNI2mm.nii.gz"

    local BRAIN_MNI_TMP="$WORK_DIR/T1_brain_MNI2mm.nii.gz"

    local MASK_MNI_TMP="$WORK_DIR/T1_brain_mask_MNI2mm.nii.gz"

    echo "Input: $input_file" > "$LOG_FILE"

    echo "Start: $(date)" >> "$LOG_FILE"

    # ------------------------------------------------------
    # E. Header QC
    # ------------------------------------------------------

    local header

    if ! header=$(fslhd "$input_file" 2>>"$LOG_FILE"); then

        touch "$STATUS_DIR/${key}.failed"

        echo "[FAILED HEADER] $relative_path"

        echo "Work directory: $WORK_DIR"

        return 0

    fi

    local qform_code
    local sform_code
    local dim0
    local dim4

    qform_code=$(awk '$1=="qform_code"{print $2}' <<< "$header")

    sform_code=$(awk '$1=="sform_code"{print $2}' <<< "$header")

    dim0=$(awk '$1=="dim0"{print $2}' <<< "$header")

    dim4=$(awk '$1=="dim4"{print $2}' <<< "$header")

    if [[ "${qform_code:-0}" -eq 0 ]] &&
       [[ "${sform_code:-0}" -eq 0 ]]; then

        echo "ERROR: Both qform and sform are invalid" >> "$LOG_FILE"

        touch "$STATUS_DIR/${key}.failed"

        echo "[FAILED ORIENTATION] $relative_path"

        return 0

    fi

    if [[ "${dim0:-0}" -ne 3 ]] ||
       [[ "${dim4:-0}" -ne 1 ]]; then

        echo "ERROR: Only 3D T1 images are accepted" >> "$LOG_FILE"

        touch "$STATUS_DIR/${key}.failed"

        echo "[FAILED DIMENSION] $relative_path"

        return 0

    fi

    # ------------------------------------------------------
    # F. Optional orientation reorientation
    # ------------------------------------------------------

    local CURRENT_IMG="$input_file"

    if [[ "$REORIENT" == "1" ]]; then

        if ! fslreorient2std \
            "$input_file" \
            "$REORIENT_IMG" >> "$LOG_FILE" 2>&1; then

            touch "$STATUS_DIR/${key}.failed"

            echo "[FAILED REORIENT] $relative_path"

            echo "Work directory: $WORK_DIR"

            return 0

        fi

        CURRENT_IMG="$REORIENT_IMG"

    fi

    # ------------------------------------------------------
    # G. Crop neck
    # ------------------------------------------------------

    if [[ "$USE_ROBUSTFOV" == "1" ]]; then

        if ! robustfov \
            -i "$CURRENT_IMG" \
            -r "$CROPPED_IMG" >> "$LOG_FILE" 2>&1; then

            touch "$STATUS_DIR/${key}.failed"

            echo "[FAILED ROBUSTFOV] $relative_path"

            echo "Work directory: $WORK_DIR"

            return 0

        fi

    else

        if ! cp "$CURRENT_IMG" "$CROPPED_IMG"; then

            touch "$STATUS_DIR/${key}.failed"

            echo "[FAILED COPY] $relative_path"

            return 0

        fi

    fi

    # ------------------------------------------------------
    # H. Brain extraction
    # ------------------------------------------------------

    if ! bet \
        "$CROPPED_IMG" \
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

        echo "Work directory: $WORK_DIR"

        return 0

    fi

    # ------------------------------------------------------
    # I. FLIRT linear registration
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
    # J. FNIRT nonlinear registration
    # ------------------------------------------------------

    if ! fnirt \
        --in="$CROPPED_IMG" \
        --aff="$AFFINE_TMP" \
        --ref="$MNI_REF" \
        --config="$FNIRT_CONFIG" \
        --cout="$WARP_TMP" >> "$LOG_FILE" 2>&1; then

        touch "$STATUS_DIR/${key}.failed"

        echo "[FAILED FNIRT] $relative_path"

        echo "Work directory: $WORK_DIR"

        return 0

    fi

    # ------------------------------------------------------
    # K. Apply nonlinear transformation
    # ------------------------------------------------------

    # Full cropped T1
    if ! applywarp \
        --in="$CROPPED_IMG" \
        --ref="$MNI_REF" \
        --warp="$WARP_TMP" \
        --out="$MNI_TMP" \
        --interp=spline >> "$LOG_FILE" 2>&1; then

        touch "$STATUS_DIR/${key}.failed"

        echo "[FAILED APPLYWARP T1] $relative_path"

        echo "Work directory: $WORK_DIR"

        return 0

    fi

    # Brain-only T1
    if ! applywarp \
        --in="$BRAIN_IMG" \
        --ref="$MNI_REF" \
        --warp="$WARP_TMP" \
        --out="$BRAIN_MNI_TMP" \
        --interp=spline >> "$LOG_FILE" 2>&1; then

        touch "$STATUS_DIR/${key}.failed"

        echo "[FAILED APPLYWARP BRAIN] $relative_path"

        echo "Work directory: $WORK_DIR"

        return 0

    fi

    # Brain mask: nearest-neighbour interpolation
    if ! applywarp \
        --in="$BRAIN_MASK" \
        --ref="$MNI_REF" \
        --warp="$WARP_TMP" \
        --out="$MASK_MNI_TMP" \
        --interp=nn >> "$LOG_FILE" 2>&1; then

        touch "$STATUS_DIR/${key}.failed"

        echo "[FAILED APPLYWARP MASK] $relative_path"

        echo "Work directory: $WORK_DIR"

        return 0

    fi

    # ------------------------------------------------------
    # L. Verify output dimensions
    # ------------------------------------------------------

    local ref_dims
    local out_dims

    ref_dims=$(fslhd "$MNI_REF" |
        awk '$1=="dim1" || $1=="dim2" || $1=="dim3" {printf "%s ",$2}')

    out_dims=$(fslhd "$MNI_TMP" |
        awk '$1=="dim1" || $1=="dim2" || $1=="dim3" {printf "%s ",$2}')

    if [[ "$ref_dims" != "$out_dims" ]]; then

        echo "ERROR: Output dimensions do not match MNI reference" >> "$LOG_FILE"

        touch "$STATUS_DIR/${key}.failed"

        echo "[FAILED DIMENSION QC] $relative_path"

        echo "Work directory: $WORK_DIR"

        return 0

    fi

    # ------------------------------------------------------
    # M. Save final outputs
    # ------------------------------------------------------

    mv "$MNI_TMP" "$MNI_OUT"

    mv "$BRAIN_MNI_TMP" "$BRAIN_OUT"

    mv "$MASK_MNI_TMP" "$MASK_OUT"

    mv "$LINEAR_TMP" "$LINEAR_OUT"

    mv "$AFFINE_TMP" "$AFFINE_OUT"

    mv "$WARP_TMP" "$WARP_OUT"

    # Clean temporary directory after success
    rm -rf "$WORK_DIR"

    echo "Finish: $(date)" >> "$LOG_FILE"

    touch "$STATUS_DIR/${key}.success"

    echo "[SUCCESS] $relative_path"

}

export -f register_one

# ==========================================================
# 6. Parallel registration
# ==========================================================

echo
echo "=========================================="
echo "PPMI T1 MNI152 2mm Registration"
echo "=========================================="

echo "Input: $INPUT_ROOT"
echo "Output: $OUTPUT_ROOT"
echo "MNI: $MNI_REF"
echo "FNIRT config: $FNIRT_CONFIG"
echo "Jobs: $JOBS"
echo "=========================================="

if [[ -n "$SERIES_FILTER" ]]; then

    find "$INPUT_ROOT" \
        -type f \
        -name '*.nii.gz' \
        -path "*${SERIES_FILTER}*" \
        -print0 |

    xargs -0 -r -n 1 -P "$JOBS" \
        bash -c 'register_one "$1"' _

else

    find "$INPUT_ROOT" \
        -type f \
        -name '*.nii.gz' \
        -print0 |

    xargs -0 -r -n 1 -P "$JOBS" \
        bash -c 'register_one "$1"' _

fi

# ==========================================================
# 7. Summary
# ==========================================================

success=$(find "$STATUS_DIR" -type f -name '*.success' | wc -l)

failed=$(find "$STATUS_DIR" -type f -name '*.failed' | wc -l)

skipped=$(find "$STATUS_DIR" -type f -name '*.skipped' | wc -l)

total=$((success + failed + skipped))

echo
echo "=========================================="
echo "Registration completed"
echo "=========================================="

echo "Total: $total"
echo "Success: $success"
echo "Skipped: $skipped"
echo "Failed: $failed"

if [[ "$total" -gt 0 ]]; then

    awk -v s="$success" -v t="$total" \
        'BEGIN {printf "Success rate: %.2f%%\n", s/t*100}'

fi

echo "Output directory: $OUTPUT_ROOT"
echo "Status directory: $STATUS_DIR"

echo "=========================================="
