#!/bin/bash
set -euo pipefail

if [[ $# -lt 3 ]]; then
    echo "Usage: bash $0 <QSM_file> <T1_file> <output_ID> [reference_anatomical_image(deprecated)]" >&2
    exit 1
fi

QSM="$1"
T1="$2"
OUT_ID="$3"
REF_ANAT="${4:-}"


if ! type module >/dev/null 2>&1; then
    for _init in /etc/profile.d/Modules.sh /usr/share/Modules/init/bash; do
        if [[ -r "$_init" ]]; then
            # shellcheck disable=SC1090
            source "$_init"
            break
        fi
    done
fi
if ! type module >/dev/null 2>&1; then
    echo "ERROR: Failed to initialize the module command; check /etc/profile.d/Modules.sh" >&2
    exit 1
fi
module load fsl/6.0.4

STANDARD="${FSLDIR}/data/standard"
MNI_BRAIN="${STANDARD}/MNI152_T1_2mm_brain.nii.gz"
MNI_BRAIN_MASK="${STANDARD}/MNI152_T1_2mm_brain_mask.nii.gz"
MNI_T1="${STANDARD}/MNI152_T1_2mm.nii.gz"
FNIRT_CFG="${FSLDIR}/etc/flirtsch/T1_2_MNI152_2mm.cnf"

OUT_DIR="${QSM_ROI_OUT_DIR:-/cwStorage/nodecw_group/FY_data/QSM_HUASHAN/roi_work}"
LOG_DIR="${OUT_DIR}/logs"
WORK_DIR="${OUT_DIR}/work/${OUT_ID}"
mkdir -p "$LOG_DIR"

STOP_AFTER="${QSM_ROI_STOP_AFTER:-mni}"
case "$STOP_AFTER" in
    mni|rigid) ;;
    *) echo "ERROR: QSM_ROI_STOP_AFTER must be mni or rigid; received "$STOP_AFTER"" >&2; exit 1 ;;
esac
KEEP_WORK="${QSM_ROI_KEEP_WORK:-0}"
KEEP_DIR="${QSM_ROI_KEEP_DIR:-${OUT_DIR}/keep}"

LOG="${LOG_DIR}/${OUT_ID}.log"

if [[ "$STOP_AFTER" == "rigid" ]]; then
    QC="${OUT_DIR}/${OUT_ID}_qc_rigid.txt"
else
    QC="${OUT_DIR}/${OUT_ID}_qc.txt"
fi
exec > >(tee -a "$LOG") 2>&1

STEP_T0=$SECONDS
elapsed() { echo "      [Elapsed $(($SECONDS - STEP_T0))s]"; STEP_T0=$SECONDS; }

echo "===== $(date '+%F %T')  node=$(hostname)  ID=${OUT_ID} ====="
echo "QSM : ${QSM}"
echo "T1  : ${T1}"
if [[ -n "$REF_ANAT" ]]; then
    echo "REF : ${REF_ANAT}  (deprecated; ignored in this run)"
else
    echo "REF : (not provided)"
fi

for f in "$QSM" "$T1" "$MNI_BRAIN" "$MNI_BRAIN_MASK" "$MNI_T1" "$FNIRT_CFG"; do
    [[ -e "$f" ]] || { echo "ERROR: Missing file $f" >&2; exit 2; }
done


FINAL="${OUT_DIR}/${OUT_ID}_QSM_MNI.nii.gz"
if [[ "$STOP_AFTER" == "rigid" ]]; then

    if [[ -s "${KEEP_DIR}/${OUT_ID}_QSM_to_T1.mat" && "${FORCE_RERUN:-0}" != "1" ]]; then
        echo "Already exists, skipping: ${KEEP_DIR}/${OUT_ID}_QSM_to_T1.mat"
        exit 0
    fi
elif [[ -s "$FINAL" && "${FORCE_RERUN:-0}" != "1" ]]; then
    echo "Already exists, skipping: $FINAL"
    exit 0
fi

QSM_DT=$(fslinfo "$QSM" | awk '/^datatype/{print $2}')
case "$QSM_DT" in
    2|4|8|16|64) ;;                
    *)
        echo "ERROR: QSM is not a scalar image (datatype=${QSM_DT}; 128 denotes RGB); cannot proceed" >&2
        {
            echo "ID=${OUT_ID}"
            echo "node=$(hostname)"
            echo "QSM=${QSM}"
            echo "T1=${T1}"
            echo "qsm_datatype=${QSM_DT}"
            echo "WARNING: QSM datatype=${QSM_DT} is unsupported (non-scalar; 128=RGB)."
            echo "      P0 may have selected a coverage map instead of a quantitative QSM image; no output for this subject."
        } > "$QC"
        exit 4
        ;;
esac

rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR"
trap 'rm -rf "$WORK_DIR"' EXIT

# ============================================================================
# 1. T1 downsampling and skull stripping
#    Do not apply BET to QSM: its reconstruction mask already restricts the brain region,
#    and its near-zero background violates BET intensity assumptions. In testing, BET
#    identified 6.94 million voxels (~3.6 L), effectively treating the whole FOV as brain.
# ============================================================================
T1_SUBSAMP="${QSM_ROI_T1_SUBSAMP:-2}"
if [[ "$T1_SUBSAMP" == "0" ]]; then
    T1_FOR_REG="$T1"
    echo "[1] T1 downsampling disabled; using native resolution"
else
    fslmaths "$T1" -subsamp2 "${WORK_DIR}/T1_sub.nii.gz"
    T1_FOR_REG="${WORK_DIR}/T1_sub.nii.gz"
    read -r S1 S2 S3 < <(fslinfo "$T1_FOR_REG" \
        | awk '/^pixdim[123]/{printf "%s ", $2} END{print ""}')
    echo "[1] T1 downsampled to ${S1}x${S2}x${S3}mm"
fi

bet "$T1_FOR_REG" "${WORK_DIR}/T1_brain" -R -f 0.5 -m

# BET can occasionally fail catastrophically on low-contrast T1 images, treating
# the whole head as brain. For subject 00455850, BET yielded 3.62 L instead of
# the expected 1.3-1.7 L. Such a mask invalidates affine registration and can
# cause FNIRT folding; otherwise the failure might only appear at step 5.
#
# Retry conservatively with -B (bias-field and neck cleanup) and a lower -f;
# accept the retry only if it yields a clearly more plausible brain mask.
# CRITICAL: Judge mask plausibility by physical volume, not voxel count.
# -subsamp2 selects every other voxel along each axis, doubling voxel spacing.
# Observed acquisition grids in 692 subjects:
#   1.0x1.0x1.0 mm -> 2x2x2 mm (8 mm3/voxel): ~170k-200k mask voxels (411 subjects)
#   1.0x0.5x0.5 mm -> 2x1x1 mm (2 mm3/voxel): ~650k-800k mask voxels (276 subjects)
#   0.5x0.5x1.0 mm -> 1x1x2 mm (2 mm3/voxel): similar counts (4 subjects)
#   0.6667 mm isotropic -> 1.3333 mm (2.37 mm3/voxel): 1 subject
# A 1.4 L brain can differ fourfold in voxel count; coarse 8 mm3 grids predominate.
# An earlier version applied a 300k-voxel threshold and falsely flagged all 2 mm3
# scans as BET failures, changing skull-stripping settings between subjects.
# fslstats -V returns "voxel_count volume_mm3" (e.g., 169224 1353792.000000).
# Divide the second field by 1000 to obtain mL; never multiply by voxel count again.
mask_volume_ml() { fslstats "$1" -V | awk '{printf "%.0f", $2/1000}'; }

BET_MASK_ML=$(mask_volume_ml "${WORK_DIR}/T1_brain_mask.nii.gz")
BET_PLAUSIBLE_MAX_ML=${QSM_ROI_BET_MAX_ML:-2000}
if (( BET_MASK_ML > BET_PLAUSIBLE_MAX_ML )); then
    echo "[1] WARNING: BET mask ${BET_MASK_ML}mL is implausibly large (expected 1.3-1.7L); retrying skull stripping"
    # The retry is a fallback and must not abort the subject if it fails.
    # In subject P190864, the second internal bet2 stage of bet -B segfaulted.
    # With set -e, an unguarded retry would abort the subject (exit 1) even
    # though the original 2311mL mask could still undergo Dice-based QC.
    set +e
    bet "$T1_FOR_REG" "${WORK_DIR}/T1_brain_retry" -R -B -f 0.3 -m
    BET_RETRY_RC=$?
    set -e
    if (( BET_RETRY_RC != 0 )) || [[ ! -s "${WORK_DIR}/T1_brain_retry_mask.nii.gz" ]]; then
        echo "[1] Skull-stripping retry failed (rc=${BET_RETRY_RC}); retaining original mask for Dice-based QC"
    else
        RETRY_ML=$(mask_volume_ml "${WORK_DIR}/T1_brain_retry_mask.nii.gz")
        echo "[1] Retried brain mask volume: ${RETRY_ML}mL"
        if (( RETRY_ML > 500 && RETRY_ML < BET_MASK_ML )); then
            mv "${WORK_DIR}/T1_brain_retry.nii.gz" "${WORK_DIR}/T1_brain.nii.gz"
            mv "${WORK_DIR}/T1_brain_retry_mask.nii.gz" "${WORK_DIR}/T1_brain_mask.nii.gz"
            BET_MASK_ML=$RETRY_ML
            echo "[1] Using the improved skull-stripping result"
        else
            echo "[1] Retry did not improve the mask; retaining original mask for Dice-based QC"
        fi
    fi
    rm -f "${WORK_DIR}/T1_brain_retry"*.nii.gz
fi

fslmaths "$QSM" -abs -bin "${WORK_DIR}/QSM_mask.nii.gz"
echo "[1] T1 skull stripping complete (brain mask ${BET_MASK_ML}mL); QSM nonzero mask volume $(mask_volume_ml "${WORK_DIR}/QSM_mask.nii.gz")mL"
elapsed

# ============================================================================
# 2. QSM -> T1_brain: six-degree-of-freedom rigid registration
# ============================================================================
SEARCH_MODE="${QSM_ROI_SEARCH_MODE:-nosearch}"
SEARCH_DEG="${QSM_ROI_SEARCH_DEG:-30}"
case "$SEARCH_MODE" in
    nosearch) SEARCH_OPTS=(-nosearch) ;;
    search)   SEARCH_OPTS=(-searchrx "-${SEARCH_DEG}" "${SEARCH_DEG}"
                           -searchry "-${SEARCH_DEG}" "${SEARCH_DEG}"
                           -searchrz "-${SEARCH_DEG}" "${SEARCH_DEG}") ;;
    wide)     SEARCH_OPTS=(-searchrx -90 90 -searchry -90 90 -searchrz -90 90) ;;
    *) echo "ERROR: Unknown QSM_ROI_SEARCH_MODE=$SEARCH_MODE (options: nosearch/search/wide)" >&2; exit 1 ;;
esac

flirt -in "$QSM" \
      -ref "${WORK_DIR}/T1_brain.nii.gz" \
      -usesqform \
      -dof 6 -cost corratio -searchcost corratio \
      "${SEARCH_OPTS[@]}" \
      -omat "${WORK_DIR}/QSM_to_T1.mat" \
      -out "${WORK_DIR}/QSM_in_T1.nii.gz"
echo "[2] QSM -> T1 rigid registration complete (search mode: ${SEARCH_MODE})"; elapsed

# 2-QC: Transform the nonzero QSM mask into T1 space with the same matrix and
# compute Dice overlap with the T1 brain mask. QSM covers ~73% of brain tissue
# because reconstruction zeros the CSF/extracerebral background; expected Dice is
flirt -in "${WORK_DIR}/QSM_mask.nii.gz" -ref "${WORK_DIR}/T1_brain.nii.gz" \
      -applyxfm -init "${WORK_DIR}/QSM_to_T1.mat" \
      -interp nearestneighbour -out "${WORK_DIR}/QSM_mask_in_T1.nii.gz"

MASK_A=$(fslstats "${WORK_DIR}/QSM_mask_in_T1.nii.gz" -V | awk '{print $1}')
MASK_B=$(fslstats "${WORK_DIR}/T1_brain_mask.nii.gz" -V | awk '{print $1}')
fslmaths "${WORK_DIR}/QSM_mask_in_T1.nii.gz" -mul "${WORK_DIR}/T1_brain_mask.nii.gz" \
         -bin "${WORK_DIR}/mask_overlap.nii.gz"
MASK_I=$(fslstats "${WORK_DIR}/mask_overlap.nii.gz" -V | awk '{print $1}')
DICE=$(awk -v a="$MASK_A" -v b="$MASK_B" -v i="$MASK_I" \
       'BEGIN{ if (a+b==0) print "0.000"; else printf "%.3f", 2*i/(a+b) }')
echo "[2-QC] QSM->T1 brain-mask Dice = ${DICE} (QSM ${MASK_A} / T1 ${MASK_B} / intersection ${MASK_I})"

# Reject poor registrations early rather than spending another 10-40 minutes on FNIRT.
# In subject 00455850, BET included the entire head (mask 3.62 L), yielding
# Dice 0.510. FNIRT on such a case folds and fails only at the later QC stage.
# Write QC and exit before FNIRT when the rigid-registration Dice gate fails.
DICE_GATE="${QSM_ROI_DICE_MIN:-0.60}"
if awk -v d="$DICE" -v g="$DICE_GATE" 'BEGIN{exit !(d < g)}'; then
    echo "WARNING: QSM->T1 brain-mask Dice=${DICE} is below ${DICE_GATE}; skipping FNIRT" >&2
    {
        echo "ID=${OUT_ID}"
        echo "node=$(hostname)"
        echo "QSM=${QSM}"
        echo "T1=${T1}"
        echo "reg_mask_dice=${DICE}"
        echo "t1_brain_mask_voxels=${MASK_B}"
        echo "t1_brain_mask_ml=$(mask_volume_ml "${WORK_DIR}/T1_brain_mask.nii.gz")"
        echo "WARNING: reg_mask_dice ${DICE} < ${DICE_GATE} (often due to failed T1 skull stripping;"
        echo "      t1_brain_mask_ml can greatly exceed 1700); stopped before FNIRT"
    } > "$QC"
    exit 5
fi

# ---------------------------------------------------------------------------
# 2b. In rigid mode, stop here; retain registration products for P8.
#     FNIRT in step 3 is unnecessary in this mode (saves ~7-50 min/subject).
# ---------------------------------------------------------------------------
keep_work() {
    mkdir -p "$KEEP_DIR"
    local f
    for f in QSM_to_T1.mat T1_brain.nii.gz T1_brain_mask.nii.gz; do
        if [[ -e "${WORK_DIR}/${f}" ]]; then
            cp -f "${WORK_DIR}/${f}" "${KEEP_DIR}/${OUT_ID}_${f}"
        else
            echo "WARNING: Cannot find ${WORK_DIR}/${f} when preserving intermediate files" >&2
        fi
    done
}

if [[ "$STOP_AFTER" == "rigid" ]]; then
    # Use an if statement rather than [[ ... ]] && keep_work: when KEEP_WORK=0,
    # the latter can return status 1 and trigger set -e.
    if [[ "$KEEP_WORK" == "1" ]]; then
        keep_work
    fi
    {
        echo "ID=${OUT_ID}"
        echo "node=$(hostname)"
        echo "QSM=${QSM}"
        echo "T1=${T1}"
        echo "stop_after=rigid"
        echo "t1_pixdim=${S1:-na}x${S2:-na}x${S3:-na}"
        echo "reg_mask_dice=${DICE}"
        echo "t1_brain_mask_voxels=${MASK_B}"
        echo "t1_brain_mask_ml=$(mask_volume_ml "${WORK_DIR}/T1_brain_mask.nii.gz")"
        echo "search_mode=${SEARCH_MODE}"
        echo "keep_dir=${KEEP_DIR}"
    } > "$QC"
    echo "[2b] Rigid mode complete (Dice=${DICE}); FNIRT skipped; intermediates in ${KEEP_DIR}"
    exit 0
fi

# ============================================================================
# 3. T1_brain -> MNI152_T1_2mm: 12-DOF affine registration, then nonlinear FNIRT
# ============================================================================
flirt -in "${WORK_DIR}/T1_brain.nii.gz" -ref "$MNI_BRAIN" \
      -dof 12 -omat "${WORK_DIR}/T1_to_MNI_lin.mat" \
      -out "${WORK_DIR}/T1_to_MNI_lin.nii.gz"

# FNIRT may take much longer in isolated cases. Subject P181633 (qsm_tra_iso1.0mm
# and 3d_sag_iso T1) converged after ~50 min (typical 7-11 min), despite recurring
# Jacobian folding warnings. The preceding 12-DOF affine Dice was 0.866.
# Using --inmask/--refmask (often advised by FSL) was worse for these data:
# T1-to-MNI Dice fell to 0.595, deformation became unstable, and runtime reached
# 63 min. Therefore, these masks are not used. The timeout is 60 min to allow
# slow valid cases without occupying a worker indefinitely. Timeouts/failures
# exit with code 6, are written to QC, and are collected in failed_ID.txt.
FNIRT_TIMEOUT=${QSM_ROI_FNIRT_TIMEOUT:-3600}
set +e
timeout "$FNIRT_TIMEOUT" fnirt --in="${WORK_DIR}/T1_brain.nii.gz" \
      --aff="${WORK_DIR}/T1_to_MNI_lin.mat" \
      --ref="$MNI_T1" \
      --config="$FNIRT_CFG" \
      --cout="${WORK_DIR}/T1_to_MNI_warp.nii.gz" \
      --iout="${WORK_DIR}/T1_in_MNI.nii.gz" \
      --logout="${LOG_DIR}/${OUT_ID}_fnirt.log"
FNIRT_RC=$?
set -e
if (( FNIRT_RC != 0 )); then
    if (( FNIRT_RC == 124 )); then
        REASON="Did not converge within ${FNIRT_TIMEOUT}s (possibly atypical T1 contrast/FOV)"
    else
        REASON="Nonzero exit code ${FNIRT_RC}"
    fi
    echo "WARNING: FNIRT ${REASON}; no output for this subject" >&2
    {
        echo "ID=${OUT_ID}"
        echo "node=$(hostname)"
        echo "QSM=${QSM}"
        echo "T1=${T1}"
        echo "reg_mask_dice=${DICE}"
        echo "t1_brain_mask_ml=$(mask_volume_ml "${WORK_DIR}/T1_brain_mask.nii.gz")"
        echo "fnirt_rc=${FNIRT_RC}"
        echo "fnirt_timeout_s=${FNIRT_TIMEOUT}"
        echo "WARNING: FNIRT incomplete (${REASON}); QSM_MNI not produced; manual review required"
    } > "$QC"
    exit 6
fi
echo "[3] T1 -> MNI affine and nonlinear registration complete"; elapsed

# 3-QC: Transform the T1 brain mask into MNI space and compute Dice against
# the template brain mask, independently assessing the T1-to-MNI registration.
applywarp -i "${WORK_DIR}/T1_brain_mask.nii.gz" -r "$MNI_T1" \
          -w "${WORK_DIR}/T1_to_MNI_warp.nii.gz" \
          --interp=nn -o "${WORK_DIR}/T1_mask_in_MNI.nii.gz" >/dev/null 2>&1
MB=$(fslstats "${WORK_DIR}/T1_mask_in_MNI.nii.gz" -V | awk '{print $1}')
MR=$(fslstats "$MNI_BRAIN_MASK" -V | awk '{print $1}')
fslmaths "${WORK_DIR}/T1_mask_in_MNI.nii.gz" -mul "$MNI_BRAIN_MASK" -bin \
         "${WORK_DIR}/T1_mask_ov.nii.gz"
MI=$(fslstats "${WORK_DIR}/T1_mask_ov.nii.gz" -V | awk '{print $1}')
T1_MNI_DICE=$(awk -v a="$MB" -v b="$MR" -v i="$MI" \
    'BEGIN{ if (a+b==0) print "0.000"; else printf "%.3f", 2*i/(a+b) }')
echo "[3-QC] T1->MNI brain-mask Dice = ${T1_MNI_DICE} (MNI-warped ${MB} / template ${MR} / intersection ${MI})"

# ============================================================================
# 4. Compose/apply: QSM --(QSM_to_T1.mat)--> T1 --(12-DOF + FNIRT)--> MNI.
#    The --premat transform was estimated between QSM and T1_brain in step 2.
# ============================================================================
applywarp -i "$QSM" -r "$MNI_T1" \
          -w "${WORK_DIR}/T1_to_MNI_warp.nii.gz" \
          --premat="${WORK_DIR}/QSM_to_T1.mat" \
          --interp=trilinear \
          -o "${WORK_DIR}/${OUT_ID}_QSM_MNI.nii.gz"
echo "[4] QSM transformed to MNI 2mm space"; elapsed

# ============================================================================
# 5. Quality control
#    voxels_in_MNI_brain is computed using MNI152_T1_2mm_brain.nii.gz, which
#    contains only ~228k mask voxels; this imposes a ~228k upper bound. Assess
#    fractional coverage of the MNI brain mask instead of an unreachable cutoff.
# ============================================================================
read -r D1 D2 D3 < <(fslinfo "${WORK_DIR}/${OUT_ID}_QSM_MNI.nii.gz" \
    | awk '/^dim[123]/{printf "%s ", $2} END{print ""}')
NONZERO=$(fslstats "${WORK_DIR}/${OUT_ID}_QSM_MNI.nii.gz" -V | awk '{print $1}')
INBRAIN=$(fslstats "${WORK_DIR}/${OUT_ID}_QSM_MNI.nii.gz" -k "$MNI_BRAIN" -V | awk '{print $1}')
RANGE=$(fslstats "${WORK_DIR}/${OUT_ID}_QSM_MNI.nii.gz" -k "$MNI_BRAIN" -R)
MNI_MASK_VOX=$(fslstats "$MNI_BRAIN_MASK" -V | awk '{print $1}')
COVER=$(awk -v a="$INBRAIN" -v b="$MNI_MASK_VOX" \
        'BEGIN{ if (b==0) print "0.000"; else printf "%.3f", a/b }')

{
    echo "ID=${OUT_ID}"
    echo "node=$(hostname)"
    echo "QSM=${QSM}"
    echo "T1=${T1}"
    echo "REF=${REF_ANAT}"
    echo "t1_pixdim=${S1:-na}x${S2:-na}x${S3:-na}"
    echo "out_dim=${D1}x${D2}x${D3}"
    echo "reg_mask_dice=${DICE}"
    echo "t1_mni_dice=${T1_MNI_DICE}"
    echo "search_mode=${SEARCH_MODE}"
    echo "search_deg=${SEARCH_DEG}"
    echo "nonzero_voxels=${NONZERO}"
    echo "voxels_in_MNI_brain=${INBRAIN}"
    echo "mni_brain_mask_voxels=${MNI_MASK_VOX}"
    echo "mni_brain_coverage=${COVER}"
    echo "mni_brain_range=${RANGE}"
} > "$QC"

mv "${WORK_DIR}/${OUT_ID}_QSM_MNI.nii.gz" "$FINAL"

if [[ "${D1}x${D2}x${D3}" != "91x109x91" ]]; then
    echo "WARNING: Output dimensions ${D1}x${D2}x${D3} are not 91x109x91; P2 extraction will fail" >&2
    echo "WARNING: Output dimensions are not 91x109x91" >> "$QC"
    exit 3
fi

DICE_GATE="${QSM_ROI_DICE_MIN:-0.60}"
if awk -v d="$DICE" -v g="$DICE_GATE" 'BEGIN{exit !(d < g)}'; then
    echo "WARNING: QSM->T1 brain-mask Dice=${DICE} is below ${DICE_GATE}" >&2
    echo "WARNING: reg_mask_dice ${DICE} < ${DICE_GATE}" >> "$QC"
    exit 5
fi
if awk -v d="$T1_MNI_DICE" -v g="$DICE_GATE" 'BEGIN{exit !(d < g)}'; then
    echo "WARNING: T1->MNI brain-mask Dice=${T1_MNI_DICE} is below ${DICE_GATE}" >&2
    echo "WARNING: t1_mni_dice ${T1_MNI_DICE} < ${DICE_GATE}" >> "$QC"
    exit 6
fi

# Coverage gate: require sufficient nonzero QSM voxels inside the MNI brain mask.
# Correct registrations typically yield 0.80-0.91 (P164054 with 180 slices: 0.914;
# P165540: 0.803). The older SWI workflow yielded ~0.43 from misregistration.
# The gate is 0.55: it rejects obvious misregistrations but accommodates scans
# with shorter QSM coverage (e.g., 118 mm for some iso1.0mm acquisitions;
# 222 of 711 subjects). Per-ROI coverage is reported by P2 roi_coverage.csv.
COVER_GATE="${QSM_COVER_MIN:-0.55}"
if awk -v c="$COVER" -v g="$COVER_GATE" 'BEGIN{exit !(c < g)}'; then
    echo "WARNING: MNI-brain QSM coverage ${COVER} is below ${COVER_GATE}; possible misregistration" >&2
    echo "WARNING: mni_brain_coverage ${COVER} < ${COVER_GATE}" >> "$QC"
    exit 7
fi

echo "[5] QC passed; output: $FINAL"
echo "===== $(date '+%F %T') ${OUT_ID} completed ====="
