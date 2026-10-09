#!/bin/bash
set -euo pipefail

if [[ $# -lt 3 ]]; then
    echo "Usage: bash $0 <QSM_FILE> <T1_FILE> <OUTPUT_ID>" >&2
    exit 1
fi

QSM="$1"
T1="$2"
OUT_ID="$3"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NATIVE_ROOT="${QSM_ROI_NATIVE_ROOT:-/cwStorage/nodecw_group/FY_data/QSM_HUASHAN/native_data}"
SYNTHSEG_THREADS="${P8_SYNTHSEG_THREADS:-16}"
PY310=/home1/fuyan/.conda/envs/py310/bin/python

LABEL_DIR="${NATIVE_ROOT}/label"
T1SEG_DIR="${NATIVE_ROOT}/t1seg"
MAT_DIR="${NATIVE_ROOT}/mat"
P1_DIR="${NATIVE_ROOT}/p1_rigid"
QC_DIR="${NATIVE_ROOT}/qc"
LOG_DIR="${NATIVE_ROOT}/logs"
mkdir -p "$LABEL_DIR" "$T1SEG_DIR" "$MAT_DIR" "$P1_DIR" "$QC_DIR" "$LOG_DIR"

QC="${QC_DIR}/${OUT_ID}_native_qc.txt"
FINAL="${LABEL_DIR}/${OUT_ID}_aparc_aseg_QSMnative.nii.gz"
T1SEG="${T1SEG_DIR}/${OUT_ID}_aparc+aseg_T1.nii.gz"

exec > >(tee -a "${LOG_DIR}/${OUT_ID}.log") 2>&1

echo "=================================================================="
echo "P8 native branch  ID=${OUT_ID}  host=$(hostname)  $(date '+%F %T')"
echo "  QSM = $QSM"
echo "  T1  = $T1"
echo "  native_root = $NATIVE_ROOT"
echo "=================================================================="

if [[ -s "$FINAL" && "${P8_FORCE_RERUN:-0}" != "1" ]]; then
    echo "Existing output; skipping: $FINAL"
    exit 0
fi

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
    echo "ERROR: Unable to initialize the module command" >&2
    exit 1
fi
module load fsl/6.0.4

for _v in QSM_ROI_SEARCH_MODE QSM_ROI_SEARCH_DEG QSM_ROI_DICE_MIN \
          QSM_ROI_T1_SUBSAMP QSM_ROI_FNIRT_TIMEOUT; do
    if [[ -n "${!_v:-}" ]]; then
        export "$_v"
    fi
done
unset _v

export FREESURFER_HOME=/public/software/apps/Freesurfer/8.2.0-1
export FSLOUTPUTTYPE=NIFTI_GZ
SYNTHSEG="${FREESURFER_HOME}/bin/mri_synthseg"
if [[ ! -x "$SYNTHSEG" ]]; then
    echo "ERROR: Cannot find $SYNTHSEG" >&2
    exit 1
fi

for f in "$QSM" "$T1"; do
    [[ -e "$f" ]] || { echo "ERROR: Missing file $f" >&2; exit 2; }
done

QSM_DT=$(fslinfo "$QSM" | awk '/^datatype/{print $2}')
if [[ "$QSM_DT" == "128" ]]; then
    echo "ERROR: QSM image is not scalar (datatype=128, RGB); cannot extract QSM" >&2
    exit 4
fi

read -r Q1 Q2 Q3 < <(fslinfo "$QSM" | awk '/^pixdim[123]/{printf "%s ", $2} END{print ""}')
QSM_PIX="${Q1}x${Q2}x${Q3}"
read -r P1D P2D P3D < <(fslinfo "$T1" | awk '/^pixdim[123]/{printf "%s ", $2} END{print ""}')
T1_PIX="${P1D}x${P2D}x${P3D}"

echo "[1] Run P1 (STOP_AFTER=rigid) to obtain QSM_to_T1.mat ..."
QSM_ROI_STOP_AFTER=rigid \
QSM_ROI_KEEP_WORK=1 \
QSM_ROI_KEEP_DIR="$MAT_DIR" \
QSM_ROI_OUT_DIR="$P1_DIR" \
    bash "${HERE}/P1_qsm_to_mni.sh" "$QSM" "$T1" "$OUT_ID"

MAT_Q2T="${MAT_DIR}/${OUT_ID}_QSM_to_T1.mat"
if [[ ! -s "$MAT_Q2T" ]]; then
    echo "ERROR: P1 did not retain $MAT_Q2T (work-file preservation may have failed)" >&2
    exit 11
fi


P1_QC="${P1_DIR}/${OUT_ID}_qc_rigid.txt"
DICE=$(awk -F= '/^reg_mask_dice/{print $2}' "$P1_QC" 2>/dev/null || echo "")

echo "[2] Run mri_synthseg --parc (--keepgeom, ${SYNTHSEG_THREADS} threads) ..."
if ! "$SYNTHSEG" --i "$T1" --o "$T1SEG" \
        --parc --keepgeom --cpu --threads "$SYNTHSEG_THREADS" \
        --vol "${QC_DIR}/${OUT_ID}_synthseg_vol.csv"; then
    echo "ERROR: mri_synthseg failed" >&2
    exit 13
fi
[[ -s "$T1SEG" ]] || { echo "ERROR: SynthSeg did not generate $T1SEG" >&2; exit 13; }


T1_TO_QSM="${MAT_DIR}/${OUT_ID}_T1_to_QSM.mat"
T1_BRAIN="${MAT_DIR}/${OUT_ID}_T1_brain.nii.gz"
if [[ -s "$T1_BRAIN" ]]; then
    "$PY310" - "$T1" "$T1_BRAIN" "$MAT_Q2T" "$T1_TO_QSM" <<'PYEOF' || exit 12
import sys
import nibabel as nib
import numpy as np

t1_p, brain_p, q2t_p, out_p = sys.argv[1:5]
full = nib.load(t1_p)
sub = nib.load(brain_p)
Df = np.diag(list(full.header.get_zooms()[:3]) + [1.0])
Ds = np.diag(list(sub.header.get_zooms()[:3]) + [1.0])
A_f, A_s = full.affine, sub.affine

# World-coordinate-preserving grid conversion between the scaled voxel coordinates of two image grids.
# p = D*v and w = A*v; therefore p_sub = D_s * inv(A_s) * A_f * inv(D_f) * p_full.
M_g = Ds @ np.linalg.inv(A_s) @ A_f @ np.linalg.inv(Df)

# The linear component of a pure grid conversion must be the identity matrix. Any rotation or scaling
# implies the grids are not resamplings of the same image (e.g., T1_brain is from another image); abort.
if not np.allclose(M_g[:3, :3], np.eye(3), atol=1e-3):
    print(f"ERROR: Grid conversion between full-resolution T1 and T1_brain is not a pure translation; "
          f"linear component=\n{np.round(M_g[:3, :3], 5)}\n"
          f"  The images are not matching resampled grids; aborting.", file=sys.stderr)
    sys.exit(12)

shift = M_g[:3, 3]
ratio = np.array(full.shape[:3], float) / np.array(sub.shape[:3], float)
print(f"  [3] Grid conversion M_g: translation {np.round(shift, 4)} (voxel ratio {np.round(ratio, 4)})")
if np.allclose(shift, 0, atol=1e-3):
    print("      WARNING: M_g is identity, which is unexpected for this cohort. Verify that T1_brain "
          "was generated with fslmaths -subsamp2.")
# The translation should be close to integer voxel units (subsamp2 samples odd-indexed voxels,
# typically yielding -1). Non-integer shifts may indicate changed subsampling behavior; record a warning.
if not np.allclose(shift, np.round(shift), atol=0.05):
    print(f"      WARNING: Translation {np.round(shift, 4)} is not an integer voxel offset; nonfatal anomaly.")

M_q2t = np.loadtxt(q2t_p)
M_f2q = np.linalg.inv(M_q2t) @ M_g
np.savetxt(out_p, M_f2q, fmt="%.10f")
print(f"  [3] Saved {out_p} (full-resolution T1 to QSM transform)")
PYEOF
else
    echo "WARNING: Cannot find $T1_BRAIN; grid conversion is unavailable and processing must stop" >&2
    exit 11
fi

echo "[4] Resample label image to native QSM grid using FLIRT (nearest-neighbor) ..."
if ! flirt -in "$T1SEG" -ref "$QSM" -applyxfm -init "$T1_TO_QSM" \
        -interp nearestneighbour -datatype short -out "$FINAL"; then
    echo "ERROR: FLIRT label resampling failed" >&2
    exit 14
fi
[[ -s "$FINAL" ]] || { echo "ERROR: FLIRT did not generate $FINAL" >&2; exit 14; }

echo "[5] QC ..."
"$PY310" - "$QSM" "$FINAL" "$T1SEG" "$QC" "${QC_DIR}/${OUT_ID}_synthseg_vol.csv" \
            "$OUT_ID" "$DICE" "$QSM_PIX" "$T1_PIX" "$QSM_DT" "$HERE" <<'PYEOF'
import csv
import json
import sys
from pathlib import Path

import nibabel as nib
import numpy as np

(qsm_p, lab_p, t1seg_p, qc_p, vol_p,
 oid, dice, qsm_pix, t1_pix, qsm_dt, here) = sys.argv[1:12]

# ---- FreeSurfer cortical labels: 34 per hemisphere with --parc, including insula; excluding 1004/2004 ----
# Previously verified against $FREESURFER_HOME/models/synthseg_parcellation_labels.npy:
# left: 1001,1002,1003,1005..1035; right: left + 1000; 34 labels per hemisphere.
CORTEX_L = [1001, 1002, 1003] + list(range(1005, 1036))
CORTEX_R = [c + 1000 for c in CORTEX_L]
# The model label list excludes corpus callosum. Keep this assertion to detect unexpected model updates.
FORBIDDEN = {1004, 2004, 1000, 2000}
WM_L, WM_R = 2, 41            # Left/right cerebral white matter in the segmentation model (outside --parc)

def fail(msg, code=15):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)

# ---------- Grid consistency checks ----------
# The label image must be voxelwise aligned with native QSM. Exit with code 12 if alignment fails,
# rather than silently using a fallback that would invalidate the native-space analysis.
qsm_img = nib.load(qsm_p)
lab_img = nib.load(lab_p)
if lab_img.shape[:3] != qsm_img.shape[:3]:
    fail(f"Label image shape {lab_img.shape[:3]} does not match native QSM {qsm_img.shape[:3]}", 12)
if not np.allclose(lab_img.affine, qsm_img.affine, atol=1e-3):
    d = float(np.abs(lab_img.affine - qsm_img.affine).max())
    fail(f"Label affine does not match native QSM (maximum difference {d:.4g})", 12)
if np.sign(np.linalg.det(lab_img.affine[:3, :3])) != \
        np.sign(np.linalg.det(qsm_img.affine[:3, :3])):
    fail("Label and QSM affine determinants have opposite signs (possible left/right flip)", 12)

qsm = np.asanyarray(qsm_img.dataobj).astype(np.float32)
lab = np.asanyarray(lab_img.dataobj).astype(np.int32)
if qsm.ndim > 3:
    qsm = qsm[..., 0]
if lab.ndim > 3:
    lab = lab[..., 0]

present = set(np.unique(lab).astype(int).tolist()) - {0}
bad = sorted(present & FORBIDDEN)
if bad:
    fmt = Path(here).parent / "Part1_Extract_QSM" / "Desikan_labels.csv"
    fail(f"Unexpected labels {bad} detected (not present in the model label list; "
         f"recheck {fmt} after any FreeSurfer update)", 14)
cortex_present = sorted(present & set(CORTEX_L + CORTEX_R))

vox_ml = float(np.abs(np.linalg.det(qsm_img.affine[:3, :3]))) / 1000.0
nz = qsm != 0

ribbon = np.isin(lab, CORTEX_L + CORTEX_R)
wm = (lab == WM_L) | (lab == WM_R)

# ---------- Volumes ----------
cortex_ml = float(ribbon.sum()) * vox_ml
n_roi_present = len(cortex_present)
counts = [int((lab == k).sum()) for k in cortex_present]
min_nvox = int(min(counts)) if counts else 0

# ---------- Anatomical placement check (native branch-specific, essential QC) ----------
# P1 Dice measures global overlap between nonzero QSM and the BET brain mask. Even a 2-mm
# misalignment may barely change Dice but can disrupt a 2-3-voxel-wide cortical ribbon.
rib_nz = ribbon & nz
wm_nz = wm & nz
rib_med = float(np.median(qsm[rib_nz])) if rib_nz.sum() else float("nan")
wm_med = float(np.median(qsm[wm_nz])) if wm_nz.sum() else float("nan")

# ---------- QSM coverage ----------
ribbon_in_qsm_frac = float(rib_nz.sum() / ribbon.sum()) if ribbon.sum() else 0.0
seg_brain = lab > 0
inter = int((seg_brain & nz).sum())
seg_qsm_dice = (2.0 * inter / float(seg_brain.sum() + nz.sum())
                if (seg_brain.sum() + nz.sum()) else 0.0)

# ---------- SynthSeg-reported volumes (lookup by column name rather than row identifier) ----------
# The row identifier is based on basename.replace(".nii.gz", ""); our T1 images end with .nii,
# so their identifiers retain ".nii" and matching filenames by row ID can fail.
icv_ml = float("nan")
ss_left_cortex_ml = float("nan")
vol_csv = Path(vol_p)
if vol_csv.exists():
    with open(vol_csv, newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh))
    if rows:
        def num(col):
            try:
                return float(rows[0].get(col, "nan")) / 1000.0
            except (TypeError, ValueError):
                return float("nan")
        icv_ml = num("total intracranial")
        ss_left_cortex_ml = num("left cerebral cortex")

# ---------- Compare Dice with archived qc_summary.csv to verify reproducibility of P1 registration ----------
dice_archive = float("nan")
dice_delta = float("nan")
qs = Path(here) / "output" / "qc_summary.csv"
if qs.exists() and dice:
    try:
        import pandas as pd
        t = pd.read_csv(qs, dtype=str).fillna("")
        col = "IID" if "IID" in t.columns else t.columns[0]
        hit = t[t[col] == oid]
        if len(hit):
            dice_archive = float(hit.iloc[0]["reg_mask_dice"])
            dice_delta = abs(float(dice) - dice_archive)
    except Exception as exc:                                       # noqa: BLE001
        print(f"  WARNING: Could not compare Dice with the archive; ignoring {type(exc).__name__}: {exc}")

# ---------- QC thresholds ----------
errors, warns = [], []
try:
    d = float(dice)
    if d < 0.60:
        errors.append(f"reg_mask_dice={d:.3f} < 0.60")
except (TypeError, ValueError):
    warns.append(f"reg_mask_dice unavailable (possibly missing P1 QC; dice={dice!r})")

if not np.isnan(icv_ml) and not (1200.0 <= icv_ml <= 1900.0):
    errors.append(f"SynthSeg ICV={icv_ml:.0f} mL is outside 1200-1900 mL")
if not (350.0 <= cortex_ml <= 700.0):
    errors.append(f"Cortical volume={cortex_ml:.0f} mL is outside 350-700 mL (reference P164054=491 mL)")
if ribbon.sum() and ribbon_in_qsm_frac < 0.50:
    errors.append(f"Only {ribbon_in_qsm_frac:.2f} of the cortical ribbon lies within QSM coverage")
if not np.isnan(rib_med) and not np.isnan(wm_med) and not (rib_med - wm_med) > 0:
    errors.append(f"gm_wm_contrast={rib_med - wm_med:+.2f} ppb <= 0 (possible cortical misalignment)")
if seg_qsm_dice < 0.60:
    errors.append(f"seg_qsm_dice={seg_qsm_dice:.3f} < 0.60 (label/QSM alignment failed)")
if n_roi_present < 68:
    warns.append(f"Only {n_roi_present}/68 cortical ROIs contain voxels")
if not np.isnan(dice_delta) and dice_delta > 0.01:
    warns.append(f"Current Dice={dice} differs from archived Dice={dice_archive:.4f} by {dice_delta:.4f}")
if not np.isnan(ss_left_cortex_ml) and not np.isnan(icv_ml):
    if abs(ss_left_cortex_ml * 2 - cortex_ml) > 0.35 * cortex_ml:
        warns.append(f"SynthSeg-reported left cortex x2={ss_left_cortex_ml * 2:.0f} mL differs from "
                     f"measured {cortex_ml:.0f} mL by >35% (some difference may arise from --keepgeom)")

with open(qc_p, "w", encoding="utf-8") as fh:
    fh.write(f"ID={oid}\n")
    fh.write(f"node={__import__('socket').gethostname()}\n")
    fh.write(f"qsm={qsm_p}\n")
    fh.write(f"qsm_dim={qsm_img.shape[0]}x{qsm_img.shape[1]}x{qsm_img.shape[2]}\n")
    fh.write(f"qsm_pixdim={qsm_pix}\n")
    fh.write(f"qsm_datatype={qsm_dt}\n")
    fh.write(f"t1seg_pixdim={t1_pix}\n")
    fh.write(f"voxel_ml={vox_ml:.6f}\n")
    fh.write(f"reg_mask_dice={dice}\n")
    fh.write(f"dice_archive={dice_archive:.4f}\n" if not np.isnan(dice_archive)
             else "dice_archive=na\n")
    fh.write(f"synthseg_icv_ml={icv_ml:.0f}\n")
    fh.write(f"synthseg_left_cortex_ml={ss_left_cortex_ml:.1f}\n")
    fh.write(f"cortex_vol_ml={cortex_ml:.1f}\n")
    fh.write(f"n_cortex_labels_present={n_roi_present}\n")
    fh.write(f"min_roi_voxels={min_nvox}\n")
    fh.write(f"ribbon_nvox={int(ribbon.sum())}\n")
    fh.write(f"ribbon_in_qsm_frac={ribbon_in_qsm_frac:.4f}\n")
    fh.write(f"ribbon_median_ppb={rib_med:.3f}\n")
    fh.write(f"wm_median_ppb={wm_med:.3f}\n")
    fh.write(f"gm_wm_contrast_ppb={rib_med - wm_med:.3f}\n")
    fh.write(f"seg_qsm_dice={seg_qsm_dice:.4f}\n")
    fh.write(f"seg_brain_nvox={int(seg_brain.sum())}\n")
    fh.write(f"qsm_nonzero_nvox={int(nz.sum())}\n")
    for w in warns:
        fh.write(f"WARNING: {w}\n")
    for e in errors:
        fh.write(f"ERROR: {e}\n")

print(f"  Native QSM grid {qsm_img.shape[:3]}; voxel volume {vox_ml:.4f} mL")
print(f"  Cortical ribbon: {int(ribbon.sum())} voxels = {cortex_ml:.0f} mL; "
       f"within QSM coverage: {ribbon_in_qsm_frac:.1%}")
print(f"  Cortical median {rib_med:+.2f} ppb / white-matter median {wm_med:+.2f} ppb "
       f"-> gm_wm_contrast {rib_med - wm_med:+.2f} ppb")
print(f"  SynthSeg ICV {icv_ml:.0f} mL; seg/QSM Dice {seg_qsm_dice:.3f}; "
       f"reg_mask_dice {dice} (archived {dice_archive})")
print(f"  Cortical ROIs containing voxels: {n_roi_present}/68; minimum {min_nvox} voxels")
for w in warns:
    print(f"  WARNING: {w}")
if errors:
    print("\n".join(f"  ERROR: {e}" for e in errors), file=sys.stderr)
    sys.exit(15)
print("  QC passed")
PYEOF

echo "[done] ${FINAL} (QC: ${QC})"
