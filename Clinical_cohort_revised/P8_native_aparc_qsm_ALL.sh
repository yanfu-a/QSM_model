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
PY310="${P8_PYTHON:-/home1/fuyan/.conda/envs/py310/bin/python}"
# Written to the QC file. Labels whose QC lacks this tag were built with the original
# T1 -> T1_brain grid conversion (no FLIRT x-flip for neurological images) and are
# regenerated; P9 accepts only tagged labels.
GRID_CONVERSION_TAG=flirt_fsl_v2

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
T1SEG_SRC="${T1SEG_DIR}/${OUT_ID}_t1seg_source.txt"
MAT_Q2T="${MAT_DIR}/${OUT_ID}_QSM_to_T1.mat"
MAT_SRC="${MAT_DIR}/${OUT_ID}_QSM_to_T1_source.txt"
T1_BRAIN="${MAT_DIR}/${OUT_ID}_T1_brain.nii.gz"
P1_QC="${P1_DIR}/${OUT_ID}_qc_rigid.txt"
# Grid-conversion and provenance fields of the current build; step [5] copies them into QC.
BUILD_INFO="${MAT_DIR}/${OUT_ID}_build_info.txt"

exec > >(tee -a "${LOG_DIR}/${OUT_ID}.log") 2>&1

echo "=================================================================="
echo "P8 native branch  ID=${OUT_ID}  host=$(hostname)  $(date '+%F %T')"
echo "  QSM = $QSM"
echo "  T1  = $T1"
echo "  native_root = $NATIVE_ROOT"
echo "=================================================================="

# One P8 run per subject at a time. Duplicate worklist rows or resubmitted jobs would
# otherwise rebuild the same files concurrently (P1 also uses one work directory per ID).
if command -v flock >/dev/null 2>&1; then
    exec 9>"${LOG_DIR}/${OUT_ID}.lock"
    if ! flock -n 9; then
        echo "ERROR: Another P8 run for ${OUT_ID} holds ${LOG_DIR}/${OUT_ID}.lock; not starting a concurrent run" >&2
        exit 16
    fi
else
    echo "WARNING: flock not found; concurrent P8 runs for ${OUT_ID} are not prevented" >&2
fi

canon_path() { readlink -f -- "$1" 2>/dev/null || printf '%s\n' "$1"; }
qc_field() { awk -v k="$1" 'index($0, k "=") == 1 {print substr($0, length(k) + 2); exit}' "$2" 2>/dev/null || true; }
file_sha1() { sha1sum -- "$1" | awk '{print $1}'; }

for f in "$QSM" "$T1"; do
    [[ -e "$f" ]] || { echo "ERROR: Missing file $f" >&2; exit 2; }
done
QSM_SHA1=$(file_sha1 "$QSM")
T1_SHA1=$(file_sha1 "$T1")

# An exported FORCE_RERUN=1 still forces P1, as it did when P1 inherited it.
P1_FORCE="${P8_FORCE_P1:-${FORCE_RERUN:-0}}"

# Registration settings that change QSM_to_T1.mat. P1 records both in its rigid QC
# (t1_pixdim=naxnaxna when the T1 is registered at full resolution), so a matrix or label
# made with other settings is never reused, e.g. by a comparison run with QSM_ROI_T1_SUBSAMP=0.
REQ_T1_RES=subsampled
if [[ "${QSM_ROI_T1_SUBSAMP:-2}" == "0" ]]; then REQ_T1_RES=native; fi
REQ_SEARCH="${QSM_ROI_SEARCH_MODE:-nosearch}"
# T1 resolution of the registration recorded in a P1 rigid QC file; empty when not recorded.
p1_t1_resolution() {
    local pix
    pix=$(qc_field t1_pixdim "$1")
    if [[ -z "$pix" ]]; then return 0; fi
    if [[ "$pix" == "naxnaxna" ]]; then echo native; else echo subsampled; fi
}

# Print why the existing label cannot be reused for these inputs; no output means it can.
# A reusable label passed QC with this grid conversion, and the QSM, T1, P1 matrix,
# segmentation and label are byte-identical to the ones recorded when it was built.
existing_output_problems() {
    if [[ "$P1_FORCE" == "1" ]]; then echo "P1 re-registration requested"; fi
    if [[ "${P8_FORCE_SYNTHSEG:-0}" == "1" ]]; then echo "SynthSeg rerun requested"; fi
    if [[ ! -s "$QC" ]]; then
        echo "no QC file"
        return
    fi
    if ! grep -qx "qc_complete=1" "$QC"; then echo "QC file incomplete"; fi
    if ! grep -qx "grid_conversion=${GRID_CONVERSION_TAG}" "$QC"; then
        echo "QC predates grid conversion ${GRID_CONVERSION_TAG}"
    fi
    if [[ "$(qc_field p1_t1_resolution "$QC")" != "$REQ_T1_RES" || \
            "$(qc_field p1_search_mode "$QC")" != "$REQ_SEARCH" ]]; then
        echo "registration settings (T1 resolution, search mode) differ or were not recorded"
    fi
    if grep -q '^ERROR:' "$QC"; then echo "QC recorded errors"; fi
    if [[ "$(canon_path "$(qc_field qsm "$QC")")" != "$(canon_path "$QSM")" ]]; then echo "QSM path differs"; fi
    if [[ "$(canon_path "$(qc_field t1 "$QC")")" != "$(canon_path "$T1")" ]]; then echo "T1 path differs"; fi
    if [[ "$(qc_field qsm_sha1 "$QC")" != "$QSM_SHA1" ]]; then echo "QSM content differs or was not recorded"; fi
    if [[ "$(qc_field t1_sha1 "$QC")" != "$T1_SHA1" ]]; then echo "T1 content differs or was not recorded"; fi
    if [[ ! -s "$MAT_Q2T" || "$(qc_field q2t_mat_sha1 "$QC")" != "$(file_sha1 "$MAT_Q2T")" ]]; then
        echo "QSM_to_T1.mat changed since the label was built"
    fi
    if [[ ! -s "$T1SEG" || "$(qc_field t1seg_sha1 "$QC")" != "$(file_sha1 "$T1SEG")" ]]; then
        echo "T1 segmentation changed since the label was built"
    fi
    if [[ "$(qc_field t1_sha1 "$T1SEG_SRC")" != "$T1_SHA1" || \
            "$(qc_field t1seg_sha1 "$T1SEG_SRC")" != "$(qc_field t1seg_sha1 "$QC")" ]]; then
        echo "T1 segmentation source not verified for this T1"
    fi
    if [[ "$(qc_field label_sha1 "$QC")" != "$(file_sha1 "$FINAL")" ]]; then
        echo "label differs from the one QC checked"
    fi
}

if [[ -s "$FINAL" && "${P8_FORCE_RERUN:-0}" != "1" ]]; then
    mapfile -t _problems < <(existing_output_problems)
    if (( ${#_problems[@]} == 0 )); then
        echo "Existing output verified for these inputs; skipping: $FINAL"
        exit 0
    fi
    # Verified P1 registrations and SynthSeg outputs are reused below, so regenerating
    # usually repeats only steps [3]-[5].
    printf -v _why '%s; ' "${_problems[@]}"
    echo "Existing output not reused (${_why%; }); regenerating"
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

export FREESURFER_HOME="${P8_FREESURFER_HOME:-/public/software/apps/Freesurfer/8.2.0-1}"
export FSLOUTPUTTYPE=NIFTI_GZ
SYNTHSEG="${FREESURFER_HOME}/bin/mri_synthseg"
if [[ ! -x "$SYNTHSEG" ]]; then
    echo "ERROR: Cannot find $SYNTHSEG" >&2
    exit 1
fi

QSM_DT=$(fslinfo "$QSM" | awk '/^datatype/{print $2}')
if [[ "$QSM_DT" == "128" ]]; then
    echo "ERROR: QSM image is not scalar (datatype=128, RGB); cannot extract QSM" >&2
    exit 4
fi

read -r Q1 Q2 Q3 < <(fslinfo "$QSM" | awk '/^pixdim[123]/{printf "%s ", $2} END{print ""}')
QSM_PIX="${Q1}x${Q2}x${Q3}"
read -r P1D P2D P3D < <(fslinfo "$T1" | awk '/^pixdim[123]/{printf "%s ", $2} END{print ""}')
T1_PIX="${P1D}x${P2D}x${P3D}"

P1_SCRIPT="${P8_P1_SCRIPT:-}"
if [[ -z "$P1_SCRIPT" ]]; then
    for _p1 in "${HERE}/P1_qsm_to_mni.sh" "${HERE}/P1_qsm_to_mni_ALL.sh"; do
        if [[ -f "$_p1" ]]; then
            P1_SCRIPT="$_p1"
            break
        fi
    done
fi
if [[ ! -f "$P1_SCRIPT" ]]; then
    echo "ERROR: Cannot find P1_qsm_to_mni[_ALL].sh in ${HERE}; set P8_P1_SCRIPT" >&2
    exit 1
fi
if [[ ! -f "${HERE}/native_grid.py" ]]; then
    echo "ERROR: ${HERE}/native_grid.py is required (grid conversion shared with check_grid_conversion.py)" >&2
    exit 1
fi

# From here on the outputs are regenerated; remove the old QC and label so that a
# failure below cannot leave a stale passing QC file next to a new or missing label.
rm -f "$QC" "$FINAL" "$BUILD_INFO"

# P1 skips rigid registration whenever QSM_to_T1.mat exists (and then writes no QC), and it
# reads FORCE_RERUN, not P8_FORCE_RERUN. Reuse the kept matrix only when it is known to come
# from these inputs: by the source record P8 wrote when it last ran P1 (input hashes, plus
# hashes of the matrix, T1_brain and P1 QC it produced) or, for a matrix from before that
# record existed, by P1's QC paths and inputs no newer than the matrix.
P1_PROVENANCE=fresh
if [[ "$P1_FORCE" != "1" && -s "$MAT_Q2T" ]]; then
    if [[ -s "$MAT_SRC" ]]; then
        if [[ "$(qc_field qsm_sha1 "$MAT_SRC")" != "$QSM_SHA1" || \
                "$(qc_field t1_sha1 "$MAT_SRC")" != "$T1_SHA1" ]]; then
            echo "  Existing $MAT_Q2T was made from different QSM/T1 content; re-running P1"
            P1_FORCE=1
        elif [[ "$(qc_field q2t_mat_sha1 "$MAT_SRC")" != "$(file_sha1 "$MAT_Q2T")" || \
                ! -s "$T1_BRAIN" || "$(qc_field t1_brain_sha1 "$MAT_SRC")" != "$(file_sha1 "$T1_BRAIN")" || \
                ! -s "$P1_QC" || "$(qc_field p1_qc_sha1 "$MAT_SRC")" != "$(file_sha1 "$P1_QC")" ]]; then
            echo "  Matrix, T1_brain or P1 QC changed since P8 recorded them; re-running P1"
            P1_FORCE=1
        else
            P1_PROVENANCE=sha1_record
        fi
    elif [[ ! -s "$P1_QC" ]]; then
        echo "  No P1 QC or source record for existing $MAT_Q2T; re-running P1"
        P1_FORCE=1
    elif [[ "$(canon_path "$(qc_field QSM "$P1_QC")")" != "$(canon_path "$QSM")" || \
            "$(canon_path "$(qc_field T1 "$P1_QC")")" != "$(canon_path "$T1")" ]]; then
        echo "  Existing $MAT_Q2T was made from different QSM/T1 files; re-running P1"
        P1_FORCE=1
    elif [[ "$QSM" -nt "$MAT_Q2T" || "$T1" -nt "$MAT_Q2T" ]]; then
        echo "  QSM/T1 modified after $MAT_Q2T was made; re-running P1"
        P1_FORCE=1
    else
        P1_PROVENANCE=legacy_path_mtime
    fi
fi
if [[ "$P1_FORCE" != "1" && -s "$MAT_Q2T" ]]; then
    _res=$(p1_t1_resolution "$P1_QC")
    _search=$(qc_field search_mode "$P1_QC")
    if [[ "$_res" != "$REQ_T1_RES" || "$_search" != "$REQ_SEARCH" ]]; then
        echo "  Existing $MAT_Q2T was registered with T1 resolution '${_res:-not recorded}' and search" \
             "mode '${_search:-not recorded}' (requested ${REQ_T1_RES}, ${REQ_SEARCH}); re-running P1"
        P1_FORCE=1
    fi
fi
P1_FRESH=0
if [[ "$P1_FORCE" == "1" || ! -s "$MAT_Q2T" ]]; then
    P1_FRESH=1
    P1_PROVENANCE=fresh
    # A failed registration below must not leave the previous record describing new files.
    rm -f "$MAT_SRC"
fi

echo "[1] Run P1 (STOP_AFTER=rigid, FORCE_RERUN=${P1_FORCE}) to obtain QSM_to_T1.mat ..."
# P1's Dice gate stops before FNIRT to save time in MNI mode. In rigid mode it would discard
# the matrix, so a low-Dice scan got no label or QC to inspect. Step [5] applies the same
# 0.60 gate, so such a scan still fails QC (exit 15 rather than 5) and is still excluded.
FORCE_RERUN="$P1_FORCE" \
QSM_ROI_DICE_MIN=0 \
QSM_ROI_STOP_AFTER=rigid \
QSM_ROI_KEEP_WORK=1 \
QSM_ROI_KEEP_DIR="$MAT_DIR" \
QSM_ROI_OUT_DIR="$P1_DIR" \
    bash "$P1_SCRIPT" "$QSM" "$T1" "$OUT_ID"

if [[ ! -s "$MAT_Q2T" ]]; then
    echo "ERROR: P1 did not retain $MAT_Q2T (work-file preservation may have failed)" >&2
    exit 11
fi
if [[ "$P1_FRESH" == "1" ]]; then
    [[ -s "$T1_BRAIN" && -s "$P1_QC" ]] || { echo "ERROR: P1 did not retain T1_brain or its QC" >&2; exit 11; }
    printf 'qsm=%s\nqsm_sha1=%s\nt1=%s\nt1_sha1=%s\nq2t_mat_sha1=%s\nt1_brain_sha1=%s\np1_qc_sha1=%s\n' \
        "$QSM" "$QSM_SHA1" "$T1" "$T1_SHA1" "$(file_sha1 "$MAT_Q2T")" \
        "$(file_sha1 "$T1_BRAIN")" "$(file_sha1 "$P1_QC")" > "${MAT_SRC}.tmp"
    mv -f "${MAT_SRC}.tmp" "$MAT_SRC"
fi
P1_T1_RES=$(p1_t1_resolution "$P1_QC")
P1_SEARCH=$(qc_field search_mode "$P1_QC")
if [[ "$P1_T1_RES" != "$REQ_T1_RES" || "$P1_SEARCH" != "$REQ_SEARCH" ]]; then
    echo "ERROR: P1 QC records T1 resolution '${P1_T1_RES:-none}' and search mode '${P1_SEARCH:-none}';" \
         "requested ${REQ_T1_RES} and ${REQ_SEARCH}" >&2
    exit 11
fi
echo "  P1 matrix provenance: ${P1_PROVENANCE}; T1 resolution: ${P1_T1_RES}; search mode: ${P1_SEARCH}"

DICE=$(qc_field reg_mask_dice "$P1_QC")

# Reuse an existing segmentation only when its source record shows it was made from a T1
# with this content, and it lies on this T1's grid (--keepgeom). Matching dimensions and
# affines alone cannot show that two T1 images are the same scan, so a segmentation without
# a source record is regenerated. P8_FORCE_SYNTHSEG=1 always reruns SynthSeg.
REUSE_T1SEG=0
if [[ -s "$T1SEG" && "${P8_FORCE_SYNTHSEG:-0}" != "1" ]]; then
    if [[ ! -s "$T1SEG_SRC" ]]; then
        echo "[2] Existing segmentation has no source record; regenerating it"
    elif [[ "$(qc_field t1_sha1 "$T1SEG_SRC")" != "$T1_SHA1" ]]; then
        echo "[2] Existing segmentation was made from a T1 with different content; regenerating it"
    elif [[ "$(qc_field t1seg_sha1 "$T1SEG_SRC")" != "$(file_sha1 "$T1SEG")" ]]; then
        echo "[2] Existing segmentation differs from the one its source record describes; regenerating it"
    elif "$PY310" - "$T1" "$T1SEG" <<'PYEOF'
import sys

import nibabel as nib
import numpy as np

t1, seg = nib.load(sys.argv[1]), nib.load(sys.argv[2])
same_grid = t1.shape[:3] == seg.shape[:3] and np.allclose(t1.affine, seg.affine, atol=1e-3)
sys.exit(0 if same_grid else 1)
PYEOF
    then
        REUSE_T1SEG=1
    else
        echo "[2] Existing segmentation is not on this T1's grid; regenerating it"
    fi
fi

if [[ "$REUSE_T1SEG" == "1" ]]; then
    echo "[2] Reusing existing SynthSeg output on the same T1 grid: $T1SEG"
else
    echo "[2] Run mri_synthseg --parc (--keepgeom, ${SYNTHSEG_THREADS} threads) ..."
    # An interrupted run must not leave the previous record describing a partial segmentation.
    rm -f "$T1SEG_SRC"
    if ! "$SYNTHSEG" --i "$T1" --o "$T1SEG" \
            --parc --keepgeom --cpu --threads "$SYNTHSEG_THREADS" \
            --vol "${QC_DIR}/${OUT_ID}_synthseg_vol.csv"; then
        echo "ERROR: mri_synthseg failed" >&2
        exit 13
    fi
    [[ -s "$T1SEG" ]] || { echo "ERROR: SynthSeg did not generate $T1SEG" >&2; exit 13; }
    printf 't1=%s\nt1_sha1=%s\nt1seg_sha1=%s\n' "$T1" "$T1_SHA1" "$(file_sha1 "$T1SEG")" > "${T1SEG_SRC}.tmp"
    mv -f "${T1SEG_SRC}.tmp" "$T1SEG_SRC"
fi


T1_TO_QSM="${MAT_DIR}/${OUT_ID}_T1_to_QSM.mat"
if [[ -s "$T1_BRAIN" ]]; then
    "$PY310" - "$HERE" "$T1" "$T1_BRAIN" "$MAT_Q2T" "$T1_TO_QSM" "$BUILD_INFO" \
            "$GRID_CONVERSION_TAG" <<'PYEOF' || exit 12
import sys

import nibabel as nib
import numpy as np

here, t1_p, brain_p, q2t_p, out_p, info_p, tag = sys.argv[1:8]
sys.path.insert(0, here)
from native_grid import flirt_grid_conversion, is_neurological, legacy_grid_conversion  # noqa: E402

full = nib.load(t1_p)
sub = nib.load(brain_p)

# Grid conversion between the FLIRT coordinates of the full-resolution T1 and of T1_brain,
# through their shared world space. FLIRT coordinates flip the first voxel axis for
# neurologically stored images; the original formula (D_s inv(A_s) A_f inv(D_f)) omitted
# that flip and, for such images, shifted labels one T1 voxel along the first axis.
M_g = flirt_grid_conversion(full, sub)

# The linear component of a pure grid conversion must be the identity matrix. Any rotation or scaling
# implies the grids are not resamplings of the same image (e.g., T1_brain is from another image); abort.
if not np.allclose(M_g[:3, :3], np.eye(3), atol=1e-3):
    print(f"ERROR: Grid conversion between full-resolution T1 and T1_brain is not a pure translation; "
          f"linear component=\n{np.round(M_g[:3, :3], 5)}\n"
          f"  The images are not matching resampled grids; aborting.", file=sys.stderr)
    sys.exit(12)

shift = M_g[:3, 3]
legacy_shift = legacy_grid_conversion(full, sub)[:3, 3]
neuro = int(is_neurological(full))
ratio = np.array(full.shape[:3], float) / np.array(sub.shape[:3], float)
print(f"  [3] Grid conversion M_g (FLIRT coordinates): translation {np.round(shift, 4) + 0.0} mm "
      f"(voxel ratio {np.round(ratio, 4)}; T1 neurological={neuro})")
if not np.allclose(legacy_shift, shift, atol=1e-3):
    print(f"      Note: the original P8 formula gave {np.round(legacy_shift, 4) + 0.0} mm here; labels built "
          "with it were offset.")
# fslmaths -subsamp2 keeps new voxels centred on old ones in FSL's voxel order, so in FLIRT
# coordinates the two grids coincide. A nonzero shift means T1_brain did not come from that.
if not np.allclose(shift, 0, atol=1e-3):
    print("      WARNING: Non-identity grid conversion; expected identity for fslmaths -subsamp2 "
          "(or QSM_ROI_T1_SUBSAMP=0). Verify how T1_brain was generated; nonfatal.")

M_q2t = np.loadtxt(q2t_p)
M_f2q = np.linalg.inv(M_q2t) @ M_g
np.savetxt(out_p, M_f2q, fmt="%.10f")

# Report-only geometry diagnostics. FSL coordinates use pixdim while world coordinates use
# the sform/qform; if their scales disagree, the two geometries differ. A 6-DOF result
# should have unit singular values in mm-to-mm FLIRT coordinates.
def pixdim_sform_rel_diff(img):
    norms = np.linalg.norm(img.affine[:3, :3], axis=0)
    zooms = np.array(img.header.get_zooms()[:3], float)
    return float(np.max(np.abs(norms - zooms) / zooms))

q2t_scale_dev = float(np.max(np.abs(np.linalg.svd(M_q2t[:3, :3], compute_uv=False) - 1.0)))
with open(info_p, "w", encoding="utf-8") as fh:
    fh.write(f"grid_conversion={tag}\n")
    fh.write(f"t1_neurological={neuro}\n")
    fh.write("grid_shift_mm=" + "x".join(f"{round(v, 4) + 0.0:.4f}" for v in shift) + "\n")
    fh.write("legacy_grid_shift_mm=" + "x".join(f"{round(v, 4) + 0.0:.4f}" for v in legacy_shift) + "\n")
    fh.write(f"t1_pixdim_sform_rel_diff={pixdim_sform_rel_diff(full):.6f}\n")
    fh.write(f"q2t_scale_dev={q2t_scale_dev:.6f}\n")
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

{
    printf 'qsm_sha1=%s\nt1_sha1=%s\n' "$QSM_SHA1" "$T1_SHA1"
    printf 'p1_matrix_provenance=%s\n' "$P1_PROVENANCE"
    printf 'p1_t1_resolution=%s\np1_search_mode=%s\n' "$P1_T1_RES" "$P1_SEARCH"
    printf 'q2t_mat_sha1=%s\n' "$(file_sha1 "$MAT_Q2T")"
    printf 't1seg_sha1=%s\n' "$(file_sha1 "$T1SEG")"
    printf 'label_sha1=%s\n' "$(file_sha1 "$FINAL")"
} >> "$BUILD_INFO"

echo "[5] QC ..."
"$PY310" - "$QSM" "$FINAL" "$T1SEG" "$QC" "${QC_DIR}/${OUT_ID}_synthseg_vol.csv" \
            "$OUT_ID" "$DICE" "$QSM_PIX" "$T1_PIX" "$QSM_DT" "$HERE" "$T1" "$BUILD_INFO" <<'PYEOF'
import csv
import os
import sys
from pathlib import Path

import nibabel as nib
import numpy as np

(qsm_p, lab_p, t1seg_p, qc_p, vol_p,
 oid, dice, qsm_pix, t1_pix, qsm_dt, here, t1_p, build_info_p) = sys.argv[1:14]

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
# NaN background (some reconstructions) is not coverage; NaN != 0 would count it as covered.
nz = np.isfinite(qsm) & (qsm != 0)

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
ss_cortex_ml = float("nan")
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
        # Whole-brain cortex of the T1 segmentation; unlike cortex_ml it does not depend on
        # how much of the brain the QSM slab covers.
        ss_cortex_ml = ss_left_cortex_ml + num("right cerebral cortex")

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
# Each failed rule is also recorded by name (qc_failed_rules) so that exclusions can be
# counted per rule (qc_exclusion_report.py) without parsing the messages.
errors, warns, failed_rules = [], [], []

def fail_rule(rule, msg):
    failed_rules.append(rule)
    errors.append(msg)

try:
    d = float(dice)
    if d < 0.60:
        fail_rule("reg_mask_dice", f"reg_mask_dice={d:.3f} < 0.60")
except (TypeError, ValueError):
    # Step [1] verifies or regenerates P1's QC, so a missing Dice means the registration is unchecked.
    fail_rule("reg_mask_dice", f"reg_mask_dice unavailable (P1 QC missing or incomplete; dice={dice!r})")

if not np.isnan(icv_ml) and not (1200.0 <= icv_ml <= 1900.0):
    fail_rule("synthseg_icv", f"SynthSeg ICV={icv_ml:.0f} mL is outside 1200-1900 mL")
if not (350.0 <= cortex_ml <= 700.0):
    fail_rule("cortex_vol_qsm_grid",
              f"Cortical volume={cortex_ml:.0f} mL is outside 350-700 mL (reference P164054=491 mL)")
if ribbon.sum() and ribbon_in_qsm_frac < 0.50:
    fail_rule("ribbon_coverage", f"Only {ribbon_in_qsm_frac:.2f} of the cortical ribbon lies within QSM coverage")
if not np.isnan(rib_med) and not np.isnan(wm_med) and not (rib_med - wm_med) > 0:
    fail_rule("gm_wm_contrast", f"gm_wm_contrast={rib_med - wm_med:+.2f} ppb <= 0 (possible cortical misalignment)")
if seg_qsm_dice < 0.60:
    fail_rule("seg_qsm_dice", f"seg_qsm_dice={seg_qsm_dice:.3f} < 0.60 (label/QSM alignment failed)")
# Evaluated but not applied (METHODOLOGICAL_DECISIONS D1): the same 350-700 mL range on the
# whole-brain T1 cortex, the candidate replacement for the slab-dependent cortex_vol_qsm_grid.
report_only_failed = []
if not np.isnan(ss_cortex_ml) and not (350.0 <= ss_cortex_ml <= 700.0):
    report_only_failed.append("cortex_vol_t1")
if n_roi_present < 68:
    warns.append(f"Only {n_roi_present}/68 cortical ROIs contain voxels")
if not np.isnan(dice_delta) and dice_delta > 0.01:
    warns.append(f"Current Dice={dice} differs from archived Dice={dice_archive:.4f} by {dice_delta:.4f}")
if not np.isnan(ss_left_cortex_ml) and not np.isnan(icv_ml):
    if abs(ss_left_cortex_ml * 2 - cortex_ml) > 0.35 * cortex_ml:
        warns.append(f"SynthSeg-reported left cortex x2={ss_left_cortex_ml * 2:.0f} mL differs from "
                     f"measured {cortex_ml:.0f} mL by >35% (some difference may arise from --keepgeom)")

# Report-only geometry diagnostics (warnings; they do not change the QC result).
build = dict(line.split("=", 1) for line in Path(build_info_p).read_text(encoding="utf-8").splitlines()
             if "=" in line)
qsm_norms = np.linalg.norm(qsm_img.affine[:3, :3], axis=0)
qsm_zooms = np.array(qsm_img.header.get_zooms()[:3], float)
qsm_pixdim_sform_rel_diff = float(np.max(np.abs(qsm_norms - qsm_zooms) / qsm_zooms))
for name, value in (("QSM", qsm_pixdim_sform_rel_diff),
                    ("T1", float(build.get("t1_pixdim_sform_rel_diff", "nan")))):
    if value > 0.01:
        warns.append(f"{name} pixdim and sform/qform voxel sizes differ by {value:.1%}; FSL and world "
                     "geometry disagree")
if float(build.get("q2t_scale_dev", "nan")) > 0.01:
    warns.append(f"QSM_to_T1.mat is not rigid (max singular-value deviation {build['q2t_scale_dev']})")

# Written to a temporary file and renamed, so an interrupted run never leaves a partial QC
# file; qc_complete=1 is the last line, and P8 and P9 require it.
tmp_qc = qc_p + ".tmp"
with open(tmp_qc, "w", encoding="utf-8") as fh:
    fh.write(f"ID={oid}\n")
    fh.write(f"node={__import__('socket').gethostname()}\n")
    fh.write(f"qsm={qsm_p}\n")
    fh.write(f"t1={t1_p}\n")
    # Grid-conversion tag and shifts, and content hashes of the inputs, P1 matrix,
    # segmentation and label; P8 checks them before reusing this label and P9 before using it.
    fh.write(Path(build_info_p).read_text(encoding="utf-8"))
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
    fh.write(f"synthseg_cortex_ml={ss_cortex_ml:.1f}\n")
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
    fh.write(f"qsm_pixdim_sform_rel_diff={qsm_pixdim_sform_rel_diff:.6f}\n")
    fh.write(f"qsm_scl_slope={qsm_img.header['scl_slope']}\n")
    fh.write(f"qsm_scl_inter={qsm_img.header['scl_inter']}\n")
    fh.write(f"qc_failed_rules={','.join(failed_rules)}\n")
    fh.write(f"qc_report_only_failed_rules={','.join(report_only_failed)}\n")
    for w in warns:
        fh.write(f"WARNING: {w}\n")
    for e in errors:
        fh.write(f"ERROR: {e}\n")
    fh.write("qc_complete=1\n")
os.replace(tmp_qc, qc_p)

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
