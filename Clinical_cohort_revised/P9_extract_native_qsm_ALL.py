#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Extract individualized cortical QSM statistics on the native QSM grid.

Required English worklist headers: image_dir_id, qsm_file, t1_file.
Optional English worklist headers: runnable (yes/no), hospital_id, model_label,
include_main_model, category, qsm_coverage_mm.
Optional ID-list headers: IID or image_dir_id.
The --which filter uses the exact English column names in that ID list.

Input files must be migrated to English headers before running this script.
This script does not translate existing input files or legacy downstream schemas.

A subject enters the primary (unreferenced) analysis set (analysis_pass in *_qc.csv,
rows of *_analysis.csv) only if all of the following hold:
  * extraction succeeded with at least --min-reportable-rois reported ROI values;
  * it has no P1 failure in a stage the native branch shares (input or QSM->T1 rigid
    registration); failures confined to T1->MNI registration do not exclude it;
  * the P8 QC file is complete and passed with a numeric registration Dice, and the
    label carries the corrected grid-conversion tag (legacy labels are never accepted);
  * the QSM and T1 (worklist qsm_file and t1_file) and the label and T1 segmentation read
    here are byte-identical (SHA-1) to the ones P8 recorded;
  * in support-mask mode, the support mask's provenance record verifies against the
    current QSM and the mask file.
CSF referencing (--csf-ref) has its own eligibility flag and never changes analysis_pass.
Every runnable worklist row is processed unless --main-model-only restricts the run to
include_main_model=yes; include_main_model is carried into every output table.
An ROI value is reported only when the ROI has enough valid voxels, valid volume (mm3)
and coverage of a finite, positive T1-space volume. Worklist IDs whose repeated rows
disagree (e.g. two hospital_ids, diagnoses or input files) are excluded; every excluded
subject is listed with its reason in *_excluded.csv, and *_run_info.json records the
command, options, software versions and input hashes of the run.

Formal extraction reads QSM support masks (--support-dir; see make_support_masks.py).
The nonzero-QSM proxy (--allow-nonzero-proxy) is for diagnostics and sensitivity
analyses only: its output file names carry _proxy and *_qc.csv records the mode.
"""

import argparse
import csv
import datetime
import hashlib
import json
import os
import platform
import re
import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from scipy.ndimage import binary_erosion

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from strata import coverage_group, si_extent_mm, voxel_size_label  # noqa: E402
from worklist_ids import resolve_duplicate_ids  # noqa: E402
ATLAS_DIR = Path("/cwStorage/home/fuyan/FY_data/QSM_STAAR_251229/Session1_Pheno_extract/"
                 "Part1_Extract_QSM")
DK_LABELS = ATLAS_DIR / "Desikan_labels.csv"
# FreeSurfer used by P8 (mri_synthseg); its LUT and --parc label list are checked here.
DEFAULT_FREESURFER_HOME = Path("/public/software/apps/Freesurfer/8.2.0-1")

DEFAULT_NATIVE_DIR = Path("/cwStorage/nodecw_group/FY_data/QSM_HUASHAN/native_data")
NATIVE_ROOT = DEFAULT_NATIVE_DIR
# P1's output directory in its default (MNI) mode; holds <IID>_qc.txt for each P1 run.
DEFAULT_P1_QC_DIR = Path("/cwStorage/nodecw_group/FY_data/QSM_HUASHAN/roi_work")

SUPPORT_TYPES = ("reconstruction", "estimated_from_qsm")

# Warnings P1 writes to <IID>_qc.txt. A failure involving the QSM input or the QSM->T1 rigid
# registration also concerns the native branch; the others are confined to T1->MNI
# registration (FNIRT) or the MNI-space output.
P1_NATIVE_WARNINGS = ("QSM datatype", "reg_mask_dice")
P1_MNI_WARNINGS = ("FNIRT", "t1_mni_dice", "mni_brain_coverage", "Output dimensions")

NON_CORTICAL = ("white_matter", "corpus")

META_COLS = ("hospital_id", "model_label", "include_main_model", "category")

# Written to the P8 QC file by the revised P8. Labels without it were built with the original
# T1 -> T1_brain grid conversion, which offsets neurologically stored T1 images.
GRID_CONVERSION_TAG = "flirt_fsl_v2"

CSF_LABELS = (4, 43)  # SynthSeg left/right lateral ventricles

# Cortical QSM spreads over tens of ppb; an interquartile range below this suggests ppm.
PPM_IQR_SUSPECT = 1.0

# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
FS_LONG_L = [
    "ctx-lh-bankssts", "ctx-lh-caudalanteriorcingulate", "ctx-lh-caudalmiddlefrontal",
    "ctx-lh-cuneus", "ctx-lh-entorhinal", "ctx-lh-fusiform", "ctx-lh-inferiorparietal",
    "ctx-lh-inferiortemporal", "ctx-lh-isthmuscingulate", "ctx-lh-lateraloccipital",
    "ctx-lh-lateralorbitofrontal", "ctx-lh-lingual", "ctx-lh-medialorbitofrontal",
    "ctx-lh-middletemporal", "ctx-lh-parahippocampal", "ctx-lh-paracentral",
    "ctx-lh-parsopercularis", "ctx-lh-parsorbitalis", "ctx-lh-parstriangularis",
    "ctx-lh-pericalcarine", "ctx-lh-postcentral", "ctx-lh-posteriorcingulate",
    "ctx-lh-precentral", "ctx-lh-precuneus", "ctx-lh-rostralanteriorcingulate",
    "ctx-lh-rostralmiddlefrontal", "ctx-lh-superiorfrontal", "ctx-lh-superiorparietal",
    "ctx-lh-superiortemporal", "ctx-lh-supramarginal", "ctx-lh-frontalpole",
    "ctx-lh-temporalpole", "ctx-lh-transversetemporal", "ctx-lh-insula",
]
# The label order strictly matches FS_LONG_L:1001,1002,1003,1005..1035
FS_L = [1001, 1002, 1003] + list(range(1005, 1036))
FS_R = [k + 1000 for k in FS_L]

INSULA_L, INSULA_R = 1035, 2035


def read_labels(path):
    """Read Desikan_labels.csv into {label ID: name}, using the same parser as P5."""
    out = {}
    with open(path, newline="", encoding="utf-8-sig") as handle:
        for row in csv.reader(handle):
            if len(row) >= 2 and row[0].strip().isdigit():
                out[int(row[0])] = row[1].strip()
    if not out:
        raise SystemExit(f"{path}  contains no parsable labels.")
    return out


def read_failed_ids(path):
    """Read the P1 batch failure list as {ID: line} (supports literal backslash-t)."""
    if path is None:
        return {}
    if not Path(path).exists():
        print(f"WARNING: Failure list {path} not found; no subject is excluded as a P1 failure")
        return {}
    ids = {}
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            text = line.strip().replace("\\t", "\t")
            tokens = text.split("\t")[0].split()
            if tokens:
                ids.setdefault(tokens[0], text)
    return ids


def classify_p1_failure(qc_path):
    """Stage of a P1 (MNI-branch) failure from P1's QC file: (stage, evidence).

    'native' when a warning concerns the QSM input or the QSM->T1 rigid registration that
    the native branch shares; 'mni_only' when every recognized warning concerns T1->MNI
    registration or the MNI output; 'unknown' without a QC file or recognized warning.
    """
    path = Path(qc_path)
    if not path.is_file():
        return "unknown", f"no P1 QC file ({path.name})"
    warnings = [line.split(":", 1)[1].strip()
                for line in path.read_text(encoding="utf-8", errors="replace").splitlines()
                if line.lstrip().upper().startswith("WARNING:")]
    native = [w for w in warnings if any(m in w for m in P1_NATIVE_WARNINGS)]
    mni = [w for w in warnings if any(m in w for m in P1_MNI_WARNINGS)]
    if native:
        return "native", "; ".join(native)
    if mni:
        return "mni_only", "; ".join(mni)
    if not warnings:
        return "unknown", "P1 QC records no failure (failure list may be outdated)"
    return "unknown", "unrecognized P1 warnings: " + "; ".join(warnings)


def read_native_qc(path):
    """Read the P8 QC summary as (status, errors, fields) without hiding extraction results.

    WARNING lines are skipped rather than parsed as key=value fields.
    """
    path = Path(path)
    if not path.is_file():
        return "missing", "", {}
    fields, errors = {}, []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        head = line.lstrip().upper()
        if head.startswith("ERROR:"):
            errors.append(line.split(":", 1)[1].strip())
        elif head.startswith("WARNING:"):
            continue
        elif "=" in line:
            key, value = line.split("=", 1)
            fields[key.strip()] = value.strip()
        elif any(ord(c) > 127 for c in line) and line.strip():
            errors.append("Unrecognized non-English P8 QC line; use ERROR: for error records")
    return ("fail" if errors else "pass"), "; ".join(errors), fields


def normalize_name(name):
    """Lower-case alphanumerics only, to compare region names across naming styles."""
    return re.sub(r"[^a-z0-9]", "", name.lower())


def build_label_map(csv_labels, labels_path, fs_home, allow_name_mismatch=False):
    """Build {FreeSurfer label ID: output column name} and validate label IDs and names.

    CSV entries are assigned to FreeSurfer labels by ID arithmetic and then
    checked by name against the FreeSurfer LUT, so a CSV listing regions in a
    different order cannot silently relabel columns.
    Return (fs_to_name, mapping_rows).
    """
    cortical = {k: v for k, v in csv_labels.items()
                if not any(s in v for s in NON_CORTICAL)}
    left = {k: v for k, v in cortical.items() if v.startswith("L_")}
    right = {k: v for k, v in cortical.items() if v.startswith("R_")}

    if len(left) != 33 or len(right) != 33:
        raise SystemExit(f"{labels_path.name}  does not contain exactly 33 cortical regions per hemisphere"
                         f" (left {len(left)} / right {len(right)}); stopping.")
    if FS_L[-1] != INSULA_L or FS_R[-1] != INSULA_R:
        raise SystemExit("The final FS_L/FS_R labels are not insula; mapping constants are inconsistent.")

    fs_lut = fs_home / "FreeSurferColorLUT.txt"
    parc_labels = fs_home / "models" / "synthseg_parcellation_labels.npy"
    lut = {}
    with open(fs_lut, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.split()
            if len(parts) >= 2 and parts[0].isdigit():
                lut[int(parts[0])] = parts[1]
    for fs, expect in zip(FS_L, FS_LONG_L):
        got = lut.get(fs)
        if got != expect:
            raise SystemExit(
                f"FreeSurfer LUT label {fs} is '{got}'; expected '{expect}'. "
                f"Label IDs or names have changed; stopping ({fs_lut}).")
    for fs in FS_R:
        if lut.get(fs) != lut[fs - 1000].replace("-lh-", "-rh-"):
            raise SystemExit(f"FreeSurfer LUT label {fs} and {fs - 1000} are not symmetric; stopping.")

    if parc_labels.exists():
        model = [int(x) for x in np.load(parc_labels)
                 if 1000 <= int(x) < 3000]
        want = sorted(FS_L + FS_R)
        if sorted(model) != want:
            extra = sorted(set(model) - set(want))
            miss = sorted(set(want) - set(model))
            raise SystemExit(
                f"synthseg --parc model labels do not match the expected labels: "
                f" (extra {extra}; missing {miss}); was FreeSurfer updated?"
                f"; model file: {parc_labels}")
    else:
        print(f"WARNING: {parc_labels} not found; SynthSeg --parc label list not checked")

    fs_to_name, rows = {}, []
    for fs_ids, side, offset in ((FS_L[:-1], left, 999), (FS_R[:-1], right, 1964)):
        for fs, csv_idx in zip(fs_ids, sorted(side)):
            if csv_idx != fs - offset:
                raise SystemExit(f"Label ID mapping mismatch: FS {fs} -> csv {csv_idx}"
                                 f" (expected {fs - offset}); stopping.")
            csv_name = side[csv_idx]
            # "L_bankssts" -> "bankssts"; "ctx-lh-bankssts" -> "bankssts"
            match = normalize_name(csv_name[2:]) == normalize_name(lut[fs].split("-", 2)[2])
            fs_to_name[fs] = csv_name
            rows.append({"fs_label": fs, "fs_name": lut[fs], "csv_index": csv_idx,
                         "csv_name": csv_name, "name_match": "yes" if match else "NO"})

    for fs, nm in ((INSULA_L, "L_insula"), (INSULA_R, "R_insula")):
        if fs in fs_to_name:
            raise SystemExit(f"Insula {fs} was already mapped to {fs_to_name[fs]}; conflicting mappings.")
        fs_to_name[fs] = nm
        rows.append({"fs_label": fs, "fs_name": lut[fs], "csv_index": "", "csv_name": nm,
                     "name_match": "added (no CSV entry)"})

    repeated = sorted(n for n, c in Counter(fs_to_name.values()).items() if c > 1)
    if repeated:
        raise SystemExit(f"Output region names are not unique: {repeated}; check {labels_path}.")

    bad = [r for r in rows if r["name_match"] == "NO"]
    if bad:
        print(f"\n{'FS label':>8}  {'FreeSurfer name':<34} {'CSV':>4}  CSV name")
        for r in bad:
            print(f"{r['fs_label']:>8}  {r['fs_name']:<34} {r['csv_index']:>4}  {r['csv_name']}")
        msg = (f"{len(bad)} region names in {labels_path.name} do not match the FreeSurfer region "
               "assigned to them by label ID")
        if not allow_name_mismatch:
            raise SystemExit(msg + "; output columns would be mislabelled. Stopping (use "
                             "--allow-label-name-mismatch only after checking the table above).")
        print(f"WARNING: {msg}; continuing because --allow-label-name-mismatch was given.")
    rows.sort(key=lambda r: r["fs_label"])
    return fs_to_name, rows


def file_digest(path, chunk=1 << 20):
    """SHA-1 of a file's contents."""
    digest = hashlib.sha1()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def same_path(a, b):
    """True when two paths resolve to the same location (says nothing about content)."""
    return bool(a) and bool(b) and os.path.realpath(a) == os.path.realpath(b)


def write_csv(df, path):
    """Write a CSV atomically: a temporary file in the same folder, then a rename."""
    tmp = path.with_name(f".{path.name}.tmp")
    df.to_csv(tmp, index=False, encoding="utf-8-sig")
    os.replace(tmp, path)


def support_record_path(mask_path):
    """<IID>_support.json next to <IID>_support.nii[.gz] (written by make_support_masks.py)."""
    name = Path(mask_path).name
    stem = name[:-7] if name.endswith(".nii.gz") else name[:-4] if name.endswith(".nii") else name
    return Path(mask_path).with_name(stem + ".json")


def read_support_record(mask_path):
    """The support mask's provenance record, {} when absent, or None when unreadable."""
    path = support_record_path(mask_path)
    if not path.is_file():
        return {}
    try:
        record = json.loads(path.read_text(encoding="utf-8"))
    except ValueError:
        return None
    return record if isinstance(record, dict) else None


def to_float(value):
    """float(value), or NaN when it is missing or not numeric."""
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def first_volume(arr):
    """3D view of a possibly 4D (x, y, z, 1) array."""
    return arr[..., 0] if arr.ndim > 3 else arr


def check_same_grid(img, ref, what):
    """Raise ValueError unless img lies on the grid of the QSM image ref."""
    if img.shape[:3] != ref.shape[:3]:
        raise ValueError(f"{what} shape {img.shape[:3]} does not match QSM {ref.shape[:3]}")
    if not np.allclose(img.affine, ref.affine, atol=1e-3):
        raise ValueError(f"{what} affine does not match QSM")


def load_native_grid(label_path, qsm_path, support_path=None):
    """Load label, QSM and valid-voxel mask on the native QSM grid.

    With a support mask, its nonzero voxels define valid QSM coverage, including
    legitimate QSM values equal to zero. Without one, nonzero QSM is used as a
    provisional coverage proxy. Non-finite QSM voxels are never valid.
    """
    lab_img, qsm_img = nib.load(str(label_path)), nib.load(str(qsm_path))
    check_same_grid(lab_img, qsm_img, "Label image")
    lab = first_volume(np.asanyarray(lab_img.dataobj)).astype(np.int32)
    qsm = first_volume(np.asanyarray(qsm_img.dataobj)).astype(np.float32)
    if support_path is not None:
        support_img = nib.load(str(support_path))
        check_same_grid(support_img, qsm_img, "Support mask")
        valid = first_volume(np.asanyarray(support_img.dataobj)) > 0
    else:
        # Provisional only: cannot distinguish a valid zero-valued voxel from
        # zero padding without an explicit reconstruction support mask.
        valid = qsm != 0
    valid &= np.isfinite(qsm)
    zooms = tuple(float(z) for z in qsm_img.header.get_zooms()[:3])
    vox_mm3 = float(abs(np.linalg.det(qsm_img.affine[:3, :3])))
    return lab, qsm, valid, zooms, vox_mm3


def t1_roi_volumes(seg_path, keys, slab=16):
    """Volume (mm3) of each label in the T1-space SynthSeg segmentation, counted slab by slab."""
    img = nib.load(str(seg_path))
    seg = first_volume(np.asanyarray(img.dataobj))
    if not np.issubdtype(seg.dtype, np.integer):
        seg = np.rint(seg)
    top = int(max(keys)) + 1  # labels above the cortical range share one overflow bin
    counts = np.zeros(top + 1, dtype=np.int64)
    for z0 in range(0, seg.shape[2], slab):
        chunk = np.clip(seg[:, :, z0:z0 + slab].astype(np.int64).ravel(), 0, top)
        counts += np.bincount(chunk, minlength=top + 1)
    return counts[np.asarray(keys)] * float(abs(np.linalg.det(img.affine[:3, :3])))


def extract_one(lab, qsm, valid, keys, vox_mm3, t1_mm3, min_voxels, min_mm3, min_coverage):
    """Per-ROI statistics, counts and coverage for one subject on the native grid.

    An ROI value is reported only with at least min_voxels valid voxels, at least
    min_mm3 of valid volume and, if min_coverage > 0, a finite, positive T1-space
    volume of which the valid volume covers at least min_coverage. Coverage
    relative to T1 also counts cortex outside the QSM field of view, which the
    label on the QSM grid cannot show. The status is OK only if at least one ROI
    value survives these gates and the label agrees with the T1 segmentation.
    """
    n_keys = len(keys)
    lut = np.full(int(max(keys)) + 1, -1, dtype=np.int32)
    lut[np.asarray(keys)] = np.arange(n_keys, dtype=np.int32)

    lab_f = lab.ravel()
    slot = lut[np.clip(lab_f, 0, len(lut) - 1)]
    slot[lab_f > len(lut) - 1] = -1
    keep = slot >= 0
    s_all = slot[keep]
    ok = valid.ravel()[keep]
    v_all = qsm.ravel()[keep]
    s, v = s_all[ok], v_all[ok]

    label_counts = np.bincount(s_all, minlength=n_keys)
    counts = np.bincount(s, minlength=n_keys)
    # Voxels with QSM exactly 0: counted as genuine values inside the support, and left out
    # as padding outside it. For quantized QSM an estimated support cannot tell the two
    # apart at the brain edge (METHODOLOGICAL_DECISIONS D3), so both are reported.
    zero_counts = np.bincount(s[v == 0], minlength=n_keys)
    unsupported_zero_counts = np.bincount(s_all[~ok & (v_all == 0)], minlength=n_keys)
    means = np.full(n_keys, np.nan)
    medians = np.full(n_keys, np.nan)
    if v.size:
        sums = np.bincount(s, weights=v.astype(np.float64), minlength=n_keys)
        order = np.argsort(s, kind="stable")
        s_sorted, v_sorted = s[order], v[order]
        bounds = np.searchsorted(s_sorted, np.arange(n_keys + 1))
        for i in np.flatnonzero(counts):
            means[i] = sums[i] / counts[i]
            medians[i] = float(np.median(v_sorted[bounds[i]:bounds[i + 1]]))

    valid_mm3 = counts * vox_mm3
    # Coverage is defined only for a finite, positive T1 volume; otherwise it is NaN
    # (never inf), and NaN fails the gate whenever the gate is on.
    t1_ok = np.isfinite(t1_mm3) & (t1_mm3 > 0)
    coverage = np.full(n_keys, np.nan)
    coverage[t1_ok] = valid_mm3[t1_ok] / t1_mm3[t1_ok]
    low_count = (counts < min_voxels) | (valid_mm3 < min_mm3)
    low_coverage = ~low_count & (min_coverage > 0) & ~(np.isfinite(coverage) & (coverage >= min_coverage))
    means[low_count | low_coverage] = np.nan
    medians[low_count | low_coverage] = np.nan
    n_reportable = int(np.isfinite(medians).sum())

    # The label is resampled from the T1 segmentation, so with that segmentation in hand
    # (any finite volume; all NaN means none was read) every ROI present on the QSM grid
    # must have a finite, positive T1 volume; otherwise they disagree.
    t1_known = bool(np.isfinite(t1_mm3).any())
    n_label_without_t1 = int(((label_counts > 0) & ~t1_ok).sum()) if t1_known else 0

    if not label_counts.any():
        status = "No cortical labels"
    elif n_label_without_t1:
        status = (f"Label and T1 segmentation disagree ({n_label_without_t1} ROIs on the QSM "
                  "grid have no T1 volume)")
    elif not counts.any():
        status = "No valid supported cortical voxels"
    elif not n_reportable:
        status = "No ROI value passes the voxel, volume and coverage gates"
    else:
        status = "OK"
    iqr = float(np.subtract(*np.percentile(v, [75, 25]))) if v.size else float("nan")
    return {"status": status, "label_counts": label_counts, "counts": counts,
            "zero_counts": zero_counts, "unsupported_zero_counts": unsupported_zero_counts,
            "means": means, "medians": medians, "coverage": coverage,
            "n_low_voxels": int(low_count.sum()), "n_low_coverage": int(low_coverage.sum()),
            "n_reportable": n_reportable, "cortex_iqr": iqr}


def native_csf_reference(lab, qsm, valid, min_voxels=30):
    """Median QSM of one-voxel-eroded bilateral lateral ventricles (aseg 4/43)."""
    csf_mask = np.zeros(lab.shape, dtype=bool)
    for k in CSF_LABELS:
        csf_mask |= binary_erosion(lab == k)
    values = qsm[csf_mask & valid]
    return (float(np.median(values)) if values.size >= min_voxels else np.nan,
            int(values.size))


def process_subject(task, keys, opts):
    """Load one subject once; return ROI statistics, T1-space coverage and the CSF reference.

    Content hashes of every input read (QSM, T1, label, T1 segmentation, support mask) are
    returned for comparison with the provenance P8 and make_support_masks.py recorded.
    Any error becomes the returned status, so one unreadable file cannot stop the run.
    """
    try:
        lab, qsm, valid, zooms, vox_mm3 = load_native_grid(task["label"], task["qsm"],
                                                           task["support"])
        t1_mm3 = (t1_roi_volumes(task["t1seg"], keys) if task["t1seg"] is not None
                  else np.full(len(keys), np.nan))
        res = extract_one(lab, qsm, valid, keys, vox_mm3, t1_mm3, opts["min_voxels"],
                          opts["min_mm3"], opts["min_coverage"])
        qsm_img = nib.load(str(task["qsm"]))
        hdr = qsm_img.header
        nonzero = np.isfinite(qsm) & (qsm != 0)
        integer_valued = bool(nonzero.any() and np.all(np.mod(qsm[nonzero], 1) == 0))
        res.update(t1_mm3=t1_mm3, vox_mm3=vox_mm3, pixdim="x".join(f"{z:.4g}" for z in zooms),
                   qsm_sha1=file_digest(task["qsm"]),
                   t1_sha1=file_digest(task["t1"]),
                   label_sha1=file_digest(task["label"]),
                   t1seg_sha1=file_digest(task["t1seg"]) if task["t1seg"] is not None else "",
                   # Reported, never used to rescale: on-disk type, NIfTI scaling, integer values.
                   qsm_dtype=str(hdr.get_data_dtype()),
                   qsm_scl_slope=float(hdr["scl_slope"]), qsm_scl_inter=float(hdr["scl_inter"]),
                   qsm_integer_valued=integer_valued,
                   # Quantized: stored as integers (any scl_slope) or integer-valued after scaling.
                   # Genuine zeros then occur, and an estimated support misses those at the edge.
                   qsm_quantized=bool(np.issubdtype(hdr.get_data_dtype(), np.integer) or integer_valued),
                   si_extent_mm=si_extent_mm(qsm_img))
        if task["support"] is not None:
            res.update(support_sha1=file_digest(task["support"]),
                       support_record=read_support_record(task["support"]),
                       qsm_nonzero_outside_support=int((nonzero & ~valid).sum()))
        if opts["csf_ref"]:
            res["csf_ref"], res["csf_voxels"] = native_csf_reference(lab, qsm, valid,
                                                                     opts["min_csf_voxels"])
        return res
    except Exception as exc:                                       # noqa: BLE001
        return {"status": f"ERROR {type(exc).__name__}: {exc}"}


def lr_concordance(out, names, label):
    """Cross-subject left/right homologous ROI correlations and their Spearman-Brown projection.

    This is anatomical left/right concordance (as in P5), not test-retest reliability:
    homologous correlations include variance shared by every ROI (QSM reference offset,
    age, global iron), so the mean correlation of non-homologous left/right pairs is
    printed as a baseline; only the excess is region-specific. test_retest() uses repeat
    scans for an actual (if small) reliability estimate.
    """
    pairs = [(c, "R_" + c[2:]) for c in names
             if c.startswith("L_") and "R_" + c[2:] in names]
    if not pairs:
        return None
    left = [lc for lc, _ in pairs]
    right = [rc for _, rc in pairs]
    data = out[left + right].apply(pd.to_numeric, errors="coerce")
    block = data.corr(min_periods=31).loc[left, right].to_numpy()
    detail = [(lc[2:], float(r)) for lc, r in zip(left, np.diag(block)) if np.isfinite(r)]
    if not detail:
        return None
    rbar = float(np.mean([r for _, r in detail]))
    sb = 2 * rbar / (1 + rbar) if rbar > -1 else float("nan")
    off = block[~np.eye(len(pairs), dtype=bool)]
    off = off[np.isfinite(off)]
    base = float(np.mean(off)) if off.size else float("nan")
    print(f"\nLeft/right concordance ({label}, {len(detail)} homologous ROI pairs; not test-retest "
          f"reliability): mean r={rbar:.3f}  Spearman-Brown projection={sb:.3f}")
    print(f"  Non-homologous left/right baseline mean r={base:.3f}; homologous excess={rbar - base:+.3f}")
    print("  Five lowest-correlated pairs:",
          ", ".join(f"{n}={r:.2f}" for n, r in sorted(detail, key=lambda t: t[1])[:5]))
    return rbar, sb, base, detail


def test_retest(out, names, label, min_pairs=5):
    """Per-ROI correlation between the first two QC-passing scans of patients scanned twice."""
    if "hospital_id" not in out.columns:
        return None
    hosp = out["hospital_id"].astype(str)
    repeated = out[hosp.ne("") & hosp.duplicated(keep=False)].sort_values(["hospital_id", "IID"])
    pairs = [(g.iloc[0], g.iloc[1]) for _, g in repeated.groupby("hospital_id") if len(g) >= 2]
    if len(pairs) < min_pairs:
        print(f"\nTest-retest ({label}): {len(pairs)} patients with two QC-passing scans; "
              f"at least {min_pairs} needed, skipped")
        return None
    first = pd.DataFrame([a for a, _ in pairs])
    second = pd.DataFrame([b for _, b in pairs])
    rs = []
    for c in names:
        x = pd.to_numeric(first[c], errors="coerce").to_numpy(float)
        y = pd.to_numeric(second[c], errors="coerce").to_numpy(float)
        m = np.isfinite(x) & np.isfinite(y)
        if m.sum() >= min_pairs and np.std(x[m]) > 0 and np.std(y[m]) > 0:
            rs.append(float(np.corrcoef(x[m], y[m])[0, 1]))
    if not rs:
        return None
    print(f"\nTest-retest ({label}, {len(pairs)} patients scanned twice, {len(rs)} ROIs): "
          f"mean r={np.mean(rs):.3f}, median r={np.median(rs):.3f} "
          "(scan intervals are not known to this script)")
    return rs


def main():
    ap = argparse.ArgumentParser(description="Native individualized cortical QSM extraction")
    ap.add_argument("--native-dir", type=Path, default=DEFAULT_NATIVE_DIR)
    ap.add_argument("--qsm-dir", type=Path, default=DEFAULT_NATIVE_DIR,
                    help="Retained option; native QSM paths are read from the worklist")
    ap.add_argument("--out-dir", type=Path, default=HERE / "output")
    ap.add_argument("--worklist", type=Path, default=HERE / "worklist.csv")
    ap.add_argument("--labels", type=Path, default=DK_LABELS)
    ap.add_argument("--freesurfer-home", type=Path, default=DEFAULT_FREESURFER_HOME,
                    help="FreeSurfer used by P8 (FreeSurferColorLUT.txt and the SynthSeg --parc label list)")
    ap.add_argument("--allow-label-name-mismatch", action="store_true",
                    help="Continue when CSV region names differ from the FreeSurfer regions assigned by ID "
                         "(only after checking the printed table)")
    ap.add_argument("--failed-list", type=Path, default=HERE / "failed_ID.txt",
                    help="P1 batch failure list; each entry is classified by --p1-qc-dir evidence")
    ap.add_argument("--p1-qc-dir", type=Path, default=DEFAULT_P1_QC_DIR,
                    help="P1 output directory holding <IID>_qc.txt; a failure confined to T1->MNI "
                         "registration does not exclude a subject from the native analysis")
    ap.add_argument("--include-failed", action="store_true",
                    help="Also extract failure-list subjects whose failure involves the native branch "
                         "or has no stage evidence (diagnostics only; they never enter the analysis set)")
    ap.add_argument("--ids-file", type=Path, default=None,
                    help="Process only subjects in this list (e.g., pilot_native.csv)")
    ap.add_argument("--main-model-only", action="store_true",
                    help="Process only worklist rows with include_main_model=yes (by default every "
                         "runnable row is processed and include_main_model is carried into the outputs)")
    ap.add_argument("--sample-root", type=Path, default=None,
                    help="Sample root; each subject directory must contain one QSM and one T1 image")
    ap.add_argument("--support-dir", type=Path, default=None,
                    help="Directory of valid QSM support masks named <IID>_support.nii[.gz]")
    ap.add_argument("--allow-nonzero-proxy", action="store_true",
                    help="Without a support mask, use nonzero QSM as a provisional coverage proxy "
                         "(diagnostics/sensitivity analysis only; every output file name carries _proxy)")
    ap.add_argument("--which", default="",
                    help="With --ids-file, filter further by 'column=value', e.g., pilot_group=primary_pair")
    ap.add_argument("--stat", choices=("both", "median", "mean"), default="both")
    ap.add_argument("--min-voxels", type=int, default=20,
                    help="Minimum valid voxels per ROI")
    ap.add_argument("--min-mm3", type=float, default=20.0,
                    help="Minimum valid volume per ROI in mm3 (voxel volume differs between QSM protocols)")
    ap.add_argument("--min-coverage", type=float, default=0.5,
                    help="Minimum fraction of the ROI's T1-space volume with valid QSM; 0 disables "
                         "(uses the P8 t1seg output)")
    ap.add_argument("--min-reportable-rois", type=int, default=1,
                    help="Minimum ROIs with a reported value for a subject to enter the analysis set "
                         "(68 keeps complete cases only, for models that need every ROI)")
    ap.add_argument("--csf-ref", action="store_true",
                    help="Also save cortical values referenced to each subject's lateral-ventricle CSF median")
    ap.add_argument("--min-csf-voxels", type=int, default=30)
    ap.add_argument("--suffix", default="")
    ap.add_argument("--jobs", type=int, default=16)
    args = ap.parse_args()
    if args.support_dir is None and not args.allow_nonzero_proxy:
        ap.error("Formal extraction requires --support-dir; use --allow-nonzero-proxy explicitly for diagnostics only")
    if args.which and not args.ids_file:
        ap.error("--which filters the --ids-file list; give --ids-file as well")
    if not 0.0 <= args.min_coverage <= 1.0:
        ap.error("--min-coverage must be between 0 and 1")
    if not 1 <= args.min_reportable_rois <= len(FS_L) + len(FS_R):
        ap.error(f"--min-reportable-rois must be between 1 and {len(FS_L) + len(FS_R)}")
    coverage_mode = "support_mask" if args.support_dir else "nonzero_proxy"
    if args.support_dir and args.allow_nonzero_proxy:
        print("Note: --support-dir given; --allow-nonzero-proxy is ignored")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    label_dir = args.native_dir / "label"
    t1seg_dir = args.native_dir / "t1seg"
    qc_dir = args.native_dir / "qc"
    # Proxy-mode (diagnostic/sensitivity) outputs never share file names with formal ones.
    sfx = ("_proxy" if coverage_mode == "nonzero_proxy" else "") + args.suffix
    if args.sample_root and not args.suffix:
        sfx += "_sample_debug"
    print(f"Coverage mode: {coverage_mode}"
          + ("" if args.support_dir else " (diagnostic/sensitivity only; output names carry _proxy)"))

    def out_path(tail):
        return args.out_dir / f"HUASHAN_NATIVE_DK{sfx}{tail}"

    csv_labels = read_labels(args.labels)
    fs_to_name, mapping_rows = build_label_map(csv_labels, args.labels, args.freesurfer_home,
                                               args.allow_label_name_mismatch)
    keys = sorted(fs_to_name)
    names = [fs_to_name[k] for k in keys]
    mapping_path = args.out_dir / f"native_label_mapping{sfx}.csv"
    write_csv(pd.DataFrame(mapping_rows), mapping_path)
    print(f"Cortical label mapping OK: {len(keys)} regions (34 per hemisphere, including insula; "
          f"{INSULA_L}/{INSULA_R} in {args.labels.name} has no corresponding CSV entry and is added in this branch); "
          f"names checked against the FreeSurfer LUT ({mapping_path.name})")

    wl = pd.read_csv(args.worklist, dtype=str).fillna("")
    # Worklist input headers must be English; no legacy-language aliases are used.
    if "image_dir_id" not in wl.columns and "IID" in wl.columns:
        wl["image_dir_id"] = wl["IID"]
    # t1_file is needed to verify that P8 segmented the T1 the worklist names.
    required = {"image_dir_id", "qsm_file", "t1_file"}
    absent = sorted(required - set(wl.columns))
    if absent:
        raise SystemExit(f"Missing worklist columns: {absent}. Rename input headers to English.")
    wl, conflicts = resolve_duplicate_ids(wl)
    excluded = list(conflicts)
    if conflicts:
        print(f"WARNING: {len(conflicts)} image_dir_id values have conflicting worklist rows and are "
              f"excluded until the worklist is corrected: {[e['IID'] for e in conflicts]}")
    if "runnable" in wl.columns:
        runnable = wl["runnable"].str.strip().str.lower() == "yes"
        excluded += [{"IID": oid, "reason": f"worklist runnable={value!r}"}
                     for oid, value in zip(wl.loc[~runnable, "image_dir_id"], wl.loc[~runnable, "runnable"])]
        wl = wl[runnable]
    if args.main_model_only:
        if "include_main_model" not in wl.columns:
            raise SystemExit("--main-model-only needs an include_main_model worklist column")
        main_model = wl["include_main_model"].str.strip().str.lower() == "yes"
        excluded += [{"IID": oid, "reason": f"worklist include_main_model={value!r} (--main-model-only)"}
                     for oid, value in zip(wl.loc[~main_model, "image_dir_id"],
                                           wl.loc[~main_model, "include_main_model"])]
        wl = wl[main_model]
        print(f"Main-model rows only (--main-model-only): {len(wl)} subjects")

    # The P1 failure list mixes stages. A failure in a stage the native branch shares (QSM
    # input, QSM->T1 rigid registration) excludes the subject, as does an entry without stage
    # evidence; a failure confined to T1->MNI registration (FNIRT) or the MNI output does not.
    failed = read_failed_ids(args.failed_list)
    p1_stage = {oid: classify_p1_failure(args.p1_qc_dir / f"{oid}_qc.txt")
                for oid in wl["image_dir_id"] if oid in failed}
    blocking = {oid for oid, (stage, _) in p1_stage.items() if stage != "mni_only"}
    if p1_stage:
        stages = Counter(stage for stage, _ in p1_stage.values())
        print(f"P1 failure-list entries in the worklist: {len(p1_stage)} "
              f"({', '.join(f'{k} {v}' for k, v in sorted(stages.items()))}); "
              f"mni_only entries stay in the native analysis")
    if args.include_failed:
        print(f"Included failure-list IDs involving the native branch or without stage evidence "
              f"(diagnostic mode; never in the analysis set): {len(blocking)} subjects")
    else:
        excluded += [{"IID": oid, "reason": f"P1 failure (stage {p1_stage[oid][0]}): {p1_stage[oid][1]}"}
                     for oid in wl["image_dir_id"] if oid in blocking]
        wl = wl[~wl["image_dir_id"].isin(blocking)]

    selection = None                    # IDs chosen by --sample-root / --ids-file, if any
    if args.sample_root:
        if not args.sample_root.is_dir():
            raise SystemExit(f"Sample directory does not exist: {args.sample_root}")
        sample_ids = sorted(p.name for p in args.sample_root.iterdir() if p.is_dir())
        if not sample_ids:
            raise SystemExit(f"No subject subdirectories in sample root: {args.sample_root}")
        wl = wl[wl["image_dir_id"].isin(sample_ids)]
        selection = set(sample_ids)
        print(f"Filtered by sample root: {len(wl)} subjects")

    if args.ids_file:
        p = pd.read_csv(args.ids_file, dtype=str).fillna("")
        idcol = next((c for c in ("IID", "image_dir_id") if c in p.columns), None)
        if idcol is None:
            raise SystemExit(f"{args.ids_file} has neither an IID nor an image_dir_id column")
        if args.which:
            k, sep, v = args.which.partition("=")
            if not sep or k not in p.columns:
                raise SystemExit(f"--which must be column=value with a column of {args.ids_file}; "
                                 f"got {args.which!r}")
            p = p[p[k] == v]
        wl = wl[wl["image_dir_id"].isin(p[idcol])]
        selection = set(p[idcol]) if selection is None else selection & set(p[idcol])
        print(f"Filtered by ID list: {len(wl)} subjects")

    qsm_of = dict(zip(wl["image_dir_id"], wl["qsm_file"]))
    t1_of = dict(zip(wl["image_dir_id"], wl["t1_file"]))
    meta = wl.set_index("image_dir_id")

    ids, tasks = [], []
    for oid in wl["image_dir_id"]:
        lab = label_dir / f"{oid}_aparc_aseg_QSMnative.nii.gz"
        qsm = Path(qsm_of.get(oid, ""))
        t1 = Path(t1_of.get(oid, ""))
        if args.sample_root:
            sample_dir = args.sample_root / oid
            qsm_candidates = [p for p in sample_dir.glob("*.nii*")
                              if "qsm" in p.name.lower() and "swi" not in p.name.lower()]
            t1_candidates = [p for p in sample_dir.glob("*.nii*")
                             if "t1" in p.name.lower()]
            if len(qsm_candidates) != 1 or len(t1_candidates) != 1:
                excluded.append({"IID": oid, "reason": f"sample directory has QSM={len(qsm_candidates)}, "
                                                       f"T1={len(t1_candidates)} candidates"})
                continue
            qsm, t1 = qsm_candidates[0], t1_candidates[0]
        support = None
        if args.support_dir:
            support_candidates = [args.support_dir / f"{oid}_support.nii.gz",
                                  args.support_dir / f"{oid}_support.nii"]
            support = next((p for p in support_candidates if p.is_file()), None)
        t1seg = t1seg_dir / f"{oid}_aparc+aseg_T1.nii.gz"
        # The first missing input in pipeline order is reported, so the reason names the cause.
        if not (lab.is_file() and lab.stat().st_size > 0):
            excluded.append({"IID": oid, "reason": "P8 native label missing"})
        elif not qsm.is_file():
            excluded.append({"IID": oid, "reason": f"QSM file missing: {qsm}"})
        elif not t1.is_file():
            excluded.append({"IID": oid, "reason": f"T1 file missing (needed to verify P8's T1): {t1}"})
        elif args.support_dir and support is None:
            excluded.append({"IID": oid, "reason": "support mask missing"})
        elif args.min_coverage > 0 and not t1seg.is_file():
            excluded.append({"IID": oid, "reason": "P8 T1 segmentation missing (needed for --min-coverage)"})
        else:
            ids.append(oid)
            tasks.append({"label": lab, "qsm": qsm, "t1": t1, "support": support,
                          "t1seg": t1seg if t1seg.is_file() else None})
    if selection is not None:
        # Worklist-level exclusions were recorded before the selection; keep the selected subjects'.
        excluded = [e for e in excluded if e["IID"] in selection]
    excluded_path = out_path("_excluded.csv")
    write_csv(pd.DataFrame(excluded, columns=["IID", "reason"]), excluded_path)
    print(f"\nSubjects with native labels: {len(ids)}; excluded or incomplete: {len(excluded)}"
          + (f" (first 10: {[e['IID'] for e in excluded[:10]]}; all in {excluded_path.name})"
             if excluded else ""))

    def write_run_info(counts):
        """Record what produced these outputs, so analysis_pass can be reproduced and audited."""
        def described(path):
            path = Path(path) if path else None
            return ({"path": str(path), "sha1": file_digest(path)} if path is not None and path.is_file()
                    else {"path": str(path) if path else None, "sha1": None})
        info = {
            "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
            "command": sys.argv,
            "coverage_mode": coverage_mode,
            "options": {k: (str(v) if isinstance(v, Path) else v) for k, v in sorted(vars(args).items())},
            "scripts": {"P9_extract_native_qsm_ALL.py": file_digest(Path(__file__))},
            "inputs": {"worklist": described(args.worklist), "labels": described(args.labels),
                       "freesurfer_lut": described(args.freesurfer_home / "FreeSurferColorLUT.txt"),
                       "failed_list": described(args.failed_list)},
            "grid_conversion_tag": GRID_CONVERSION_TAG,
            "versions": {"python": platform.python_version(), "numpy": np.__version__,
                         "nibabel": nib.__version__, "pandas": pd.__version__},
            "counts": counts,
        }
        path = out_path("_run_info.json")
        tmp = path.with_name(f".{path.name}.tmp")
        tmp.write_text(json.dumps(info, indent=1) + "\n", encoding="utf-8")
        os.replace(tmp, path)

    if not ids:
        write_run_info({"extracted": 0, "excluded": len(excluded), "analysis_pass": 0})
        raise SystemExit("No eligible subjects. Complete P8 first.")

    opts = {"min_voxels": args.min_voxels, "min_mm3": args.min_mm3,
            "min_coverage": args.min_coverage, "csf_ref": args.csf_ref,
            "min_csf_voxels": args.min_csf_voxels}
    print(f"Extracting with {args.jobs} threads...")
    with ThreadPoolExecutor(args.jobs) as ex:
        results = list(ex.map(lambda t: process_subject(t, keys, opts), tasks))

    n_sub, n_keys = len(ids), len(keys)
    counts = np.zeros((n_sub, n_keys), dtype=np.int64)
    label_counts = np.zeros((n_sub, n_keys), dtype=np.int64)
    zero_counts = np.zeros((n_sub, n_keys), dtype=np.int64)
    unsupported_zero = np.zeros((n_sub, n_keys), dtype=np.int64)
    means = np.full((n_sub, n_keys), np.nan)
    medians = np.full((n_sub, n_keys), np.nan)
    coverage = np.full((n_sub, n_keys), np.nan)
    t1_mm3 = np.full((n_sub, n_keys), np.nan)
    vox_mm3 = np.full(n_sub, np.nan)
    status = [r["status"] for r in results]
    for i, r in enumerate(results):
        if "counts" not in r:
            continue
        counts[i], label_counts[i] = r["counts"], r["label_counts"]
        zero_counts[i], unsupported_zero[i] = r["zero_counts"], r["unsupported_zero_counts"]
        means[i], medians[i] = r["means"], r["medians"]
        coverage[i], t1_mm3[i], vox_mm3[i] = r["coverage"], r["t1_mm3"], r["vox_mm3"]

    # Median/mean columns stay adjacent for each ROI.
    stats = {"both": ("med", "mean"), "median": ("med",), "mean": ("mean",)}[args.stat]
    primary = stats[0]
    value_cols = {}
    for j, nm in enumerate(names):
        if "med" in stats:
            value_cols[f"{nm}_med"] = medians[:, j]
        if "mean" in stats:
            value_cols[f"{nm}_mean"] = means[:, j]
    roi_out = pd.DataFrame(value_cols)

    # Metadata first; build in one operation to avoid a fragmented DataFrame.
    sub = meta.reindex(ids)
    meta_cols = [c for c in META_COLS if c in sub.columns]
    meta_out = sub[meta_cols].reset_index(drop=True)
    meta_out.insert(0, "IID", ids)
    out = pd.concat([meta_out, roi_out], axis=1)

    csv_path = out_path(".csv")
    write_csv(out, csv_path)

    csf_values = np.array([r.get("csf_ref", np.nan) for r in results], dtype=float)
    csf_counts = np.array([r.get("csf_voxels", 0) for r in results], dtype=np.int64)

    # ---------- Analysis set ----------
    # Each file read in this run is compared by content (SHA-1) with the provenance P8 and
    # make_support_masks.py recorded; an unchanged path alone is never accepted.
    p8_qc = [read_native_qc(qc_dir / f"{oid}_native_qc.txt") for oid in ids]
    p1_failed = np.asarray([oid in failed for oid in ids], dtype=bool)
    p8_dice = np.array([to_float(f.get("reg_mask_dice")) for _, _, f in p8_qc])
    grid_ok = np.array([f.get("grid_conversion") == GRID_CONVERSION_TAG for _, _, f in p8_qc])

    def matches_p8(key, field):
        """True/False when both the file hash and P8's record exist, else None."""
        return [None if not (r.get(key) and f.get(field)) else r[key] == f[field]
                for r, (_, _, f) in zip(results, p8_qc)]

    qsm_match = matches_p8("qsm_sha1", "qsm_sha1")
    t1_match = matches_p8("t1_sha1", "t1_sha1")
    label_match = matches_p8("label_sha1", "label_sha1")
    t1seg_match = matches_p8("t1seg_sha1", "t1seg_sha1")
    qsm_path_match = [same_path(t["qsm"], f.get("qsm", "")) for t, (_, _, f) in zip(tasks, p8_qc)]
    t1_path_match = [same_path(t["t1"], f.get("t1", "")) for t, (_, _, f) in zip(tasks, p8_qc)]
    n_reportable = np.array([r.get("n_reportable", 0) for r in results], dtype=int)

    def support_problem(res):
        """Why the support mask cannot be used with this QSM, or '' when its record verifies."""
        record = res.get("support_record")
        if record is None:
            return "support mask provenance record unreadable"
        if not record:
            return "support mask has no provenance record (rebuild it with make_support_masks.py)"
        if record.get("support_type") not in SUPPORT_TYPES:
            return f"support mask has unknown support_type {record.get('support_type')!r}"
        if record.get("qsm_sha1") != res.get("qsm_sha1"):
            return "support mask was built from a different QSM"
        if record.get("mask_sha1") != res.get("support_sha1"):
            return "support mask file differs from its provenance record"
        return ""

    fail_reasons, support_ok = [], []
    for i in range(n_sub):
        why = []
        res, (qc_status, _, fields) = results[i], p8_qc[i]
        if status[i] != "OK":
            why.append(f"extraction: {status[i]}")
        if ids[i] in blocking:          # present only with --include-failed
            stage, evidence = p1_stage[ids[i]]
            why.append(f"P1 failure (stage {stage}): {evidence}")
        if qc_status == "missing":
            why.append("P8 QC missing")
        else:
            if fields.get("qc_complete") != "1":
                why.append("P8 QC incomplete")
            if qc_status != "pass":
                why.append("P8 QC failed")
            if not np.isfinite(p8_dice[i]):
                why.append("P8 reg_mask_dice missing")
            if not grid_ok[i]:
                why.append("P8 label predates grid-conversion fix (rerun P8)")
            checks = [("qsm_sha1", qsm_match[i], "QSM content differs from the file P8 used", "QSM"),
                      ("t1_sha1", t1_match[i], "T1 content differs from the file P8 segmented", "T1"),
                      ("label_sha1", label_match[i], "label differs from the one P8 QC checked", "label")]
            if args.min_coverage > 0:
                checks.append(("t1seg_sha1", t1seg_match[i], "T1 segmentation differs from the one P8 used",
                               "segmentation"))
            for field, match, mismatch, what in checks:
                if match is False:
                    why.append(mismatch)
                elif field not in fields:
                    why.append(f"P8 QC has no {what} provenance")
        problem = ""
        if tasks[i]["support"] is not None and "qsm_sha1" in res:
            problem = support_problem(res)
            if problem:
                why.append(problem)
        support_ok.append(None if tasks[i]["support"] is None or "qsm_sha1" not in res else not problem)
        if status[i] == "OK" and n_reportable[i] < args.min_reportable_rois:
            why.append(f"{n_reportable[i]} reportable ROIs < --min-reportable-rois {args.min_reportable_rois}")
        fail_reasons.append("; ".join(why))
    analysis_pass = np.array([not w for w in fail_reasons], dtype=bool)
    # CSF referencing has its own eligibility flag; it never changes the primary analysis set.
    csfref_pass = analysis_pass & np.isfinite(csf_values) if args.csf_ref else np.zeros(n_sub, dtype=bool)

    analysis_path = out_path("_analysis.csv")
    write_csv(out.loc[analysis_pass], analysis_path)

    if args.csf_ref:
        ref_roi = roi_out.subtract(csf_values, axis=0)
        ref_out = pd.concat([meta_out, ref_roi], axis=1)
        ref_path = out_path("_csfref.csv")
        write_csv(ref_out, ref_path)
        write_csv(ref_out.loc[csfref_pass], out_path("_csfref_analysis.csv"))
        print(f"CSF-referenced results: {ref_path} ({int(csfref_pass.sum())} of {int(analysis_pass.sum())} "
              "analysis-set subjects have a CSF reference)")

    qsm_label_mm3 = label_counts * vox_mm3[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        fov_fraction = qsm_label_mm3 / t1_mm3
    finite_fov = np.where(np.isfinite(fov_fraction), fov_fraction, -np.inf).max(axis=1)
    hosp = (meta_out["hospital_id"] if "hospital_id" in meta_out.columns
            else pd.Series([""] * n_sub))
    repeat_patient = (hosp.ne("") & hosp.duplicated(keep=False)).to_numpy()
    cortex_iqr = np.array([r.get("cortex_iqr", np.nan) for r in results], dtype=float)
    records = [r.get("support_record") or {} for r in results]
    # Protocol stratum: QSM voxel size and axial coverage (worklist qsm_coverage_mm, else the
    # superior-inferior extent of the QSM field of view).
    fov_si_mm = np.array([r.get("si_extent_mm", np.nan) for r in results], dtype=float)
    wl_coverage = (pd.to_numeric(sub["qsm_coverage_mm"], errors="coerce").to_numpy(float)
                   if "qsm_coverage_mm" in sub.columns else np.full(n_sub, np.nan))
    coverage_mm = np.where(np.isfinite(wl_coverage), wl_coverage, fov_si_mm)
    protocol = [voxel_size_label(r["pixdim"]) if r.get("pixdim") else "unknown" for r in results]
    coverage_grp = [coverage_group(c) for c in coverage_mm]
    stratum = np.array([f"{p} | {c}" for p, c in zip(protocol, coverage_grp)], dtype=object)
    with np.errstate(divide="ignore", invalid="ignore"):
        cortex_zero_frac = zero_counts.sum(axis=1) / counts.sum(axis=1)
    qc_out = pd.DataFrame({
        "IID": ids,
        "analysis_pass": analysis_pass,
        "fail_reasons": fail_reasons,
        "csfref_analysis_pass": csfref_pass if args.csf_ref else [""] * n_sub,
        "coverage_mode": coverage_mode,
        "support_mask": [str(t["support"] or "") for t in tasks],
        "support_type": [rec.get("support_type", "") for rec in records],
        "support_method": [rec.get("method", "") for rec in records],
        "support_max_hole_mm3": [rec.get("max_hole_mm3", "") for rec in records],
        "support_provenance_ok": support_ok,
        "qsm_nonzero_outside_support": [r.get("qsm_nonzero_outside_support", "") for r in results],
        "extract_status": status,
        "extract_detail": ["" if s == "OK" else s for s in status],
        "P1_failed_list": p1_failed,
        "P1_failure_stage": [p1_stage.get(oid, ("", ""))[0] for oid in ids],
        "P1_failure_evidence": [p1_stage.get(oid, ("", ""))[1] for oid in ids],
        "P8_native_qc_status": [x[0] for x in p8_qc],
        "P8_native_qc_errors": [x[1] for x in p8_qc],
        "P8_qc_complete": [x[2].get("qc_complete", "") == "1" for x in p8_qc],
        "P8_reg_mask_dice": p8_dice,
        "P8_seg_qsm_dice": [x[2].get("seg_qsm_dice", "") for x in p8_qc],
        "P8_gm_wm_contrast_ppb": [x[2].get("gm_wm_contrast_ppb", "") for x in p8_qc],
        "P8_ribbon_in_qsm_frac": [x[2].get("ribbon_in_qsm_frac", "") for x in p8_qc],
        "P8_cortex_vol_ml": [x[2].get("cortex_vol_ml", "") for x in p8_qc],
        "P8_synthseg_icv_ml": [x[2].get("synthseg_icv_ml", "") for x in p8_qc],
        "P8_synthseg_cortex_ml": [x[2].get("synthseg_cortex_ml", "") for x in p8_qc],
        "P8_qc_failed_rules": [x[2].get("qc_failed_rules", "") for x in p8_qc],
        "P8_qc_report_only_failed_rules": [x[2].get("qc_report_only_failed_rules", "") for x in p8_qc],
        "P8_p1_t1_resolution": [x[2].get("p1_t1_resolution", "") for x in p8_qc],
        "P8_p1_search_mode": [x[2].get("p1_search_mode", "") for x in p8_qc],
        "P8_grid_conversion": [x[2].get("grid_conversion", "") for x in p8_qc],
        "P8_t1_neurological": [x[2].get("t1_neurological", "") for x in p8_qc],
        "P8_p1_matrix_provenance": [x[2].get("p1_matrix_provenance", "") for x in p8_qc],
        "P8_q2t_scale_dev": [x[2].get("q2t_scale_dev", "") for x in p8_qc],
        "P8_qsm_pixdim_sform_rel_diff": [x[2].get("qsm_pixdim_sform_rel_diff", "") for x in p8_qc],
        "P8_t1_pixdim_sform_rel_diff": [x[2].get("t1_pixdim_sform_rel_diff", "") for x in p8_qc],
        "qsm_matches_P8": qsm_match,
        "qsm_path_matches_P8": qsm_path_match,
        "t1_matches_P8": t1_match,
        "t1_path_matches_P8": t1_path_match,
        "label_matches_P8": label_match,
        "t1seg_matches_P8": t1seg_match,
        "qsm_pixdim": [r.get("pixdim", "") for r in results],
        "qsm_voxel_mm3": vox_mm3,
        "qsm_dtype": [r.get("qsm_dtype", "") for r in results],
        "qsm_scl_slope": [r.get("qsm_scl_slope", "") for r in results],
        "qsm_scl_inter": [r.get("qsm_scl_inter", "") for r in results],
        "qsm_integer_valued": [r.get("qsm_integer_valued", "") for r in results],
        "qsm_quantized": [r.get("qsm_quantized", "") for r in results],
        "qsm_protocol": protocol,
        "qsm_fov_si_mm": fov_si_mm.round(1),
        "qsm_coverage_group": coverage_grp,
        "qsm_coverage_mm": (sub["qsm_coverage_mm"].to_numpy() if "qsm_coverage_mm" in sub.columns
                            else [""] * n_sub),
        "t1_cortex_ml": np.where(np.isfinite(t1_mm3).any(axis=1), np.nansum(t1_mm3, axis=1) / 1000.0, np.nan),
        "max_fov_fraction": np.where(np.isfinite(finite_fov), finite_fov, np.nan),
        "repeat_patient": repeat_patient,
        "csf_ref_ppb": csf_values,
        "csf_voxels": csf_counts,
        "label_voxels_without_qsm_support": (label_counts - counts).sum(axis=1),
        "cortex_valid_zero_fraction": np.round(cortex_zero_frac, 5),
        "cortex_unsupported_zero_voxels": unsupported_zero.sum(axis=1),
        "n_reportable_rois": n_reportable,
        "n_nan_DK": out[[f"{n}_{primary}" for n in names]].isna().sum(axis=1).to_numpy(),
        "n_low_voxels": [r.get("n_low_voxels", n_keys) for r in results],
        "n_low_coverage": [r.get("n_low_coverage", 0) for r in results],
        "cortex_valid_iqr": cortex_iqr,
        "possible_ppm_units": cortex_iqr < PPM_IQR_SUSPECT,
    })
    write_csv(qc_out, out_path("_qc.csv"))

    used = counts.sum(axis=0)
    atlas_voxels = label_counts.sum(axis=0)
    write_csv(pd.DataFrame({
        "ROI": names,
        "atlas_label_voxels": atlas_voxels,
        "used_qsm_voxels": used,
        "excluded_cerebellar_voxel_fraction": 0.0,
    }), args.out_dir / f"native_mask_voxel_counts{sfx}.csv")

    sizes = label_counts.astype(float)
    cov = np.divide(counts, sizes, out=np.zeros(counts.shape), where=sizes > 0)
    reported = medians if primary == "med" else means
    write_csv(pd.DataFrame({
        "IID": np.repeat(ids, n_keys),
        "ROI": np.tile(names, n_sub),
        "atlas_label_voxels": label_counts.ravel(),
        "valid_qsm_voxels": counts.ravel(),
        "valid_zero_voxels": zero_counts.ravel(),
        "unsupported_zero_voxels": unsupported_zero.ravel(),
        "coverage_fraction": cov.ravel().round(4),
        "coverage_definition": ("explicit_support_mask" if args.support_dir
                                else "nonzero_QSM_proxy; valid zeros may be excluded"),
        "t1_label_mm3": t1_mm3.ravel().round(1),
        "qsm_grid_label_mm3": qsm_label_mm3.ravel().round(1),
        "valid_qsm_mm3": (counts * vox_mm3[:, None]).ravel().round(1),
        "fov_fraction": fov_fraction.ravel().round(4),
        "coverage_vs_t1": coverage.ravel().round(4),
        "roi_value_reported": np.isfinite(reported).ravel(),
    }), args.out_dir / f"native_roi_coverage{sfx}.csv")

    # ROI missingness in the analysis set, overall and by protocol stratum. Missing values
    # depend on slab coverage, so they are not missing at random across protocols (D2).
    present = np.isfinite(reported)
    miss_rows = []
    for level in ["all"] + sorted(set(stratum[analysis_pass])):
        rows_in = analysis_pass & (np.ones(n_sub, dtype=bool) if level == "all" else stratum == level)
        n = int(rows_in.sum())
        for j, nm in enumerate(names):
            k = int(present[rows_in, j].sum())
            cov_j = coverage[rows_in, j]
            miss_rows.append({"stratum": level, "ROI": nm, "n_subjects": n, "n_reported": k,
                              "missing_fraction": round(1.0 - k / n, 4) if n else np.nan,
                              "median_coverage_vs_t1": (round(float(np.median(cov_j[np.isfinite(cov_j)])), 4)
                                                        if np.isfinite(cov_j).any() else np.nan)})
    write_csv(pd.DataFrame(miss_rows, columns=["stratum", "ROI", "n_subjects", "n_reported",
                                               "missing_fraction", "median_coverage_vs_t1"]),
              args.out_dir / f"native_roi_missingness{sfx}.csv")
    write_run_info({"extracted": n_sub, "excluded": len(excluded),
                    "extraction_ok": sum(1 for s in status if s == "OK"),
                    "analysis_pass": int(analysis_pass.sum()),
                    "csfref_analysis_pass": int(csfref_pass.sum()) if args.csf_ref else None})

    n_ok = sum(1 for s in status if s == "OK")
    print(f"QC-passing subjects for downstream analysis: {int(analysis_pass.sum())}/{n_sub} "
          f"(rows of {analysis_path.name}); other subjects remain in {csv_path.name} but are not QC-passing")
    tally = Counter("extraction failed" if w.startswith("extraction:") else w
                    for reasons in fail_reasons for w in reasons.split("; ") if w)
    for reason, k in tally.most_common():
        print(f"  {k:5d}  {reason}")
    print(f"\nCompleted: {csv_path}  ({out.shape[0]} subjects x {len(names)} regions, "
          f"successful extractions: {n_ok}/{n_sub})")
    if n_sub - n_ok:
        bad = [(ids[i], status[i]) for i in range(n_sub) if status[i] != "OK"]
        print(f"  Failed: {len(bad)} subjects (first 5): {bad[:5]}")
    ok_rows = [i for i in range(n_sub) if status[i] == "OK"]
    print(f"ROI values withheld in successful extractions: "
          f"{sum(results[i]['n_low_voxels'] for i in ok_rows)} below --min-voxels {args.min_voxels}/"
          f"--min-mm3 {args.min_mm3:g}, {sum(results[i]['n_low_coverage'] for i in ok_rows)} below "
          f"--min-coverage {args.min_coverage:g}")
    if ok_rows:
        print(f"Reportable ROIs per successful extraction: median {np.median(n_reportable[ok_rows]):.0f}, "
              f"minimum {n_reportable[ok_rows].min()} (--min-reportable-rois {args.min_reportable_rois})")
    if (cortex_iqr < PPM_IQR_SUSPECT).any():
        print(f"  WARNING: {int((cortex_iqr < PPM_IQR_SUSPECT).sum())} subjects have cortical QSM IQR < "
              f"{PPM_IQR_SUSPECT} ppb (possible ppm units; see possible_ppm_units)")
    if repeat_patient.any():
        print(f"  Note: {int(repeat_patient.sum())} rows belong to patients with more than one scan "
              "(repeat_patient); they are not independent observations")
    if analysis_pass.any():
        complete = analysis_pass & (n_reportable == n_keys)
        by_stratum = Counter(stratum[analysis_pass])
        full = Counter(stratum[complete])
        print(f"Analysis-set subjects with all {n_keys} ROI values: {int(complete.sum())}/{int(analysis_pass.sum())} ("
              + "; ".join(f"{k}: {full.get(k, 0)}/{v}" for k, v in sorted(by_stratum.items()))
              + f"); missingness by ROI and stratum in native_roi_missingness{sfx}.csv")
    if not args.main_model_only and "include_main_model" in meta_out.columns:
        outside = analysis_pass & meta_out["include_main_model"].str.strip().str.lower().ne("yes").to_numpy()
        if outside.any():
            print(f"  Note: {int(outside.sum())} analysis-set subjects have include_main_model other than "
                  "'yes'; use --main-model-only, or filter on include_main_model, for the main model")

    if not analysis_pass.any():
        print("\nNo QC-passing subjects; per-ROI summaries skipped.")
    else:
        per_roi = pd.DataFrame(counts[analysis_pass], columns=names).mean(axis=0)
        print(f"\nMean voxel count by ROI: median {per_roi.median():.0f},"
              f"minimum {per_roi.min():.0f}({per_roi.idxmin()})")
        low = [n for n in names if per_roi[n] < args.min_voxels * 2]
        if low:
            print(f"  WARNING Regions with low voxel counts: {len(low)} regions: {low[:8]}")

        print("\nCross-subject median QSM by ROI (ppb; P1/P8 QC-passing subjects only):")
        for c in names:
            v = pd.to_numeric(out.loc[analysis_pass, f"{c}_{primary}"], errors="coerce").dropna()
            if v.size:
                print(f"  {c:<34} {v.median():>7.2f}  SD={v.std():>6.2f}  n={v.size}")

        analysis_out = out.loc[analysis_pass].copy()
        for st, label in (("med", "native"), ("mean", "native (mean statistic)")):
            if st in stats:
                renamed = analysis_out.rename(columns={f"{n}_{st}": n for n in names})
                lr_concordance(renamed, names, label)
                test_retest(renamed, names, label)
        if args.csf_ref:
            lr_concordance(ref_out.loc[csfref_pass].rename(columns={f"{n}_{primary}": n for n in names}),
                           names, "native, CSF-referenced")

    if excluded:
        print(f"\nExcluded or incomplete: {len(excluded)} subjects (first 10): "
              f"{[e['IID'] for e in excluded[:10]]}; reasons in {excluded_path.name}")


if __name__ == "__main__":
    main()
