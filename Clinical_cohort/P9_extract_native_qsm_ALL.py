#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Extract individualized cortical QSM statistics on the native QSM grid.

Required English worklist headers: image_dir_id, qsm_file.
Optional English worklist headers: runnable (yes/no), hospital_id, model_label,
include_main_model, category.
Optional ID-list headers: IID or image_dir_id.
The --which filter uses the exact English column names in that ID list.

Input files must be migrated to English headers before running this script.
This script does not translate existing input files or legacy downstream schemas.
"""

import argparse
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from scipy.ndimage import binary_erosion

HERE = Path(__file__).resolve().parent
ATLAS_DIR = Path("/cwStorage/home/fuyan/FY_data/QSM_STAAR_251229/Session1_Pheno_extract/"
                 "Part1_Extract_QSM")
DK_LABELS = ATLAS_DIR / "Desikan_labels.csv"
FS_LUT = Path("/public/software/apps/Freesurfer/8.2.0-1/FreeSurferColorLUT.txt")
FS_PARC_LABELS = Path("/public/software/apps/Freesurfer/8.2.0-1/models/"
                      "synthseg_parcellation_labels.npy")

DEFAULT_NATIVE_DIR = Path("/cwStorage/nodecw_group/FY_data/QSM_HUASHAN/native_data")
NATIVE_ROOT = DEFAULT_NATIVE_DIR

NON_CORTICAL = ("white_matter", "corpus")

META_COLS = ("hospital_id", "model_label", "include_main_model", "category")

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
    import csv
    out = {}
    with open(path, newline="", encoding="utf-8-sig") as handle:
        for row in csv.reader(handle):
            if len(row) >= 2 and row[0].strip().isdigit():
                out[int(row[0])] = row[1].strip()
    if not out:
        raise SystemExit(f"{path}  contains no parsable labels.")
    return out


def read_failed_ids(path):
    """Read failed subject IDs from the P1 batch list (supports literal backslash-t)."""
    if path is None or not Path(path).exists():
        return set()
    ids = set()
    for line in open(path, encoding="utf-8", errors="replace"):
        line = line.strip().replace("\\t", "\t")
        if line:
            ids.add(line.split("\t")[0].split()[0])
    return ids


def read_native_qc(path):
    """Read the P8 QC summary without hiding extraction results for review."""
    path = Path(path)
    if not path.is_file():
        return "missing", "", {}
    fields, errors = {}, []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.lstrip().upper().startswith("ERROR:"):
            errors.append(line.split(":", 1)[1].strip())
        elif "=" in line:
            key, value = line.split("=", 1)
            fields[key.strip()] = value.strip()
        elif any(ord(c) > 127 for c in line) and line.strip():
            errors.append("Unrecognized non-English P8 QC line; use ERROR: for error records")
    return ("fail" if errors else "pass"), "; ".join(errors), fields


def build_label_map(csv_labels):
    """Build {FreeSurfer label ID: output column name} and validate label IDs and names.

    Return (fs_to_name, cortex_names_L, cortex_names_R).
    """
    cortical = {k: v for k, v in csv_labels.items()
                if not any(s in v for s in NON_CORTICAL)}
    left = {k: v for k, v in cortical.items() if v.startswith("L_")}
    right = {k: v for k, v in cortical.items() if v.startswith("R_")}

    if len(left) != 33 or len(right) != 33:
        raise SystemExit(f"{DK_LABELS.name}  does not contain exactly 33 cortical regions per hemisphere"
                         f" (left {len(left)} / right {len(right)}); stopping.")

    fs_to_name = {}
    for fs, csv_idx in zip(FS_L[:-1], sorted(left)):
        if csv_idx != fs - 999:
            raise SystemExit(f"Label ID mapping mismatch: FS {fs} -> csv {csv_idx}"
                             f" (expected {fs - 999}); stopping.")
        fs_to_name[fs] = left[csv_idx]
    for fs, csv_idx in zip(FS_R[:-1], sorted(right)):
        if csv_idx != fs - 1964:
            raise SystemExit(f"Label ID mapping mismatch: FS {fs} -> csv {csv_idx}"
                             f" (expected {fs - 1964}); stopping.")
        fs_to_name[fs] = right[csv_idx]
    if FS_L[-1] != INSULA_L or FS_R[-1] != INSULA_R:
        raise SystemExit("The final FS_L/FS_R labels are not insula; mapping constants are inconsistent.")

    lut = {}
    for line in open(FS_LUT, encoding="utf-8", errors="replace"):
        parts = line.split()
        if len(parts) >= 2 and parts[0].isdigit():
            lut[int(parts[0])] = parts[1]
    for fs, expect in zip(FS_L, FS_LONG_L):
        got = lut.get(fs)
        if got != expect:
            raise SystemExit(
                f"FreeSurfer LUT label {fs} is '{got}'; expected '{expect}'. "
                f"Label IDs or names have changed; stopping ({FS_LUT}).")
    for fs in FS_R:
        if lut.get(fs) != lut[fs - 1000].replace("-lh-", "-rh-"):
            raise SystemExit(f"FreeSurfer LUT label {fs} and {fs - 1000} are not symmetric; stopping.")

    if FS_PARC_LABELS.exists():
        model = [int(x) for x in np.load(FS_PARC_LABELS)
                 if 1000 <= int(x) < 3000]
        want = sorted(FS_L + FS_R)
        if sorted(model) != want:
            extra = sorted(set(model) - set(want))
            miss = sorted(set(want) - set(model))
            raise SystemExit(
                f"synthseg --parc model labels do not match the expected labels: "
                f" (extra {extra}; missing {miss}); was FreeSurfer updated?"
                f"; model file: {FS_PARC_LABELS}")

    for fs, nm in ((INSULA_L, "L_insula"), (INSULA_R, "R_insula")):
        if fs in fs_to_name:
            raise SystemExit(f"Insula {fs} was already mapped to {fs_to_name[fs]}; conflicting mappings.")
        fs_to_name[fs] = nm

    cortex_names_L = [fs_to_name[k] for k in FS_L]
    cortex_names_R = [fs_to_name[k] for k in FS_R]
    return fs_to_name, cortex_names_L, cortex_names_R


def reliability(out, names, label):
    """Cross-subject left/right homologous ROI correlations and Spearman-Brown coefficient (as in P5)."""
    lr = {c: "R_" + c[2:] for c in names
          if c.startswith("L_") and "R_" + c[2:] in names}
    rs, detail = [], []
    for lc, rc in lr.items():
        x = pd.to_numeric(out[lc], errors="coerce").to_numpy()
        y = pd.to_numeric(out[rc], errors="coerce").to_numpy()
        m = np.isfinite(x) & np.isfinite(y)
        if m.sum() > 30:
            r = float(np.corrcoef(x[m], y[m])[0, 1])
            rs.append(r)
            detail.append((lc[2:], r))
    if not rs:
        return None
    rbar = float(np.mean(rs))
    sb = 2 * rbar / (1 + rbar) if rbar > -1 else float("nan")
    print(f"\nReliability check ({label}, {len(rs)} left/right homologous ROI pairs): "
          f"mean r={rbar:.3f}  Spearman-Brown coefficient={sb:.3f}")
    print("  Five lowest-correlated pairs:",
          ", ".join(f"{n}={r:.2f}" for n, r in sorted(detail, key=lambda t: t[1])[:5]))
    return rbar, sb, detail


def extract_one(label_path, qsm_path, fs_to_name, min_voxels, support_path=None):
    """Return ROI statistics, valid counts, label counts, and status for one subject.

    When a support mask is provided, its nonzero voxels define valid QSM coverage,
    including legitimate QSM values equal to zero. Without one, nonzero QSM is
    used as a provisional coverage proxy and this limitation is recorded.
    """
    lab_img = nib.load(str(label_path))
    qsm_img = nib.load(str(qsm_path))
    if lab_img.shape[:3] != qsm_img.shape[:3]:
        return None, None, None, (f"Shape mismatch {lab_img.shape[:3]} vs {qsm_img.shape[:3]}")
    if not np.allclose(lab_img.affine, qsm_img.affine, atol=1e-3):
        return None, None, None, ("Label image and QSM affine mismatch")
    lab = np.asanyarray(lab_img.dataobj).astype(np.int32)
    qsm = np.asanyarray(qsm_img.dataobj).astype(np.float32)
    if lab.ndim > 3:
        lab = lab[..., 0]
    if qsm.ndim > 3:
        qsm = qsm[..., 0]

    keys = sorted(fs_to_name)

    lut = np.full(int(max(keys)) + 1, -1, dtype=np.int32)
    lut[np.asarray(keys)] = np.arange(len(keys), dtype=np.int32)

    lab_f = lab.ravel()
    slot = lut[np.clip(lab_f, 0, len(lut) - 1)]
    slot[lab_f > len(lut) - 1] = -1
    q = qsm.ravel()

    keep = slot >= 0
    s_all, q_all = slot[keep], q[keep]
    if support_path is not None:
        support_img = nib.load(str(support_path))
        if support_img.shape[:3] != qsm_img.shape[:3]:
            return None, None, None, ("Support mask and QSM shape mismatch")
        if not np.allclose(support_img.affine, qsm_img.affine, atol=1e-3):
            return None, None, None, ("Support mask and QSM affine mismatch")
        support = np.asanyarray(support_img.dataobj)
        if support.ndim > 3:
            support = support[..., 0]
        supported = support.ravel()[keep] > 0
        supported &= np.isfinite(q_all)
    else:
        # Provisional only: cannot distinguish a valid zero-valued voxel from
        # zero padding without an explicit reconstruction support mask.
        supported = np.isfinite(q_all) & (q_all != 0)

    label_counts = {}
    for i, k in enumerate(keys):
        label_counts[fs_to_name[k]] = int(np.count_nonzero(s_all == i))
    s, v = s_all[supported], q_all[supported]
    if not v.size:
        return None, None, None, "No cortical labels or valid supported cortical voxels"

    order = np.argsort(s, kind="stable")
    s, v = s[order], v[order]
    bounds = np.searchsorted(s, np.arange(len(keys) + 1))

    vals, counts = {}, {}
    for i, k in enumerate(keys):
        seg = v[bounds[i]:bounds[i + 1]]
        counts[fs_to_name[k]] = int(seg.size)
        if seg.size >= min_voxels:
            vals[fs_to_name[k]] = (float(np.mean(seg)), float(np.median(seg)))
        else:
            vals[fs_to_name[k]] = (np.nan, np.nan)
    return vals, counts, label_counts, "OK"


def native_csf_reference(label_path, qsm_path, support_path=None, min_voxels=30):
    """Median QSM of one-voxel-eroded bilateral lateral ventricles (aseg 4/43)."""
    lab_img, qsm_img = nib.load(str(label_path)), nib.load(str(qsm_path))
    if lab_img.shape != qsm_img.shape or not np.allclose(
            lab_img.affine, qsm_img.affine, atol=1e-3):
        return np.nan, 0
    lab = np.asanyarray(lab_img.dataobj)
    qsm = np.asanyarray(qsm_img.dataobj)
    csf_mask = binary_erosion(lab == 4) | binary_erosion(lab == 43)
    if support_path is not None:
        sup_img = nib.load(str(support_path))
        if sup_img.shape != qsm_img.shape or not np.allclose(
                sup_img.affine, qsm_img.affine, atol=1e-3):
            return np.nan, 0
        csf_mask &= np.asanyarray(sup_img.dataobj) > 0
    else:
        # Diagnostic-only proxy, as in extract_one; legitimate zeros are lost.
        csf_mask &= qsm != 0
    values = qsm[csf_mask & np.isfinite(qsm)]
    return (float(np.median(values)) if len(values) >= min_voxels else np.nan,
            int(len(values)))


def main():
    ap = argparse.ArgumentParser(description="Native individualized cortical QSM extraction")
    ap.add_argument("--native-dir", type=Path, default=DEFAULT_NATIVE_DIR)
    ap.add_argument("--qsm-dir", type=Path, default=DEFAULT_NATIVE_DIR,
                    help="Retained option; native QSM paths are read from the worklist")
    ap.add_argument("--out-dir", type=Path, default=HERE / "output")
    ap.add_argument("--worklist", type=Path, default=HERE / "worklist.csv")
    ap.add_argument("--labels", type=Path, default=DK_LABELS)
    ap.add_argument("--failed-list", type=Path, default=HERE / "failed_ID.txt")
    ap.add_argument("--include-failed", action="store_true",
                    help="Include subjects from the failure list (diagnostics only)")
    ap.add_argument("--ids-file", type=Path, default=None,
                    help="Process only subjects in this list (e.g., pilot_native.csv)")
    ap.add_argument("--sample-root", type=Path, default=None,
                    help="Sample root; each subject directory must contain one QSM and one T1 image")
    ap.add_argument("--support-dir", type=Path, default=None,
                    help="Directory of valid QSM support masks named <IID>_support.nii[.gz]")
    ap.add_argument("--allow-nonzero-proxy", action="store_true",
                    help="Without a support mask, use nonzero QSM as a provisional coverage proxy (diagnostics/sensitivity analysis only)")
    ap.add_argument("--which", default="",
                    help="With --ids-file, filter further by 'column=value', e.g., pilot_group=primary_pair")
    ap.add_argument("--stat", choices=("both", "median", "mean"), default="both")
    ap.add_argument("--min-voxels", type=int, default=20)
    ap.add_argument("--csf-ref", action="store_true",
                    help="Also save cortical values referenced to each subject's lateral-ventricle CSF median")
    ap.add_argument("--min-csf-voxels", type=int, default=30)
    ap.add_argument("--suffix", default="")
    ap.add_argument("--jobs", type=int, default=16)
    args = ap.parse_args()
    if args.support_dir is None and not args.allow_nonzero_proxy:
        ap.error("Formal extraction requires --support-dir; use --allow-nonzero-proxy explicitly for diagnostics only")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    label_dir = args.native_dir / "label"

    csv_labels = read_labels(args.labels)
    fs_to_name, _, _ = build_label_map(csv_labels)
    keys = sorted(fs_to_name)
    print(f"Cortical label mapping OK: {len(keys)} regions (34 per hemisphere, including insula; "
          f"{INSULA_L}/{INSULA_R} in {args.labels.name} has no corresponding CSV entry and is added in this branch)")

    wl = pd.read_csv(args.worklist, dtype=str).fillna("")
    # Worklist input headers must be English; no legacy-language aliases are used.
    if "image_dir_id" not in wl.columns and "IID" in wl.columns:
        wl["image_dir_id"] = wl["IID"]
    required = {"image_dir_id", "qsm_file"}
    absent = sorted(required - set(wl.columns))
    if absent:
        raise SystemExit(f"Missing worklist columns: {absent}. Rename input headers to English.")
    if "runnable" in wl.columns:
        wl = wl[wl["runnable"].str.lower() == "yes"]
    wl = wl.drop_duplicates("image_dir_id")
    failed = read_failed_ids(args.failed_list)
    if args.include_failed:
        print(f"Included failed-list IDs (diagnostic mode): {len(set(wl['image_dir_id']) & failed)} subjects")
    else:
        wl = wl[~wl["image_dir_id"].isin(failed)]

    if args.sample_root:
        if not args.sample_root.is_dir():
            raise SystemExit(f"Sample directory does not exist: {args.sample_root}")
        sample_ids = sorted(p.name for p in args.sample_root.iterdir() if p.is_dir())
        wl = wl[wl["image_dir_id"].isin(sample_ids)]
        if not sample_ids:
            raise SystemExit(f"No subject subdirectories in sample root: {args.sample_root}")
        print(f"Filtered by sample root: {len(wl)} subjects")

    if args.ids_file:
        p = pd.read_csv(args.ids_file, dtype=str).fillna("")
        idcol = "IID" if "IID" in p.columns else "image_dir_id"
        ids = p[idcol].tolist()
        if args.which:
            k, v = args.which.split("=", 1)
            ids = p.loc[p[k] == v, idcol].tolist()
        wl = wl[wl["image_dir_id"].isin(ids)]
        print(f"Filtered by ID list: {len(wl)} subjects")

    qsm_of = dict(zip(wl["image_dir_id"], wl["qsm_file"]))
    meta = wl.set_index("image_dir_id")

    ids, paths, missing = [], [], []
    for oid in wl["image_dir_id"]:
        lab = label_dir / f"{oid}_aparc_aseg_QSMnative.nii.gz"
        qsm = Path(qsm_of.get(oid, ""))
        if args.sample_root:
            sample_dir = args.sample_root / oid
            qsm_candidates = [p for p in sample_dir.glob("*.nii*")
                              if "qsm" in p.name.lower() and "swi" not in p.name.lower()]
            t1_candidates = [p for p in sample_dir.glob("*.nii*")
                              if "t1" in p.name.lower()]
            if len(qsm_candidates) != 1 or len(t1_candidates) != 1:
                missing.append(f"{oid} (QSM={len(qsm_candidates)}, T1={len(t1_candidates)})")
                continue
            qsm = qsm_candidates[0]
        support = None
        if args.support_dir:
            support_candidates = [args.support_dir / f"{oid}_support.nii.gz",
                                  args.support_dir / f"{oid}_support.nii"]
            support = next((p for p in support_candidates if p.is_file()), None)
            if support is None:
                missing.append(f"{oid} (support mask missing)")
                continue
        if lab.is_file() and lab.stat().st_size > 0 and qsm.is_file():
            ids.append(oid)
            paths.append((lab, qsm, support))
        else:
            missing.append(oid)
    print(f"\nSubjects with native labels: {len(ids)}; missing outputs: {len(missing)}"
          + (f" (first 10: {missing[:10]})" if missing else ""))
    if not ids:
        raise SystemExit("No eligible subjects. Complete P8 first.")

    print(f"Extracting with {args.jobs} threads...")
    with ThreadPoolExecutor(args.jobs) as ex:
        results = list(ex.map(
            lambda t: extract_one(t[0], t[1], fs_to_name, args.min_voxels, t[2]), paths))

    counts = np.zeros((len(ids), len(keys)), dtype=np.int32)
    label_counts = np.zeros((len(ids), len(keys)), dtype=np.int32)
    status, detail = [], []
    for i, (vals, cnt, denom, st) in enumerate(results):
        status.append(st)
        detail.append("" if st == "OK" else st)
        if vals is None:
            continue
        for j, k in enumerate(keys):
            nm = fs_to_name[k]
            counts[i, j] = cnt.get(nm, 0)
            label_counts[i, j] = denom.get(nm, 0)

    names_med = [fs_to_name[k] for k in keys]
    value_cols = {f"{nm}_med": np.full(len(ids), np.nan) for nm in names_med}
    if args.stat in ("both", "mean"):
        for nm in names_med:
            value_cols[f"{nm}_mean"] = np.full(len(ids), np.nan)
    for i, (vals, _cnt, _denom, _st) in enumerate(results):
        if vals is None:
            continue
        for nm in names_med:
            value_cols[f"{nm}_med"][i] = vals[nm][1]
            if f"{nm}_mean" in value_cols:
                value_cols[f"{nm}_mean"][i] = vals[nm][0]

    roi_out = pd.DataFrame(value_cols)
    if args.stat in ("both", "mean"):
        # Keep median/mean columns adjacent for each ROI.
        ordered = []
        for nm in names_med:
            ordered.append(f"{nm}_med")
            ordered.append(f"{nm}_mean")
        roi_out = roi_out[ordered]

    # Metadata first; build in one operation to avoid a fragmented DataFrame.
    sub = meta.reindex(ids)
    meta_cols = [c for c in META_COLS if c in sub.columns]
    meta_out = sub[meta_cols].reset_index(drop=True)
    meta_out.insert(0, "IID", ids)
    out = pd.concat([meta_out, roi_out], axis=1)

    sfx = args.suffix
    if args.sample_root and not sfx:
        sfx = "_sample_debug"
    csv_path = args.out_dir / f"HUASHAN_NATIVE_DK{sfx}.csv"
    out.to_csv(csv_path, index=False, encoding="utf-8-sig")

    csf_values = np.full(len(ids), np.nan)
    csf_counts = np.zeros(len(ids), dtype=np.int32)
    if args.csf_ref:
        with ThreadPoolExecutor(args.jobs) as ex:
            references = list(ex.map(
                lambda t: native_csf_reference(t[0], t[1], t[2],
                                               args.min_csf_voxels), paths))
        csf_values[:] = [x[0] for x in references]
        csf_counts[:] = [x[1] for x in references]
        ref_roi = roi_out.subtract(csf_values, axis=0)
        ref_out = pd.concat([meta_out, ref_roi], axis=1)
        ref_path = args.out_dir / f"HUASHAN_NATIVE_DK{sfx}_csfref.csv"
        ref_out.to_csv(ref_path, index=False, encoding="utf-8-sig")
        print(f"CSF-referenced results: {ref_path}")

    p8_qc = [read_native_qc(args.native_dir / "qc" / f"{oid}_native_qc.txt")
             for oid in ids]
    p1_failed = np.asarray([oid in failed for oid in ids], dtype=bool)
    pd.DataFrame({
        "IID": ids,
        "extract_status": status,
        "extract_detail": detail,
        "P1_failed_list": p1_failed,
        "P8_native_qc_status": [x[0] for x in p8_qc],
        "P8_native_qc_errors": [x[1] for x in p8_qc],
        "P8_ribbon_in_qsm_frac": [x[2].get("ribbon_in_qsm_frac", "") for x in p8_qc],
        "csf_ref_ppb": csf_values,
        "csf_voxels": csf_counts,
        "label_voxels_without_qsm_support": (label_counts - counts).sum(axis=1),
        "n_nan_DK": out[[f"{n}_med" for n in names_med]].isna().sum(axis=1).to_numpy(),
        "n_low_voxels": (counts < args.min_voxels).sum(axis=1),
    }).to_csv(args.out_dir / f"HUASHAN_NATIVE_DK{sfx}_qc.csv", index=False,
              encoding="utf-8-sig")

    used = counts.sum(axis=0)
    atlas_voxels = label_counts.sum(axis=0)
    pd.DataFrame({
        "ROI": names_med,
        "atlas_label_voxels": atlas_voxels,
        "used_qsm_voxels": used,
        "excluded_cerebellar_voxel_fraction": 0.0,
    }).to_csv(args.out_dir / f"native_mask_voxel_counts{sfx}.csv", index=False,
              encoding="utf-8-sig")

    sizes = label_counts.astype(float)
    cov = np.divide(counts, sizes, out=np.zeros(counts.shape), where=sizes > 0)
    pd.DataFrame({
        "IID": np.repeat(ids, len(keys)),
        "ROI": np.tile(names_med, len(ids)),
        "atlas_label_voxels": label_counts.ravel(),
        "valid_qsm_voxels": counts.ravel(),
        "coverage_fraction": cov.ravel().round(4),
        "coverage_definition": ("explicit_support_mask" if args.support_dir
                                else "nonzero_QSM_proxy; valid zeros may be excluded"),
    }).to_csv(args.out_dir / f"native_roi_coverage{sfx}.csv", index=False,
              encoding="utf-8-sig")

    n_ok = sum(1 for s in status if s == "OK")
    analysis_pass = np.asarray([(args.include_failed or not p1_failed[i])
                                and p8_qc[i][0] == "pass"
                                and (not args.csf_ref or np.isfinite(csf_values[i]))
                                for i in range(len(ids))])
    print(f"QC-passing subjects for downstream analysis: {int(analysis_pass.sum())}/{len(ids)};"
          f"other subjects remain in the output but are not QC-passing")
    print(f"\nCompleted: {csv_path}  ({out.shape[0]} subjects x {len(names_med)} regions, "
          f"successful extractions: {n_ok}/{len(ids)})")
    if len(ids) - n_ok:
        bad = [(ids[i], status[i]) for i in range(len(ids)) if status[i] != "OK"]
        print(f"  Failed: {len(bad)} subjects (first 5): {bad[:5]}")

    per_roi = pd.DataFrame(counts[analysis_pass], columns=names_med).mean(axis=0)
    print(f"\nMean voxel count by ROI: median {per_roi.median():.0f},"
          f"minimum {per_roi.min():.0f}({per_roi.idxmin()})")
    low = [n for n in names_med if per_roi[n] < args.min_voxels * 2]
    if low:
        print(f"  WARNING Regions with low voxel counts: {len(low)} regions: {low[:8]}")

    print(f"\nCross-subject median QSM by ROI (ppb; P1/P8 QC-passing subjects only):")
    for c in names_med:
        v = pd.to_numeric(out.loc[analysis_pass, f"{c}_med"], errors="coerce").dropna()
        if v.size:
            print(f"  {c:<34} {v.median():>7.2f}  SD={v.std():>6.2f}  n={v.size}")

    analysis_out = out.loc[analysis_pass].copy()
    reliability(analysis_out.rename(columns={f"{n}_med": n for n in names_med}),
                names_med, "native")
    if f"{names_med[0]}_mean" in out.columns:
        reliability(analysis_out.rename(columns={f"{n}_mean": n for n in names_med}),
                    names_med, "native (mean statistic)")

    if missing:
        print(f"\nIncomplete: {len(missing)} subjects (first 10): {missing[:10]}")


if __name__ == "__main__":
    main()
