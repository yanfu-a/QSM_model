#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Audit existing P8 native labels for the T1 -> T1_brain grid-conversion error.

The original P8 converted between the full-resolution T1 grid and the kept
T1_brain grid without FLIRT's x-flip for neurologically stored images, so for
affected subjects T1_to_QSM.mat (and the label resampled with it) is offset by
one T1 voxel along the first voxel axis. For each worklist subject with a kept
T1_brain, this script compares the original and corrected conversions (image
headers only) and, when T1_to_QSM.mat exists, identifies which one produced it.

The revised P8 regenerates untagged labels automatically; this audit tells you
how many existing labels, and any results derived from them, were affected.
"""

import argparse
import sys
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from native_grid import flirt_grid_conversion, is_neurological, legacy_grid_conversion  # noqa: E402

DEFAULT_NATIVE_DIR = Path("/cwStorage/nodecw_group/FY_data/QSM_HUASHAN/native_data")


def fmt(vec):
    return "x".join(f"{round(v, 4) + 0.0:.4f}" for v in vec)


def audit_subject(oid, t1_path, mat_dir):
    """One audit row for a subject."""
    row = {"IID": oid, "t1_file": t1_path}
    brain = mat_dir / f"{oid}_T1_brain.nii.gz"
    if not brain.is_file():
        return {**row, "status": "no kept T1_brain"}
    if not Path(t1_path).is_file():
        return {**row, "status": "T1 file missing"}
    full, sub = nib.load(t1_path), nib.load(str(brain))
    fixed = flirt_grid_conversion(full, sub)
    legacy = legacy_grid_conversion(full, sub)
    offset_mm = float(np.abs(fixed - legacy).max())
    row.update(status="ok", t1_neurological=is_neurological(full),
               t1_shape="x".join(str(n) for n in full.shape[:3]),
               legacy_shift_mm=fmt(legacy[:3, 3]), fixed_shift_mm=fmt(fixed[:3, 3]),
               legacy_offset_mm=round(offset_mm, 4))

    saved, q2t = mat_dir / f"{oid}_T1_to_QSM.mat", mat_dir / f"{oid}_QSM_to_T1.mat"
    if saved.is_file() and q2t.is_file():
        m_saved, m_q2t_inv = np.loadtxt(saved), np.linalg.inv(np.loadtxt(q2t))
        by_legacy = np.allclose(m_saved, m_q2t_inv @ legacy, atol=1e-5)
        by_fixed = np.allclose(m_saved, m_q2t_inv @ fixed, atol=1e-5)
        built = ("both (identical here)" if by_legacy and by_fixed else
                 "legacy" if by_legacy else "fixed" if by_fixed else "unknown")
    else:
        built = "no T1_to_QSM.mat"
    row["saved_matrix_formula"] = built
    # An offset matters only for a label that was actually built with the legacy matrix.
    row["label_offset"] = bool(offset_mm > 1e-3 and built in ("legacy", "unknown"))
    return row


def main():
    ap = argparse.ArgumentParser(description="Audit P8 labels for the T1 -> T1_brain grid-conversion error")
    ap.add_argument("--worklist", type=Path, default=HERE / "worklist.csv")
    ap.add_argument("--native-dir", type=Path, default=DEFAULT_NATIVE_DIR)
    ap.add_argument("--out", type=Path, default=HERE / "output" / "native_grid_audit.csv")
    args = ap.parse_args()

    wl = pd.read_csv(args.worklist, dtype=str).fillna("")
    if "runnable" in wl.columns:
        wl = wl[wl["runnable"].str.strip().str.lower() == "yes"]
    wl = wl.drop_duplicates("image_dir_id")
    mat_dir = args.native_dir / "mat"

    rows = []
    for oid, t1 in zip(wl["image_dir_id"], wl["t1_file"]):
        try:
            rows.append(audit_subject(oid, t1, mat_dir))
        except Exception as exc:                                   # noqa: BLE001
            rows.append({"IID": oid, "t1_file": t1, "status": f"ERROR {type(exc).__name__}: {exc}"})
    audit = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    audit.to_csv(args.out, index=False, encoding="utf-8-sig")

    ok = audit[audit["status"] == "ok"]
    print(f"Audited {len(audit)} subjects; {len(ok)} with a kept T1_brain ({args.out})")
    if len(ok):
        print(f"  neurologically stored T1: {int(ok['t1_neurological'].sum())}")
        print(f"  legacy conversion differs from FLIRT's: {int((ok['legacy_offset_mm'] > 1e-3).sum())} "
              f"(max {ok['legacy_offset_mm'].max():.3f} mm)")
        print("  saved T1_to_QSM.mat built by:", ok["saved_matrix_formula"].value_counts().to_dict())
        print(f"  existing labels offset (rerun P8): {int(ok['label_offset'].sum())}")
    other = audit["status"].value_counts().drop("ok", errors="ignore")
    if len(other):
        print("  not audited:", other.to_dict())


if __name__ == "__main__":
    main()
