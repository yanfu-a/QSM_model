#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Write QSM support masks (<IID>_support.nii.gz) for P9's --support-dir.

The support mask is the region where the QSM reconstruction produced values. When
the reconstruction pipeline exports its own brain/ROI mask, use that mask instead of
this script. The cohort's vendor QSM maps come without one, so this script recovers
the support from each QSM image: the nonzero region plus zero-valued holes enclosed
by it that are no larger than --max-hole-mm3. Small enclosed zeros are genuine values
(common when QSM is stored as quantized integers); large or open zero regions (masked
CSF, outside the brain) remain unsupported. A summary CSV records, per subject, how
many voxels were added, so the effect of the hole size can be reviewed.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from scipy import ndimage

DEFAULT_NATIVE_DIR = Path("/cwStorage/nodecw_group/FY_data/QSM_HUASHAN/native_data")
HERE = Path(__file__).resolve().parent


def support_from_qsm(qsm, vox_mm3, max_hole_mm3):
    """Nonzero finite QSM plus enclosed zero holes of at most max_hole_mm3.

    Return (support, n_holes_filled, filled_voxels, largest_filled_hole_mm3).
    """
    nonzero = np.isfinite(qsm) & (qsm != 0)
    holes, n_comp = ndimage.label(~nonzero)       # face-connected zero regions
    if n_comp == 0:
        return nonzero, 0, 0, 0.0
    sizes = np.bincount(holes.ravel(), minlength=n_comp + 1)
    # Regions touching the array border are outside the reconstruction, never holes.
    border = np.unique(np.concatenate([
        holes[0].ravel(), holes[-1].ravel(), holes[:, 0].ravel(), holes[:, -1].ravel(),
        holes[:, :, 0].ravel(), holes[:, :, -1].ravel()]))
    fill = (sizes * vox_mm3 <= max_hole_mm3)
    fill[0] = False
    fill[border] = False
    support = nonzero | fill[holes]
    filled = int(sizes[fill].sum())
    largest = float(sizes[fill].max() * vox_mm3) if fill.any() else 0.0
    return support, int(fill.sum()), filled, largest


def make_one(oid, qsm_path, out_dir, max_hole_mm3, overwrite):
    """Write one subject's support mask; return its summary row."""
    out = out_dir / f"{oid}_support.nii.gz"
    row = {"IID": oid, "qsm_file": qsm_path, "support_mask": str(out)}
    if out.is_file() and not overwrite:
        return {**row, "status": "exists (use --overwrite to rebuild)"}
    try:
        img = nib.load(qsm_path)
        qsm = np.asanyarray(img.dataobj)
        if qsm.ndim > 3:
            qsm = qsm[..., 0]
        vox_mm3 = float(abs(np.linalg.det(img.affine[:3, :3])))
        support, n_holes, filled, largest = support_from_qsm(qsm, vox_mm3, max_hole_mm3)
        # Same grid as the QSM, so P9's shape/affine check accepts it.
        mask = nib.Nifti1Image(support.astype(np.uint8), img.affine)
        mask.set_sform(img.affine, int(img.header["sform_code"]) or 1)
        mask.set_qform(img.affine, int(img.header["qform_code"]) or 1)
        nib.save(mask, str(out))
        nonzero = int((np.isfinite(qsm) & (qsm != 0)).sum())
        return {**row, "status": "ok", "qsm_dtype": str(img.get_data_dtype()),
                "nonzero_voxels": nonzero, "holes_filled": n_holes, "filled_voxels": filled,
                "largest_filled_hole_mm3": round(largest, 2), "support_voxels": int(support.sum()),
                "filled_fraction": round(filled / max(nonzero, 1), 6)}
    except Exception as exc:                                       # noqa: BLE001
        return {**row, "status": f"ERROR {type(exc).__name__}: {exc}"}


def main():
    ap = argparse.ArgumentParser(description="Write QSM support masks for P9 --support-dir")
    ap.add_argument("--worklist", type=Path, default=HERE / "worklist.csv")
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_NATIVE_DIR / "support")
    ap.add_argument("--ids-file", type=Path, default=None,
                    help="Only these subjects (IID or image_dir_id column)")
    ap.add_argument("--max-hole-mm3", type=float, default=10.0,
                    help="Largest enclosed zero region treated as genuine QSM values")
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--jobs", type=int, default=8)
    args = ap.parse_args()

    wl = pd.read_csv(args.worklist, dtype=str).fillna("")
    if "runnable" in wl.columns:
        wl = wl[wl["runnable"].str.strip().str.lower() == "yes"]
    wl = wl.drop_duplicates("image_dir_id")
    if args.ids_file:
        p = pd.read_csv(args.ids_file, dtype=str).fillna("")
        idcol = next((c for c in ("IID", "image_dir_id") if c in p.columns), None)
        if idcol is None:
            raise SystemExit(f"{args.ids_file} has neither an IID nor an image_dir_id column")
        wl = wl[wl["image_dir_id"].isin(p[idcol])]
    args.out_dir.mkdir(parents=True, exist_ok=True)

    with ThreadPoolExecutor(args.jobs) as ex:
        rows = list(ex.map(lambda t: make_one(t[0], t[1], args.out_dir, args.max_hole_mm3,
                                              args.overwrite),
                           zip(wl["image_dir_id"], wl["qsm_file"])))
    summary = pd.DataFrame(rows)
    summary["max_hole_mm3"] = args.max_hole_mm3
    summary_path = args.out_dir / "support_masks_summary.csv"
    if summary_path.is_file() and not args.overwrite:
        # Keep earlier rows for subjects not rebuilt in this run.
        old = pd.read_csv(summary_path, dtype=str)
        rebuilt = set(summary.loc[summary["status"] == "ok", "IID"])
        summary = pd.concat([old[~old["IID"].isin(rebuilt)], summary[summary["IID"].isin(rebuilt)]])
    summary.to_csv(summary_path, index=False, encoding="utf-8-sig")

    status = pd.Series([r["status"] for r in rows])
    print(f"Support masks in {args.out_dir}: {int((status == 'ok').sum())} written, "
          f"{int(status.str.startswith('exists').sum())} kept, "
          f"{int(status.str.startswith('ERROR').sum())} failed (summary: {summary_path.name})")
    done = pd.DataFrame([r for r in rows if r["status"] == "ok"])
    if len(done):
        print(f"  enclosed zero voxels added: median fraction {done['filled_fraction'].median():.2e}, "
              f"max {done['filled_fraction'].max():.2e}; largest hole filled "
              f"{done['largest_filled_hole_mm3'].max():.1f} mm3 (limit {args.max_hole_mm3:g})")


if __name__ == "__main__":
    main()
