#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Pick a pilot set of scans that spans the cohort's image geometries.

Strata combine T1 storage orientation (neurological/radiological), T1 voxel size,
the parity of the T1 first dimension (which decided whether the original P8 grid
conversion was offset), QSM voxel size, and QSM axial coverage (< or >= the
--short-coverage-mm cut). Up to --per-stratum subjects are drawn from each stratum
(deterministically, by ID). Only image headers are read. The output CSV feeds the
--ids-file options of make_support_masks.py and P9 and lists the IDs to run through
P8 first.
"""

import argparse
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd

from strata import SHORT_COVERAGE_MM, coverage_group, si_extent_mm, voxel_size_label
from worklist_ids import resolve_duplicate_ids

HERE = Path(__file__).resolve().parent


def geometry(path):
    """Header-only geometry summary of one image."""
    img = nib.load(path)
    det = np.linalg.det(img.affine[:3, :3])
    return {"neurological": bool(det > 0), "pixdim": voxel_size_label(img.header.get_zooms()),
            "first_dim_even": img.shape[0] % 2 == 0, "si_extent_mm": si_extent_mm(img)}


def main():
    ap = argparse.ArgumentParser(description="Select representative pilot scans")
    ap.add_argument("--worklist", type=Path, default=HERE / "worklist.csv")
    ap.add_argument("--out", type=Path, default=HERE / "pilot_ids.csv")
    ap.add_argument("--per-stratum", type=int, default=2)
    ap.add_argument("--short-coverage-mm", type=float, default=SHORT_COVERAGE_MM)
    args = ap.parse_args()

    wl = pd.read_csv(args.worklist, dtype=str).fillna("")
    wl, conflicts = resolve_duplicate_ids(wl)
    if conflicts:
        print(f"Skipped {len(conflicts)} image_dir_id values with conflicting worklist rows: "
              f"{[e['IID'] for e in conflicts]}")
    if "runnable" in wl.columns:
        wl = wl[wl["runnable"].str.strip().str.lower() == "yes"]

    rows = []
    for _, r in wl.iterrows():
        try:
            t1, q = geometry(r["t1_file"]), geometry(r["qsm_file"])
        except Exception as exc:                                   # noqa: BLE001
            print(f"  skipped {r['image_dir_id']}: {type(exc).__name__}: {exc}")
            continue
        coverage = pd.to_numeric(r.get("qsm_coverage_mm", ""), errors="coerce")
        if not np.isfinite(coverage):
            coverage = q["si_extent_mm"]
        rows.append({"IID": r["image_dir_id"],
                     "t1_orientation": "neurological" if t1["neurological"] else "radiological",
                     "t1_pixdim": t1["pixdim"], "t1_first_dim": "even" if t1["first_dim_even"] else "odd",
                     "qsm_pixdim": q["pixdim"],
                     "qsm_coverage": coverage_group(coverage, args.short_coverage_mm),
                     "qsm_coverage_mm": round(float(coverage), 1)})
    geo = pd.DataFrame(rows)
    if geo.empty:
        raise SystemExit("No readable scans in the worklist")
    keys = ["t1_orientation", "t1_pixdim", "t1_first_dim", "qsm_pixdim", "qsm_coverage"]
    geo["stratum"] = geo[keys].agg(" | ".join, axis=1)
    pilot = (geo.sort_values("IID").groupby("stratum", sort=True).head(args.per_stratum)
             .assign(pilot_group=lambda d: d["stratum"]))
    pilot.to_csv(args.out, index=False, encoding="utf-8-sig")

    sizes = geo["stratum"].value_counts()
    print(f"{len(geo)} scans in {len(sizes)} strata; pilot of {len(pilot)} written to {args.out}")
    for stratum, n in sizes.items():
        print(f"  {n:5d}  {stratum}")


if __name__ == "__main__":
    main()
