#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Paired comparison of two P9 runs on the same scans (sensitivity analyses).

Typical uses (README, "Sensitivity analyses"):
  * support-mask hole limit: masks rebuilt with make_support_masks.py --max-hole-mm3 0 (or 50)
    into another directory, P9 run on them with --suffix, and the two *_analysis.csv compared;
  * registration resolution: P8 run with QSM_ROI_T1_SUBSAMP=0 into a second native root, P9
    run on it, and the result compared with the default (downsampled-T1) registration.

Per ROI, overall and per protocol stratum (QSM voxel size | coverage, from run A's *_qc.csv
when given): scans with a value in both runs or in one only (changes in missingness), mean
and SD of the difference B - A, mean and maximum absolute difference, and Pearson r. Scans
present in only one of the two tables are counted as well (changes in the analysis set).
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def read_table(path):
    df = pd.read_csv(path, dtype={"IID": str})
    if "IID" not in df.columns:
        raise SystemExit(f"{path} has no IID column (expected a P9 output table)")
    if df["IID"].duplicated().any():
        raise SystemExit(f"{path} has repeated IIDs")
    return df.set_index("IID")


def compare(a, b, rois, stat, strata=None):
    """Long table of per-ROI paired statistics, overall and per stratum."""
    common = a.index.intersection(b.index)
    levels = [("all", common)]
    if strata is not None:
        s = strata.reindex(common).fillna("unknown")
        levels += [(level, s.index[s == level]) for level in sorted(s.unique())]
    rows = []
    for level, ids in levels:
        for roi in rois:
            col = f"{roi}_{stat}"
            x = pd.to_numeric(a.loc[ids, col], errors="coerce").to_numpy(float)
            y = pd.to_numeric(b.loc[ids, col], errors="coerce").to_numpy(float)
            fx, fy = np.isfinite(x), np.isfinite(y)
            both = fx & fy
            d = y[both] - x[both]
            r = (float(np.corrcoef(x[both], y[both])[0, 1])
                 if both.sum() >= 3 and np.std(x[both]) > 0 and np.std(y[both]) > 0 else np.nan)
            rows.append({"stratum": level, "ROI": roi, "n_scans": len(ids), "n_both": int(both.sum()),
                         "n_only_a": int((fx & ~fy).sum()), "n_only_b": int((fy & ~fx).sum()),
                         "mean_diff_b_minus_a": float(d.mean()) if d.size else np.nan,
                         "sd_diff": float(d.std(ddof=1)) if d.size > 1 else np.nan,
                         "mean_abs_diff": float(np.abs(d).mean()) if d.size else np.nan,
                         "max_abs_diff": float(np.abs(d).max()) if d.size else np.nan,
                         "pearson_r": r})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description="Compare ROI values of two P9 runs scan by scan")
    ap.add_argument("--a", type=Path, required=True, help="Reference P9 table (e.g. HUASHAN_NATIVE_DK_analysis.csv)")
    ap.add_argument("--b", type=Path, required=True, help="Alternative P9 table from the sensitivity run")
    ap.add_argument("--qc", type=Path, default=None,
                    help="Run A's *_qc.csv, for protocol strata (qsm_protocol, qsm_coverage_group)")
    ap.add_argument("--stat", choices=("med", "mean"), default="med")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    a, b = read_table(args.a), read_table(args.b)
    suffix = f"_{args.stat}"
    rois = [c[:-len(suffix)] for c in a.columns if c.endswith(suffix) and c in b.columns]
    if not rois:
        raise SystemExit(f"No common *{suffix} ROI columns in {args.a.name} and {args.b.name}")
    strata = None
    if args.qc:
        qc = pd.read_csv(args.qc, dtype=str).fillna("").set_index("IID")
        need = {"qsm_protocol", "qsm_coverage_group"}
        if not need <= set(qc.columns):
            raise SystemExit(f"{args.qc} lacks {sorted(need - set(qc.columns))} (rerun the current P9)")
        strata = qc["qsm_protocol"] + " | " + qc["qsm_coverage_group"]

    table = compare(a, b, rois, args.stat, strata)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.out, index=False, encoding="utf-8-sig")

    only_a, only_b = a.index.difference(b.index), b.index.difference(a.index)
    print(f"Scans in both tables: {len(a.index.intersection(b.index))}; only in A: {len(only_a)}; "
          f"only in B: {len(only_b)}")
    for name, ids in (("A", only_a), ("B", only_b)):
        if len(ids):
            print(f"  only in {name} (first 10): {list(ids[:10])}")
    overall = table[table["stratum"] == "all"]
    print(f"Per-ROI {args.stat} over {len(rois)} ROIs: median |B-A| {overall['mean_abs_diff'].median():.3g}, "
          f"median r {overall['pearson_r'].median():.3f}; ROI values gained {int(overall['n_only_b'].sum())}, "
          f"lost {int(overall['n_only_a'].sum())}")
    worst = overall.sort_values("mean_abs_diff", ascending=False).head(5)
    print("  Largest mean |B-A|: " + ", ".join(f"{r.ROI}={r.mean_abs_diff:.3g}" for r in worst.itertuples()))
    print(f"Written: {args.out}")


if __name__ == "__main__":
    main()
