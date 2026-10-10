#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Count P8 QC exclusions by rule, overall and by protocol, coverage, diagnosis and sex.

Report only: thresholds are applied by P8, which records in each QC file the rules a scan
failed (qc_failed_rules) and the failures of candidate rules that are evaluated but not
applied (qc_report_only_failed_rules). For scans where an earlier P8 stopped at P1's Dice
gate, P1's rigid QC is read instead. Use the tables to see whether a rule removes scans of
one protocol, diagnosis or sex before deciding on thresholds (METHODOLOGICAL_DECISIONS D1).

Outputs in --out-dir:
  qc_rules_by_subject.csv     one row per scan: strata, stage reached, failed rules, metrics
  qc_exclusions_by_rule.csv   per stratifier level and rule: scans evaluated, failed, and
                              failed by that rule alone (for a report-only rule: failed while
                              passing every applied rule, i.e. exclusions it would add); scans
                              that stopped at P1's Dice gate were evaluated on that rule only
  qc_cortex_rule_comparison.csv  cortex-volume rule on the QSM grid (applied) vs on the
                              whole-brain T1 segmentation (candidate), per stratifier level

Strata: protocol (QSM voxel size), coverage (short/full axial slab; strata.py), diagnosis
(worklist model_label) and any column given with --by, from the worklist or a --covariates
file keyed by IID/image_dir_id or hospital_id (e.g. sex, which the worklist lacks).
"""

import argparse
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd

from strata import SHORT_COVERAGE_MM, coverage_group, si_extent_mm, voxel_size_label
from worklist_ids import resolve_duplicate_ids

HERE = Path(__file__).resolve().parent
DEFAULT_NATIVE_DIR = Path("/cwStorage/nodecw_group/FY_data/QSM_HUASHAN/native_data")

# P8 rule name -> (metric written to the P8 QC file, applied by P8)
RULES = {
    "reg_mask_dice": ("reg_mask_dice", True),
    "synthseg_icv": ("synthseg_icv_ml", True),
    "cortex_vol_qsm_grid": ("cortex_vol_ml", True),
    "ribbon_coverage": ("ribbon_in_qsm_frac", True),
    "gm_wm_contrast": ("gm_wm_contrast_ppb", True),
    "seg_qsm_dice": ("seg_qsm_dice", True),
    "cortex_vol_t1": ("synthseg_cortex_ml", False),
}
APPLIED = [r for r, (_, applied) in RULES.items() if applied]
METRICS = [m for m, _ in RULES.values()]


def read_fields(path):
    """key=value lines of a QC file, and its WARNING/ERROR lines."""
    fields, flagged = {}, []
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("WARNING:", "ERROR:")):
            flagged.append(line)
        elif "=" in line:
            key, value = line.split("=", 1)
            fields[key.strip()] = value.strip()
    return fields, flagged


def rule_list(text):
    return [r for r in (text or "").split(",") if r]


def finite(value):
    try:
        return bool(np.isfinite(float(value)))
    except (TypeError, ValueError):
        return False


def subject_record(oid, native_dir):
    """Stage reached and rule outcomes of one scan."""
    qc = native_dir / "qc" / f"{oid}_native_qc.txt"
    p1 = native_dir / "p1_rigid" / f"{oid}_qc_rigid.txt"
    rec = {"IID": oid, "stage": "no_p8_qc", "failed_rules": "", "report_only_failed_rules": ""}
    if qc.is_file():
        fields, _ = read_fields(qc)
        rec.update({m: fields.get(m, "") for m in METRICS})
        rec["p8_qsm_pixdim"] = fields.get("qsm_pixdim", "")
        if fields.get("qc_complete") != "1":
            rec["stage"] = "p8_qc_incomplete"
        elif "qc_failed_rules" not in fields:
            rec["stage"] = "p8_qc_without_rule_record"          # built before rules were recorded
        else:
            rec.update(stage="p8_qc", failed_rules=fields["qc_failed_rules"],
                       report_only_failed_rules=fields.get("qc_report_only_failed_rules", ""))
    elif p1.is_file():
        fields, flagged = read_fields(p1)
        if any("reg_mask_dice" in line for line in flagged):
            # Earlier P8 versions stopped where P1's Dice gate failed; only that rule was evaluated.
            rec.update(stage="p1_dice_stop", failed_rules="reg_mask_dice",
                       reg_mask_dice=fields.get("reg_mask_dice", ""))
    return rec


def qsm_strata(qsm_file, p8_pixdim, worklist_coverage):
    """(protocol, coverage group, coverage mm) from the QSM header, else recorded values."""
    protocol, extent = "unknown", np.nan
    try:
        img = nib.load(qsm_file)
        protocol, extent = voxel_size_label(img.header.get_zooms()), si_extent_mm(img)
    except Exception:                                              # noqa: BLE001
        if p8_pixdim:
            protocol = voxel_size_label(p8_pixdim)
    cov = pd.to_numeric(worklist_coverage, errors="coerce")
    cov = float(cov) if np.isfinite(cov) else extent
    return protocol, cov


def evaluate(df):
    """Per scan: evaluated/failed flags for every rule (columns ev_<rule>, fail_<rule>)."""
    failed = df["failed_rules"].map(rule_list)
    report_only = df["report_only_failed_rules"].map(rule_list)
    with_qc = df["stage"].eq("p8_qc")
    for rule, (metric, applied) in RULES.items():
        fails = (failed if applied else report_only).map(lambda rs, r=rule: r in rs)
        if rule == "reg_mask_dice":
            # P8 fails a missing Dice; the stopped P1 runs evaluated only this rule.
            evaluated = with_qc | df["stage"].eq("p1_dice_stop")
        else:
            evaluated = with_qc & (df.get(metric, pd.Series("", index=df.index)).map(finite) | fails)
        df[f"ev_{rule}"] = evaluated
        df[f"fail_{rule}"] = evaluated & fails
    n_applied_failed = df[[f"fail_{r}" for r in APPLIED]].sum(axis=1)
    df["qc_evaluated"] = with_qc | df["stage"].eq("p1_dice_stop")
    df["qc_pass"] = with_qc & n_applied_failed.eq(0)
    df["n_applied_rules_failed"] = n_applied_failed
    return df


def tabulate(df, by):
    rule_rows, cortex_rows = [], []
    levels = [("all", "all", df)]
    for col in by:
        for level, grp in df.groupby(col, sort=True):
            levels.append((col, level, grp))
    for stratifier, level, g in levels:
        base = {"stratifier": stratifier, "level": level, "n_scans": len(g)}
        for stage, n in g["stage"].value_counts().sort_index().items():
            rule_rows.append({**base, "rule": f"stage: {stage}", "applied": "stage", "n_evaluated": len(g),
                              "n_failed": int(n), "pct_failed": round(100.0 * n / len(g), 1),
                              "n_failed_only_this_rule": np.nan})
        n_ev = int(g["qc_evaluated"].sum())
        n_fail = int((g["qc_evaluated"] & ~g["qc_pass"]).sum())
        rule_rows.append({**base, "rule": "any applied rule", "applied": "summary", "n_evaluated": n_ev,
                          "n_failed": n_fail, "pct_failed": round(100.0 * n_fail / n_ev, 1) if n_ev else np.nan,
                          "n_failed_only_this_rule": np.nan})
        for rule, (_, applied) in RULES.items():
            ev, fail = g[f"ev_{rule}"], g[f"fail_{rule}"]
            if applied:
                alone = fail & g["n_applied_rules_failed"].eq(1)
            else:
                alone = fail & g["qc_pass"]
            n_ev_r = int(ev.sum())
            rule_rows.append({**base, "rule": rule, "applied": "yes" if applied else "report-only",
                              "n_evaluated": n_ev_r, "n_failed": int(fail.sum()),
                              "pct_failed": round(100.0 * fail.sum() / n_ev_r, 1) if n_ev_r else np.nan,
                              "n_failed_only_this_rule": int(alone.sum())})
        grid, t1 = g["fail_cortex_vol_qsm_grid"], g["fail_cortex_vol_t1"]
        grid_alone = grid & g["n_applied_rules_failed"].eq(1)
        cortex_rows.append({**base,
                            "n_both_evaluated": int((g["ev_cortex_vol_qsm_grid"] & g["ev_cortex_vol_t1"]).sum()),
                            "fail_qsm_grid": int(grid.sum()), "fail_t1": int(t1.sum()),
                            "excluded_only_by_qsm_grid_rule": int(grid_alone.sum()),
                            "of_which_pass_t1_rule": int((grid_alone & ~t1).sum()),
                            "pass_qc_but_fail_t1_rule": int((g["qc_pass"] & t1).sum())})
    return pd.DataFrame(rule_rows), pd.DataFrame(cortex_rows)


def main():
    ap = argparse.ArgumentParser(description="Count P8 QC exclusions by rule and stratum (report only)")
    ap.add_argument("--worklist", type=Path, default=HERE / "worklist.csv")
    ap.add_argument("--native-dir", type=Path, default=DEFAULT_NATIVE_DIR)
    ap.add_argument("--out-dir", type=Path, default=HERE / "output")
    ap.add_argument("--covariates", type=Path, default=None,
                    help="CSV with IID/image_dir_id or hospital_id and extra columns (e.g. sex)")
    ap.add_argument("--by", nargs="+", default=None,
                    help="Stratifiers (default: protocol coverage model_label, plus sex when available)")
    ap.add_argument("--ids-file", type=Path, default=None, help="Only these scans (IID or image_dir_id column)")
    ap.add_argument("--main-model-only", action="store_true", help="Only worklist rows with include_main_model=yes")
    ap.add_argument("--short-coverage-mm", type=float, default=SHORT_COVERAGE_MM)
    args = ap.parse_args()

    wl = pd.read_csv(args.worklist, dtype=str).fillna("")
    wl, conflicts = resolve_duplicate_ids(wl)
    if conflicts:
        print(f"Skipped {len(conflicts)} image_dir_id values with conflicting worklist rows: "
              f"{[e['IID'] for e in conflicts]}")
    if "runnable" in wl.columns:
        wl = wl[wl["runnable"].str.strip().str.lower() == "yes"]
    if args.main_model_only:
        wl = wl[wl["include_main_model"].str.strip().str.lower() == "yes"]
    if args.ids_file:
        ids = pd.read_csv(args.ids_file, dtype=str).fillna("")
        col = next((c for c in ("IID", "image_dir_id") if c in ids.columns), None)
        if col is None:
            raise SystemExit(f"{args.ids_file} has neither an IID nor an image_dir_id column")
        wl = wl[wl["image_dir_id"].isin(ids[col])]
    if wl.empty:
        raise SystemExit("No scans selected")

    df = pd.DataFrame([subject_record(oid, args.native_dir) for oid in wl["image_dir_id"]])
    df = df.merge(wl.rename(columns={"image_dir_id": "IID"}), on="IID", how="left")
    strata = [qsm_strata(q, p, c) for q, p, c in
              zip(df["qsm_file"], df.get("p8_qsm_pixdim", pd.Series("", index=df.index)).fillna(""),
                  df.get("qsm_coverage_mm", pd.Series("", index=df.index)))]
    df["protocol"] = [p for p, _ in strata]
    df["coverage_mm"] = [c for _, c in strata]
    df["coverage"] = [coverage_group(c, args.short_coverage_mm) for c in df["coverage_mm"]]
    if "model_label" in df.columns:
        df["model_label"] = df["model_label"].replace("", "(no model_label)")

    if args.covariates:
        cov = pd.read_csv(args.covariates, dtype=str).fillna("")
        key = next((c for c in ("IID", "image_dir_id", "hospital_id") if c in cov.columns), None)
        if key is None:
            raise SystemExit(f"{args.covariates} needs an IID, image_dir_id or hospital_id column")
        if key == "image_dir_id":
            cov, key = cov.rename(columns={"image_dir_id": "IID"}), "IID"
        if cov[key].duplicated().any():
            raise SystemExit(f"{args.covariates}: {key} values are not unique")
        extra = [c for c in cov.columns if c != key and c not in df.columns]
        df = df.merge(cov[[key] + extra], on=key, how="left")
        for c in extra:
            df[c] = df[c].fillna("").replace("", "(missing)")

    by = args.by or (["protocol", "coverage", "model_label"] + (["sex"] if "sex" in df.columns else []))
    missing = [c for c in by if c not in df.columns]
    if missing:
        raise SystemExit(f"Unknown stratifier(s) {missing}; available: {sorted(df.columns)}")

    df = evaluate(df)
    rules, cortex = tabulate(df, by)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    keep = (["IID", "hospital_id"] + by + ["coverage_mm", "stage", "qc_pass", "failed_rules",
            "report_only_failed_rules"] + METRICS)
    df[[c for c in dict.fromkeys(keep) if c in df.columns]].to_csv(
        args.out_dir / "qc_rules_by_subject.csv", index=False, encoding="utf-8-sig")
    rules.to_csv(args.out_dir / "qc_exclusions_by_rule.csv", index=False, encoding="utf-8-sig")
    cortex.to_csv(args.out_dir / "qc_cortex_rule_comparison.csv", index=False, encoding="utf-8-sig")

    print(f"{len(df)} scans; stages: {df['stage'].value_counts().to_dict()}")
    if df["stage"].isin(["p8_qc_without_rule_record", "p8_qc_incomplete", "no_p8_qc"]).any():
        print("  Scans without a complete, rule-recording P8 QC are not counted per rule; rerun P8 for them.")
    overall = rules[(rules["stratifier"] == "all") & rules["applied"].isin(["yes", "report-only"])]
    print("Rule (all scans): evaluated, failed (%), failed by this rule alone")
    for _, r in overall.iterrows():
        tag = "" if r["applied"] == "yes" else "  [report-only; 'alone' = would add exclusions]"
        print(f"  {r['rule']:<22} {r['n_evaluated']:5d} {r['n_failed']:5d} ({r['pct_failed']}%) "
              f"{int(r['n_failed_only_this_rule']):5d}{tag}")
    print(f"Tables in {args.out_dir}: qc_rules_by_subject.csv, qc_exclusions_by_rule.csv, "
          "qc_cortex_rule_comparison.csv")


if __name__ == "__main__":
    main()
