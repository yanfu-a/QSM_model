#!/usr/bin/env python3
"""Match each ADNI MRI scan to the diagnosis record nearest to its acquisition date.

Supplementary match, kept separate from any time-of-scan labelling: the nearest
record may come from a visit after the scan.

Scans:      ADNI_T1_registration_QC.xlsx, sheet of QC-passed scans (one per participant).
Diagnoses:  All_Subjects_DXSUM (ADNI diagnostic summary), one record per visit:
            DIAGNOSIS (1 = CN, 2 = MCI, 3 = Dementia) dated by EXAMDATE, with the
            dementia etiology fields DXDDUE (1 = due to AD, 2 = other) and DXODES kept.
Predictions (optional): adni_1.csv, only to flag which scans have a model prediction.

Dates in both sources are exact days, so the interval is in days:
interval_days = diagnosis date - scan date (negative = before the scan).
"""

import argparse
from collections import Counter
from pathlib import Path

import pandas as pd

DIAGNOSIS = {"1": "CN", "2": "MCI", "3": "DEMENTIA"}
NA = {"", "nan", "-4", "-1"}


def clean(value):
    value = "" if pd.isna(value) else str(value).strip()
    if value.endswith(".0"):
        value = value[:-2]
    return "" if value.lower() in NA else value


def label_of(code, ddue):
    if code == "3":
        return {"1": "DEMENTIA_AD", "2": "DEMENTIA_OTHER_ETIOLOGY"}.get(ddue, "DEMENTIA_ETIOLOGY_NOT_RECORDED")
    return DIAGNOSIS.get(code, "")


def nearest(records, scan_date):
    """records: list of dicts with keys date, code, ... Returns (record, status, tie, n_nearest)."""
    if pd.isna(scan_date):
        return None, "missing_scan_date", False, 0
    if not records:
        return None, "no_dated_diagnosis", False, 0
    for r in records:
        r["delta"] = (r["date"] - scan_date).days
    best = min(abs(r["delta"]) for r in records)
    near = [r for r in records if abs(r["delta"]) == best]
    tie = best > 0 and len({r["delta"] > 0 for r in near}) > 1
    if len({r["label"] for r in near}) > 1:
        return None, ("tie_before_after_disagree" if tie else "conflicting_records_same_date"), tie, len(near)
    chosen = max(near, key=lambda r: (r["delta"] <= 0, r["visit"]))  # prefer the pre-scan side of a tie
    status = "matched" if chosen["label"] else "unmapped_diagnosis_code"
    return chosen, status, tie, len(near)


def interval_bin(days):
    d = abs(days)
    return ("same_day" if d == 0 else "<=30d" if d <= 30 else "31-90d" if d <= 90 else "91-180d" if d <= 180
            else "181-365d" if d <= 365 else ">365d")


def build(qc_path, dxsum_path, pred_path=None):
    scans = pd.read_excel(qc_path, sheet_name="最早合格检查", dtype=str)
    scans = scans.rename(columns={"受试者ID": "PTID", "检查时间": "scan_time"})[["PTID", "scan_time"]]
    if scans.PTID.duplicated().any():
        raise SystemExit("QC scan list has more than one scan per participant.")
    predicted = set()
    if pred_path is not None:
        predicted = set(pd.read_csv(pred_path, dtype=str)["受试者ID"].str.strip())
    dx = pd.read_csv(dxsum_path, dtype=str, low_memory=False)
    dx["date"] = pd.to_datetime(dx.EXAMDATE, errors="coerce")
    undated = Counter(dx.loc[dx.date.isna(), "PTID"])
    records = {}
    for row in dx[dx.date.notna()].itertuples(index=False):
        code, ddue = clean(row.DIAGNOSIS), clean(row.DXDDUE)
        if not code:
            continue
        records.setdefault(row.PTID.strip(), []).append(dict(
            date=row.date, code=code, ddue=ddue, odes=clean(row.DXODES), label=label_of(code, ddue),
            phase=clean(row.PHASE), visit=clean(row.VISCODE2) or clean(row.VISCODE), record_id=clean(row.ID),
            qc_error=clean(row.HAS_QC_ERROR)))
    out = []
    for s in scans.itertuples(index=False):
        scan_dt = pd.to_datetime(s.scan_time, errors="coerce")
        scan_date = scan_dt.normalize() if pd.notna(scan_dt) else scan_dt
        recs = records.get(s.PTID.strip(), [])
        rec, status, tie, n_near = nearest(recs, scan_date)
        delta = rec["delta"] if rec else None
        out.append(dict(
            PTID=s.PTID, scan_time=s.scan_time, scan_date=scan_date.date() if pd.notna(scan_date) else "",
            has_model_prediction=int(s.PTID in predicted) if pred_path is not None else "",
            nearest_label=rec["label"] if rec and status == "matched" else "",
            nearest_status=status,
            diagnosis_code=rec["code"] if rec else "",
            dementia_due_to=rec["ddue"] if rec else "",
            other_dementia_etiology_code=rec["odes"] if rec else "",
            diagnosis_date=rec["date"].date() if rec else "",
            source_file="All_Subjects_DXSUM" if rec else "",
            source_phase=rec["phase"] if rec else "",
            source_visit=rec["visit"] if rec else "",
            source_record_id=rec["record_id"] if rec else "",
            source_has_qc_error=rec["qc_error"] if rec else "",
            timing="" if rec is None else "before_scan" if delta < 0 else "same_day" if delta == 0 else "after_scan",
            interval_days=delta if rec else "",
            interval_years=round(delta / 365.25, 2) if rec else "",
            tie_before_after=int(tie),
            records_at_nearest_distance=n_near,
            dated_diagnosis_records=len(recs),
            undated_diagnosis_records=undated.get(s.PTID, 0)))
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--qc-scans", type=Path, required=True, help="ADNI_T1_registration_QC.xlsx")
    ap.add_argument("--dxsum", type=Path, required=True, help="All_Subjects_DXSUM export")
    ap.add_argument("--predictions", type=Path, help="adni_1.csv (optional, to flag predicted scans)")
    ap.add_argument("--out-dir", type=Path, required=True)
    a = ap.parse_args()
    df = build(a.qc_scans, a.dxsum, a.predictions)
    df["abs_interval"] = [interval_bin(d) if d != "" else "NA" for d in df.interval_days]
    summary = (df.assign(nearest_label=df.nearest_label.replace("", "UNASSIGNED"), timing=df.timing.replace("", "NA"))
                 .groupby(["nearest_status", "nearest_label", "timing", "abs_interval"]).size()
                 .rename("scan_count").reset_index())
    a.out_dir.mkdir(parents=True, exist_ok=True)
    df.drop(columns="abs_interval").to_csv(a.out_dir / "ADNI_scan_nearest_diagnosis_PRIVATE.csv", index=False)
    summary.to_csv(a.out_dir / "ADNI_scan_nearest_diagnosis_summary.csv", index=False)
    print(f"Wrote nearest-diagnosis match for {len(df)} scans and {len(summary)} summary rows.")


if __name__ == "__main__":
    main()
