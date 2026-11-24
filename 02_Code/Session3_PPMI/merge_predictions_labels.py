#!/usr/bin/env python3
"""Merge PPMI model predictions with the MRI-matched disease labels, one row per scan.

Inputs (all in --result-dir unless given explicitly):
  PD_1.csv                               model probabilities (subject_id, exam_key, 7 classes)
  PD_scan_diagnosis_PRIVATE.csv          time-of-scan labels (build_scan_diagnosis.py)
  PD_scan_nearest_diagnosis_PRIVATE.csv  nearest-diagnosis labels (build_nearest_diagnosis.py)
and Participant_Status (for age at enrollment).

Rows are joined on subject_id = PATNO and exam_key. Scans present in only one input are
kept and marked in merge_status. Image paths from PD_1.csv are not copied.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

CLASSES = {"REF": "REF", "阿尔茨海默病": "AD", "多系统萎缩": "MSA", "进行性核上性麻痹": "PSP",
           "帕金森病": "PD", "路易体痴呆": "DLB", "额颞叶痴呆": "FTD"}
NO_DISORDER = "NO_PD_OR_OTHER_NEUROLOGICAL_DISORDER"


def analysis_group(label, cohort):
    if label == "PD":
        return "PD"
    if label == NO_DISORDER:
        return {"Healthy Control": "HC", "Prodromal": "At-risk, no diagnosis",
                "SWEDD": "SWEDD, no PD"}.get(cohort, "No disorder, other cohort")
    if label == "PRODROMAL":
        return "Prodromal"
    if isinstance(label, str) and label:
        return "Other neurological"
    return "Unassigned"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--result-dir", type=Path, required=True)
    ap.add_argument("--participant-status", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    pred = pd.read_csv(a.result_dir / "PD_1.csv", dtype={"subject_id": str})
    if pred.duplicated(["subject_id", "exam_key"]).any():
        raise SystemExit("PD_1.csv has repeated (subject_id, exam_key).")
    pred = pred.rename(columns={k: f"prob_{v}" for k, v in CLASSES.items()})
    pred["pred_disease_en"] = pred["pred_disease"].map(CLASSES)
    pred = pred.drop(columns=["selected_mni2mm_nifti"], errors="ignore").rename(columns={"subject_id": "PATNO"})
    tos = pd.read_csv(a.result_dir / "PD_scan_diagnosis_PRIVATE.csv", dtype={"PATNO": str})
    tos = tos[["PATNO", "exam_key", "exam_date", "label_at_scan", "label_at_scan_month", "diagnosis_status",
               "diagnosis_source", "diagnosis_code", "diagnosis_visit", "diagnosis_confidence_code"]].rename(
        columns={c: f"scan_{c}" for c in ["diagnosis_status", "diagnosis_source", "diagnosis_code",
                                          "diagnosis_visit", "diagnosis_confidence_code"]})
    near = pd.read_csv(a.result_dir / "PD_scan_nearest_diagnosis_PRIVATE.csv", dtype={"PATNO": str})
    near = near[["PATNO", "exam_key", "nearest_label", "nearest_status", "source_file", "source_form",
                 "source_record_id", "diagnosis_code", "diagnosis_month", "timing", "interval_months",
                 "enroll_cohort", "enroll_date"]].rename(
        columns={"source_file": "nearest_source_file", "source_form": "nearest_source_form",
                 "source_record_id": "nearest_source_record_id", "diagnosis_code": "nearest_diagnosis_code",
                 "diagnosis_month": "nearest_diagnosis_month", "timing": "nearest_timing",
                 "interval_months": "nearest_interval_months"})
    labels = tos.merge(near, on=["PATNO", "exam_key"], how="outer", validate="one_to_one")
    m = pred.merge(labels, on=["PATNO", "exam_key"], how="outer", suffixes=("", "_label"),
                   validate="one_to_one", indicator=True)
    m["merge_status"] = m.pop("_merge").map({"both": "matched", "left_only": "prediction_only",
                                             "right_only": "label_only"})
    m["exam_date"] = m["exam_date"].fillna(m.pop("exam_date_label"))
    status = pd.read_csv(a.participant_status, dtype={"PATNO": str})[["PATNO", "ENROLL_AGE", "ENROLL_DATE"]]
    m = m.merge(status, on="PATNO", how="left")
    enrolled = pd.to_datetime(m.pop("ENROLL_DATE"), format="%m/%Y", errors="coerce")
    m["age_at_scan"] = (m.pop("ENROLL_AGE") + (pd.to_datetime(m["exam_date"], errors="coerce") - enrolled).dt.days
                        / 365.25).round(2)
    m["score_pairwise_PD_REF"] = m["prob_PD"] / (m["prob_PD"] + m["prob_REF"])
    m["score_log_PD_REF"] = np.log(m["prob_PD"] / m["prob_REF"])
    m["group_time_of_scan"] = [analysis_group(l, c) for l, c in zip(m["label_at_scan_month"], m["enroll_cohort"])]
    m["group_nearest"] = [analysis_group(l, c) for l, c in zip(m["nearest_label"], m["enroll_cohort"])]
    m["group_nearest_within_3m"] = np.where(m["nearest_interval_months"].abs() <= 3, m["group_nearest"],
                                            "Unassigned")
    first = ["PATNO", "exam_key", "exam_date", "merge_status", "age_at_scan", "enroll_cohort",
             "group_time_of_scan", "group_nearest", "group_nearest_within_3m"]
    m = m[first + [c for c in m.columns if c not in first]].sort_values(["PATNO", "exam_key"])
    a.out.parent.mkdir(parents=True, exist_ok=True)
    m.to_csv(a.out, index=False)
    print(f"Wrote {len(m)} rows: " + ", ".join(f"{k}={v}" for k, v in m["merge_status"].value_counts().items()))


if __name__ == "__main__":
    main()
