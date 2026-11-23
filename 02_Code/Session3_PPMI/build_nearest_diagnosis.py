#!/usr/bin/env python3
"""Match each PPMI MRI scan to the diagnosis record nearest to its acquisition date.

This is a supplementary match, kept separate from build_scan_diagnosis.py. It may
use a diagnosis recorded after the scan, so it is not a time-of-scan label.

How the three PPMI files contribute (highest precedence first):
  1. Clinical_Diagnosis (NEWCLINDX form): NEWDIAG dated by INFODT.
  2. Primary_Research_Diagnosis, PRIMDXPD form: PRIMDIAG dated by INFODT.
  3. Primary_Research_Diagnosis, CLINDX, DIAGQUES and PRODDIAG forms: PRIMDIAG
     dated by INFODT (same tier; undated rows are counted, not used).
  4. Participant_Status: enrollment COHORT dated by ENROLL_DATE. Used only when a
     participant has no dated record from tiers 1-3, and always reported in the
     enroll_* columns for comparison.
The nearest record is chosen by distance alone. Precedence is used only when
records at the same nearest distance disagree: the highest tier present wins if it
is unanimous, and the overruled records are listed. Disagreement within the top
tier, or between a record before and a record after the scan at equal distance,
leaves the scan unassigned.

PPMI dates are month/year, while scan dates are exact days. interval_months is
the record month minus the scan month (negative = before the scan, 0 = same month,
positive = after); interval_days_min/max give the possible range in days.
"""

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from build_scan_diagnosis import KNOWN_CODES, LABEL, DateRange, date_range, rows, token, write_csv

RESEARCH_FORMS = ("PRIMDXPD", "CLINDX", "DIAGQUES", "PRODDIAG")
COHORT_LABEL = {"1": "PD", "2": "NO_PD_OR_OTHER_NEUROLOGICAL_DISORDER", "3": "SWEDD", "4": "PRODROMAL",
                "7": "PD", "8": "UNAFFECTED_CARRIER"}
CONFIDENCE_RANK = {"1": 4, "2": 3, "3": 1, "4": 0}


@dataclass(frozen=True)
class Record:
    code: str            # diagnosis code, or "COHORT:<n>" for enrollment records
    label: str
    when: DateRange
    file: str            # Clinical_Diagnosis / Primary_Research_Diagnosis / Participant_Status
    form: str            # NEWCLINDX, PRIMDXPD, CLINDX, DIAGQUES, PRODDIAG or PATIENT_STATUS
    tier: int            # 1 = highest precedence
    record_id: str
    visit: str
    confidence: str
    reported_date: Optional[DateRange]

    @property
    def source(self):
        return f"{self.file}:{self.form}"


def tier_of(form):
    return 1 if form == "NEWCLINDX" else 2 if form == "PRIMDXPD" else 4 if form == "PATIENT_STATUS" else 3


def diagnosis_records(path, field, file, forms, undated, excluded):
    for row in rows(path, ("PATNO", "REC_ID", "EVENT_ID", "PAG_NAME", "INFODT", field)):
        pid, code = token(row["PATNO"]), token(row[field])
        form = row["PAG_NAME"].strip().upper()
        if not pid or not code:
            continue
        if forms is not None and form not in forms:
            excluded[pid] += 1
            continue
        when = date_range(row["INFODT"])
        if when is None:
            undated[pid] += 1
            continue
        label = LABEL.get(code, "OTHER" if code in KNOWN_CODES else "")
        yield pid, Record(code, label, when, file, form, tier_of(form), row["REC_ID"].strip(),
                          row["EVENT_ID"].strip().upper(), token(row.get("DXLVL", "")),
                          date_range(row.get("DIAGDT", "")))


def cohort_records(path):
    for row in rows(path, ("PATNO", "COHORT", "ENROLL_DATE")):
        pid, cohort = token(row["PATNO"]), token(row["COHORT"])
        if not pid:
            continue
        when = date_range(row["ENROLL_DATE"])
        record = None
        if cohort and when:
            record = Record(f"COHORT:{cohort}", COHORT_LABEL.get(cohort, ""), when, "Participant_Status",
                            "PATIENT_STATUS", 4, f"PATNO {pid}", "ENROLLMENT", "", None)
        yield pid, cohort, row.get("COHORT_DEFINITION", "").strip(), row["ENROLL_DATE"].strip(), record


def month_index(day):
    return day.year * 12 + day.month - 1


def nearest_record(records, scan):
    """Return (record, status, tie_before_after, nearest_records, overruled_records)."""
    if scan is None:
        return None, "missing_scan_date", False, [], []
    if not records:
        return None, "no_dated_record", False, [], []
    scan_month = month_index(scan.lo)
    distance = {id(r): month_index(r.when.lo) - scan_month for r in records}
    best = min(abs(d) for d in distance.values())
    nearest = [r for r in records if abs(distance[id(r)]) == best]
    tie = best > 0 and len({distance[id(r)] > 0 for r in nearest}) > 1
    if len({r.code for r in nearest}) > 1:
        if tie:
            return None, "tie_before_after_disagree", tie, nearest, []
        top = min(r.tier for r in nearest)
        winners = [r for r in nearest if r.tier == top]
        if len({r.code for r in winners}) > 1:
            return None, "conflicting_records_in_top_source", tie, nearest, []
        overruled = [r for r in nearest if r.code != winners[0].code]
        status, pool = "resolved_by_source_precedence", winners
    else:
        overruled, status, pool = [], "matched", nearest
    # Same code: prefer the pre-scan side, higher precedence, higher confidence, stable order.
    chosen = max(pool, key=lambda r: (distance[id(r)] <= 0, -r.tier, CONFIDENCE_RANK.get(r.confidence, 2),
                                      r.record_id))
    if not chosen.label:
        status = "unexpected_diagnosis_code"
    return chosen, status, tie, nearest, overruled


def interval_bin(months):
    m = abs(months)
    return "0" if m == 0 else "1-3" if m <= 3 else "4-6" if m <= 6 else "7-12" if m <= 12 else ">12"


def describe(records):
    return "|".join(sorted({f"{r.source}={r.code}" for r in records}))


def build(args):
    diagnoses, undated, excluded = defaultdict(list), Counter(), Counter()
    for path, field, file in ((args.clinical_dx, "NEWDIAG", "Clinical_Diagnosis"),
                              (args.research_dx, "PRIMDIAG", "Primary_Research_Diagnosis")):
        if path is None:
            continue
        forms = None if file == "Clinical_Diagnosis" else args.research_forms
        for pid, record in diagnosis_records(path, field, file, forms, undated, excluded):
            diagnoses[pid].append(record)
    cohort = {}
    if args.participant_status is not None:
        for pid, code, definition, enrolled, record in cohort_records(args.participant_status):
            cohort[pid] = (code, definition, enrolled, record)
    output, seen = [], set()
    for row in rows(args.scan_list, ("SUBJECT_ID", "EXAM_DATE", "EXAM_KEY")):
        pid, key = token(row["SUBJECT_ID"]), row["EXAM_KEY"].strip()
        if not pid or not key or (pid, key) in seen:
            raise ValueError("Scan list has a missing or repeated (subject_id, exam_key).")
        seen.add((pid, key))
        scan = date_range(row["EXAM_DATE"])
        enroll_code, enroll_def, enroll_date, enroll_record = cohort.get(pid, ("", "", "", None))
        records = diagnoses[pid]
        used_fallback = not records and enroll_record is not None
        rec, status, tie, nearest, overruled = nearest_record(
            records or ([enroll_record] if enroll_record else []), scan)
        if used_fallback and status in {"matched", "unexpected_diagnosis_code"}:
            status = "enrollment_cohort_only" if rec.label else status
        months = month_index(rec.when.lo) - month_index(scan.lo) if rec else ""
        assigned = status in {"matched", "resolved_by_source_precedence", "enrollment_cohort_only"}
        label = rec.label if rec and assigned else ""
        enroll_label = COHORT_LABEL.get(enroll_code, "")
        output.append(dict(
            PATNO=pid, exam_key=key, exam_date=row["EXAM_DATE"],
            nearest_label=label,
            nearest_status=status,
            source_file=rec.file if rec else "",
            source_form=rec.form if rec else "",
            source_precedence_tier=rec.tier if rec else "",
            source_record_id=rec.record_id if rec else "",
            source_visit=rec.visit if rec else "",
            diagnosis_code=rec.code if rec else "",
            diagnosis_month=f"{rec.when.lo:%Y-%m}" if rec else "",
            reported_diagnosis_date=(f"{rec.reported_date.lo:%Y-%m}" if rec and rec.reported_date else ""),
            diagnosis_confidence_code=rec.confidence if rec else "",
            timing=("" if rec is None else "before_scan" if months < 0
                    else "same_month" if months == 0 else "after_scan"),
            interval_months=months,
            interval_days_min=(rec.when.lo - scan.lo).days if rec else "",
            interval_days_max=(rec.when.hi - scan.lo).days if rec else "",
            tie_before_after=int(tie),
            records_at_nearest_distance=len(nearest),
            sources_at_nearest_distance=describe(nearest),
            overruled_records=describe(overruled),
            enroll_cohort_code=enroll_code,
            enroll_cohort=enroll_def,
            enroll_date=enroll_date,
            label_differs_from_enroll_cohort=(int(label != enroll_label) if label and enroll_label else ""),
            dated_diagnosis_records=len(records),
            undated_diagnosis_records=undated[pid],
            excluded_form_records=excluded[pid]))
    return output


def summarize(records):
    counts = Counter((r["nearest_status"], r["nearest_label"] or "UNASSIGNED", r["source_form"] or "NA",
                      r["timing"] or "NA",
                      interval_bin(r["interval_months"]) if r["interval_months"] != "" else "NA")
                     for r in records)
    return [dict(nearest_status=s, nearest_label=l, source_form=f, timing=t, abs_interval_months=b, scan_count=n)
            for (s, l, f, t, b), n in sorted(counts.items())]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for opt in ("scan-list", "out-dir"):
        parser.add_argument(f"--{opt}", type=Path, required=True)
    parser.add_argument("--research-dx", type=Path, help="Primary_Research_Diagnosis export.")
    parser.add_argument("--clinical-dx", type=Path, help="Clinical_Diagnosis export (NEWCLINDX).")
    parser.add_argument("--participant-status", type=Path, help="Participant_Status export (enrollment cohort).")
    parser.add_argument("--research-page", action="append", metavar="PAG_NAME",
                        help=f"Restrict Primary_Research_Diagnosis forms; repeat (default: {', '.join(RESEARCH_FORMS)}).")
    args = parser.parse_args()
    if args.research_dx is None and args.clinical_dx is None:
        parser.error("Provide at least one diagnosis export: --research-dx or --clinical-dx.")
    args.research_forms = {page.strip().upper() for page in (args.research_page or RESEARCH_FORMS)}
    try:
        matched = build(args)
        if not matched:
            raise ValueError("No selected scans in scan list.")
        summary = summarize(matched)
        args.out_dir.mkdir(mode=0o700, parents=True, exist_ok=True)
        write_csv(args.out_dir / "PPMI_scan_nearest_diagnosis_PRIVATE.csv", matched)
        write_csv(args.out_dir / "PPMI_scan_nearest_diagnosis_summary.csv", summary)
    except (OSError, ValueError) as exc:
        parser.exit(2, f"Cannot build nearest-diagnosis match: {exc}\n")
    print(f"Wrote nearest-diagnosis match for {len(matched)} scan rows and {len(summary)} summary rows.")


if __name__ == "__main__":
    main()
