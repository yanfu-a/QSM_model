#!/usr/bin/env python3
"""Match each PPMI MRI scan to the diagnosis record nearest to its acquisition date.

This is a supplementary match, kept separate from build_scan_diagnosis.py. It may
use a diagnosis recorded after the scan, so it is not a time-of-scan label.

PPMI diagnosis dates (INFODT) are month/year, while scan dates are exact days.
Distance is therefore measured in calendar months: interval_months is the
diagnosis month minus the scan month (negative = before the scan, 0 = same
month, positive = after). The possible range in days is also reported.
"""

import argparse
from collections import Counter, defaultdict
from pathlib import Path

from build_scan_diagnosis import KNOWN_CODES, date_range, diagnosis_label, read_dx, rows, token, write_csv

CONFIDENCE_RANK = {"1": 4, "2": 3, "3": 1, "4": 0}


def month_index(day):
    return day.year * 12 + day.month - 1


def nearest_record(records, scan):
    """Return (record, status, tie_before_after, n_at_nearest)."""
    if scan is None:
        return None, "missing_scan_date", False, 0
    if not records:
        return None, "no_dated_diagnosis", False, 0
    scan_month = month_index(scan.lo)
    distance = {id(r): month_index(r.when.lo) - scan_month for r in records}
    best = min(abs(d) for d in distance.values())
    nearest = [r for r in records if abs(distance[id(r)]) == best]
    sides = {distance[id(r)] > 0 for r in nearest} if best else {False}
    tie = len(sides) > 1
    if len({r.code for r in nearest}) > 1:
        status = "tie_before_after_disagree" if tie else "conflicting_codes_at_nearest"
        return None, status, tie, len(nearest)
    # Same code everywhere: prefer the pre-scan side of a tie, then higher confidence,
    # then the Clinical Diagnosis form, then the visit ID for a stable choice.
    chosen = max(nearest, key=lambda r: (distance[id(r)] <= 0, CONFIDENCE_RANK.get(r.confidence, 2),
                                         r.source == "NEWCLINDX", r.event))
    if chosen.code not in KNOWN_CODES:
        return chosen, "unexpected_diagnosis_code", tie, len(nearest)
    return chosen, "matched", tie, len(nearest)


def interval_bin(months):
    m = abs(months)
    return "0" if m == 0 else "1-3" if m <= 3 else "4-6" if m <= 6 else "7-12" if m <= 12 else ">12"


def build(args):
    diagnoses, undated = defaultdict(list), Counter()
    for path, field, source in ((args.research_dx, "PRIMDIAG", "PRIMDXPD"),
                                (args.clinical_dx, "NEWDIAG", "NEWCLINDX")):
        if path is None:
            continue
        for pid, evidence, flag in read_dx(path, field, source, args.research_page):
            if flag == "undated_diagnosis":
                undated[pid] += 1
            elif evidence is not None:
                diagnoses[pid].append(evidence)
    output, seen = [], set()
    for row in rows(args.scan_list, ("SUBJECT_ID", "EXAM_DATE", "EXAM_KEY")):
        pid, key = token(row["SUBJECT_ID"]), row["EXAM_KEY"].strip()
        if not pid or not key or (pid, key) in seen:
            raise ValueError("Scan list has a missing or repeated (subject_id, exam_key).")
        seen.add((pid, key))
        scan = date_range(row["EXAM_DATE"])
        rec, status, tie, n_nearest = nearest_record(diagnoses[pid], scan)
        months = month_index(rec.when.lo) - month_index(scan.lo) if rec else ""
        output.append(dict(
            PATNO=pid, exam_key=key, exam_date=row["EXAM_DATE"],
            nearest_label=diagnosis_label(rec) if status == "matched" else "",
            nearest_status=status,
            diagnosis_code=rec.code if rec else "",
            diagnosis_source=rec.source if rec else "",
            diagnosis_visit=rec.event if rec else "",
            diagnosis_month=f"{rec.when.lo:%Y-%m}" if rec else "",
            reported_diagnosis_date=(f"{rec.reported_date.lo:%Y-%m}" if rec and rec.reported_date else ""),
            diagnosis_confidence_code=rec.confidence if rec else "",
            timing=("" if rec is None else "before_scan" if months < 0
                    else "same_month" if months == 0 else "after_scan"),
            interval_months=months,
            interval_days_min=(rec.when.lo - scan.lo).days if rec else "",
            interval_days_max=(rec.when.hi - scan.lo).days if rec else "",
            tie_before_after=int(tie),
            records_at_nearest_distance=n_nearest,
            dated_diagnosis_records=len(diagnoses[pid]),
            undated_diagnosis_records=undated[pid]))
    return output


def summarize(records):
    counts = Counter((r["nearest_status"], r["nearest_label"] or "UNASSIGNED", r["timing"] or "NA",
                      interval_bin(r["interval_months"]) if r["interval_months"] != "" else "NA")
                     for r in records)
    return [dict(nearest_status=s, nearest_label=l, timing=t, abs_interval_months=b, scan_count=n)
            for (s, l, t, b), n in sorted(counts.items())]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for opt in ("scan-list", "out-dir"):
        parser.add_argument(f"--{opt}", type=Path, required=True)
    parser.add_argument("--research-dx", type=Path, help="Primary Research Diagnosis (PRIMDXPD).")
    parser.add_argument("--clinical-dx", type=Path, help="Clinical Diagnosis (NEWCLINDX).")
    parser.add_argument("--research-page", action="append", metavar="PAG_NAME",
                        help="Accepted research form page; repeat if needed (default: PRIMDXPD only).")
    args = parser.parse_args()
    if args.research_dx is None and args.clinical_dx is None:
        parser.error("Provide at least one diagnosis export: --research-dx or --clinical-dx.")
    args.research_page = {page.strip().upper() for page in (args.research_page or ["PRIMDXPD"])}
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
