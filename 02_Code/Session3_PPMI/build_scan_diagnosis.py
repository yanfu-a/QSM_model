#!/usr/bin/env python3
"""Build PPMI scan-time clinical labels offline from approved local CSV exports."""

import argparse
import calendar
import csv
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import date, datetime
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Optional


LABEL = {"1": "PD", "2": "AD", "3": "FTD", "4": "CBS", "5": "DLB",
         "11": "MSA", "14": "PSP", "17": "NO_PD_OR_OTHER_NEUROLOGICAL_DISORDER", "23": "PRODROMAL",
         "24": "PRODROMAL", "25": "PRODROMAL"}
KNOWN_CODES = set(LABEL) | {"6", "7", "8", "9", "10", "12", "13", "15", "16", "18", "97"}
NA = {"", "na", "nan", "null", "-4", "-9"}
CONFIDENCE_RANK = {"1": 4, "2": 3, "3": 1, "4": 0}


@dataclass(frozen=True)
class DateRange:
    lo: date
    hi: date

    def overlaps(self, other):
        return self.lo <= other.hi and other.lo <= self.hi


def date_range(value):
    value = value.strip()
    if value.lower() in NA:
        return None
    if re.fullmatch(r"\d{4}-\d{1,2}", value):
        year, month = map(int, value.split("-"))
    elif re.fullmatch(r"\d{1,2}/\d{4}", value):
        month, year = map(int, value.split("/"))
    else:
        if re.fullmatch(r"\d{4}-\d\d-\d\d(?:[T _].+)?", value):
            day = date.fromisoformat(value[:10])
        elif re.fullmatch(r"\d{4}/\d{1,2}/\d{1,2}", value):
            day = datetime.strptime(value, "%Y/%m/%d").date()
        else:
            day = datetime.strptime(value, "%m/%d/%Y").date()
        return DateRange(day, day)
    return DateRange(date(year, month, 1), date(year, month, calendar.monthrange(year, month)[1]))


def token(value):
    value = value.strip()
    if value.lower() in NA:
        return ""
    return str(int(value.split(".")[0])) if re.fullmatch(r"\d+(?:\.0+)?", value) else value


def rows(path, required):
    with open(path, encoding="utf-8-sig", newline="") as f:
        source = csv.reader(f)
        header = [name.strip().upper() for name in next(source, [])]
        if len(header) != len(set(header)) or not set(required) <= set(header):
            raise ValueError(f"{Path(path).name}: invalid header; expected {', '.join(required)}")
        for row in source:
            if not row:
                continue
            if len(row) != len(header):
                raise ValueError(f"{Path(path).name}: malformed CSV near line {source.line_num}")
            yield dict(zip(header, row))


@dataclass(frozen=True)
class Evidence:
    code: str
    when: DateRange
    event: str
    source: str
    confidence: str
    reported_date: Optional[DateRange]
    between_visits: bool


def read_dx(path, field, source, research_pages):
    required = ("PATNO", "INFODT", "EVENT_ID", field)
    if source == "PRIMDXPD":
        required += ("PAG_NAME",)
    for row in rows(path, required):
        pid, code = token(row["PATNO"]), token(row[field])
        if not pid:
            continue
        if source == "PRIMDXPD" and row["PAG_NAME"].strip().upper() not in research_pages:
            yield pid, None, "excluded_research_form"
            continue
        if not code:
            continue
        when = date_range(row["INFODT"])
        if when is None:
            yield pid, None, "undated_diagnosis"
        else:
            yield pid, Evidence(code, when, row["EVENT_ID"].strip().upper(), source,
                                token(row.get("DXLVL", "")), date_range(row.get("DIAGDT", "")),
                                token(row.get("DIAGSINLV", "")) == "1"), None


def latest(records, date_key=lambda r: r.when):
    if not records:
        return None, False
    boundary = max(date_key(r).lo for r in records)
    tied = [r for r in records if date_key(r).hi >= boundary]
    if len({r.code for r in tied}) > 1:
        return None, True
    return max(tied, key=lambda r: (CONFIDENCE_RANK.get(r.confidence, 2),
                                   r.source == "NEWCLINDX", date_key(r).hi, date_key(r).lo,
                                   r.when.hi, r.event)), False


def evidence_status(observed, scan):
    if observed.reported_date and observed.reported_date.lo > scan.hi:
        return "recorded_diagnosis_date_after_scan"
    if observed.code not in KNOWN_CODES:
        return "unexpected_diagnosis_code"
    if observed.confidence in {"3", "4"}:
        return "low_research_diagnosis_confidence"
    return "matched"


def scan_evidence(records, scan):
    if scan is None:
        return None, None, "missing_scan_date"
    known = [r for r in records if r.when.hi <= scan.lo]
    observed, conflict = latest(known)
    if conflict:
        return None, None, "conflicting_latest_diagnoses"
    uncertain = [r for r in records if r.when.hi > scan.lo and r.when.overlaps(scan)]
    # A change within an unresolved date interval cannot be placed before/after the MRI.
    relevant_codes = {r.code for r in uncertain} | ({observed.code} if observed else set())
    if len(relevant_codes) > 1:
        return observed, None, "uncertain_within_month_order"
    if observed:
        status = evidence_status(observed, scan)
        return observed, observed if status == "matched" else None, status
    if uncertain:
        concurrent, _ = latest(uncertain)
        status = evidence_status(concurrent, scan)
        return None, concurrent, "concurrent" if status == "matched" else status
    return None, None, "no_pre_scan_diagnosis"


def retrospective_evidence(records, scan, prior):
    """Return only supplemental evidence recorded at a definitely later visit."""
    if scan is None:
        return None, "not_applicable"
    later = [r for r in records if r.source == "NEWCLINDX" and r.when.lo > scan.hi
             and r.code in KNOWN_CODES and (prior is None or r.code != prior.code)]
    effective = [r for r in later if r.reported_date and r.reported_date.hi <= scan.lo]
    dated, conflict = latest(effective, lambda r: r.reported_date)
    if conflict:
        return None, "conflicting_reported_changes"
    if dated:
        if any(r.code != dated.code and r.when.lo >= dated.reported_date.hi and r.when.hi <= scan.lo
               for r in records):
            return dated, "reported_change_conflicts_with_pre_scan_assessment"
        return dated, "dated_report_at_later_visit"
    if any(r.reported_date and r.reported_date.overlaps(scan) for r in later):
        return None, "reported_change_order_uncertain"
    next_records = [r for r in records if r.when.lo > scan.hi]
    if not prior or not next_records:
        return None, "not_applicable"
    first_possible_visit = min(r.when.hi for r in next_records)
    next_records = [r for r in next_records if r.when.lo <= first_possible_visit]
    if len({r.code for r in next_records}) != 1:
        return None, "not_applicable"
    if any(r.source == "NEWCLINDX" and r.between_visits and not r.reported_date
           and r.code in KNOWN_CODES and r.code != prior.code for r in next_records):
        return None, "possible_transition_between_assessments"
    return None, "not_applicable"


def diagnosis_label(record):
    if record is None:
        return ""
    return LABEL.get(record.code, "OTHER" if record.code in KNOWN_CODES else "")


def build(args):
    diagnoses, input_flags = defaultdict(list), defaultdict(Counter)
    for path, field, source in ((args.research_dx, "PRIMDIAG", "PRIMDXPD"),
                                (args.clinical_dx, "NEWDIAG", "NEWCLINDX")):
        if path is None:
            continue
        for pid, evidence, flag in read_dx(path, field, source, args.research_page):
            if flag:
                input_flags[pid][flag] += 1
            else:
                diagnoses[pid].append(evidence)
    output, seen = [], set()
    for row in rows(args.scan_list, ("SUBJECT_ID", "EXAM_DATE", "EXAM_KEY")):
        pid, key = token(row["SUBJECT_ID"]), row["EXAM_KEY"].strip()
        if not pid or not key or (pid, key) in seen:
            raise ValueError("Scan list has a missing or repeated (subject_id, exam_key).")
        seen.add((pid, key))
        scan = date_range(row["EXAM_DATE"])
        observed, month_evidence, status = scan_evidence(diagnoses[pid], scan)
        reported, retrospective_status = retrospective_evidence(diagnoses[pid], scan, observed)
        evidence = month_evidence or observed
        prior_label = diagnosis_label(observed)
        label = prior_label if status == "matched" else ""
        month_label = diagnosis_label(month_evidence) if status in {"matched", "concurrent"} else ""
        output.append(dict(PATNO=pid, exam_key=key, exam_date=row["EXAM_DATE"], label_at_scan=label,
                           label_at_scan_month=month_label,
                           last_pre_scan_label=prior_label,
                           diagnosis_code=evidence.code if evidence else "",
                           diagnosis_source=evidence.source if evidence else "",
                           diagnosis_visit=evidence.event if evidence else "",
                           diagnosis_assessment_start=evidence.when.lo if evidence else "",
                           diagnosis_assessment_end=evidence.when.hi if evidence else "",
                           diagnosis_lag_days_min=(scan.lo - evidence.when.hi).days if evidence else "",
                           diagnosis_lag_days_max=(scan.hi - evidence.when.lo).days if evidence else "",
                           diagnosis_confidence_code=evidence.confidence if evidence else "",
                           diagnosis_status=status,
                           has_undated_diagnosis=int(input_flags[pid]["undated_diagnosis"] > 0),
                           undated_diagnosis_records=input_flags[pid]["undated_diagnosis"],
                           excluded_research_form_records=input_flags[pid]["excluded_research_form"],
                           retrospective_label=(diagnosis_label(reported)
                                                if retrospective_status == "dated_report_at_later_visit" else ""),
                           retrospective_status=retrospective_status,
                           reported_diagnosis_start=reported.reported_date.lo if reported else "",
                           reported_diagnosis_end=reported.reported_date.hi if reported else "",
                           retrospective_assessment_start=reported.when.lo if reported else "",
                           retrospective_assessment_end=reported.when.hi if reported else ""))
    return output


def write_csv(path, records):
    with NamedTemporaryFile("w", dir=path.parent, encoding="utf-8", newline="", delete=False) as f:
        tmp = Path(f.name)
        try:
            writer = csv.DictWriter(f, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
            f.close()
            tmp.replace(path)
        finally:
            tmp.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for opt in ("scan-list", "out-dir"):
        parser.add_argument(f"--{opt}", type=Path, required=True)
    parser.add_argument("--research-dx", type=Path, help="Primary Research Diagnosis (PRIMDXPD).")
    parser.add_argument("--clinical-dx", type=Path, help="Clinical Diagnosis (NEWCLINDX).")
    parser.add_argument("--research-page", action="append", metavar="PAG_NAME",
                        help="Accepted research form page; repeat if needed (default: PRIMDXPD only).")
    args = parser.parse_args()
    if args.research_dx is None and args.clinical_dx is None:
        parser.error("Provide at least one visit-level diagnosis export: --research-dx or --clinical-dx.")
    args.research_page = {page.strip().upper() for page in (args.research_page or ["PRIMDXPD"])}
    try:
        matched = build(args)
        if not matched:
            raise ValueError("No selected scans in scan list.")
        counts = Counter((r["diagnosis_status"], r["label_at_scan"] or "UNASSIGNED",
                          r["label_at_scan_month"] or "UNASSIGNED", r["retrospective_status"])
                         for r in matched)
        summary = [dict(diagnosis_status=d, label=l, label_at_scan_month=m, retrospective_status=s, scan_count=n)
                   for (d, l, m, s), n in sorted(counts.items())]
        args.out_dir.mkdir(mode=0o700, parents=True, exist_ok=True)
        write_csv(args.out_dir / "PPMI_scan_diagnosis_PRIVATE.csv", matched)
        write_csv(args.out_dir / "PPMI_scan_diagnosis_summary.csv", summary)
    except (OSError, ValueError) as exc:
        parser.exit(2, f"Cannot build PPMI labels: {exc}\n")
    print(f"Wrote labels for {len(matched)} scan rows and {len(summary)} summary rows.")


if __name__ == "__main__":
    main()
