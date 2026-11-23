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


def read_dx(path, field, source):
    for row in rows(path, ("PATNO", "INFODT", "EVENT_ID", field)):
        pid, code = token(row["PATNO"]), token(row[field])
        if not pid or not code:
            continue
        when = date_range(row["INFODT"])
        if when is None:
            yield pid, None
        else:
            yield pid, Evidence(code, when, row["EVENT_ID"].strip().upper(), source,
                                token(row.get("DXLVL", "")), date_range(row.get("DIAGDT", "")))


def latest(records):
    if not records:
        return None, False
    boundary = max(r.when.lo for r in records)
    tied = [r for r in records if r.when.hi >= boundary]
    if len({r.code for r in tied}) > 1:
        return None, True
    return max(tied, key=lambda r: (r.when.hi, r.when.lo, r.source == "PRIMDXPD",
                                   r.confidence in {"3", "4"}, r.event)), False


def scan_evidence(records, scan):
    if scan is None:
        return None, "missing_scan_date"
    known = [r for r in records if r.when.hi <= scan.lo]
    observed, conflict = latest(known)
    if conflict:
        return None, "conflicting_latest_diagnoses"
    uncertain = [r for r in records if r.when.hi > scan.lo and r.when.overlaps(scan)]
    if observed is None or any(r.code != observed.code for r in uncertain):
        return observed, "uncertain_within_month_order" if uncertain else "no_pre_scan_diagnosis"
    if observed.reported_date and observed.reported_date.lo > scan.hi:
        return observed, "recorded_diagnosis_date_after_scan"
    if observed.code not in KNOWN_CODES:
        return observed, "unexpected_diagnosis_code"
    if observed.confidence in {"3", "4"}:
        return observed, "low_research_diagnosis_confidence"
    return observed, "matched"


def build(args):
    diagnoses, missing_date = defaultdict(list), Counter()
    for path, field, source in ((args.research_dx, "PRIMDIAG", "PRIMDXPD"),
                                (args.clinical_dx, "NEWDIAG", "NEWCLINDX")):
        if path is None:
            continue
        for pid, evidence in read_dx(path, field, source):
            if evidence is None:
                missing_date[pid] += 1
            else:
                diagnoses[pid].append(evidence)
    output, seen = [], set()
    for row in rows(args.scan_list, ("SUBJECT_ID", "EXAM_DATE", "EXAM_KEY")):
        pid, key = token(row["SUBJECT_ID"]), row["EXAM_KEY"].strip()
        if not pid or not key or (pid, key) in seen:
            raise ValueError("Scan list has a missing or repeated (subject_id, exam_key).")
        seen.add((pid, key))
        scan = date_range(row["EXAM_DATE"])
        observed, status = scan_evidence(diagnoses[pid], scan)
        code = observed.code if observed else ""
        prior_label = LABEL.get(code, "OTHER" if code in KNOWN_CODES else "")
        label = prior_label if status == "matched" else ""
        output.append(dict(PATNO=pid, exam_key=key, exam_date=row["EXAM_DATE"], label_at_scan=label,
                           last_pre_scan_label=prior_label,
                           diagnosis_code=code, diagnosis_source=observed.source if observed else "",
                           diagnosis_visit=observed.event if observed else "",
                           diagnosis_assessment_start=observed.when.lo if observed else "",
                           diagnosis_assessment_end=observed.when.hi if observed else "",
                           diagnosis_lag_days_min=(scan.lo - observed.when.hi).days if observed else "",
                           diagnosis_lag_days_max=(scan.hi - observed.when.lo).days if observed else "",
                           diagnosis_confidence_code=observed.confidence if observed else "",
                           diagnosis_status=status, undated_diagnosis_records=missing_date[pid]))
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
    args = parser.parse_args()
    if args.research_dx is None and args.clinical_dx is None:
        parser.error("Provide at least one visit-level diagnosis export: --research-dx or --clinical-dx.")
    try:
        matched = build(args)
        if not matched:
            raise ValueError("No selected scans in scan list.")
        counts = Counter((r["diagnosis_status"], r["label_at_scan"] or "UNASSIGNED") for r in matched)
        summary = [dict(diagnosis_status=d, label=l, scan_count=n) for (d, l), n in sorted(counts.items())]
        args.out_dir.mkdir(mode=0o700, parents=True, exist_ok=True)
        write_csv(args.out_dir / "PPMI_scan_diagnosis_PRIVATE.csv", matched)
        write_csv(args.out_dir / "PPMI_scan_diagnosis_summary.csv", summary)
    except (OSError, ValueError) as exc:
        parser.exit(2, f"Cannot build PPMI labels: {exc}\n")
    print(f"Wrote labels for {len(matched)} scan rows and {len(summary)} summary rows.")


if __name__ == "__main__":
    main()
