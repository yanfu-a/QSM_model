import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from build_nearest_diagnosis import Record, nearest_record, tier_of
from build_scan_diagnosis import LABEL, date_range

HERE = Path(__file__).resolve().parent


def dx(code, month, *, form="PRIMDXPD", confidence="1", rid="r"):
    file = "Clinical_Diagnosis" if form == "NEWCLINDX" else "Primary_Research_Diagnosis"
    return Record(code, LABEL.get(code, "OTHER"), date_range(month), file, form, tier_of(form),
                  rid, "BL", confidence, None)


class NearestDiagnosisTests(unittest.TestCase):
    scan = date_range("2020/6/15")

    def test_picks_closest_on_either_side(self):
        rec, status, tie, _, _ = nearest_record([dx("1", "01/2020"), dx("17", "08/2020")], self.scan)
        self.assertEqual((rec.code, status, tie), ("17", "matched", False))

    def test_same_month_wins(self):
        rec, status, _, _, _ = nearest_record([dx("1", "06/2020"), dx("17", "05/2020")], self.scan)
        self.assertEqual((rec.code, status), ("1", "matched"))

    def test_tie_same_code_prefers_pre_scan(self):
        rec, status, tie, _, _ = nearest_record([dx("1", "04/2020"), dx("1", "08/2020")], self.scan)
        self.assertEqual((rec.when.lo.month, status, tie), (4, "matched", True))

    def test_tie_different_codes_is_unassigned(self):
        rec, status, tie, _, _ = nearest_record([dx("1", "04/2020"), dx("11", "08/2020")], self.scan)
        self.assertEqual((rec, status, tie), (None, "tie_before_after_disagree", True))

    def test_precedence_resolves_cross_source_disagreement(self):
        records = [dx("1", "07/2020", form="NEWCLINDX"), dx("14", "07/2020", form="PRIMDXPD"),
                   dx("17", "07/2020", form="PRODDIAG")]
        rec, status, _, nearest, overruled = nearest_record(records, self.scan)
        self.assertEqual((rec.form, rec.code, status), ("NEWCLINDX", "1", "resolved_by_source_precedence"))
        self.assertEqual((len(nearest), sorted(r.code for r in overruled)), (3, ["14", "17"]))

    def test_disagreement_within_top_source_is_unassigned(self):
        records = [dx("1", "07/2020", rid="a"), dx("14", "07/2020", rid="b"), dx("17", "07/2020", form="DIAGQUES")]
        rec, status, _, _, _ = nearest_record(records, self.scan)
        self.assertEqual((rec, status), (None, "conflicting_records_in_top_source"))

    def test_prefers_higher_confidence(self):
        rec, _, _, _, _ = nearest_record([dx("1", "06/2020", confidence="3", rid="a"), dx("1", "06/2020", rid="b")],
                                         self.scan)
        self.assertEqual(rec.confidence, "1")

    def test_missing_inputs(self):
        self.assertEqual(nearest_record([], self.scan)[1], "no_dated_record")
        self.assertEqual(nearest_record([dx("1", "06/2020")], None)[1], "missing_scan_date")

    def test_cli_uses_cohort_only_as_fallback(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            files = {
                "scans.csv": [["subject_id", "exam_date", "exam_key"], ["1", "2020/6/15", "1_a"],
                              ["2", "2020/6/15", "2_a"]],
                "research.csv": [["REC_ID", "PATNO", "EVENT_ID", "PAG_NAME", "INFODT", "PRIMDIAG", "DXLVL"],
                                 ["r1", "1", "BL", "PRODDIAG", "05/2020", "23", "1"]],
                "clinical.csv": [["REC_ID", "PATNO", "EVENT_ID", "PAG_NAME", "INFODT", "NEWDIAG", "DIAGSINLV",
                                  "DIAGDT"]],
                "status.csv": [["PATNO", "COHORT", "COHORT_DEFINITION", "ENROLL_DATE"],
                               ["1", "4", "Prodromal", "05/2020"], ["2", "2", "Healthy Control", "06/2020"]],
            }
            for name, data in files.items():
                with open(tmp / name, "w", newline="") as f:
                    csv.writer(f).writerows(data)
            subprocess.run([sys.executable, str(HERE / "build_nearest_diagnosis.py"),
                            "--scan-list", str(tmp / "scans.csv"), "--research-dx", str(tmp / "research.csv"),
                            "--clinical-dx", str(tmp / "clinical.csv"),
                            "--participant-status", str(tmp / "status.csv"), "--out-dir", str(tmp / "out")],
                           check=True, capture_output=True)
            with open(tmp / "out" / "PPMI_scan_nearest_diagnosis_PRIVATE.csv") as f:
                out = {r["PATNO"]: r for r in csv.DictReader(f)}
        self.assertEqual((out["1"]["nearest_label"], out["1"]["source_form"], out["1"]["timing"]),
                         ("PRODROMAL", "PRODDIAG", "before_scan"))
        self.assertEqual(out["1"]["label_differs_from_enroll_cohort"], "0")
        self.assertEqual((out["2"]["nearest_status"], out["2"]["nearest_label"], out["2"]["source_file"]),
                         ("enrollment_cohort_only", "NO_PD_OR_OTHER_NEUROLOGICAL_DISORDER", "Participant_Status"))


if __name__ == "__main__":
    unittest.main()
