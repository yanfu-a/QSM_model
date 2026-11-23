import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from build_scan_diagnosis import Evidence, date_range, retrospective_evidence, scan_evidence


def diagnosis(code, assessment, *, source="PRIMDXPD", confidence="1", effective="", change=False):
    return Evidence(code, date_range(assessment), "BL", source, confidence, date_range(effective), change)


class ScanDiagnosisTests(unittest.TestCase):
    def test_same_month_agreement_and_disagreement(self):
        scan = date_range("2024/6/15")
        current = diagnosis("1", "06/2024")
        earlier = diagnosis("1", "05/2024")
        other = diagnosis("11", "06/2024")
        self.assertEqual(scan_evidence([current], scan), (None, current, "concurrent"))
        self.assertEqual(scan_evidence([earlier, current], scan), (earlier, earlier, "matched"))
        self.assertEqual(scan_evidence([earlier, other], scan), (earlier, None, "uncertain_within_month_order"))
        self.assertEqual(scan_evidence([current, other], scan), (None, None, "uncertain_within_month_order"))

    def test_confidence_and_exact_future_date(self):
        scan = date_range("2024/6/15")
        high = diagnosis("1", "05/2024", confidence="1")
        low = diagnosis("1", "05/2024", confidence="3")
        future = diagnosis("11", "06/28/2024")
        self.assertEqual(scan_evidence([low, future, high], scan), (high, high, "matched"))
        self.assertEqual(scan_evidence([future], scan), (None, None, "no_pre_scan_diagnosis"))

    def test_retrospective_date_is_separate(self):
        scan = date_range("2024/6/15")
        earlier = diagnosis("1", "02/2024")
        later = diagnosis("11", "12/2024", source="NEWCLINDX", effective="03/2024", change=True)
        self.assertEqual(scan_evidence([earlier, later], scan), (earlier, earlier, "matched"))
        self.assertEqual(retrospective_evidence([earlier, later], scan, earlier), (later, "dated_report_at_later_visit"))
        intervening = diagnosis("1", "05/2024")
        self.assertEqual(retrospective_evidence([earlier, later, intervening], scan, intervening),
                         (later, "reported_change_conflicts_with_pre_scan_assessment"))
        undated = diagnosis("11", "12/2024", source="NEWCLINDX", change=True)
        self.assertEqual(retrospective_evidence([earlier, undated], scan, earlier),
                         (None, "possible_transition_between_assessments"))

    def test_cli_reports_concurrent_and_flags_research_forms(self):
        with tempfile.TemporaryDirectory() as root:
            base = Path(root)
            scan = base / "scan.csv"
            dx = base / "dx.csv"
            scan.write_text("subject_id,exam_date,exam_key\n1,2024/6/15,study\n")
            dx.write_text("PATNO,INFODT,EVENT_ID,PRIMDIAG,DXLVL,PAG_NAME\n"
                          "1,06/2024,BL,1,1,PRIMDXPD\n1,06/2024,BL,11,1,PRODDIAG\n"
                          "1,,BL,1,1,PRIMDXPD\n")
            command = [sys.executable, str(Path(__file__).with_name("build_scan_diagnosis.py")),
                       "--scan-list", str(scan), "--research-dx", str(dx), "--out-dir", str(base / "output")]
            subprocess.run(command, check=True, capture_output=True)
            with (base / "output/PPMI_scan_diagnosis_PRIVATE.csv").open() as f:
                row, = csv.DictReader(f)
            self.assertEqual(row["label_at_scan"], "")
            self.assertEqual(row["label_at_scan_month"], "PD")
            self.assertEqual(row["diagnosis_status"], "concurrent")
            self.assertEqual(row["has_undated_diagnosis"], "1")
            self.assertEqual(row["excluded_research_form_records"], "1")


if __name__ == "__main__":
    unittest.main()
