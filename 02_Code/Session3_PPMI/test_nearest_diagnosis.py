import unittest

from build_nearest_diagnosis import nearest_record
from build_scan_diagnosis import Evidence, date_range


def dx(code, month, *, source="PRIMDXPD", confidence="1"):
    return Evidence(code, date_range(month), "BL", source, confidence, None, False)


class NearestDiagnosisTests(unittest.TestCase):
    scan = date_range("2020/6/15")

    def test_picks_closest_on_either_side(self):
        rec, status, tie, _ = nearest_record([dx("1", "01/2020"), dx("17", "08/2020")], self.scan)
        self.assertEqual((rec.code, status, tie), ("17", "matched", False))

    def test_same_month_wins(self):
        rec, status, _, _ = nearest_record([dx("1", "06/2020"), dx("17", "05/2020")], self.scan)
        self.assertEqual((rec.code, status), ("1", "matched"))

    def test_tie_same_code_prefers_pre_scan(self):
        rec, status, tie, _ = nearest_record([dx("1", "04/2020"), dx("1", "08/2020")], self.scan)
        self.assertEqual((rec.when.lo.month, status, tie), (4, "matched", True))

    def test_tie_different_codes_is_unassigned(self):
        rec, status, tie, _ = nearest_record([dx("1", "04/2020"), dx("11", "08/2020")], self.scan)
        self.assertEqual((rec, status, tie), (None, "tie_before_after_disagree", True))

    def test_conflict_in_nearest_month(self):
        rec, status, _, n = nearest_record([dx("1", "07/2020"), dx("14", "07/2020")], self.scan)
        self.assertEqual((rec, status, n), (None, "conflicting_codes_at_nearest", 2))

    def test_prefers_higher_confidence(self):
        rec, _, _, _ = nearest_record([dx("1", "06/2020", confidence="3"), dx("1", "06/2020")], self.scan)
        self.assertEqual(rec.confidence, "1")

    def test_missing_inputs(self):
        self.assertEqual(nearest_record([], self.scan)[1], "no_dated_diagnosis")
        self.assertEqual(nearest_record([dx("1", "06/2020")], None)[1], "missing_scan_date")


if __name__ == "__main__":
    unittest.main()
