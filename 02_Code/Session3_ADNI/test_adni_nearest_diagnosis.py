import unittest

import pandas as pd

from build_adni_nearest_diagnosis import label_of, nearest


def rec(date, code, ddue="", visit="v"):
    return dict(date=pd.Timestamp(date), code=code, ddue=ddue, odes="", label=label_of(code, ddue),
                phase="ADNI2", visit=visit, record_id="1", qc_error="")


class AdniNearestTests(unittest.TestCase):
    scan = pd.Timestamp("2015-06-15")

    def test_closest_after(self):
        r, status, tie, _ = nearest([rec("2015-01-01", "1"), rec("2015-07-01", "2")], self.scan)
        self.assertEqual((r["label"], r["delta"], status, tie), ("MCI", 16, "matched", False))

    def test_tie_same_label_prefers_before(self):
        r, status, tie, _ = nearest([rec("2015-06-05", "1", visit="a"), rec("2015-06-25", "1", visit="b")], self.scan)
        self.assertEqual((r["delta"], status, tie), (-10, "matched", True))

    def test_tie_different_labels_unassigned(self):
        r, status, _, _ = nearest([rec("2015-06-05", "1"), rec("2015-06-25", "3", "1")], self.scan)
        self.assertEqual((r, status), (None, "tie_before_after_disagree"))

    def test_same_date_conflict(self):
        r, status, _, n = nearest([rec("2015-06-15", "2"), rec("2015-06-15", "3", "1")], self.scan)
        self.assertEqual((r, status, n), (None, "conflicting_records_same_date", 2))

    def test_dementia_etiology(self):
        self.assertEqual(label_of("3", "1"), "DEMENTIA_AD")
        self.assertEqual(label_of("3", "2"), "DEMENTIA_OTHER_ETIOLOGY")
        self.assertEqual(label_of("3", ""), "DEMENTIA_ETIOLOGY_NOT_RECORDED")
        self.assertEqual(label_of("10", ""), "")

    def test_unmapped_and_missing(self):
        self.assertEqual(nearest([rec("2015-06-15", "10")], self.scan)[1], "unmapped_diagnosis_code")
        self.assertEqual(nearest([], self.scan)[1], "no_dated_diagnosis")
        self.assertEqual(nearest([rec("2015-06-15", "1")], pd.NaT)[1], "missing_scan_date")


if __name__ == "__main__":
    unittest.main()
