"""qc_exclusion_report.py and compare_p9_runs.py on synthetic inputs."""

import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

from conftest import REVISED, run_p8, write_worklist


def _qc(native, oid, rules="", report_only="", pixdim="1.000000x1.000000x1.000000", **metrics):
    base = {"reg_mask_dice": "0.85", "synthseg_icv_ml": "1500", "cortex_vol_ml": "450", "ribbon_in_qsm_frac": "0.9",
            "gm_wm_contrast_ppb": "20", "seg_qsm_dice": "0.8", "synthseg_cortex_ml": "480", **metrics}
    lines = [f"ID={oid}", f"qsm_pixdim={pixdim}"] + [f"{k}={v}" for k, v in base.items()]
    if rules is not None:
        lines += [f"qc_failed_rules={rules}", f"qc_report_only_failed_rules={report_only}"]
    lines += [f"ERROR: {r} failed" for r in (rules or "").split(",") if r] + ["qc_complete=1"]
    (native / "qc").mkdir(parents=True, exist_ok=True)
    (native / "qc" / f"{oid}_native_qc.txt").write_text("\n".join(lines) + "\n")


def _report(tmp_path, wl, native, *extra):
    out = tmp_path / "report"
    proc = subprocess.run([sys.executable, str(REVISED / "qc_exclusion_report.py"), "--worklist", str(wl),
                           "--native-dir", str(native), "--out-dir", str(out), *map(str, extra)],
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    return out, proc.stdout


def test_exclusions_are_counted_per_rule_and_stratum(tmp_path):
    native = tmp_path / "native"
    _qc(native, "PASS")
    _qc(native, "GRIDONLY", "cortex_vol_qsm_grid", cortex_vol_ml="300")       # rescued by the T1 rule
    _qc(native, "TWO", "reg_mask_dice,seg_qsm_dice", reg_mask_dice="0.5", seg_qsm_dice="0.5")
    _qc(native, "T1ONLY", "", "cortex_vol_t1", synthseg_cortex_ml="320", pixdim="0.8x0.8x0.8")
    _qc(native, "OLDQC", None)                                                # no rule record
    _qc(native, "GRIDNOT1", "cortex_vol_qsm_grid", cortex_vol_ml="300", synthseg_cortex_ml="nan")
    _qc(native, "NOVOL", synthseg_icv_ml="nan", synthseg_cortex_ml="nan")     # earlier P8: ICV rule skipped
    (native / "p1_rigid").mkdir()
    (native / "p1_rigid" / "P1STOP_qc_rigid.txt").write_text("reg_mask_dice=0.41\nWARNING: reg_mask_dice 0.41 < 0.60\n")
    rows = [{"image_dir_id": oid, "hospital_id": f"H{i}", "qsm_file": str(tmp_path / "absent.nii"), "t1_file": "",
             "qsm_coverage_mm": cov, "model_label": dx}
            for i, (oid, cov, dx) in enumerate([("PASS", "110", "AD"), ("GRIDONLY", "110", "AD"),
                                                ("TWO", "144", "PSP"), ("T1ONLY", "144", "PSP"),
                                                ("OLDQC", "144", "AD"), ("P1STOP", "110", "AD"),
                                                ("NONE", "110", ""), ("GRIDNOT1", "144", "AD"),
                                                ("NOVOL", "144", "AD")])]
    wl = write_worklist(tmp_path / "wl.csv", rows)
    pd.DataFrame({"hospital_id": ["H0", "H1", "H2", "H3", "H4", "H5"], "sex": list("FMFMFM")}).to_csv(
        tmp_path / "cov.csv", index=False)
    out, log = _report(tmp_path, wl, native, "--covariates", tmp_path / "cov.csv")

    subj = pd.read_csv(out / "qc_rules_by_subject.csv", keep_default_na=False).set_index("IID")
    assert subj.loc["PASS", "stage"] == "p8_qc" and str(subj.loc["PASS", "qc_pass"]) == "True"
    assert subj.loc["OLDQC", "stage"] == "p8_qc_without_rule_record"
    assert subj.loc["P1STOP", "stage"] == "p1_dice_stop" and subj.loc["NONE", "stage"] == "no_p8_qc"
    assert subj.loc["T1ONLY", "protocol"] == "0.8x0.8x0.8" and subj.loc["PASS", "coverage"] == "short"
    assert subj.loc["NONE", "sex"] == "(missing)" and subj.loc["NONE", "model_label"] == "(no model_label)"

    rules = pd.read_csv(out / "qc_exclusions_by_rule.csv")
    get = lambda strat, level, rule: rules[(rules.stratifier == strat) & (rules.level == level)
                                           & (rules.rule == rule)].iloc[0]
    # "Alone" needs every other applied rule measured and passed: TWO fails two rules, and P1STOP
    # was evaluated on reg_mask_dice only.
    r = get("all", "all", "reg_mask_dice")
    assert (r.n_evaluated, r.n_not_evaluated, r.n_failed, r.n_failed_only_this_rule) == (7, 0, 2, 0)
    r = get("all", "all", "synthseg_icv")
    assert (r.n_evaluated, r.n_not_evaluated, r.n_failed) == (5, 1, 0)             # NOVOL not measured
    r = get("all", "all", "cortex_vol_qsm_grid")
    assert (r.n_evaluated, r.n_failed, r.n_failed_only_this_rule) == (6, 2, 2)
    r = get("all", "all", "cortex_vol_t1")
    assert r.applied == "report-only" and (r.n_evaluated, r.n_not_evaluated) == (4, 2)
    assert (r.n_failed, r.n_failed_only_this_rule) == (1, 1)
    r = get("all", "all", "any applied rule")
    assert (r.n_evaluated, r.n_failed) == (7, 4)
    assert get("coverage", "short", "reg_mask_dice").n_failed == 1
    assert get("model_label", "PSP", "seg_qsm_dice").n_failed == 1
    assert get("sex", "F", "any applied rule").n_failed == 1                        # TWO (H2)
    assert get("all", "all", "stage: p8_qc_without_rule_record").n_failed == 1

    cortex = pd.read_csv(out / "qc_cortex_rule_comparison.csv").set_index(["stratifier", "level"])
    c = cortex.loc[("all", "all")]
    # GRIDNOT1 has no T1 cortex volume: counted as not measured, never as passing the T1 rule.
    assert (c.excluded_only_by_qsm_grid_rule, c.of_which_pass_t1_rule, c.of_which_t1_not_measured) == (2, 1, 1)
    assert (c.t1_not_measured, c.pass_qc_but_fail_t1_rule) == (2, 1)
    assert "rerun P8" in log


def test_report_keeps_identical_rows_and_excludes_conflicts_before_main_model_filter(tmp_path):
    native = tmp_path / "native"
    for oid in ("SAME", "DIFFERENT", "UNIQUE"):
        _qc(native, oid)
    base = lambda oid: {"image_dir_id": oid, "qsm_file": str(tmp_path / "absent.nii"),
                        "t1_file": "t1", "include_main_model": "yes"}
    rows = [base("SAME"), base("SAME"), base("DIFFERENT"),
            {**base("DIFFERENT"), "include_main_model": "no"}, base("UNIQUE")]
    wl = write_worklist(tmp_path / "wl.csv", rows)
    out, log = _report(tmp_path, wl, native, "--main-model-only")
    subj = pd.read_csv(out / "qc_rules_by_subject.csv", dtype=str)
    assert list(subj["IID"]) == ["SAME", "UNIQUE"]
    assert "DIFFERENT" in log and "Skipped 1" in log


@pytest.mark.slow
def test_report_reads_p8_qc_files(tools, subject, tmp_path):
    native = tmp_path / "native"
    assert run_p8(tools, subject, native, "OK")[0] == 0
    code, out, _ = run_p8(tools, subject, native, "LOWDICE", env={"FAKE_P1_DICE": "0.500"})
    assert code == 15, out
    rows = [{"image_dir_id": oid, "qsm_file": str(subject["qsm"]), "t1_file": str(subject["t1"])}
            for oid in ("OK", "LOWDICE")]
    out, _ = _report(tmp_path, write_worklist(tmp_path / "wl.csv", rows), native)
    subj = pd.read_csv(out / "qc_rules_by_subject.csv", keep_default_na=False).set_index("IID")
    assert subj.loc["LOWDICE", "failed_rules"] == "reg_mask_dice" and str(subj.loc["OK", "qc_pass"]) == "True"
    assert subj.loc["OK", "protocol"] == "1x1x1"


def _p9_table(path, values):
    df = pd.DataFrame({"IID": list(values), "model_label": "AD"})
    for roi in ("L_a", "R_a"):
        df[f"{roi}_med"] = [v[roi] for v in values.values()]
        df[f"{roi}_mean"] = df[f"{roi}_med"] + 1
    df.to_csv(path, index=False)
    return path


def test_compare_p9_runs(tmp_path):
    a = _p9_table(tmp_path / "a.csv", {"S1": {"L_a": 10, "R_a": 5}, "S2": {"L_a": 20, "R_a": 6},
                                       "S3": {"L_a": 30, "R_a": np.nan}, "S4": {"L_a": 1, "R_a": 1}})
    b = _p9_table(tmp_path / "b.csv", {"S1": {"L_a": 11, "R_a": 5}, "S2": {"L_a": 21, "R_a": 6},
                                       "S3": {"L_a": 33, "R_a": 7}, "S5": {"L_a": 9, "R_a": 9}})
    pd.DataFrame({"IID": ["S1", "S2", "S3", "S4"], "qsm_protocol": ["1x1x1", "1x1x1", "0.8x0.8x0.8", "1x1x1"],
                  "qsm_coverage_group": ["short", "short", "full", "short"]}).to_csv(tmp_path / "qc.csv", index=False)
    out = tmp_path / "cmp.csv"
    proc = subprocess.run([sys.executable, str(REVISED / "compare_p9_runs.py"), "--a", str(a), "--b", str(b),
                           "--qc", str(tmp_path / "qc.csv"), "--out", str(out)], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "only in A: 1; only in B: 1" in proc.stdout
    t = pd.read_csv(out).set_index(["stratum", "ROI"])
    la = t.loc[("all", "L_a")]
    assert la.n_both == 3 and np.isclose(la.mean_diff_b_minus_a, 5 / 3) and la.max_abs_diff == 3
    assert np.isclose(la.pearson_r, np.corrcoef([10, 20, 30], [11, 21, 33])[0, 1])
    ra = t.loc[("all", "R_a")]
    assert (ra.n_both, ra.n_only_b, ra.mean_abs_diff) == (2, 1, 0)
    assert t.loc[("1x1x1 | short", "L_a")].n_both == 2 and t.loc[("0.8x0.8x0.8 | full", "L_a")].n_both == 1
