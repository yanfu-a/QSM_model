"""P9 on a synthetic cohort built with P8 (stub tools): eligibility, provenance, exclusions.

Each subject ID exercises one rule; see COHORT below. Outputs of a formal (support-mask) run
are checked subject by subject, and one subject's ROI values against an independent
calculation from the image files.
"""

import hashlib
import json
import shutil
import subprocess
import sys

import numpy as np
import nibabel as nib
import pandas as pd
import pytest

from conftest import REVISED, clone_outputs, run_p8, run_p9, write_worklist
from synth import make_subject

pytestmark = pytest.mark.slow

# ID -> expected outcome. True: in the analysis set; str: substring of fail_reasons or of the
# excluded.csv reason (prefixed "excluded:").
COHORT = {
    "OK1": True,
    "REP": True,                     # second scan of OK1's patient
    "FAILMNI": True,                 # P1 failure confined to FNIRT: stays in the native analysis
    "QSM4D": True,                   # singleton-4D QSM
    "NOCSF": True,                   # no CSF reference: primary set unchanged, CSF set excludes it
    "QREP": "QSM content differs from the file P8 used",
    "T1SAME": "T1 content differs from the file P8 segmented",
    "T1MIS": "T1 content differs from the file P8 segmented",
    "LEG_AFF": "P8 label predates grid-conversion fix",
    "LEG_UNAFF": "P8 label predates grid-conversion fix",
    "QCMISS": "P8 QC missing",
    "QCSTALE": "label differs from the one P8 QC checked",
    "QCINC": "P8 QC incomplete",
    "CORRUPT": "extraction: ERROR",
    "ZEROSUP": "extraction: No valid supported cortical voxels",
    "STALESUP": "support mask was built from a different QSM",
    "NOREC": "support mask has no provenance record",
    "DUP": "excluded:worklist rows disagree: t1_file=",
    "NONRUN": "excluded:worklist runnable='no'",
    "FAILRIG": "excluded:P1 failure (stage native)",
    "FAILUNK": "excluded:P1 failure (stage unknown)",
    "NOLABEL": "excluded:P8 native label missing",
    "NOQSM": "excluded:QSM file missing",
    "NOT1": "excluded:T1 file missing",
}


def _modify(path, delta=1.0):
    img = nib.load(path)
    data = np.asanyarray(img.dataobj).astype(np.float32)
    data[data != 0] += delta
    nib.save(nib.Nifti1Image(data, img.affine, img.header), path)


def _sha1(path):
    return hashlib.sha1(path.read_bytes()).hexdigest()


@pytest.fixture(scope="module")
def cohort(tools, subject, tmp_path_factory):
    root = tmp_path_factory.mktemp("cohort")
    native, files = root / "native", root / "files"
    files.mkdir()
    s = subject
    row = lambda oid, qsm=None, t1=None, **kw: {"image_dir_id": oid, "hospital_id": kw.pop("hosp", oid),
                                                 "qsm_file": str(qsm or s["qsm"]), "t1_file": str(t1 or s["t1"]), **kw}

    assert run_p8(tools, s, native, "OK1")[0] == 0
    for oid in ("REP", "QCMISS", "QCSTALE", "QCINC", "CORRUPT", "ZEROSUP", "STALESUP", "NOREC",
                "FAILMNI", "FAILRIG", "FAILUNK", "T1MIS", "DUP", "NONRUN", "NOQSM", "NOT1"):
        clone_outputs(native, "OK1", oid)
    (native / "qc" / "QCMISS_native_qc.txt").unlink()
    _modify(native / "label" / "QCSTALE_aparc_aseg_QSMnative.nii.gz", 0)  # rewritten: new bytes
    lab = native / "label" / "QCSTALE_aparc_aseg_QSMnative.nii.gz"
    img = nib.load(lab)
    a = np.asanyarray(img.dataobj).astype(np.int16)
    a[a == 1005] = 0
    nib.save(nib.Nifti1Image(a, img.affine, img.header), lab)
    qc = native / "qc" / "QCINC_native_qc.txt"
    qc.write_text("".join(l for l in qc.read_text().splitlines(True) if l.strip() != "qc_complete=1"))
    corrupt = native / "label" / "CORRUPT_aparc_aseg_QSMnative.nii.gz"
    corrupt.write_bytes(corrupt.read_bytes()[:20000])

    own = {}
    for oid in ("QREP", "T1SAME", "QSM4D", "NOCSF"):
        d = files / oid
        d.mkdir()
        own[oid] = {**s, "qsm": d / "QSM.nii", "t1": d / "T1.nii"}
        shutil.copy2(s["qsm"], own[oid]["qsm"])
        shutil.copy2(s["t1"], own[oid]["t1"])
    q = nib.load(s["qsm"])
    nib.save(nib.Nifti1Image(np.asanyarray(q.dataobj)[..., None], q.affine, q.header), own["QSM4D"]["qsm"])
    truth = np.asanyarray(nib.load(s["truth"]).dataobj)
    data = np.asanyarray(q.dataobj).copy()
    data[np.isin(truth, (4, 43))] = 0                       # reconstruction masked the ventricles
    nib.save(nib.Nifti1Image(data, q.affine, q.header), own["NOCSF"]["qsm"])
    for oid, o in own.items():
        code, out, _ = run_p8(tools, o, native, oid, qsm=o["qsm"], t1=o["t1"])
        assert code == 0, out
    _modify(own["QREP"]["qsm"])                             # replaced at the same path after P8
    _modify(own["T1SAME"]["t1"])
    t1_other = files / "T1_other.nii"
    shutil.copy2(s["t1"], t1_other)
    _modify(t1_other)

    radio = make_subject(root / "radio", seed=8, t1_axes="-x,+y,+z")
    if tools["original_p8"] is not None:
        assert run_p8(tools, s, native, "LEG_AFF", original=True)[0] == 0
        assert run_p8(tools, radio, native, "LEG_UNAFF", original=True)[0] == 0

    rows = [row("OK1", hosp="H1"), row("REP", hosp="H1"), row("FAILMNI"), row("FAILRIG"), row("FAILUNK"),
            row("QSM4D", own["QSM4D"]["qsm"], own["QSM4D"]["t1"]),
            row("NOCSF", own["NOCSF"]["qsm"], own["NOCSF"]["t1"]),
            row("QREP", own["QREP"]["qsm"], own["QREP"]["t1"]),
            row("T1SAME", own["T1SAME"]["qsm"], own["T1SAME"]["t1"]),
            row("T1MIS", t1=t1_other), row("LEG_AFF"), row("LEG_UNAFF", radio["qsm"], radio["t1"]),
            row("QCMISS"), row("QCSTALE"), row("QCINC"), row("CORRUPT"), row("ZEROSUP"), row("STALESUP"),
            row("NOREC"), row("DUP"), row("DUP", t1=t1_other), row("NONRUN", runnable="no"),
            row("NOLABEL"), row("NOQSM", qsm=files / "absent.nii"), row("NOT1", t1=files / "absent_T1.nii")]
    wl = write_worklist(root / "worklist.csv", rows)

    support = native / "support"
    make = [sys.executable, str(REVISED / "make_support_masks.py"), "--out-dir", str(support), "--jobs", "2"]
    assert subprocess.run(make + ["--worklist", str(wl)], capture_output=True).returncode == 0
    # STALESUP: mask rebuilt from another QSM; ZEROSUP: an all-zero reconstruction mask; NOREC: no record.
    other_wl = write_worklist(root / "wl_other.csv", [row("STALESUP", own["QREP"]["qsm"])])
    assert subprocess.run(make + ["--worklist", str(other_wl), "--overwrite"], capture_output=True).returncode == 0
    nib.save(nib.Nifti1Image(np.zeros(q.shape, np.uint8), q.affine), files / "ZEROSUP_recon.nii.gz")
    zero_wl = write_worklist(root / "wl_zero.csv", [row("ZEROSUP")])
    assert subprocess.run(make + ["--worklist", str(zero_wl), "--recon-mask-template",
                                  str(files / "{IID}_recon.nii.gz")], capture_output=True).returncode == 0
    (support / "NOREC_support.json").unlink()

    p1qc = root / "roi_work"
    p1qc.mkdir()
    (p1qc / "FAILMNI_qc.txt").write_text("ID=FAILMNI\nreg_mask_dice=0.85\nfnirt_rc=124\n"
                                         "WARNING: FNIRT incomplete (Did not converge); QSM_MNI not produced\n")
    (p1qc / "FAILRIG_qc.txt").write_text("ID=FAILRIG\nreg_mask_dice=0.51\nWARNING: reg_mask_dice 0.51 < 0.60\n")
    failed = root / "failed_ID.txt"
    failed.write_text("FAILMNI\tP1 exit 6\nFAILRIG\tP1 exit 5\nFAILUNK\tP1 exit 1\n")
    common = ["--failed-list", failed, "--p1-qc-dir", p1qc]
    return {"root": root, "native": native, "worklist": wl, "support": support, "common": common,
            "own": own, "subject": s}


@pytest.fixture(scope="module")
def formal(cohort, fs_fixtures):
    out = cohort["root"] / "out"
    code, log = run_p9(fs_fixtures, cohort["native"], out, cohort["worklist"], "--support-dir", cohort["support"],
                       "--csf-ref", *cohort["common"])
    assert code == 0, log
    return out, log


def test_every_subject_has_the_expected_outcome(formal, tools):
    out, _ = formal
    qc = pd.read_csv(out / "HUASHAN_NATIVE_DK_qc.csv", dtype={"IID": str}).set_index("IID")
    excl = pd.read_csv(out / "HUASHAN_NATIVE_DK_excluded.csv").groupby("IID")["reason"].apply(" | ".join)
    for oid, expected in COHORT.items():
        if oid.startswith("LEG") and tools["original_p8"] is None:
            continue
        if expected is True:
            assert qc.loc[oid, "analysis_pass"], (oid, qc.loc[oid, "fail_reasons"])
        elif expected.startswith("excluded:"):
            assert oid not in qc.index and expected[9:] in excl.get(oid, ""), (oid, excl.get(oid))
        else:
            assert not qc.loc[oid, "analysis_pass"] and expected in str(qc.loc[oid, "fail_reasons"]), \
                (oid, qc.loc[oid, "fail_reasons"])
    assert set(pd.read_csv(out / "HUASHAN_NATIVE_DK_analysis.csv")["IID"]) == set(qc.index[qc["analysis_pass"]])


def test_provenance_columns(formal):
    out, _ = formal
    qc = pd.read_csv(out / "HUASHAN_NATIVE_DK_qc.csv").set_index("IID")
    assert qc.loc["QREP", "qsm_path_matches_P8"] and not qc.loc["QREP", "qsm_matches_P8"]   # same path, new bytes
    assert not qc.loc["T1MIS", "t1_path_matches_P8"] and not qc.loc["T1MIS", "t1_matches_P8"]
    assert qc.loc["OK1", "support_type"] == "estimated_from_qsm" and qc.loc["ZEROSUP", "support_type"] == "reconstruction"
    assert qc.loc["OK1", "coverage_mode"] == "support_mask" and qc.loc["OK1", "support_provenance_ok"]
    assert qc.loc["FAILMNI", "P1_failure_stage"] == "mni_only" and qc.loc["OK1", "P8_p1_matrix_provenance"] == "fresh"
    assert qc.loc["OK1", "repeat_patient"] and qc.loc["REP", "repeat_patient"]
    assert not qc.loc["FAILMNI", "repeat_patient"]


def test_csf_reference_does_not_change_the_primary_set(formal, cohort, fs_fixtures):
    out, _ = formal
    qc = pd.read_csv(out / "HUASHAN_NATIVE_DK_qc.csv").set_index("IID")
    assert qc.loc["NOCSF", "analysis_pass"] and not qc.loc["NOCSF", "csfref_analysis_pass"]
    assert qc.loc["OK1", "csfref_analysis_pass"]
    assert "NOCSF" not in set(pd.read_csv(out / "HUASHAN_NATIVE_DK_csfref_analysis.csv")["IID"])
    out2 = cohort["root"] / "out_nocsf"
    code, log = run_p9(fs_fixtures, cohort["native"], out2, cohort["worklist"], "--support-dir", cohort["support"],
                       *cohort["common"])
    assert code == 0, log
    a1 = pd.read_csv(out / "HUASHAN_NATIVE_DK_analysis.csv")
    a2 = pd.read_csv(out2 / "HUASHAN_NATIVE_DK_analysis.csv")
    pd.testing.assert_frame_equal(a1, a2)                    # identical primary cohort and values


def test_roi_values_match_independent_calculation(formal, cohort):
    out, _ = formal
    native, s = cohort["native"], cohort["subject"]
    lab = np.asanyarray(nib.load(native / "label" / "OK1_aparc_aseg_QSMnative.nii.gz").dataobj).astype(int)
    qimg = nib.load(s["qsm"])
    q = np.asanyarray(qimg.dataobj).astype(np.float64)
    sup = np.asanyarray(nib.load(cohort["support"] / "OK1_support.nii.gz").dataobj) > 0
    seg_img = nib.load(native / "t1seg" / "OK1_aparc+aseg_T1.nii.gz")
    seg = np.asanyarray(seg_img.dataobj).astype(int)
    vox_q = abs(np.linalg.det(qimg.affine[:3, :3]))
    vox_t = abs(np.linalg.det(seg_img.affine[:3, :3]))
    mapping = pd.read_csv(out / "native_label_mapping.csv").set_index("fs_label")["csv_name"]
    vals = pd.read_csv(out / "HUASHAN_NATIVE_DK.csv").set_index("IID").loc["OK1"]
    cov = pd.read_csv(out / "native_roi_coverage.csv")
    cov = cov[cov["IID"] == "OK1"].set_index("ROI")
    checked = 0
    for fs_label, name in mapping.items():
        m = (lab == fs_label) & sup & np.isfinite(q)
        n, t1 = int(m.sum()), (seg == fs_label).sum() * vox_t
        assert cov.loc[name, "valid_qsm_voxels"] == n and cov.loc[name, "atlas_label_voxels"] == (lab == fs_label).sum()
        reported = n >= 20 and n * vox_q >= 20 and t1 > 0 and n * vox_q / t1 >= 0.5
        if reported:
            assert np.isclose(vals[f"{name}_mean"], q[m].mean(), atol=1e-6)
            assert np.isclose(vals[f"{name}_med"], np.median(q[m]), atol=1e-5)
            checked += 1
        else:
            assert np.isnan(vals[f"{name}_med"])
    assert checked == 68


def test_proxy_outputs_never_touch_formal_outputs(formal, cohort, fs_fixtures):
    out, _ = formal
    before = {p.name: _sha1(p) for p in out.iterdir()}
    code, log = run_p9(fs_fixtures, cohort["native"], out, cohort["worklist"], "--allow-nonzero-proxy",
                       *cohort["common"])
    assert code == 0, log
    after = {p.name: _sha1(p) for p in out.iterdir() if p.name in before}
    assert after == before
    proxy = pd.read_csv(out / "HUASHAN_NATIVE_DK_proxy_qc.csv")
    assert (proxy["coverage_mode"] == "nonzero_proxy").all()
    assert (out / "HUASHAN_NATIVE_DK_proxy_analysis.csv").is_file()


def test_run_info_records_inputs_and_counts(formal):
    out, _ = formal
    info = json.loads((out / "HUASHAN_NATIVE_DK_run_info.json").read_text())
    assert len(info["inputs"]["worklist"]["sha1"]) == 40 and info["coverage_mode"] == "support_mask"
    assert info["counts"]["analysis_pass"] == sum(v is True for v in COHORT.values())
    assert len(info["scripts"]["P9_extract_native_qsm_ALL.py"]) == 40


def test_include_failed_extracts_but_never_passes(cohort, fs_fixtures):
    out = cohort["root"] / "out_include_failed"
    code, log = run_p9(fs_fixtures, cohort["native"], out, cohort["worklist"], "--support-dir", cohort["support"],
                       "--include-failed", *cohort["common"])
    assert code == 0, log
    qc = pd.read_csv(out / "HUASHAN_NATIVE_DK_qc.csv").set_index("IID")
    assert not qc.loc["FAILRIG", "analysis_pass"] and "P1 failure (stage native)" in qc.loc["FAILRIG", "fail_reasons"]
    assert not qc.loc["FAILUNK", "analysis_pass"]


def test_run_with_no_passing_subject(cohort, fs_fixtures, tmp_path):
    ids = tmp_path / "ids.csv"
    pd.DataFrame({"IID": ["QCMISS", "QREP", "CORRUPT"]}).to_csv(ids, index=False)
    out = tmp_path / "out"
    code, log = run_p9(fs_fixtures, cohort["native"], out, cohort["worklist"], "--support-dir", cohort["support"],
                       "--ids-file", ids, *cohort["common"])
    assert code == 0, log
    assert "No QC-passing subjects" in log
    assert len(pd.read_csv(out / "HUASHAN_NATIVE_DK_analysis.csv")) == 0
    assert len(pd.read_csv(out / "HUASHAN_NATIVE_DK_qc.csv")) == 3


def test_run_with_no_extractable_subject_still_records_why(cohort, fs_fixtures, tmp_path):
    ids = tmp_path / "ids.csv"
    pd.DataFrame({"IID": ["NOLABEL", "NOQSM"]}).to_csv(ids, index=False)
    out = tmp_path / "out"
    code, log = run_p9(fs_fixtures, cohort["native"], out, cohort["worklist"], "--support-dir", cohort["support"],
                       "--ids-file", ids, *cohort["common"])
    assert code != 0 and "No eligible subjects" in log
    assert len(pd.read_csv(out / "HUASHAN_NATIVE_DK_excluded.csv")) == 2
    assert (out / "HUASHAN_NATIVE_DK_run_info.json").is_file()
