"""P8 end-to-end with stub tools: label correctness, verified reuse, invalidation, failures."""

import fcntl
import hashlib
import shutil

import numpy as np
import nibabel as nib
import pytest

from conftest import run_p8

pytestmark = pytest.mark.slow


def qc_fields(native, out_id):
    text = (native / "qc" / f"{out_id}_native_qc.txt").read_text()
    return dict(line.split("=", 1) for line in text.splitlines() if "=" in line and not line.startswith(("WARNING", "ERROR")))


def agreement(native, out_id, subject):
    truth = np.asanyarray(nib.load(subject["truth"]).dataobj)
    lab = np.asanyarray(nib.load(native / "label" / f"{out_id}_aparc_aseg_QSMnative.nii.gz").dataobj)
    brain, cortex = truth > 0, truth >= 1000
    return float((lab[brain] == truth[brain]).mean()), float((lab[cortex] == truth[cortex]).mean())


@pytest.fixture(scope="session")
def baseline(tools, subject, tmp_path_factory):
    native = tmp_path_factory.mktemp("baseline") / "native"
    code, out, calls = run_p8(tools, subject, native, "S")
    assert code == 0, out
    return native


@pytest.fixture
def native(baseline, tmp_path):
    dst = tmp_path / "native"
    shutil.copytree(baseline, dst)
    return dst


def test_fresh_run_reproduces_ground_truth_and_records_provenance(baseline, subject):
    brain, cortex = agreement(baseline, "S", subject)
    assert brain > 0.9999 and cortex > 0.9999
    f = qc_fields(baseline, "S")
    assert f["qc_complete"] == "1" and f["grid_conversion"] == "flirt_fsl_v2"
    assert f["p1_matrix_provenance"] == "fresh" and f["grid_shift_mm"] == "0.0000x0.0000x0.0000"
    for key in ("qsm_sha1", "t1_sha1", "q2t_mat_sha1", "t1seg_sha1", "label_sha1"):
        assert len(f[key]) == 40
    assert float(f["q2t_scale_dev"]) < 1e-6 and float(f["qsm_pixdim_sform_rel_diff"]) < 1e-6
    assert f["p1_t1_resolution"] == "subsampled" and f["p1_search_mode"] == "nosearch"
    assert f["qc_failed_rules"] == "" and f["qc_report_only_failed_rules"] == ""
    assert float(f["synthseg_cortex_ml"]) == 433.0          # left + right cortex of the T1 segmentation
    assert (baseline / "qc" / "S_native_qc.txt").read_text().splitlines()[-1] == "qc_complete=1"


def test_original_p8_offsets_neurological_even_t1(tools, subject, tmp_path):
    if tools["original_p8"] is None:
        pytest.skip("original P8 not found")
    code, out, _ = run_p8(tools, subject, tmp_path / "native", "S", original=True)
    assert code == 0, out
    brain, cortex = agreement(tmp_path / "native", "S", subject)
    assert cortex < 0.95                                  # one-voxel shift of the cortical ribbon


def test_verified_output_is_skipped(tools, subject, native):
    code, out, calls = run_p8(tools, subject, native, "S")
    assert code == 0 and "verified for these inputs; skipping" in out and calls == []


def _append_error(native):
    with open(native / "qc" / "S_native_qc.txt", "a") as fh:
        fh.write("ERROR: seg_qsm_dice=0.41 < 0.60\n")


def _drop_complete(native):
    qc = native / "qc" / "S_native_qc.txt"
    qc.write_text("".join(l for l in qc.read_text().splitlines(True) if l.strip() != "qc_complete=1"))


def _edit_label(native):
    p = native / "label" / "S_aparc_aseg_QSMnative.nii.gz"
    img = nib.load(p)
    a = np.asanyarray(img.dataobj).copy()
    a[0, 0, 0] = 2
    nib.save(nib.Nifti1Image(a, img.affine, img.header), p)


def _drop_seg_record(native):
    (native / "t1seg" / "S_t1seg_source.txt").unlink()


def _drop_reg_settings(native):
    qc = native / "qc" / "S_native_qc.txt"
    qc.write_text("".join(l for l in qc.read_text().splitlines(True) if not l.startswith("p1_t1_resolution=")))


def _edit_matrix(native):
    p = native / "mat" / "S_QSM_to_T1.mat"
    m = np.loadtxt(p)
    m[0, 3] += 0.5
    np.savetxt(p, m, fmt="%.10f")


@pytest.mark.parametrize("change,env,reason,p1_forced,synthseg", [
    (_append_error, {}, "QC recorded errors", False, False),
    (_drop_complete, {}, "QC file incomplete", False, False),
    (_edit_label, {}, "label differs from the one QC checked", False, False),
    (_drop_seg_record, {}, "T1 segmentation source not verified", False, True),
    (_drop_reg_settings, {}, "registration settings (T1 resolution, search mode) differ", False, False),
    (_edit_matrix, {}, "QSM_to_T1.mat changed since the label was built", True, False),
    (None, {"P8_FORCE_P1": "1"}, "P1 re-registration requested", True, False),
    (None, {"FORCE_RERUN": "1"}, "P1 re-registration requested", True, False),
    (None, {"P8_FORCE_SYNTHSEG": "1"}, "SynthSeg rerun requested", False, True),
])
def test_invalid_existing_output_is_rebuilt(tools, subject, native, change, env, reason, p1_forced, synthseg):
    if change:
        change(native)
    code, out, calls = run_p8(tools, subject, native, "S", env=env)
    assert code == 0, out
    assert "not reused" in out and reason in out
    assert ("P1 FORCE_RERUN=1" in calls) == p1_forced
    assert any(c.startswith("synthseg") for c in calls) == synthseg
    assert agreement(native, "S", subject)[1] > 0.9999
    assert qc_fields(native, "S")["qc_complete"] == "1"


def test_legacy_matrix_is_reused_with_recorded_provenance(tools, subject, native):
    (native / "mat" / "S_QSM_to_T1_source.txt").unlink()
    code, out, calls = run_p8(tools, subject, native, "S", env={"P8_FORCE_RERUN": "1"})
    assert code == 0 and calls == ["P1 FORCE_RERUN=0"]
    assert qc_fields(native, "S")["p1_matrix_provenance"] == "legacy_path_mtime"


@pytest.mark.parametrize("which", ["qsm", "t1"])
def test_input_content_change_at_same_path_is_detected(tools, subject, tmp_path, which):
    own = dict(subject)
    for key in ("qsm", "t1"):
        own[key] = tmp_path / subject[key].name
        shutil.copy2(subject[key], own[key])
    native = tmp_path / "native"
    assert run_p8(tools, own, native, "S", qsm=own["qsm"], t1=own["t1"])[0] == 0
    img = nib.load(own[which])
    data = np.asanyarray(img.dataobj).copy()
    data[data != 0] += 1
    stat = own[which].stat()
    nib.save(nib.Nifti1Image(data, img.affine, img.header), own[which])
    import os
    os.utime(own[which], (stat.st_atime, stat.st_mtime - 3600))   # older mtime: content check, not timestamps
    code, out, calls = run_p8(tools, own, native, "S", qsm=own["qsm"], t1=own["t1"])
    assert code == 0, out
    assert f"{which.upper()} content differs" in out
    assert "P1 FORCE_RERUN=1" in calls
    assert any(c.startswith("synthseg") for c in calls) == (which == "t1")


def test_failure_after_cleanup_leaves_no_stale_outputs(tools, subject, native):
    code, out, _ = run_p8(tools, subject, native, "S", env={"P8_FORCE_RERUN": "1", "FAKE_FLIRT_FAIL": "1"})
    assert code == 14
    assert not (native / "qc" / "S_native_qc.txt").exists()
    assert not (native / "label" / "S_aparc_aseg_QSMnative.nii.gz").exists()


def test_missing_registration_dice_fails_qc(tools, subject, native):
    code, out, _ = run_p8(tools, subject, native, "S", env={"P8_FORCE_P1": "1", "FAKE_NO_DICE": "1"})
    assert code == 15
    text = (native / "qc" / "S_native_qc.txt").read_text()
    assert "ERROR: reg_mask_dice unavailable" in text and text.endswith("qc_complete=1\n")
    code, out, _ = run_p8(tools, subject, native, "S")
    assert code == 15 and "QC recorded errors" in out          # not silently skipped


def test_failed_registration_removes_matrix_record(tools, subject, native):
    code, out, _ = run_p8(tools, subject, native, "S", env={"P8_FORCE_P1": "1", "FAKE_P1_CRASH": "1"})
    assert code == 3
    assert not (native / "mat" / "S_QSM_to_T1_source.txt").exists()
    assert not (native / "qc" / "S_native_qc.txt").exists()


def test_low_registration_dice_keeps_label_and_fails_qc(tools, subject, native):
    """P1 no longer stops at its Dice gate in rigid mode; P8's own 0.60 gate still fails the scan."""
    env = {"P8_FORCE_P1": "1", "FAKE_P1_DICE": "0.410"}
    code, out, _ = run_p8(tools, subject, native, "S", env=env)
    assert code == 15, out
    assert (native / "label" / "S_aparc_aseg_QSMnative.nii.gz").is_file()
    assert (native / "mat" / "S_QSM_to_T1_source.txt").is_file()
    f = qc_fields(native, "S")
    text = (native / "qc" / "S_native_qc.txt").read_text()
    assert "ERROR: reg_mask_dice=0.410 < 0.60" in text and f["qc_failed_rules"] == "reg_mask_dice"
    assert f["qc_complete"] == "1" and float(f["reg_mask_dice"]) == 0.41
    code, out, calls = run_p8(tools, subject, native, "S", env={"FAKE_P1_DICE": "0.410"})
    assert code == 15 and "QC recorded errors" in out and calls == ["P1 FORCE_RERUN=0"]


def test_registration_settings_are_recorded_and_enforced(tools, subject, native):
    code, out, calls = run_p8(tools, subject, native, "S", env={"QSM_ROI_T1_SUBSAMP": "0"})
    assert code == 0, out
    assert "registration settings (T1 resolution, search mode) differ" in out
    assert "registered with T1 resolution 'subsampled'" in out and "P1 FORCE_RERUN=1" in calls
    assert qc_fields(native, "S")["p1_t1_resolution"] == "native"
    code, out, calls = run_p8(tools, subject, native, "S", env={"QSM_ROI_T1_SUBSAMP": "0"})
    assert code == 0 and "skipping" in out and calls == []
    code, out, calls = run_p8(tools, subject, native, "S", env={"QSM_ROI_SEARCH_MODE": "search"})
    assert code == 0 and "P1 FORCE_RERUN=1" in calls
    f = qc_fields(native, "S")
    assert f["p1_t1_resolution"] == "subsampled" and f["p1_search_mode"] == "search"


def test_legacy_p1_qc_without_settings_is_re_registered(tools, subject, native):
    (native / "mat" / "S_QSM_to_T1_source.txt").unlink()
    p1qc = native / "p1_rigid" / "S_qc_rigid.txt"
    p1qc.write_text("".join(l for l in p1qc.read_text().splitlines(True)
                            if not l.startswith(("t1_pixdim=", "search_mode="))))
    code, out, calls = run_p8(tools, subject, native, "S", env={"P8_FORCE_RERUN": "1"})
    assert code == 0, out
    assert "T1 resolution 'not recorded'" in out and calls == ["P1 FORCE_RERUN=1"]
    assert qc_fields(native, "S")["p1_matrix_provenance"] == "fresh"



@pytest.mark.parametrize("change", ["delete", "edit"])
def test_segmentation_reuse_requires_its_volume_table(tools, subject, native, change):
    """ICV and the T1 cortex volume come from SynthSeg's volume table, so it is part of the segmentation."""
    vol = native / "qc" / "S_synthseg_vol.csv"
    if change == "delete":
        vol.unlink()
    else:
        vol.write_text(vol.read_text().replace("1500000", "1100000"))
    code, out, calls = run_p8(tools, subject, native, "S")
    assert code == 0, out
    assert "SynthSeg volume table changed since the label was built" in out
    assert any(c.startswith("synthseg") for c in calls)                     # segmentation regenerated
    f = qc_fields(native, "S")
    assert float(f["synthseg_icv_ml"]) == 1500 and float(f["synthseg_cortex_ml"]) == 433.0
    assert f["synthseg_vol_sha1"] == hashlib.sha1(vol.read_bytes()).hexdigest()


def test_missing_volume_measures_fail_or_stop(tools, subject, native):
    code, out, _ = run_p8(tools, subject, native, "S", env={"P8_FORCE_SYNTHSEG": "1", "FAKE_SYNTHSEG_NO_VOL": "1"})
    assert code == 13 and "did not write" in out
    code, out, _ = run_p8(tools, subject, native, "S", env={"P8_FORCE_SYNTHSEG": "1", "FAKE_SYNTHSEG_NO_ICV": "1"})
    assert code == 15, out
    f = qc_fields(native, "S")
    text = (native / "qc" / "S_native_qc.txt").read_text()
    assert "ERROR: SynthSeg ICV unavailable" in text and f["qc_failed_rules"] == "synthseg_icv"

def test_concurrent_run_is_refused(tools, subject, native):
    if shutil.which("flock") is None:
        pytest.skip("flock not available")
    with open(native / "logs" / "S.lock", "w") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
        code, out, calls = run_p8(tools, subject, native, "S", env={"P8_FORCE_RERUN": "1"})
    assert code == 16 and "Another P8 run" in out and calls == []
