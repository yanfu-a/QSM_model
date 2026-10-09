"""P8 end-to-end with stub tools: label correctness, verified reuse, invalidation, failures."""

import fcntl
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
    code, out, _ = run_p8(tools, subject, native, "S", env={"P8_FORCE_P1": "1", "FAKE_P1_FAIL": "1"})
    assert code == 5
    assert not (native / "mat" / "S_QSM_to_T1_source.txt").exists()
    assert not (native / "qc" / "S_native_qc.txt").exists()


def test_concurrent_run_is_refused(tools, subject, native):
    if shutil.which("flock") is None:
        pytest.skip("flock not available")
    with open(native / "logs" / "S.lock", "w") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
        code, out, calls = run_p8(tools, subject, native, "S", env={"P8_FORCE_RERUN": "1"})
    assert code == 16 and "Another P8 run" in out and calls == []
