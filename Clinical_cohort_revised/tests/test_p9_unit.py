"""Unit tests of P9's extraction, gating, parsing and provenance helpers."""

import importlib.util
import json
import sys

import numpy as np
import nibabel as nib
import pandas as pd
import pytest

from conftest import REVISED

sys.dont_write_bytecode = True
_spec = importlib.util.spec_from_file_location("p9", REVISED / "P9_extract_native_qsm_ALL.py")
p9 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(p9)
KEYS = sorted(p9.FS_L + p9.FS_R)


def reference_stats(lab, qsm, valid, keys):
    """Independent reference: pandas groupby over the valid labelled voxels."""
    df = pd.DataFrame({"label": lab.ravel(), "q": qsm.ravel().astype(np.float64), "ok": valid.ravel()})
    labelled = df[df["label"].isin(keys)]
    total = labelled.groupby("label").size().reindex(keys, fill_value=0)
    g = labelled[labelled["ok"]].groupby("label")["q"]
    return (total.to_numpy(), g.size().reindex(keys, fill_value=0).to_numpy(),
            g.mean().reindex(keys).to_numpy(), g.median().reindex(keys).to_numpy())


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_roi_statistics_match_independent_reference(seed):
    rng = np.random.default_rng(seed)
    shape = (40, 44, 30)
    pool = np.array(KEYS + [0, 2, 41, 4, 43, 24, 1004, 2004, 3000, -5])
    lab = rng.choice(pool, size=shape).astype(np.int32)
    qsm = rng.normal(10, 30, size=shape).astype(np.float32)
    qsm[rng.random(shape) < 0.1] = 0
    qsm[rng.random(shape) < 0.01] = np.nan
    valid = (rng.random(shape) < 0.85) & np.isfinite(qsm)
    res = p9.extract_one(lab, qsm, valid, KEYS, 1.0, np.full(len(KEYS), np.nan), 1, 0.0, 0.0)
    total, n, mean, median = reference_stats(lab, qsm, valid, KEYS)
    assert np.array_equal(res["label_counts"], total)
    assert np.array_equal(res["counts"], n)
    assert np.allclose(res["means"], mean, rtol=0, atol=1e-9, equal_nan=True)
    assert np.allclose(res["medians"], median, rtol=0, atol=1e-6, equal_nan=True)


def _one_roi(n_valid, n_label, vox_mm3, t1_mm3, **gates):
    lab = np.zeros((n_label,), np.int32) + 1001
    qsm = np.arange(n_label, dtype=np.float32) + 1
    valid = np.zeros(n_label, bool)
    valid[:n_valid] = True
    keys = [1001]
    g = {"min_voxels": 20, "min_mm3": 20.0, "min_coverage": 0.5, **gates}
    return p9.extract_one(lab.reshape(-1, 1, 1), qsm.reshape(-1, 1, 1), valid.reshape(-1, 1, 1), keys,
                          vox_mm3, np.array([t1_mm3], float), g["min_voxels"], g["min_mm3"], g["min_coverage"])


def test_voxel_volume_and_coverage_thresholds_are_inclusive():
    assert np.isfinite(_one_roi(20, 40, 1.0, 40.0)["medians"][0])          # 20 vox, 20 mm3, 0.5 coverage
    assert np.isnan(_one_roi(19, 40, 1.0, 38.0)["medians"][0])             # 19 voxels
    r = _one_roi(30, 40, 0.512, 40.0)                                      # 30 vox = 15.4 mm3 < 20 mm3
    assert np.isnan(r["medians"][0]) and r["n_low_voxels"] == 1
    r = _one_roi(39, 60, 0.512, 40.0)                                      # 19.97 mm3, just below
    assert np.isnan(r["medians"][0])
    r = _one_roi(40, 60, 0.512, 40.0)                                      # 20.48 mm3, coverage 0.512
    assert np.isfinite(r["medians"][0])
    r = _one_roi(25, 40, 1.0, 60.0)                                        # coverage 0.417
    assert np.isnan(r["medians"][0]) and r["n_low_coverage"] == 1
    assert r["status"] == "No ROI value passes the voxel, volume and coverage gates"


@pytest.mark.parametrize("t1", [0.0, -1.0, np.inf])
def test_roi_without_finite_positive_t1_volume_fails_subject(t1):
    lab = np.array([1001] * 50 + [1002] * 50).reshape(10, 10, 1)
    q = np.full(lab.shape, 12.0, np.float32)
    r = p9.extract_one(lab, q, np.ones(lab.shape, bool), [1001, 1002], 1.0, np.array([t1, 100.0]),
                       20, 20.0, 0.5)
    assert np.isnan(r["coverage"][0]) and np.isnan(r["medians"][0])
    assert r["status"].startswith("Label and T1 segmentation disagree")


def test_coverage_gate_needs_t1_but_gate_off_does_not():
    lab = np.full((10, 10, 1), 1001)
    q = np.full(lab.shape, 12.0, np.float32)
    ok = np.ones(lab.shape, bool)
    no_t1 = np.array([np.nan])
    assert p9.extract_one(lab, q, ok, [1001], 1.0, no_t1, 20, 20.0, 0.5)["status"].startswith("No ROI value")
    assert p9.extract_one(lab, q, ok, [1001], 1.0, no_t1, 20, 20.0, 0.0)["status"] == "OK"


def test_zero_supported_voxels_keeps_label_counts():
    lab = np.full((10, 10, 1), 1001)
    r = p9.extract_one(lab, np.ones(lab.shape, np.float32), np.zeros(lab.shape, bool), [1001], 1.0,
                       np.array([100.0]), 20, 20.0, 0.5)
    assert r["status"] == "No valid supported cortical voxels"
    assert r["label_counts"][0] == 100 and r["counts"][0] == 0


def test_singleton_4d_qsm_and_support_are_read_as_3d(tmp_path):
    aff = np.diag([-1.0, 1.0, 1.0, 1.0])
    lab = np.full((6, 6, 6), 1001, np.int16)
    q = np.full((6, 6, 6), 5.0, np.float32)
    nib.save(nib.Nifti1Image(lab, aff), tmp_path / "lab.nii.gz")
    nib.save(nib.Nifti1Image(q[..., None], aff), tmp_path / "q4.nii.gz")
    nib.save(nib.Nifti1Image(np.ones((6, 6, 6, 1), np.uint8), aff), tmp_path / "s4.nii.gz")
    lab3, q3, valid, _, vox = p9.load_native_grid(tmp_path / "lab.nii.gz", tmp_path / "q4.nii.gz",
                                                  tmp_path / "s4.nii.gz")
    assert q3.shape == (6, 6, 6) and valid.shape == (6, 6, 6) and valid.all() and vox == 1.0


def test_grid_mismatch_raises(tmp_path):
    nib.save(nib.Nifti1Image(np.zeros((6, 6, 6), np.int16), np.eye(4)), tmp_path / "lab.nii.gz")
    nib.save(nib.Nifti1Image(np.zeros((6, 6, 6), np.float32), np.diag([2, 1, 1, 1.0])), tmp_path / "q.nii.gz")
    with pytest.raises(ValueError, match="affine"):
        p9.load_native_grid(tmp_path / "lab.nii.gz", tmp_path / "q.nii.gz")


def test_csf_reference_uses_eroded_ventricles_and_support():
    lab = np.zeros((20, 20, 20), np.int32)
    lab[2:8, 2:8, 2:8] = 4
    lab[12:18, 12:18, 12:18] = 43
    q = np.full(lab.shape, 3.0, np.float32)
    q[lab == 4] = 7.0
    med, n = p9.native_csf_reference(lab, q, np.ones(lab.shape, bool), 30)
    assert n == 2 * 4 ** 3 and med == 5.0          # 64 voxels at 7 and 64 at 3 after 1-voxel erosion
    assert np.isnan(p9.native_csf_reference(lab, q, np.zeros(lab.shape, bool), 30)[0])


@pytest.mark.parametrize("text,stage", [
    ("WARNING: QSM datatype=128 is unsupported (non-scalar; 128=RGB).", "native"),
    ("WARNING: reg_mask_dice 0.510 < 0.60 (often due to failed T1 skull stripping;", "native"),
    ("WARNING: FNIRT incomplete (Did not converge within 3600s); QSM_MNI not produced", "mni_only"),
    ("WARNING: t1_mni_dice 0.55 < 0.60", "mni_only"),
    ("WARNING: mni_brain_coverage 0.40 < 0.55", "mni_only"),
    ("WARNING: Output dimensions are not 91x109x91", "mni_only"),
    ("WARNING: reg_mask_dice 0.5 < 0.60\nWARNING: t1_mni_dice 0.5 < 0.60", "native"),
    ("reg_mask_dice=0.85\nmni_brain_coverage=0.90", "unknown"),
    ("WARNING: something else", "unknown"),
])
def test_p1_failure_stage_classification(tmp_path, text, stage):
    qc = tmp_path / "X_qc.txt"
    qc.write_text("ID=X\n" + text + "\n")
    assert p9.classify_p1_failure(qc)[0] == stage
    assert p9.classify_p1_failure(tmp_path / "missing_qc.txt")[0] == "unknown"


def test_failed_list_parsing(tmp_path):
    f = tmp_path / "failed.txt"
    f.write_text("A\tP1 exit 6\nB\\tliteral tab\n  C extra words\n\n")
    assert set(p9.read_failed_ids(f)) == {"A", "B", "C"}
    assert p9.read_failed_ids(tmp_path / "absent.txt") == {}


def test_duplicate_ids_conflicting_t1_are_excluded():
    wl = pd.DataFrame({"image_dir_id": ["A", "A", "B", "B", "C"],
                       "qsm_file": ["q", "q", "qb", "qb", "qc"],
                       "t1_file": ["t1", "t1_other", "tb", "tb", "tc"],
                       "hospital_id": ["1", "1", "2", "2", "3"]})
    kept, excluded = p9.resolve_duplicate_ids(wl)
    assert list(kept["image_dir_id"]) == ["B", "C"]
    assert excluded[0]["IID"] == "A" and "t1_file=t1 | t1_other" in excluded[0]["reason"]


def test_support_record_reading(tmp_path):
    mask = tmp_path / "S_support.nii.gz"
    assert p9.read_support_record(mask) == {}
    p9.support_record_path(mask).write_text(json.dumps({"support_type": "estimated_from_qsm"}))
    assert p9.read_support_record(mask)["support_type"] == "estimated_from_qsm"
    p9.support_record_path(mask).write_text("{not json")
    assert p9.read_support_record(mask) is None


def test_label_names_are_checked_against_freesurfer(fs_fixtures):
    labels = p9.read_labels(fs_fixtures / "Desikan_labels.csv")
    fs_to_name, rows = p9.build_label_map(labels, fs_fixtures / "Desikan_labels.csv", fs_fixtures)
    assert len(fs_to_name) == 68 and len(set(fs_to_name.values())) == 68
    assert fs_to_name[1001] == "L_bankssts" and fs_to_name[2035] == "R_insula"
    bad = p9.read_labels(fs_fixtures / "Desikan_labels_alphabetical.csv")
    with pytest.raises(SystemExit, match="do not match"):
        p9.build_label_map(bad, fs_fixtures / "Desikan_labels_alphabetical.csv", fs_fixtures)


def test_concordance_and_test_retest_run(capsys):
    rng = np.random.default_rng(0)
    names = [f"L_r{i}" for i in range(5)] + [f"R_r{i}" for i in range(5)]
    n = 60
    offset = rng.normal(0, 10, n)
    df = pd.DataFrame({c: offset + rng.normal(0, 3, n) for c in names})
    df.insert(0, "IID", [f"P{i:03d}" for i in range(n)])
    df.insert(1, "hospital_id", [str(i) for i in range(n)])
    for i in range(6):
        df.loc[50 + i, "hospital_id"] = str(i)
    assert p9.lr_concordance(df, names, "t") is not None
    assert p9.test_retest(df, names, "t") is not None
    out = capsys.readouterr().out
    assert "not test-retest reliability" in out and "Test-retest" in out
