"""make_support_masks.py: hole filling, provenance records and cache invalidation."""

import json
import subprocess
import sys

import numpy as np
import nibabel as nib
import pandas as pd

from conftest import REVISED, write_worklist
from make_support_masks import provenance_path, support_from_qsm


def _brain_qsm(dtype=np.float32):
    q = np.zeros((60, 60, 60), np.float32)
    zz, yy, xx = np.indices(q.shape)
    r2 = (xx - 30) ** 2 + (yy - 30) ** 2 + (zz - 30) ** 2
    q[r2 < 25 ** 2] = 7.0
    q[30, 30, 10] = 0                          # single enclosed zero
    q[20:22, 20:22, 30] = 0                    # 4-voxel enclosed zero cluster
    q[r2 < 6 ** 2] = 0                         # masked "ventricle" (~900 voxels)
    q[30, 30, 54:56] = 0                       # zero notch open to the outside
    return q.astype(dtype)


def test_hole_filling_keeps_only_small_enclosed_zeros():
    q = _brain_qsm()
    zz, yy, xx = np.indices(q.shape)
    vent = (xx - 30) ** 2 + (yy - 30) ** 2 + (zz - 30) ** 2 < 6 ** 2
    sup, n, filled, largest = support_from_qsm(q, 1.0, 10.0)
    assert (n, filled, largest) == (2, 5, 4.0)
    assert sup[30, 30, 10] and sup[20:22, 20:22, 30].all()
    assert not sup[vent].any() and not sup[30, 30, 55] and not sup[0, 0, 0]
    # Sensitivity to the hole-size limit: 0 mm3 = nonzero proxy; a huge limit also fills the ventricle.
    assert support_from_qsm(q, 1.0, 0.0)[2] == 0
    assert support_from_qsm(q, 1.0, 1e6)[0][vent].all()


def test_quantized_integer_qsm_limits_of_estimated_masks():
    """Documents a limitation (METHODOLOGICAL_DECISIONS.md): with quantized QSM, zero-valued
    voxels touching the brain edge cannot be told apart from background, so an estimated mask
    misses part of the outermost shell, where the cortex lies; deeper zeros are recovered."""
    rng = np.random.default_rng(0)
    q = np.zeros((60, 60, 60), np.float32)
    zz, yy, xx = np.indices(q.shape)
    r = np.sqrt((xx - 30) ** 2 + (yy - 30) ** 2 + (zz - 30) ** 2)
    inside = r < 25
    q[inside] = rng.normal(0, 3, inside.sum())
    qi = np.rint(q).astype(np.int16)           # 1 ppb quantization, SD 3 ppb: ~13% exact zeros
    sup = support_from_qsm(qi, 1.0, 10.0)[0]
    outer, deep = inside & (r >= 22), inside & (r < 22)
    assert np.mean(qi[inside] == 0) > 0.10
    assert np.mean(~sup[deep]) < 0.01                    # interior zeros recovered
    assert 0.02 < np.mean(~sup[outer]) < 0.10            # edge zeros partly lost (measured 5.6%)
    assert np.mean(~support_from_qsm(qi, 1.0, 0.0)[0][deep]) > 0.10   # nonzero proxy loses all zeros


def _run(worklist, out_dir, *extra):
    proc = subprocess.run([sys.executable, str(REVISED / "make_support_masks.py"), "--worklist", str(worklist),
                           "--out-dir", str(out_dir), "--jobs", "1", *map(str, extra)],
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    return pd.read_csv(out_dir / "support_masks_summary.csv", dtype=str).set_index("IID")["status"]


def test_provenance_reuse_and_invalidation(tmp_path):
    aff = np.diag([-1.0, 1.0, 1.0, 1.0])
    qsm = tmp_path / "Q.nii"
    nib.save(nib.Nifti1Image(_brain_qsm(), aff), qsm)
    wl = write_worklist(tmp_path / "wl.csv", [{"image_dir_id": "S", "qsm_file": str(qsm), "t1_file": "t1"}])
    out = tmp_path / "support"
    mask = out / "S_support.nii.gz"

    assert _run(wl, out)["S"] == "written"
    record = json.loads(provenance_path(mask).read_text())
    assert record["support_type"] == "estimated_from_qsm" and record["max_hole_mm3"] == 10.0
    assert record["filled_voxels"] == 5 and record["qsm_nonzero_outside_support"] == 0
    assert _run(wl, out)["S"] == "kept (provenance verified)"
    assert _run(wl, out, "--max-hole-mm3", "20")["S"] == "stale (provenance changed)"   # parameter change
    q2 = _brain_qsm()
    q2[q2 != 0] += 1
    nib.save(nib.Nifti1Image(q2, aff), qsm)                                             # same path, new content
    assert _run(wl, out, "--max-hole-mm3", "20")["S"] == "stale (provenance changed)"
    assert json.loads(provenance_path(mask).read_text())["qsm_sha1"] != record["qsm_sha1"]
    m = nib.load(mask)
    nib.save(nib.Nifti1Image(np.zeros(m.shape, np.uint8), m.affine), mask)              # tampered mask
    assert _run(wl, out, "--max-hole-mm3", "20")["S"] == "stale (provenance changed)"
    provenance_path(mask).unlink()                                                      # record lost
    assert _run(wl, out, "--max-hole-mm3", "20")["S"] == "written"


def test_support_keeps_identical_rows_and_removes_old_conflicting_summary(tmp_path):
    qsm = tmp_path / "Q.nii"
    nib.save(nib.Nifti1Image(_brain_qsm(), np.eye(4)), qsm)
    base = lambda oid: {"image_dir_id": oid, "qsm_file": str(qsm), "t1_file": "t1"}
    out = tmp_path / "support"
    first = write_worklist(tmp_path / "first.csv", [base("SAME"), base("DIFFERENT")])
    assert set(_run(first, out).index) == {"SAME", "DIFFERENT"}
    changed = write_worklist(tmp_path / "changed.csv", [
        base("SAME"), base("SAME"), base("DIFFERENT"),
        {**base("DIFFERENT"), "qsm_coverage_mm": "120"}])
    status = _run(changed, out)
    assert set(status.index) == {"SAME"}
    assert status["SAME"] == "kept (provenance verified)"


def test_reconstruction_masks_are_registered_and_grid_checked(tmp_path):
    aff = np.diag([-1.0, 1.0, 1.0, 1.0])
    q = _brain_qsm()
    nib.save(nib.Nifti1Image(q, aff), tmp_path / "Q.nii")
    recon = (q != 0).astype(np.uint8)
    recon[q.shape[0] // 2:] = 0                         # a mask smaller than the nonzero region
    nib.save(nib.Nifti1Image(recon, aff), tmp_path / "S_recon.nii.gz")
    nib.save(nib.Nifti1Image(recon, np.diag([-2.0, 1, 1, 1])), tmp_path / "T_recon.nii.gz")
    wl = write_worklist(tmp_path / "wl.csv", [
        {"image_dir_id": "S", "qsm_file": str(tmp_path / "Q.nii"), "t1_file": "t1"},
        {"image_dir_id": "T", "qsm_file": str(tmp_path / "Q.nii"), "t1_file": "t1"},
        {"image_dir_id": "U", "qsm_file": str(tmp_path / "Q.nii"), "t1_file": "t1"}])
    out = tmp_path / "support"
    status = _run(wl, out, "--recon-mask-template", str(tmp_path / "{IID}_recon.nii.gz"))
    assert status["S"] == "written"
    assert status["T"].startswith("ERROR reconstruction mask is not on the QSM grid")
    assert status["U"].startswith("ERROR reconstruction mask missing")
    record = json.loads(provenance_path(out / "S_support.nii.gz").read_text())
    assert record["support_type"] == "reconstruction" and record["max_hole_mm3"] is None
    assert record["qsm_nonzero_outside_support"] > 0    # reported for review, not hidden


def test_quantization_is_recorded_and_float_hole_filling_is_flagged(tmp_path):
    aff = np.diag([-1.0, 1.0, 1.0, 1.0])
    q = _brain_qsm() * 1.37                            # non-integer float QSM with small enclosed zero holes
    nib.save(nib.Nifti1Image(q, aff), tmp_path / "F.nii")
    scaled = nib.Nifti1Image(np.rint(q * 10).astype(np.int16), aff)
    scaled.header.set_slope_inter(0.001, 0.0)          # stored integers, not integer-valued after scaling
    nib.save(scaled, tmp_path / "I.nii")
    wl = write_worklist(tmp_path / "wl.csv", [{"image_dir_id": oid, "qsm_file": str(tmp_path / f"{oid}.nii"),
                                               "t1_file": "t1"} for oid in ("F", "I")])
    proc = subprocess.run([sys.executable, str(REVISED / "make_support_masks.py"), "--worklist", str(wl),
                           "--out-dir", str(tmp_path / "support"), "--jobs", "1"], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    rec = {oid: json.loads(provenance_path(tmp_path / "support" / f"{oid}_support.nii.gz").read_text())
           for oid in ("F", "I")}
    assert rec["F"]["qsm_quantized"] is False and rec["F"]["filled_voxels"] > 0
    assert rec["I"]["qsm_quantized"] is True and rec["I"]["qsm_integer_valued"] is False
    assert "1 non-quantized QSM maps had enclosed zero holes filled" in proc.stdout
