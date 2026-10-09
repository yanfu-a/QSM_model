"""Checks with a real FSL installation (run on the cluster; skipped where FSL is absent).

    module load fsl/6.0.4
    python -m pytest Clinical_cohort_revised/tests -m fsl -v

To also run every P8 test with the real flirt/fslinfo instead of the stubs:

    QSM_TEST_REAL_FSL=1 python -m pytest Clinical_cohort_revised/tests -v
"""

import shutil
import subprocess

import numpy as np
import nibabel as nib
import pytest

from native_grid import flirt_grid_conversion, vox_to_fsl

pytestmark = [pytest.mark.fsl,
              pytest.mark.skipif(not (shutil.which("flirt") and shutil.which("fslmaths")),
                                 reason="needs FSL (flirt, fslmaths) on PATH")]


def _fsl(*args):
    proc = subprocess.run([str(a) for a in args], capture_output=True, text=True,
                          env={**__import__("os").environ, "FSLOUTPUTTYPE": "NIFTI_GZ"})
    assert proc.returncode == 0, proc.stdout + proc.stderr


@pytest.mark.parametrize("axes,shape", [("+x,+y,+z", (64, 72, 60)), ("+x,+y,+z", (63, 72, 60)),
                                        ("-x,+y,+z", (64, 72, 60)), ("+y,+z,+x", (72, 60, 64))])
def test_fslmaths_subsamp2_grid_conversion_is_identity(tmp_path, axes, shape):
    """The assumption behind P8's 'expected identity' message (the fix itself uses headers)."""
    perm = np.zeros((3, 3))
    for col, a in enumerate(axes.split(",")):
        perm["xyz".index(a[1]), col] = 1.0 if a[0] == "+" else -1.0
    aff = np.eye(4)
    aff[:3, :3] = perm @ np.diag([1.0, 0.5, 0.5])
    aff[:3, 3] = [-30, -20, -15]
    img = nib.Nifti1Image(np.random.default_rng(0).random(shape).astype(np.float32), aff)
    img.set_sform(aff, 1)
    img.set_qform(aff, 1)
    nib.save(img, tmp_path / "full.nii.gz")
    _fsl("fslmaths", tmp_path / "full.nii.gz", "-subsamp2", tmp_path / "sub.nii.gz")
    m = flirt_grid_conversion(nib.load(tmp_path / "full.nii.gz"), nib.load(tmp_path / "sub.nii.gz"))
    assert np.allclose(m, np.eye(4), atol=1e-4), m


def test_flirt_applyxfm_with_p8_matrix_reproduces_truth(tmp_path):
    pytest.importorskip("fsl.data.image")
    from synth import make_subject
    s = make_subject(tmp_path / "s", seed=21, scale=0.5)
    m_g = flirt_grid_conversion(nib.load(s["t1"]), nib.load(s["t1_brain"]))
    np.savetxt(tmp_path / "T1_to_QSM.mat", np.linalg.inv(np.loadtxt(s["q2t"])) @ m_g, fmt="%.10f")
    _fsl("flirt", "-in", s["seg"], "-ref", s["qsm"], "-applyxfm", "-init", tmp_path / "T1_to_QSM.mat",
         "-interp", "nearestneighbour", "-datatype", "short", "-out", tmp_path / "label.nii.gz")
    lab = np.asanyarray(nib.load(tmp_path / "label.nii.gz").dataobj)
    truth = np.asanyarray(nib.load(s["truth"]).dataobj)
    cortex = truth >= 1000
    assert (lab[cortex] == truth[cortex]).mean() > 0.999


def test_rigid_registration_recovers_known_misalignment(tmp_path):
    """P1's call: flirt -usesqform -dof 6 -cost corratio -nosearch, QSM-like onto T1-like images."""
    pytest.importorskip("fsl.data.image")
    from synth import make_subject
    s = make_subject(tmp_path / "s", seed=22, scale=0.8)
    t1 = nib.load(s["t1"])
    brain = nib.Nifti1Image(np.asanyarray(t1.dataobj), t1.affine)
    nib.save(brain, tmp_path / "T1_full_ref.nii.gz")
    q = nib.load(s["qsm"])
    t = np.deg2rad(2.0)
    misalign = np.array([[1, 0, 0, 2.0], [0, np.cos(t), -np.sin(t), -1.5], [0, np.sin(t), np.cos(t), 1.0], [0, 0, 0, 1]])
    wrong = misalign @ q.affine                       # header now misplaces the QSM by a known rigid motion
    moved = nib.Nifti1Image(np.asanyarray(q.dataobj), wrong)
    moved.set_sform(wrong, 1)
    moved.set_qform(wrong, 1)
    nib.save(moved, tmp_path / "qsm_moved.nii.gz")
    ref = nib.load(tmp_path / "T1_full_ref.nii.gz")
    _fsl("flirt", "-in", tmp_path / "qsm_moved.nii.gz", "-ref", tmp_path / "T1_full_ref.nii.gz", "-usesqform",
         "-dof", "6", "-cost", "corratio", "-searchcost", "corratio", "-nosearch",
         "-omat", tmp_path / "q2t.mat", "-out", tmp_path / "q_in_t1.nii.gz")
    got = np.loadtxt(tmp_path / "q2t.mat")
    moved_img = nib.load(tmp_path / "qsm_moved.nii.gz")
    expected = vox_to_fsl(ref) @ np.linalg.inv(ref.affine) @ q.affine @ np.linalg.inv(vox_to_fsl(moved_img))
    sv = np.linalg.svd(got[:3, :3], compute_uv=False)
    assert np.allclose(sv, 1.0, atol=1e-3)                          # rigid
    rng = np.random.default_rng(0)
    pts = np.vstack([rng.uniform(20, 120, (3, 200)), np.ones((1, 200))])  # points in QSM FLIRT space
    assert np.abs((got - expected) @ pts)[:3].max() < 1.0           # < 1 mm everywhere
