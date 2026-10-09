"""FLIRT coordinate conventions of native_grid.py (the T1 -> T1_brain grid conversion in P8)."""

import types

import numpy as np
import nibabel as nib
import pytest

from native_grid import flirt_grid_conversion, is_neurological, legacy_grid_conversion, vox_to_fsl

fslpy = pytest.importorskip("fsl.data.image", reason="fslpy provides FSL's reference conventions")


def _random_header(rng, mode):
    shape = tuple(int(x) for x in rng.integers(20, 300, 3))
    zooms = rng.uniform(0.4, 2.5, 3)
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    q *= np.sign(np.linalg.det(q))                           # proper rotation
    a = np.eye(4)
    a[:3, :3] = q @ np.diag(zooms) @ np.diag([rng.choice([-1, 1]), 1, 1])
    a[:3, 3] = rng.uniform(-120, 120, 3)
    hdr = nib.Nifti1Header()
    hdr.set_data_shape(shape)
    hdr.set_zooms(tuple(zooms))
    if mode == 0:
        hdr.set_sform(a, 1)
        hdr.set_qform(a, 1)
    elif mode == 1:
        hdr.set_sform(a, 1)
        hdr.set_qform(None, 0)
    else:
        hdr.set_sform(None, 0)
        hdr.set_qform(a, 1)
    # What nib.load() gives: affine = the header's best affine (sform, then qform).
    ours = types.SimpleNamespace(affine=hdr.get_best_affine(), shape=shape, header=hdr)
    ref = fslpy.Image(nib.Nifti1Image(np.zeros(shape, np.uint8), None, hdr))
    return ours, ref


def test_matches_fslpy_on_random_oblique_grids():
    rng = np.random.default_rng(1)
    worst, n_neuro = 0.0, 0
    for trial in range(300):
        (s, s_ref), (d, d_ref) = _random_header(rng, trial % 3), _random_header(rng, (trial + 1) % 3)
        assert np.allclose(vox_to_fsl(s), s_ref.getAffine("voxel", "fsl"), atol=1e-6)
        assert is_neurological(s) == s_ref.isNeurological()
        n_neuro += is_neurological(s)
        ref = d_ref.getAffine("world", "fsl") @ s_ref.getAffine("fsl", "world")
        worst = max(worst, float(np.abs(flirt_grid_conversion(s, d) - ref).max()))
    assert 50 < n_neuro < 250                                 # both storage orientations covered
    assert worst < 1e-4                                        # float32 header precision


def _subsamp2_pair(axes, zooms, shape):
    """Full grid and the grid fslmaths -subsamp2 gives (centred in FSL's radiological voxel order)."""
    perm = np.zeros((3, 3))
    for col, a in enumerate(axes.split(",")):
        perm["xyz".index(a[1]), col] = 1.0 if a[0] == "+" else -1.0
    af = np.eye(4)
    af[:3, :3] = perm @ np.diag(zooms)
    af[:3, 3] = [-90, -110, -80]
    neuro = np.linalg.det(af[:3, :3]) > 0
    c0 = 1.0 if (neuro and shape[0] % 2 == 0) else 0.0
    asub = af @ np.array([[2, 0, 0, c0], [0, 2, 0, 0], [0, 0, 2, 0], [0, 0, 0, 1.0]])
    full = nib.Nifti1Image(np.zeros(shape, np.uint8), af)
    sub = nib.Nifti1Image(np.zeros(tuple((n + 1) // 2 for n in shape), np.uint8), asub)
    for img, a in ((full, af), (sub, asub)):
        img.set_sform(a, 1)
        img.set_qform(a, 1)
    return full, sub, neuro


@pytest.mark.parametrize("axes,zooms,shape", [
    ("+x,+y,+z", (1, 1, 1), (180, 216, 180)),     # neurological, even first dimension
    ("+x,+y,+z", (1, 1, 1), (179, 216, 180)),     # neurological, odd
    ("-x,+y,+z", (1, 1, 1), (180, 216, 180)),     # radiological, even
    ("-x,+y,+z", (1, 1, 1), (179, 216, 181)),     # radiological, odd
    ("+x,+y,+z", (1, 0.5, 0.5), (176, 512, 512)),  # anisotropic sagittal grid
    ("-x,+y,+z", (1, 0.5, 0.5), (176, 512, 512)),
    ("+y,+z,+x", (0.5, 0.5, 1), (360, 260, 150)),  # permuted axes, neurological
    ("+y,-z,+x", (0.5, 0.5, 1), (360, 260, 150)),  # permuted axes, radiological
])
def test_subsamp2_grid_conversion_is_identity(axes, zooms, shape):
    full, sub, neuro = _subsamp2_pair(axes, zooms, shape)
    fixed = flirt_grid_conversion(full, sub)
    assert np.allclose(fixed, np.eye(4), atol=1e-6)
    # Against FSL's own library, not only against the expectation of identity.
    ref = fslpy.Image(sub).getAffine("world", "fsl") @ fslpy.Image(full).getAffine("fsl", "world")
    assert np.allclose(fixed, ref, atol=1e-6)
    legacy = legacy_grid_conversion(full, sub)
    affected = neuro and shape[0] % 2 == 0
    if affected:
        # The original P8 shifted by one full-resolution voxel along the first voxel axis.
        assert np.isclose(abs(legacy[0, 3] - fixed[0, 3]), zooms[0])
        assert np.allclose(legacy[1:3, 3], fixed[1:3, 3])
    else:
        assert np.allclose(legacy, fixed, atol=1e-6)
