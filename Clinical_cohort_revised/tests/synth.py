"""Synthetic subjects with known ground truth for the regression tests.

A subject is an ellipsoidal "brain" in world space with SynthSeg-like labels: cerebral
white matter (2/41), lateral ventricles (4/43) and a 7 mm cortical shell split into the 34
Desikan-Killiany parcels per hemisphere. It provides:

* T1.nii and seg_T1.nii.gz on a configurable T1 grid (storage orientation, voxel size,
  dimensions, axis permutation);
* T1_brain.nii.gz on the grid fslmaths -subsamp2 produces (centred subsampling in FSL's
  radiological voxel order: stored offset (1, 0, 0) only for neurological images with an
  even first dimension);
* QSM.nii, a radiological axial slab of configurable voxel size, coverage, data type and
  dimensionality, with tissue-dependent values and an eroded reconstruction mask;
* Q2T.mat, the true QSM -> T1_brain FLIRT matrix (computed with fslpy, FSL's own library);
* truth_label_QSM.nii.gz, the T1 labels sampled at QSM voxel centres straight through world
  space (no FLIRT conventions involved), against which P8's label is compared.
"""

import numpy as np
import nibabel as nib

FS_L = [1001, 1002, 1003] + list(range(1005, 1036))
FS_R = [k + 1000 for k in FS_L]
RADII = np.array([70.0, 85.0, 60.0])


def _rotation(deg_xyz):
    out = np.eye(3)
    for axis, deg in enumerate(deg_xyz):
        t = np.deg2rad(deg)
        c, s = np.cos(t), np.sin(t)
        r = np.eye(3)
        i, j = [k for k in range(3) if k != axis]
        r[i, i], r[i, j], r[j, i], r[j, j] = c, -s, s, c
        out = out @ r
    return out


def tissue(w, radii):
    """Labels, depth below the surface (mm) and parcel index for world points w (3, N)."""
    x, y, z = w
    rx, ry, rz = radii
    re = np.sqrt((x / rx) ** 2 + (y / ry) ** 2 + (z / rz) ** 2)
    depth = (1 - re) * np.sqrt(x * x + y * y + z * z) / np.maximum(re, 1e-6)
    lab = np.zeros(x.shape, np.int16)
    brain = re <= 1
    lab[brain] = np.where(x[brain] < 0, 2, 41)
    s = rx / 70.0
    vent = ((np.abs(x) - 12 * s) / (6 * s)) ** 2 + (y / (20 * s)) ** 2 + ((z - 10 * s) / (8 * s)) ** 2 <= 1
    lab[vent & brain] = np.where(x[vent & brain] < 0, 4, 43)
    ctx = brain & (depth < 7)
    az = np.arctan2(y, np.abs(x))
    el = np.arctan2(z, np.hypot(x, y))
    idx = np.clip(((az + np.pi / 2) / np.pi * 17).astype(int), 0, 16) * 2 + (el >= 0)
    lab[ctx] = np.where(x[ctx] < 0, np.array(FS_L)[idx[ctx]], np.array(FS_R)[idx[ctx]])
    return lab, depth, brain, idx


def _sample(affine, shape, fn, slab=8):
    outs = None
    for k0 in range(0, shape[2], slab):
        k1 = min(k0 + slab, shape[2])
        ii = np.indices((shape[0], shape[1], k1 - k0)).reshape(3, -1).astype(float)
        ii[2] += k0
        res = fn(affine[:3, :3] @ ii + affine[:3, 3:])
        if outs is None:
            outs = [np.zeros(shape, r.dtype) for r in res]
        for o, r in zip(outs, res):
            o[:, :, k0:k1] = r.reshape(shape[0], shape[1], k1 - k0)
    return outs


def _save(arr, affine, path):
    img = nib.Nifti1Image(arr, affine)
    img.set_qform(affine, 1)
    img.set_sform(affine, 1)
    nib.save(img, str(path))
    return img


def make_subject(out, *, seed=0, t1_axes="+x,+y,+z", t1_zooms=(1.0, 1.0, 1.0), t1_shape=None,
                 qsm_zoom=1.0, qsm_cov_mm=144.0, scale=1.0, qsm_dtype="float32", qsm_4d=False,
                 qsm_zero_csf=False):
    """Write one synthetic subject into directory out; return a dict of its file paths."""
    from fsl.data.image import Image   # FSL's own conventions for the true FLIRT matrix

    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    radii = RADII * scale
    zooms = np.asarray(t1_zooms, float)
    perm = np.zeros((3, 3))
    for col, a in enumerate(t1_axes.split(",")):
        perm["xyz".index(a[1]), col] = 1.0 if a[0] == "+" else -1.0
    if t1_shape is None:
        world_extent = np.abs(perm.T) @ (2.6 * radii)
        t1_shape = tuple(int(np.ceil(e / z / 2) * 2) for e, z in zip(world_extent, zooms))
    at = np.eye(4)
    at[:3, :3] = _rotation([2, -3, 5]) @ perm @ np.diag(zooms)
    at[:3, 3] = -at[:3, :3] @ ((np.array(t1_shape) - 1) / 2.0)
    (seg,) = _sample(at, t1_shape, lambda w: (tissue(w, radii)[0],))
    neuro = np.linalg.det(at[:3, :3]) > 0
    _save((seg > 0).astype(np.float32) * 100 + 1, at, out / "T1.nii")
    _save(seg, at, out / "seg_T1.nii.gz")

    c0 = 1.0 if (neuro and t1_shape[0] % 2 == 0) else 0.0
    ab = at @ np.array([[2, 0, 0, c0], [0, 2, 0, 0], [0, 0, 2, 0], [0, 0, 0, 1.0]])
    brain = _save(np.zeros(tuple((n + 1) // 2 for n in t1_shape), np.float32), ab, out / "T1_brain.nii.gz")

    rq = _rotation([3, 0, 0])
    q_shape = (int(np.ceil(2.57 * radii[0] / qsm_zoom)), int(np.ceil(2.54 * radii[1] / qsm_zoom)),
               int(round(qsm_cov_mm * scale / qsm_zoom)))
    aq = np.eye(4)
    aq[:3, :3] = rq @ np.diag([-qsm_zoom, qsm_zoom, qsm_zoom])
    aq[:3, 3] = rq @ np.array([1.285 * radii[0], -1.27 * radii[1], -1.13 * radii[2]])

    def qfn(w):
        lab, depth, inb, idx = tissue(w, radii)
        inside = inb & (depth > 1.5)
        q = np.zeros(lab.shape, np.float32)
        ctx = (lab >= 1000) & inside
        q[ctx] = 10 + idx[ctx] * 0.5 + rng.normal(0, 4, ctx.sum())
        wm = np.isin(lab, (2, 41)) & inside
        q[wm] = -20 + rng.normal(0, 4, wm.sum())
        ve = np.isin(lab, (4, 43)) & inside
        q[ve] = 0.0 if qsm_zero_csf else 3 + rng.normal(0, 2, ve.sum())
        v = np.linalg.inv(at) @ np.vstack([w, np.ones(w.shape[1])])
        cc = np.rint(v[:3]).astype(int)
        ok = np.all((cc >= 0) & (cc < np.array(t1_shape)[:, None]), axis=0)
        truth = np.zeros(w.shape[1], np.int16)
        truth[ok] = seg[cc[0, ok], cc[1, ok], cc[2, ok]]
        return q, truth

    qsm, truth = _sample(aq, q_shape, qfn)
    if qsm_dtype != "float32":
        qsm = np.rint(qsm).astype(qsm_dtype)
    data = qsm[..., None] if qsm_4d else qsm
    q_img = _save(data, aq, out / "QSM.nii")
    _save(truth, aq, out / "truth_label_QSM.nii.gz")
    q2t = Image(brain).getAffine("world", "fsl") @ Image(q_img).getAffine("fsl", "world")
    np.savetxt(out / "Q2T.mat", q2t, fmt="%.10f")
    return {"dir": out, "t1": out / "T1.nii", "seg": out / "seg_T1.nii.gz",
            "t1_brain": out / "T1_brain.nii.gz", "qsm": out / "QSM.nii", "q2t": out / "Q2T.mat",
            "truth": out / "truth_label_QSM.nii.gz", "neurological": bool(neuro), "t1_shape": t1_shape}
