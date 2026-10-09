#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""FLIRT coordinate helpers shared by P8 and check_grid_conversion.py.

FLIRT matrices map between "FSL coordinates": voxel indices scaled by pixdim,
with the first voxel axis flipped when the image's sform/qform has a positive
determinant (neurological storage). This matches fslpy's
Image.getAffine('voxel', 'fsl'). The original P8 omitted the flip, which is
harmless for radiological images but offsets neurological ones by one voxel
along the first axis when the subsampled grid is not symmetric about the
full grid (e.g. fslmaths -subsamp2 on an even-sized first axis).
"""

import numpy as np


def vox_to_fsl(img):
    """4x4 matrix from voxel indices to FLIRT scaled-voxel coordinates for a nibabel image."""
    zooms = [float(z) for z in img.header.get_zooms()[:3]]
    scale = np.diag(zooms + [1.0])
    # nibabel's affine follows FSL's precedence (sform, then qform, then pixdim scaling).
    if np.linalg.det(img.affine[:3, :3]) > 0:
        flip = np.eye(4)
        flip[0, 0] = -1.0
        flip[0, 3] = img.shape[0] - 1
        scale = scale @ flip
    return scale


def flirt_grid_conversion(src, dst):
    """FLIRT matrix mapping src FSL coordinates to dst FSL coordinates through world space."""
    return (vox_to_fsl(dst) @ np.linalg.inv(dst.affine) @ src.affine
            @ np.linalg.inv(vox_to_fsl(src)))


def legacy_grid_conversion(src, dst):
    """Grid conversion of the original P8 (no neurological x-flip); for audits only."""
    d_src = np.diag([float(z) for z in src.header.get_zooms()[:3]] + [1.0])
    d_dst = np.diag([float(z) for z in dst.header.get_zooms()[:3]] + [1.0])
    return d_dst @ np.linalg.inv(dst.affine) @ src.affine @ np.linalg.inv(d_src)


def is_neurological(img):
    """True when FSL treats the image as neurologically stored (positive determinant)."""
    return bool(np.linalg.det(img.affine[:3, :3]) > 0)
