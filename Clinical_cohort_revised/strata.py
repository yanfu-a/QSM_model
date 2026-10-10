"""Image-protocol strata shared by select_pilot.py, P9 and the QC/sensitivity reports.

A scan's protocol stratum is its QSM voxel size and its axial (superior-inferior) coverage,
short or full. In this cohort the two are confounded: most iso1.0mm scans are short slabs.
"""

import numpy as np

# Coverage below this is a short slab (the QSM_HUASHAN worklist: 250 of 711 scans).
SHORT_COVERAGE_MM = 125.0


def voxel_size_label(zooms):
    """Voxel-size label, e.g. '0.8x0.8x0.8', from voxel sizes (numbers or an 'AxBxC' string)."""
    if isinstance(zooms, str):
        try:
            zooms = [float(z) for z in zooms.split("x")]
        except ValueError:
            return "unknown"
    zooms = list(zooms)[:3]
    if len(zooms) != 3 or not all(np.isfinite(zooms)):
        return "unknown"
    return "x".join(f"{float(z):.2g}" for z in zooms)


def si_extent_mm(img):
    """Field-of-view extent along the voxel axis closest to the world superior-inferior axis."""
    zooms = img.header.get_zooms()[:3]
    axis = int(np.argmax(np.abs(img.affine[2, :3])))
    return float(img.shape[axis] * zooms[axis])


def coverage_group(coverage_mm, cut=SHORT_COVERAGE_MM):
    """'short' or 'full' axial coverage, or 'unknown' when the coverage is not a number."""
    try:
        value = float(coverage_mm)
    except (TypeError, ValueError):
        return "unknown"
    if not np.isfinite(value):
        return "unknown"
    return "short" if value < cut else "full"
