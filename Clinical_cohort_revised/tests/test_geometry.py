"""P8 label geometry across T1/QSM orientations, voxel sizes, dimensions and coverage.

Brains are scaled to half size for speed, so P8's volume QC may fail (exit 15); the label is
written before QC, and only its agreement with the ground truth is tested here.
"""

import subprocess
import sys

import numpy as np
import nibabel as nib
import pandas as pd
import pytest

from conftest import REVISED, run_p8, write_worklist
from synth import make_subject

pytestmark = pytest.mark.slow

GEOMETRIES = {
    # name: (T1 voxel axes, T1 voxel size, QSM voxel size, QSM coverage mm, QSM data type, 4D QSM)
    "neuro_even_1mm": ("+x,+y,+z", (1, 1, 1), 0.8, 144, "float32", False),
    "neuro_odd_1mm": ("+x,+y,+z", (1, 1, 1), 1.0, 118, "float32", False),
    "radio_even_1mm": ("-x,+y,+z", (1, 1, 1), 0.8, 118, "float32", False),
    "neuro_aniso_1x05x05": ("+x,+y,+z", (1, 0.5, 0.5), 1.0, 118, "int16", False),
    "radio_aniso_1x05x05": ("-x,+y,+z", (1, 0.5, 0.5), 0.8, 144, "float32", True),
    "neuro_permuted_05x05x1": ("+y,+z,+x", (0.5, 0.5, 1), 0.8, 118, "float32", False),
    "radio_permuted_05x05x1": ("+y,-z,+x", (0.5, 0.5, 1), 1.0, 144, "float32", True),
}


def _agreement(label_path, truth_path):
    truth = np.asanyarray(nib.load(truth_path).dataobj)
    lab = np.asanyarray(nib.load(label_path).dataobj)
    cortex = truth >= 1000
    return float((lab[truth > 0] == truth[truth > 0]).mean()), float((lab[cortex] == truth[cortex]).mean())


@pytest.fixture(scope="module")
def runs(tools, tmp_path_factory):
    root = tmp_path_factory.mktemp("geometry")
    out = {}
    for i, (name, (axes, zooms, qz, cov, dtype, four_d)) in enumerate(GEOMETRIES.items()):
        shape = None
        if name == "neuro_odd_1mm":
            shape = (93, 110, 78)                        # odd first dimension
        subj = make_subject(root / name, seed=10 + i, t1_axes=axes, t1_zooms=zooms, t1_shape=shape,
                            qsm_zoom=qz, qsm_cov_mm=cov, scale=0.5, qsm_dtype=dtype, qsm_4d=four_d)
        codes = {}
        for variant in ("revised", "original"):
            if variant == "original" and tools["original_p8"] is None:
                continue
            code, log, _ = run_p8(tools, subj, root / f"native_{variant}", name, original=(variant == "original"))
            codes[variant] = (code, log)
        out[name] = (subj, codes)
    return root, out


@pytest.mark.parametrize("name", list(GEOMETRIES))
def test_revised_p8_label_matches_ground_truth(runs, name):
    root, out = runs
    subj, codes = out[name]
    code, log = codes["revised"]
    assert code in (0, 15), log
    brain, cortex = _agreement(root / "native_revised" / "label" / f"{name}_aparc_aseg_QSMnative.nii.gz",
                               subj["truth"])
    assert brain > 0.9995 and cortex > 0.999, (brain, cortex)


@pytest.mark.parametrize("name", list(GEOMETRIES))
def test_original_p8_error_only_for_neurological_even_first_dimension(runs, name):
    root, out = runs
    subj, codes = out[name]
    if "original" not in codes:
        pytest.skip("original P8 not found")
    _, cortex = _agreement(root / "native_original" / "label" / f"{name}_aparc_aseg_QSMnative.nii.gz",
                           subj["truth"])
    affected = subj["neurological"] and subj["t1_shape"][0] % 2 == 0
    assert (cortex < 0.99) == affected, (name, cortex)


def test_grid_audit_classifies_legacy_labels(runs, tmp_path):
    root, out = runs
    if "original" not in next(iter(out.values()))[1]:
        pytest.skip("original P8 not found")
    wl = write_worklist(tmp_path / "wl.csv", [
        {"image_dir_id": name, "qsm_file": str(s["qsm"]), "t1_file": str(s["t1"])} for name, (s, _) in out.items()])
    audit = tmp_path / "audit.csv"
    proc = subprocess.run([sys.executable, str(REVISED / "check_grid_conversion.py"), "--worklist", str(wl),
                           "--native-dir", str(root / "native_original"), "--out", str(audit)],
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    a = pd.read_csv(audit).set_index("IID")
    for name, (s, _) in out.items():
        affected = s["neurological"] and s["t1_shape"][0] % 2 == 0
        assert bool(a.loc[name, "label_offset"]) == affected
        assert a.loc[name, "saved_matrix_formula"] == ("legacy" if affected else "both (identical here)")
        assert pd.isna(a.loc[name, "p8_qc_grid_conversion"])             # legacy QC has no tag
