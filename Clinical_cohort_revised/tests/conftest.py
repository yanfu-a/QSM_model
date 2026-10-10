"""Shared fixtures: stub FSL/FreeSurfer/P1 tools, P8/P9 runners and FreeSurfer label files.

FSL, FreeSurfer and real MRI data are not required. The stubs are:
* module (no-op) and fslinfo (header fields via nibabel);
* flirt -applyxfm -interp nearestneighbour, which maps reference voxels to input voxels
  through FLIRT coordinates using fslpy (FSL's own library), independently of native_grid;
* mri_synthseg, which copies the subject's prepared segmentation ($FAKE_SEG);
* a rigid-mode P1 with the real P1's skip rule (matrix kept and FORCE_RERUN != 1 -> exit 0),
  Dice gate (QSM_ROI_DICE_MIN; below it, exit 5 without keeping files) and QC fields.
Set QSM_TEST_REAL_FSL=1 to use the real flirt and fslinfo on PATH instead of the stubs.
"""

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

TESTS = Path(__file__).resolve().parent
REVISED = TESTS.parent
REPO = REVISED.parent
ORIGINAL_P8 = REPO / "Clinical_cohort" / "P8_native_aparc_qsm_ALL.sh"
sys.path.insert(0, str(TESTS))
sys.path.insert(0, str(REVISED))

NAMES = ["bankssts", "caudalanteriorcingulate", "caudalmiddlefrontal", "corpuscallosum", "cuneus",
         "entorhinal", "fusiform", "inferiorparietal", "inferiortemporal", "isthmuscingulate",
         "lateraloccipital", "lateralorbitofrontal", "lingual", "medialorbitofrontal", "middletemporal",
         "parahippocampal", "paracentral", "parsopercularis", "parsorbitalis", "parstriangularis",
         "pericalcarine", "postcentral", "posteriorcingulate", "precentral", "precuneus",
         "rostralanteriorcingulate", "rostralmiddlefrontal", "superiorfrontal", "superiorparietal",
         "superiortemporal", "supramarginal", "frontalpole", "temporalpole", "transversetemporal", "insula"]

STUB_FLIRT = '''#!{py}
import os, sys
import numpy as np, nibabel as nib
from fsl.data.image import Image
if os.environ.get("FAKE_FLIRT_FAIL") == "1":
    sys.exit(1)
a = sys.argv[1:]
opt = lambda k: a[a.index(k) + 1]
assert "-applyxfm" in a and opt("-interp") == "nearestneighbour", a
src, ref = nib.load(opt("-in")), nib.load(opt("-ref"))
m = np.loadtxt(opt("-init"))
vox2vox = Image(src).getAffine("fsl", "voxel") @ np.linalg.inv(m) @ Image(ref).getAffine("voxel", "fsl")
ii = np.indices(ref.shape[:3]).reshape(3, -1)
c = np.rint(vox2vox[:3, :3] @ ii + vox2vox[:3, 3:]).astype(np.int64)
data = np.asanyarray(src.dataobj)
ok = np.all((c >= 0) & (c < np.array(src.shape[:3])[:, None]), axis=0)
out = np.zeros(ii.shape[1], np.int16)
out[ok] = data[c[0, ok], c[1, ok], c[2, ok]]
hdr = ref.header.copy()
img = nib.Nifti1Image(out.reshape(ref.shape[:3]), ref.affine, hdr)
img.set_data_dtype(np.int16)
path = opt("-out")
nib.save(img, path if path.endswith((".nii", ".nii.gz")) else path + ".nii.gz")
'''

STUB_FSLINFO = '''#!{py}
import sys, nibabel as nib
h = nib.load(sys.argv[1]).header
print(f"datatype\\t{{int(h['datatype'])}}")
for i, z in enumerate(h.get_zooms()[:3], 1):
    print(f"pixdim{{i}}\\t{{z:.6f}}")
'''

STUB_SYNTHSEG = '''#!{py}
import os, shutil, sys
if os.environ.get("FAKE_SYNTHSEG_FAIL") == "1":
    sys.exit(1)
a = sys.argv[1:]
opt = lambda k: a[a.index(k) + 1]
shutil.copy(os.environ["FAKE_SEG"], opt("--o"))
if os.environ.get("FAKE_SYNTHSEG_NO_ICV") == "1":
    open(opt("--vol"), "w").write("subject,left cerebral cortex,right cerebral cortex\\nx,215000,218000\\n")
elif os.environ.get("FAKE_SYNTHSEG_NO_VOL") != "1":
    open(opt("--vol"), "w").write("subject,total intracranial,left cerebral cortex,right cerebral cortex\\n"
                                  "x,1500000,215000,218000\\n")
open(os.environ["STUB_LOG"], "a").write("synthseg " + opt("--i") + "\\n")
'''

FAKE_P1 = '''#!/bin/bash
set -euo pipefail
QSM="$1"; T1="$2"; ID="$3"
echo "P1 FORCE_RERUN=${FORCE_RERUN:-unset}" >> "$STUB_LOG"
if [[ -s "${QSM_ROI_KEEP_DIR}/${ID}_QSM_to_T1.mat" && "${FORCE_RERUN:-0}" != "1" ]]; then exit 0; fi
mkdir -p "$QSM_ROI_KEEP_DIR" "$QSM_ROI_OUT_DIR"
QC="${QSM_ROI_OUT_DIR}/${ID}_qc_rigid.txt"
if [[ "${FAKE_P1_CRASH:-0}" == "1" ]]; then exit 3; fi
DICE="${FAKE_P1_DICE:-0.850}"
# Like the real P1: below its Dice gate it writes a warning QC and exits before keeping files.
if awk -v d="$DICE" -v g="${QSM_ROI_DICE_MIN:-0.60}" 'BEGIN{exit !(d < g)}'; then
    printf 'ID=%s\\nQSM=%s\\nT1=%s\\nreg_mask_dice=%s\\nWARNING: reg_mask_dice %s < %s\\n' \\
        "$ID" "$QSM" "$T1" "$DICE" "$DICE" "${QSM_ROI_DICE_MIN:-0.60}" > "$QC"
    exit 5
fi
cp "$FAKE_Q2T" "${QSM_ROI_KEEP_DIR}/${ID}_QSM_to_T1.mat"
cp "$FAKE_T1BRAIN" "${QSM_ROI_KEEP_DIR}/${ID}_T1_brain.nii.gz"
if [[ "${QSM_ROI_T1_SUBSAMP:-2}" == "0" ]]; then PIX=naxnaxna; else PIX=2.000000x2.000000x2.000000; fi
{
    printf 'ID=%s\\nQSM=%s\\nT1=%s\\nstop_after=rigid\\nt1_pixdim=%s\\n' "$ID" "$QSM" "$T1" "$PIX"
    if [[ "${FAKE_NO_DICE:-0}" != "1" ]]; then printf 'reg_mask_dice=%s\\n' "$DICE"; fi
    printf 'search_mode=%s\\n' "${QSM_ROI_SEARCH_MODE:-nosearch}"
} > "$QC"
'''


@pytest.fixture(scope="session")
def tools(tmp_path_factory):
    """Directory with the stub tools; returns a dict of paths."""
    root = tmp_path_factory.mktemp("tools")
    bin_dir, fs_bin = root / "bin", root / "freesurfer" / "bin"
    bin_dir.mkdir()
    fs_bin.mkdir(parents=True)
    py = sys.executable
    (bin_dir / "module").write_text("#!/bin/bash\nexit 0\n")
    real_fsl = os.environ.get("QSM_TEST_REAL_FSL") == "1" and shutil.which("flirt")
    if not real_fsl:
        (bin_dir / "flirt").write_text(STUB_FLIRT.format(py=py))
        (bin_dir / "fslinfo").write_text(STUB_FSLINFO.format(py=py))
    (fs_bin / "mri_synthseg").write_text(STUB_SYNTHSEG.format(py=py))
    p1 = root / "fake_P1.sh"
    p1.write_text(FAKE_P1)
    for f in list(bin_dir.iterdir()) + [fs_bin / "mri_synthseg", p1]:
        f.chmod(0o755)
    original = None
    if ORIGINAL_P8.is_file():
        text = ORIGINAL_P8.read_text()
        text = re.sub(r"^PY310=.*$", f"PY310={py}", text, flags=re.M)
        text = re.sub(r"^export FREESURFER_HOME=.*$", f"export FREESURFER_HOME={root / 'freesurfer'}", text, flags=re.M)
        text = text.replace('bash "${HERE}/P1_qsm_to_mni.sh"', f'bash "{p1}"')
        original = root / "P8_original_stubbed.sh"
        original.write_text(text)
    return {"bin": bin_dir, "fs_home": root / "freesurfer", "p1": p1, "original_p8": original,
            "log": root / "stub.log", "real_fsl": bool(real_fsl)}


def run_p8(tools, subject, native, out_id, *, qsm=None, t1=None, original=False, env=None):
    """Run P8 (revised, or the original with stubbed paths) for one subject.

    Returns (exit code, combined output, list of stub calls made during the run).
    """
    tools["log"].write_text("")
    script = tools["original_p8"] if original else REVISED / "P8_native_aparc_qsm_ALL.sh"
    e = dict(os.environ)
    e.update({"PATH": f"{tools['bin']}:{e['PATH']}", "P8_PYTHON": sys.executable,
              "P8_FREESURFER_HOME": str(tools["fs_home"]), "P8_P1_SCRIPT": str(tools["p1"]),
              "STUB_LOG": str(tools["log"]), "P8_SYNTHSEG_THREADS": "1",
              "FAKE_SEG": str(subject["seg"]), "FAKE_Q2T": str(subject["q2t"]),
              "FAKE_T1BRAIN": str(subject["t1_brain"]), "QSM_ROI_NATIVE_ROOT": str(native)})
    for k in ("P8_FORCE_RERUN", "P8_FORCE_P1", "P8_FORCE_SYNTHSEG", "FORCE_RERUN"):
        e.pop(k, None)
    e.update(env or {})
    proc = subprocess.run(["bash", str(script), str(qsm or subject["qsm"]), str(t1 or subject["t1"]), out_id],
                          env=e, capture_output=True, text=True)
    calls = [line for line in tools["log"].read_text().splitlines() if line]
    return proc.returncode, proc.stdout + proc.stderr, calls


@pytest.fixture(scope="session")
def fs_fixtures(tmp_path_factory):
    """FreeSurferColorLUT.txt, SynthSeg --parc label list and Desikan CSVs (FreeSurfer order)."""
    root = tmp_path_factory.mktemp("fs")
    (root / "models").mkdir()
    lines = ["#No. Label Name: R G B A", "0 Unknown 0 0 0 0", "2 Left-Cerebral-White-Matter 245 245 245 0"]
    for base, hemi in ((1000, "lh"), (2000, "rh")):
        lines.append(f"{base} ctx-{hemi}-unknown 25 5 25 0")
        lines += [f"{base + i} ctx-{hemi}-{n} 25 100 40 0" for i, n in enumerate(NAMES, 1)]
    (root / "FreeSurferColorLUT.txt").write_text("\n".join(lines) + "\n")
    fs_l = [1001, 1002, 1003] + list(range(1005, 1036))
    np.save(root / "models" / "synthseg_parcellation_labels.npy",
            np.array([0, 2, 3, 4, 41, 42, 43] + fs_l + [k + 1000 for k in fs_l]))

    def csv_rows(order):
        rows = ["index,name", "1,L_white_matter"] + [f"{i + 1},L_{n}" for i, n in enumerate(order, 1)]
        rows += ["36,R_white_matter"] + [f"{i + 36},R_{n}" for i, n in enumerate(order, 1)]
        return "\n".join(rows) + "\n"
    (root / "Desikan_labels.csv").write_text(csv_rows(NAMES[:34]))
    (root / "Desikan_labels_alphabetical.csv").write_text(csv_rows(sorted(NAMES[:34])))
    return root


def run_p9(fs_fixtures, native, out_dir, worklist, *extra):
    """Run P9; returns (exit code, combined output)."""
    cmd = [sys.executable, str(REVISED / "P9_extract_native_qsm_ALL.py"),
           "--native-dir", str(native), "--out-dir", str(out_dir), "--worklist", str(worklist),
           "--labels", str(fs_fixtures / "Desikan_labels.csv"), "--freesurfer-home", str(fs_fixtures),
           "--jobs", "2", *map(str, extra)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    return proc.returncode, proc.stdout + proc.stderr


def clone_outputs(native, src_id, dst_id):
    """Copy one subject's P8 outputs to another ID (contents, and so recorded hashes, unchanged)."""
    for sub, pattern in (("label", "{}_aparc_aseg_QSMnative.nii.gz"), ("qc", "{}_native_qc.txt"),
                         ("t1seg", "{}_aparc+aseg_T1.nii.gz"), ("t1seg", "{}_t1seg_source.txt")):
        src = native / sub / pattern.format(src_id)
        if src.is_file():
            shutil.copy(src, native / sub / pattern.format(dst_id))


def write_worklist(path, rows):
    """rows: dicts with image_dir_id, qsm_file, t1_file and optional metadata."""
    import pandas as pd
    cols = ["hospital_id", "image_dir_id", "qsm_file", "t1_file", "qsm_coverage_mm", "model_label",
            "include_main_model", "category", "runnable"]
    full = [{**{"hospital_id": "", "qsm_coverage_mm": "144", "model_label": "Alzheimer's disease",
                "include_main_model": "yes", "category": "Neurodegenerative disease", "runnable": "yes"},
             **r} for r in rows]
    pd.DataFrame(full, columns=cols).to_csv(path, index=False)
    return path


@pytest.fixture(scope="session")
def subject(tmp_path_factory):
    """Realistic-size subject: neurological 1 mm T1 with even dimensions (the case the original
    P8 got wrong) and a radiological 1 mm QSM slab; passes P8 QC."""
    from synth import make_subject
    return make_subject(tmp_path_factory.mktemp("subject"), seed=3)
