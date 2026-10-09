# Clinical_cohort_revised

Revised native-space cortical QSM pipeline (P8 labels, P9 extraction). The original
scripts in `../Clinical_cohort/` are unchanged.

| File | Status |
|---|---|
| `P8_native_aparc_qsm_ALL.sh` | revised |
| `P9_extract_native_qsm_ALL.py` | revised |
| `native_grid.py` | new; FLIRT-coordinate grid conversion shared by P8 and the audit |
| `make_support_masks.py` | new; writes the `<IID>_support.nii.gz` masks P9 needs |
| `select_pilot.py` | new; picks pilot scans that span the cohort's image geometries |
| `check_grid_conversion.py` | new; reports which existing P8 labels were built with the faulty conversion |
| `P1_qsm_to_mni_ALL.sh` | unchanged copy, so P8 finds P1 in this folder |

## P8 changes

1. **Grid conversion fix.** The T1 → T1_brain conversion now uses FLIRT's coordinate
   convention, which flips the first voxel axis for neurologically stored images (positive
   determinant). The original formula left the flip out. For such T1s with an even first
   dimension it produced the "−1 mm" shift the old comments called typical. That offset
   every label by one T1 voxel along the first axis. The correct conversion for
   `fslmaths -subsamp2` is the identity. The QC file records `grid_conversion=flirt_fsl_v2`,
   the T1 orientation, and the corrected and legacy shifts.
2. **An existing label is reused only when verified.** P8 skips a subject only if all of
   these hold:
   - the QC carries the grid-conversion tag and has no `ERROR` line;
   - the QSM and T1 paths match the QC;
   - the QSM, T1, P1 matrix, segmentation and label are byte-identical to the ones whose
     SHA-1 hashes the QC recorded;
   - the segmentation's source record verifies (item 4);
   - no rebuild was requested (`P8_FORCE_RERUN`, `P8_FORCE_P1`, `P8_FORCE_SYNTHSEG`, or an
     exported `FORCE_RERUN=1`).

   Otherwise it rebuilds and prints the reasons.
3. **P1 matrix reuse.**
   - **Matrices P8 registered:** P8 now writes a source record with the QSM and T1 hashes,
     and the matrix is reused only if those still match.
   - **Older matrices:** reused only if P1's QC names the same QSM and T1 paths and neither
     file is newer than the matrix.
   - Otherwise P1 is forced to register again. P1 reads `FORCE_RERUN`, not
     `P8_FORCE_RERUN`, and skips silently whenever a matrix exists.
4. **SynthSeg reuse.** A segmentation is reused only if its source record names a T1 with
   this content (SHA-1) and it lies on the T1 grid. A segmentation without a record is
   regenerated: matching dimensions and affines cannot show that two T1s are the same scan.
5. **No stale outputs.** The old QC and label are deleted before a rebuild, so a crash
   cannot leave a passing QC next to a new or missing label.
6. **Stricter QC.** A missing P1 Dice is a QC error. NaN QSM background is not coverage.
7. **P1 script lookup.** P8 runs `P1_qsm_to_mni.sh`, or `P1_qsm_to_mni_ALL.sh` if the first
   is absent.

Environment variables:
- `P8_PYTHON`, `P8_FREESURFER_HOME`, `P8_P1_SCRIPT`: path overrides.
- `P8_FORCE_RERUN=1`: rebuild the label; verified P1 and SynthSeg outputs are reused.
- `P8_FORCE_P1=1`: also re-register. An exported `FORCE_RERUN=1` does the same.
- `P8_FORCE_SYNTHSEG=1`: also rerun SynthSeg.

## P9 changes

- **Label names:** each CSV region name is checked against the FreeSurfer LUT name assigned
  to it by ID; mismatches stop the run unless `--allow-label-name-mismatch` is given. The
  mapping is written to `native_label_mapping*.csv`.
- **Duplicate worklist IDs:** an ID whose repeated rows disagree is excluded instead of
  taking the first row.
- **ROI values:** a value is reported only if the ROI has:
  - at least `--min-voxels` (20) valid voxels;
  - at least `--min-mm3` (20) of valid volume;
  - a finite, positive T1-space volume (from P8's segmentation), of which valid QSM covers
    at least `--min-coverage` (0.5).

  Coverage relative to T1 also counts cortex cut off by the QSM slab. If a label ROI has no
  volume in the T1 segmentation, the two disagree and extraction fails.
- **Analysis set** (`analysis_pass`, `fail_reasons` in `*_qc.csv`; rows of
  `*_analysis.csv`). A subject needs all of these:
  - extraction succeeded with at least `--min-reportable-rois` reported values (default 1);
  - it is not on the P1 failure list;
  - P8 QC passed with a numeric Dice and the grid-conversion tag;
  - the QSM is the file P8 used;
  - the label and T1 segmentation are byte-identical to the ones P8 QC recorded;
  - with `--csf-ref`, the CSF reference is finite.

  `--include-failed` subjects are extracted but never pass.
- **ROI completeness:** choose `--min-reportable-rois` to suit the downstream model.
  - Per-ROI (mass-univariate) models can use subjects with some ROIs missing: keep 1 and
    rely on each ROI column's own non-missing values.
  - Models needing every ROI per subject (multivariate, complete case): set
    `--min-reportable-rois 68`.

  `n_reportable_rois` is in `*_qc.csv` either way.
- **Proxy mode is kept apart.** Formal extraction requires `--support-dir`.
  `--allow-nonzero-proxy` (diagnostics/sensitivity only) adds `_proxy` to every output file
  name, so it cannot overwrite formal outputs. `*_qc.csv` records `coverage_mode`
  (`support_mask` or `nonzero_proxy`) and each subject's `support_mask` path.
- **Robustness and bookkeeping:**
  - Per-subject errors don't stop the run.
  - The CSF reference works with 4D QSM.
  - `*_excluded.csv` lists every dropped subject.
  - The QC CSV holds the P8 QC fields, protocol columns, `repeat_patient` and
    `possible_ppm_units`.
  - Reliability output adds a non-homologous baseline and a test–retest check.

All original output columns are kept. With the coverage and mm³ gates off, ROI medians
equal the original's and means differ by less than 1e-5 ppb.

## Support masks

`<IID>_support.nii[.gz]` marks the voxels where the QSM reconstruction produced values. It
separates genuine zero susceptibility from zero padding.

- **If your reconstruction pipeline exports its brain/ROI mask** (e.g. the mask used for
  background-field removal), save it under that name and use it. It must be on the QSM
  grid: P9 rejects a mask whose shape or affine differs from the QSM's.
- **The cohort's vendor QSM maps (`qsm_tra_*_QSM.nii`) come without one.** For these,
  `make_support_masks.py` recovers the support from each QSM image:
  - the nonzero, finite region;
  - plus zero-valued holes fully enclosed by it of at most `--max-hole-mm3` (default 10).

  Enclosed small zeros are genuine values, common when QSM is stored as quantized integers.
  Larger or open zero regions (masked CSF, outside the brain) stay unsupported.

`support_masks_summary.csv` lists, per subject, the QSM data type and how many voxels were
added. Check it on the pilot before relying on the masks.

## Running on the cluster

Set once (adjust paths):

```bash
cd /path/to/Clinical_cohort_revised          # P8, P1 and native_grid.py must be together
NATIVE=/cwStorage/nodecw_group/FY_data/QSM_HUASHAN/native_data
WL=$PWD/worklist.csv                         # copy of Clinical_cohort/worklist.csv
PY=/home1/fuyan/.conda/envs/py310/bin/python
```

**1. Fix the worklist.** These IDs have conflicting rows (two hospital_ids, some with two
diagnoses) and stay excluded by P9 until corrected: P191936, P194010, P195951, P196693,
P196837, P201485.

**2. Pilot on representative scans.** Pick up to 2 scans per geometry stratum (T1
orientation, T1 voxel size, T1 first-dimension parity, QSM voxel size, QSM coverage), then
run the whole chain on them:

```bash
$PY select_pilot.py --worklist $WL --out pilot_ids.csv --per-stratum 2

$PY - "$WL" pilot_ids.csv > pilot_jobs.tsv <<'EOF'
import sys
import pandas as pd
w = pd.read_csv(sys.argv[1], dtype=str).fillna("")
ids = set(pd.read_csv(sys.argv[2], dtype=str)["IID"])
w = w[w["image_dir_id"].isin(ids)].drop_duplicates("image_dir_id")
w[["qsm_file", "t1_file", "image_dir_id"]].to_csv(sys.stdout, sep="\t", index=False, header=False)
EOF
while IFS=$'\t' read -r qsm t1 iid; do
    QSM_ROI_NATIVE_ROOT=$NATIVE bash P8_native_aparc_qsm_ALL.sh "$qsm" "$t1" "$iid"
done < pilot_jobs.tsv

$PY make_support_masks.py --worklist $WL --ids-file pilot_ids.csv --out-dir $NATIVE/support
$PY check_grid_conversion.py --worklist $WL --native-dir $NATIVE --out output/native_grid_audit.csv
$PY P9_extract_native_qsm_ALL.py --worklist $WL --native-dir $NATIVE \
    --support-dir $NATIVE/support --failed-list failed_ID.txt \
    --ids-file pilot_ids.csv --suffix _pilot --out-dir output --jobs 8
```

Pilot checks, per stratum:
- **P8:** exit code 0 and no `ERROR:` lines in `$NATIVE/qc/<IID>_native_qc.txt`. Expect
  `grid_shift_mm=0.0000x0.0000x0.0000`; `legacy_grid_shift_mm` shows what the old code
  used.
- **Labels:** open each QSM in FSLeyes with `$NATIVE/label/<IID>_aparc_aseg_QSMnative.nii.gz`
  as a label overlay. Check the ribbon sits on cortex in both hemispheres, at the top and
  bottom of the slab.
- **Support masks:** `support_masks_summary.csv` should show small `filled_fraction` values.
  Integer-stored QSM will show more filled voxels than float.
- **P9:** check `output/HUASHAN_NATIVE_DK_pilot_qc.csv` (`analysis_pass`, `fail_reasons`,
  `n_reportable_rois`). In `output/native_roi_coverage_pilot.csv`, short-coverage strata
  should show lower `fov_fraction` for superior ROIs. Use this to settle `--min-coverage`
  and `--min-reportable-rois` before the full run.

**3. Full cohort.** The same commands without the ID filter:

```bash
$PY - "$WL" > p8_jobs.tsv <<'EOF'
import sys
import pandas as pd
w = pd.read_csv(sys.argv[1], dtype=str).fillna("")
w = w[w["runnable"].str.strip().str.lower() == "yes"].drop_duplicates("image_dir_id")
w[["qsm_file", "t1_file", "image_dir_id"]].to_csv(sys.stdout, sep="\t", index=False, header=False)
EOF
while IFS=$'\t' read -r qsm t1 iid; do        # or submit each line through your scheduler
    QSM_ROI_NATIVE_ROOT=$NATIVE bash P8_native_aparc_qsm_ALL.sh "$qsm" "$t1" "$iid"
done < p8_jobs.tsv

$PY make_support_masks.py --worklist $WL --out-dir $NATIVE/support --jobs 8
$PY P9_extract_native_qsm_ALL.py --worklist $WL --native-dir $NATIVE \
    --support-dir $NATIVE/support --failed-list failed_ID.txt --out-dir output --jobs 8
```

Add `--min-reportable-rois 68` for complete-case models and `--csf-ref` for CSF-referenced
values. A nonzero-proxy sensitivity run uses `--allow-nonzero-proxy` instead of
`--support-dir` and writes `*_proxy*` files beside the formal ones.

**First full P8 run after this update:** no existing output has the provenance records, so
every label is rebuilt once.
- SynthSeg reruns for every subject, because older segmentations have no T1 hash.
- P1 matrices are reused where P1's QC paths match and the inputs are older than the
  matrix; otherwise P1 re-registers.

Later runs skip verified subjects in seconds. P9 reads each subject's T1-space segmentation
(about 200 MB for a 176×512×512 T1), so peak memory is roughly `--jobs` × 0.4 GB.

## Testing and limits

FSL and FreeSurfer were not available where these changes were made. Testing used stub
tools:
- a fake `flirt -applyxfm` built on fslpy's FLIRT conventions;
- fake `fslinfo`, `mri_synthseg` and rigid-mode P1;
- synthetic brains, each with a ground-truth label on the QSM grid.

**Geometry matrix:** seven T1/QSM geometries were run through both versions of P8:

| Case | T1 storage | T1 shape, voxel (mm) | QSM | Revised P8, cortex | Original P8, cortex |
|---|---|---|---|---|---|
| A | neurological | 180×216×180, 1×1×1 | 0.8 mm, 144 mm slab | 99.9999% | 85.4% |
| B | neurological | 179×216×180, 1×1×1 | 1.0 mm, 118 mm slab | 100% | 100% |
| C | radiological | 180×216×180, 1×1×1 | 0.8 mm, 118 mm slab | 100% | 100% |
| D | neurological | 150×360×260, 1×0.5×0.5 | 1.0 mm, 118 mm slab | 99.9997% | 85.9% |
| E | radiological | 150×360×260, 1×0.5×0.5 | 0.8 mm, 144 mm slab | 99.9999% | 99.9999% |
| F | radiological, permuted axes | 360×260×150, 0.5×0.5×1 | 1.0 mm, 144 mm slab | 100% | 100% |
| G | neurological, permuted axes | 360×260×150, 0.5×0.5×1 | 0.8 mm, 118 mm slab | 99.9999% | 93.6% |

All seven passed P8 QC and P9's analysis set with 68 reportable ROIs. The minimum
`fov_fraction` was about 0.80 for 118 mm slabs and 0.98 for 144 mm slabs.

**Other checks:**
- Every P8 skip, reuse and provenance case: QC errors; changed QSM or T1 content, including
  a same-path replacement with an older timestamp; a changed P1 matrix, segmentation or
  label; a missing or mismatched segmentation source; each force flag.
- P9 cases: zero, negative, NaN or infinite T1 volume; all ROIs withheld; the complete-case
  threshold; label or segmentation altered after P8; formal and proxy runs in one folder.

**Assumption:** the synthetic T1_brain header follows `fslmaths -subsamp2` (centred
subsampling in FSL's radiological voxel order). The fix itself does not depend on it,
because it computes the conversion from the real headers.

**Before the full cohort,** run the pilot in step 2 on real scans. These tests use synthetic
images, not your data.

Unchanged on purpose: P1 registration, the `QSM_ROI_T1_SUBSAMP` default, the worklist,
and the choice of CSF reference region.
