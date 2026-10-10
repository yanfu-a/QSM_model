# Clinical_cohort_revised: native-space individualized cortical QSM

Each subject's 68 Desikan-Killiany cortical ROIs come from SynthSeg on that subject's T1.
They are mapped onto the subject's native QSM grid, and cortical QSM is summarised there.
No MNI registration or FNIRT is needed. The original scripts in `../Clinical_cohort/` are
unchanged.

```
P1 (rigid QSM->T1, run by P8) -> P8 (SynthSeg + label on the QSM grid + QC)
    -> make_support_masks.py (QSM support + provenance) -> P9 (ROI extraction, eligibility, QC)
```

| File | Role |
|---|---|
| `P1_qsm_to_mni_ALL.sh` | Unchanged copy. P8 runs it in rigid mode (QSM->T1_brain, 6 DOF) |
| `P8_native_aparc_qsm_ALL.sh` | Per subject: P1 rigid, SynthSeg `--parc`, T1->QSM label resampling, QC, provenance |
| `native_grid.py` | FLIRT-coordinate grid conversion (shared by P8 and the audit) |
| `make_support_masks.py` | Support masks `<IID>_support.nii.gz` with provenance records `<IID>_support.json` |
| `P9_extract_native_qsm_ALL.py` | Cohort ROI extraction, analysis eligibility, QC tables, run provenance |
| `select_pilot.py` | Picks representative pilot scans by image geometry |
| `strata.py` | Protocol strata (QSM voxel size, short/full coverage) shared by the scripts |
| `qc_exclusion_report.py` | Counts P8 QC exclusions per rule by protocol, coverage, diagnosis and sex (report only) |
| `compare_p9_runs.py` | Compares two P9 runs ROI by ROI (sensitivity analyses) |
| `check_grid_conversion.py` | Audits existing (historical) labels for the former grid-conversion error |
| `tests/` | Regression tests (pytest); see `TEST_REPORT.md` |
| `CODE_REVIEW.md`, `TEST_REPORT.md`, `METHODOLOGICAL_DECISIONS.md`, `CHANGELOG.md` | Review, test results, decisions awaiting approval, change log |

**Requirements:**
- FSL 6.0.4 and FreeSurfer 8.2.0 (SynthSeg), used by P8;
- Python with numpy, nibabel, pandas and scipy;
- `flock` (util-linux) and `sha1sum` (coreutils);
- for the tests, also pytest and fslpy.

## What makes a subject eligible

A subject is in the primary (unreferenced) analysis set (`analysis_pass` in `*_qc.csv`, rows
of `*_analysis.csv`) only if all of the following hold:
- **P8 QC:** the QC file is complete (`qc_complete=1`), has no `ERROR` line, records a
  numeric registration Dice, and carries the corrected grid-conversion tag. Labels from
  before the fix are never accepted.
- **Provenance:** the QSM and T1 that the worklist names, and the label and T1 segmentation
  that P9 reads, are byte-identical (SHA-1) to the ones P8 recorded.
- **Support mask:** its provenance record matches the current QSM and the mask file.
- **Extraction:** extraction succeeded with at least `--min-reportable-rois` ROI values.
- **P1 failure list:** the subject has no P1 failure in a stage the native branch shares
  (QSM input or rigid registration). Failures confined to T1->MNI/FNIRT do not exclude it.

Every excluded subject is listed with its reason in `*_excluded.csv`. `*_run_info.json`
records the command, options, software versions and input hashes, so the result can be
reproduced and audited.

## Setup (cluster)

```bash
cd /path/to/Clinical_cohort_revised            # P8, P1 and native_grid.py must stay together
module load fsl/6.0.4
NATIVE=/cwStorage/nodecw_group/FY_data/QSM_HUASHAN/native_data
P1QC=/cwStorage/nodecw_group/FY_data/QSM_HUASHAN/roi_work   # P1's MNI-mode QC files (<IID>_qc.txt)
WL=$PWD/worklist.csv                           # copy of Clinical_cohort/worklist.csv
PY=/home1/fuyan/.conda/envs/py310/bin/python
```

**Fix the worklist first.** Six image_dir_ids have conflicting rows: two hospital_ids, and
some have two diagnoses. They stay excluded until corrected: P191936, P194010, P195951,
P196693, P196837, P201485.

## 1. Representative pilot

Run the whole chain on a few scans from each image-geometry stratum: T1 storage orientation,
T1 voxel size, T1 first-dimension parity, QSM voxel size, and QSM coverage (short or full).

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

# P8. P8_FORCE_P1=1 re-registers so every matrix gets content provenance (decision D7).
while IFS=$'\t' read -r qsm t1 iid; do
    QSM_ROI_NATIVE_ROOT=$NATIVE P8_FORCE_P1=1 bash P8_native_aparc_qsm_ALL.sh "$qsm" "$t1" "$iid"
    echo "$iid exit=$?"
done < pilot_jobs.tsv

# Support masks: register reconstruction masks if they exist (preferred) ...
#   $PY make_support_masks.py --worklist $WL --ids-file pilot_ids.csv --out-dir $NATIVE/support \
#       --recon-mask-template '/path/to/{IID}_qsm_mask.nii.gz'
# ... otherwise estimate them from the QSM (see "Support masks" below):
$PY make_support_masks.py --worklist $WL --ids-file pilot_ids.csv --out-dir $NATIVE/support

$PY check_grid_conversion.py --worklist $WL --native-dir $NATIVE --out output/native_grid_audit_pilot.csv

$PY P9_extract_native_qsm_ALL.py --worklist $WL --native-dir $NATIVE --support-dir $NATIVE/support \
    --failed-list failed_ID.txt --p1-qc-dir $P1QC --ids-file pilot_ids.csv --suffix _pilot \
    --out-dir output --jobs 8

# QC exclusions per rule and stratum; add --covariates demographics.csv (IID or hospital_id, sex)
$PY qc_exclusion_report.py --worklist $WL --native-dir $NATIVE --ids-file pilot_ids.csv \
    --out-dir output/qc_report_pilot

# Optional, on the cluster: run the real-FSL checks
python -m pytest tests -m fsl -v
```

**Pilot checklist (per stratum):**
- **P8 exit code and QC file:** exit 0, and `$NATIVE/qc/<IID>_native_qc.txt` has no
  `ERROR:` lines. Check:
  - `grid_shift_mm=0.0000x0.0000x0.0000`;
  - `q2t_scale_dev` near 0 (a rigid registration);
  - `qsm_pixdim_sform_rel_diff` and `t1_pixdim_sform_rel_diff` near 0.
- **Visual QC:** open the QSM in FSLeyes with `$NATIVE/label/<IID>_aparc_aseg_QSMnative.nii.gz`
  as a label overlay. The ribbon should lie on cortex in both hemispheres, at the top and
  bottom of the slab.
- **Support masks** (`$NATIVE/support/support_masks_summary.csv`): note `support_type`,
  `qsm_dtype` and `qsm_quantized`, and check that `filled_voxels` is small. Pay particular
  attention to non-quantized maps, for which `make_support_masks.py` warns.
- **P9 QC** (`output/HUASHAN_NATIVE_DK_pilot_qc.csv`): review `analysis_pass`,
  `fail_reasons`, `n_reportable_rois`, `P1_failure_stage`, `P8_cortex_vol_ml` vs
  `t1_cortex_ml`, and `max_fov_fraction`.
- **Coverage** (`output/native_roi_coverage_pilot.csv`): `coverage_vs_t1` and
  `fov_fraction` by ROI and stratum. Short slabs lose superior and inferior ROIs. Also look
  at `valid_zero_voxels` and `unsupported_zero_voxels`, the zeros counted as values and
  those left out as padding.
- **Missingness** (`output/native_roi_missingness_pilot.csv`): missing fraction per ROI and
  protocol stratum.
- **QC rules** (`output/qc_report_pilot/`): exclusions per rule and stratum, and the
  QSM-grid versus T1 cortex-volume rule (decision D1).
- **Decisions:** settle D1-D10 in `METHODOLOGICAL_DECISIONS.md` before the full run.

## 2. Native-space extraction (after the pilot and decisions)

```bash
$PY - "$WL" > p8_jobs.tsv <<'EOF'
import sys
import pandas as pd
w = pd.read_csv(sys.argv[1], dtype=str).fillna("")
w = w[w["runnable"].str.strip().str.lower() == "yes"].drop_duplicates("image_dir_id")
w[["qsm_file", "t1_file", "image_dir_id"]].to_csv(sys.stdout, sep="\t", index=False, header=False)
EOF
while IFS=$'\t' read -r qsm t1 iid; do          # or submit each line through your scheduler
    QSM_ROI_NATIVE_ROOT=$NATIVE bash P8_native_aparc_qsm_ALL.sh "$qsm" "$t1" "$iid"
done < p8_jobs.tsv

$PY make_support_masks.py --worklist $WL --out-dir $NATIVE/support --jobs 8
$PY P9_extract_native_qsm_ALL.py --worklist $WL --native-dir $NATIVE --support-dir $NATIVE/support \
    --failed-list failed_ID.txt --p1-qc-dir $P1QC --out-dir output --jobs 8
```

**Options:**
- `--main-model-only`: only `include_main_model=yes` rows (D11). Without it, every runnable
  row is processed and `include_main_model` is a column of every table.
- `--min-reportable-rois 68`: complete cases only (D2).
- `--csf-ref`: adds CSF-referenced tables. They have their own `csfref_analysis_pass` and
  never change `analysis_pass`.
- `--allow-nonzero-proxy` (instead of `--support-dir`): a diagnostic or sensitivity run.
  Every output name carries `_proxy`, so formal outputs are never overwritten.

**Output files:**
- `HUASHAN_NATIVE_DK.csv` (all extracted subjects) and `HUASHAN_NATIVE_DK_analysis.csv`
  (analysis set);
- `*_qc.csv` (per-subject QC, provenance and eligibility);
- `*_excluded.csv`;
- `*_run_info.json`;
- `native_roi_coverage.csv` (per subject and ROI);
- `native_roi_missingness.csv` (analysis set, per ROI and protocol stratum);
- `native_label_mapping.csv`;
- `native_mask_voxel_counts.csv`;
- with `--csf-ref`: `*_csfref.csv` and `*_csfref_analysis.csv`.

**First full run.** Labels without the current QC fields (all labels made before this
revision) are rebuilt once. SynthSeg reruns where the segmentation's source record cannot
be verified. Afterwards, verified subjects are skipped in seconds.

Peak P9 memory is roughly `--jobs` x 0.4 GB, because each subject's T1-space segmentation is
read for coverage.

## Sensitivity analyses

Each variant goes into its own directory and suffix, so primary outputs are never
overwritten. `compare_p9_runs.py` reports, per ROI and protocol stratum:
- values gained or lost;
- the mean and SD of the difference B - A;
- the mean and maximum absolute difference, and Pearson r.

```bash
# Support-mask hole limit (D3): 0 = no hole filling
$PY make_support_masks.py --worklist $WL --out-dir $NATIVE/support_h0 --max-hole-mm3 0
$PY P9_extract_native_qsm_ALL.py --worklist $WL --native-dir $NATIVE --support-dir $NATIVE/support_h0 \
    --failed-list failed_ID.txt --p1-qc-dir $P1QC --suffix _h0 --out-dir output --jobs 8
$PY compare_p9_runs.py --a output/HUASHAN_NATIVE_DK_analysis.csv --b output/HUASHAN_NATIVE_DK_h0_analysis.csv \
    --qc output/HUASHAN_NATIVE_DK_qc.csv --out output/compare_hole0.csv

# Registration at full T1 resolution (D12), in a second native root
NATIVE_FULL=/cwStorage/nodecw_group/FY_data/QSM_HUASHAN/native_data_t1full
while IFS=$'\t' read -r qsm t1 iid; do
    QSM_ROI_NATIVE_ROOT=$NATIVE_FULL QSM_ROI_T1_SUBSAMP=0 bash P8_native_aparc_qsm_ALL.sh "$qsm" "$t1" "$iid"
done < pilot_jobs.tsv
$PY make_support_masks.py --worklist $WL --ids-file pilot_ids.csv --out-dir $NATIVE_FULL/support
$PY P9_extract_native_qsm_ALL.py --worklist $WL --native-dir $NATIVE_FULL --support-dir $NATIVE_FULL/support \
    --failed-list failed_ID.txt --p1-qc-dir $P1QC --ids-file pilot_ids.csv --suffix _pilot_t1full \
    --out-dir output --jobs 8
$PY compare_p9_runs.py --a output/HUASHAN_NATIVE_DK_pilot_analysis.csv \
    --b output/HUASHAN_NATIVE_DK_pilot_t1full_analysis.csv --qc output/HUASHAN_NATIVE_DK_pilot_qc.csv \
    --out output/compare_t1full_pilot.csv
```

P8 records the registration's T1 resolution and search mode (`p1_t1_resolution`,
`p1_search_mode`). It never reuses a matrix or label made with other settings, even in the
same native root.

## Support masks

A support mask marks the voxels where the QSM reconstruction produced values. It separates
genuine zero susceptibility from padding. Each mask's record states its `support_type`:
- **`reconstruction`:** the reconstruction pipeline's own brain or ROI mask, registered with
  `--recon-mask-template`. It must lie on the QSM grid. This is the correct definition, and
  the preferred one.
- **`estimated_from_qsm`:** for QSM exported without its mask, as the cohort's vendor QSM
  appears to be. It is the nonzero region plus enclosed zero holes of at most
  `--max-hole-mm3` (default 10). It is an estimate:
  - for float QSM it equals the nonzero region;
  - for integer-quantized QSM it misses genuine zeros at the brain edge (measured in
    METHODOLOGICAL_DECISIONS D3).

A mask is rebuilt automatically when its source QSM, parameters or file change. P9 rejects
masks whose record does not verify.

## P8 reference

**Environment variables:**
- `QSM_ROI_NATIVE_ROOT`: output root.
- `P8_FORCE_RERUN=1`: rebuild the label, reusing verified P1 and SynthSeg outputs.
- `P8_FORCE_P1=1`, or an exported `FORCE_RERUN=1`: re-register.
- `P8_FORCE_SYNTHSEG=1`: rerun SynthSeg.
- `P8_PYTHON`, `P8_FREESURFER_HOME`, `P8_P1_SCRIPT`: path overrides.
- `P8_SYNTHSEG_THREADS`: SynthSeg thread count.
- `QSM_ROI_T1_SUBSAMP=0`: register to the full-resolution T1 (default: downsampled).
- `QSM_ROI_SEARCH_MODE`: P1 search mode (default `nosearch`).

These two settings are recorded and checked before any reuse. `QSM_ROI_SEARCH_DEG` is not
recorded, so use `P8_FORCE_P1=1` after changing it. P1's own Dice gate (`QSM_ROI_DICE_MIN`)
is disabled inside P8, which applies its 0.60 gate in the QC step instead.

**QC record:** `qc_failed_rules` names the failed rules, and `qc_report_only_failed_rules`
names the failures of rules that are evaluated but not applied. See
`qc_exclusion_report.py`.

**Exit codes:**

| Code | Meaning |
|---|---|
| 0 | Done, or verified existing output skipped |
| 1 | Setup error |
| 2 | Missing input |
| 4 | Non-scalar QSM |
| 5 | P1 rigid Dice below its gate (no longer raised: P8 disables that gate, so such a scan exits 15) |
| 11 | P1 outputs missing |
| 12 | Grid-conversion or label-grid error |
| 13 | SynthSeg failure |
| 14 | FLIRT failure or unexpected labels |
| 15 | QC errors (label and QC are written) |
| 16 | Another P8 run holds the subject's lock |

## Tests

```bash
python -m pytest tests -v                      # about 3 min; FSL/FreeSurfer not needed (stubs)
python -m pytest tests -m "not slow" -v        # unit tests only
python -m pytest tests -m fsl -v               # with real FSL (cluster)
QSM_TEST_REAL_FSL=1 python -m pytest tests -v  # all P8 tests with the real flirt/fslinfo
```

The tests use synthetic images only; see `TEST_REPORT.md` for results and remaining
validation gaps. They do not establish clinical validity.
