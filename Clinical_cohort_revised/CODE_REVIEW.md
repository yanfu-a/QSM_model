# Code review: native-space individualized cortical QSM pipeline

**Scope.** The native branch only: P1 rigid registration (`P1_qsm_to_mni_ALL.sh`, rigid
mode), P8 segmentation and native-space label (`P8_native_aparc_qsm_ALL.sh`, with
`native_grid.py`), support-mask generation (`make_support_masks.py`), and P9 extraction and
QC (`P9_extract_native_qsm_ALL.py`). It also covers the helpers `check_grid_conversion.py`
and `select_pilot.py` and their compatibility with `worklist.csv`.

**Method.**
- Every file was read in full. Every suspected issue was reproduced on synthetic data
  before being fixed.
- Coordinate conventions were checked against fslpy (FSL's own library).
- ROI labels were checked against an independent copy of `FreeSurferColorLUT.txt` (the
  copy shipped with MNE-Python 1.13.2).
- Line numbers refer to the audited commit `915e78b` unless stated otherwise. Locations in
  the fixed code are given by function name.

**Classification.** Each finding is one of: confirmed bug, potential risk, methodological
decision (needs your approval; see `METHODOLOGICAL_DECISIONS.md`), or non-issue.

**Severity.**
- **High:** can change QSM values, which subjects are analysed, or the validity of a
  result.
- **Medium:** can do so under realistic but less common conditions.
- **Low:** auditability or robustness only.

## Summary

| ID | Finding | Class | Severity | Status |
|---|---|---|---|---|
| CR-01 | P9 accepted a replaced QSM when its path was unchanged | Confirmed bug | High | Fixed |
| CR-02 | Support masks reused without checking source QSM or parameters | Confirmed bug | High | Fixed |
| CR-03 | Worklist `t1_file` never compared with the T1 P8 segmented; `t1_file` ignored in duplicate conflicts | Confirmed bug | High | Fixed |
| CR-04 | `--allow-legacy-labels` admitted labels possibly built with the faulty grid conversion | Confirmed bug | High | Fixed (option removed) |
| CR-05 | Any P1 failure-list entry excluded the subject, including MNI/FNIRT-only failures | Confirmed bug | High (sample size) | Fixed; unknown-stage rule is decision D5 |
| CR-06 | Registration conventions (`native_grid.py`, P8 matrix composition, P1 `-usesqform -dof 6`) | Verified; residual risk | n/a | Verified on synthetic data; real-FSL tests added |
| CR-07 | P8 QC thresholds can select by protocol, sex/head size and disease severity | Methodological | High (bias) | Decision D1; report columns added |
| CR-08 | ROI coverage measures: `fov_fraction` and `coverage_vs_t1` can exceed 1; QSM-grid coverage blind to slab truncation | Methodological + reporting gap | Medium | Decision D2; `max_fov_fraction` added |
| CR-09 | Estimated support masks are not reconstruction masks; they lose edge zeros in quantized QSM | Methodological | Medium | Decision D3; `support_type` recorded |
| CR-10 | QSM units, scaling, integer data, reference offsets | Non-issue (no rescaling) + reporting gap | Low | dtype/scaling/integer columns added; decision D4 |
| CR-11 | `--csf-ref` changed the primary (unreferenced) analysis set; L/R concordance labelled "reliability" | Confirmed bug | High / Low | Fixed |
| CR-12 | Non-atomic QC and outputs; reuse records without artifact hashes; no protection against concurrent runs | Confirmed bugs / potential risk | Medium | Fixed; legacy-matrix rule is decision D7 |
| CR-13 | ROI label IDs, LUT names, 68 unique names | Non-issue | n/a | Verified against independent LUT |
| CR-14 | Exclusion reasons incomplete or misleading; no run provenance | Confirmed bug | Low | Fixed |
| CR-15 | English-only text; compatibility of P1/P8/P9/worklist | Verified | n/a | One schema change (`t1_file` required) |
| CR-16 | Header pixdim vs sform/qform inconsistency not detected (e.g. 00413156, reported 1.83 mm slices) | Potential risk | Medium | Reported, not gating (decision D10) |
| CR-17 | Both scans of a re-scanned patient enter the analysis set | Methodological | Medium | Flagged (`repeat_patient`); decision D6 |

## Priority 0: bugs that could affect scientific results

### CR-01 QSM provenance accepted an unchanged path (confirmed bug, High; fixed)

**Where:** `P9 915e78b:265-276` (`same_file()`) and `:731-732`. `same_file()` returned True
as soon as the two paths resolved to the same location. A QSM re-exported or re-selected
into the same path therefore passed as the file P8 had used.

**Impact:** ROI values and their P8 QC verdict (registration Dice, GM-WM contrast, coverage)
could refer to different images.

**Reproduced:** subject QREP (QSM replaced at the same path after P8) was
`analysis_pass=True` under `915e78b` (`TEST_REPORT.md`, section 4).

**Fix:**
- P9 hashes the QSM it reads (`process_subject`) and compares it with P8's `qsm_sha1`.
  A mismatch, or a QC without `qsm_sha1`, prevents eligibility.
- `same_file()` is replaced by `same_path()`, which is reported only as
  `qsm_path_matches_P8`.
- `qsm_matches_P8` now means "content matches".

### CR-02 Stale support masks (confirmed bug, High; fixed)

**Where:**
- `make_support_masks.py 915e78b:56-57` reused any existing mask
  ("exists (use --overwrite to rebuild)").
- P9 (`:653-658`) checked only that a mask existed and lay on the QSM grid.

**Impact:** a mask computed from a superseded QSM (or with different parameters) would
define valid voxels for the new QSM. In the demonstration, an 8,000-voxel mask from the old
QSM stayed in use when the current QSM had 4,000 nonzero voxels.

**Fix:**
- Each mask now has a provenance record `<IID>_support.json` containing:
  - `support_type` (`reconstruction` or `estimated_from_qsm`), `method` and `max_hole_mm3`;
  - the SHA-1 of the source QSM (and of the reconstruction mask, if one was used);
  - the SHA-1 of the written mask.
- A mask is reused only if the record matches the current inputs and parameters and the
  mask file is unchanged. Otherwise it is rebuilt automatically.
- Writes are atomic, and the old record is deleted before a new mask is written.
- P9 verifies the record against the QSM it reads and the mask file. A missing, unreadable
  or mismatched record prevents eligibility.
- Genuine reconstruction masks can be registered with `--recon-mask-template`. They are
  grid-checked against the QSM, and nonzero QSM outside the mask is reported.

### CR-03 T1 provenance (confirmed bug, High; fixed)

**Where:** `P9 915e78b:593` required only `image_dir_id, qsm_file` and never read
`t1_file`. Duplicate-conflict detection (`:242`) compared `qsm_file` and metadata but not
`t1_file`.

**Impact:** the segmentation and labels could come from a different T1 than the one the
worklist now names, and nothing would report it. Duplicate rows naming different T1s were
silently collapsed to the first row.

**Fix:**
- `t1_file` is a required worklist column.
- P9 hashes the worklist T1 and compares it with P8's `t1_sha1`, reported as
  `t1_matches_P8` and `t1_path_matches_P8`.
- A missing T1 file is an explicit exclusion.
- `t1_file` is part of duplicate-conflict detection.

All six currently duplicated IDs in `Clinical_cohort/worklist.csv` name the same T1 in both
rows; their hospital_id and diagnosis conflicts were already excluded.

### CR-04 Legacy labels (confirmed bug, High; fixed)

**Where:** `P9 915e78b:543-545, 756, 763, 769, 870-871`. `--allow-legacy-labels` accepted
any label without the grid-conversion tag and waived the label and segmentation provenance
checks. The help text asked the user to restrict it to "unaffected" subjects, but nothing
enforced that.

**Reproduced:** LEG_AFF, a label built by the original P8 from a neurologically stored T1
with an even first dimension (offset by one voxel), was `analysis_pass=True`.

**Fix:**
- The option is removed, and legacy labels are never eligible ("P8 label predates
  grid-conversion fix (rerun P8)").
- They also fail the content-provenance checks, which legacy QC files cannot satisfy.
- `check_grid_conversion.py` remains for assessing historical outputs. It reports, per
  subject, the legacy and corrected conversions, which formula built the saved matrix,
  whether the label is offset, and the P8 QC tag.
- The revised P8 rebuilds legacy labels automatically.

### CR-05 Stage-blind use of the P1 failure list (confirmed bug, High for sample size; fixed)

**Where:** `P9 915e78b:600-608, 748`. Every ID in `failed_ID.txt` was excluded. The list is
filled by P1's MNI-branch batch, as `P1:296` notes ("Timeouts/failures exit with code 6,
are written to QC, and are collected in failed_ID.txt"). P1's exit codes fall into two
groups:

| Stage | Exit codes |
|---|---|
| Shared with the native branch | 2 (missing input), 4 (non-scalar QSM), 5 (QSM->T1 rigid Dice < 0.60; `P1:240, 401`) |
| MNI branch only | 3 (output not 91x109x91; `P1:394`), 6 (FNIRT failure or T1->MNI Dice; `P1:326, 406`), 7 (MNI coverage; `P1:419`) |

**Impact:** subjects whose native-space analysis is valid were excluded because FNIRT
failed. FNIRT failures are more likely with atypical anatomy, such as atrophy, so this bias
can also be disease-related.

**Fix:**
- P9 classifies each failure-list entry from P1's QC file (`--p1-qc-dir`, default P1's
  `roi_work`) into `native`, `mni_only` or `unknown`.
- `mni_only` entries stay in the native analysis. `native` and `unknown` entries are
  excluded, with the stage and evidence in `*_excluded.csv`.
- `P1_failure_stage` and `P1_failure_evidence` are reported.
- The native branch also re-checks rigid registration independently: P8 runs P1 in rigid
  mode and enforces Dice >= 0.60.

**Unverified:** the real `failed_ID.txt` is not in the repository, so the classification
of the actual entries could not be checked. Excluding `unknown` entries is the
conservative choice, but it is a decision (D5).

### CR-06 Registration correctness (verified on synthetic data; residual risk)

`native_grid.py`:
- `vox_to_fsl` and `flirt_grid_conversion` match fslpy to within 1e-4 mm (float32 header
  precision) on 300 random oblique grid pairs.
- The pairs cover both storage orientations and sform+qform, sform-only and qform-only
  headers.
- For grids that `fslmaths -subsamp2` would produce, the conversion is the identity for
  neurological and radiological storage, odd and even dimensions, anisotropic and permuted
  axes. The legacy formula is off by exactly one voxel only for neurological images with an
  even first dimension.

**P8 composition:** `M_T1->QSM = inv(M_QSM->T1_brain) @ M_T1->T1_brain` (P8 step 3) is
correct. With FLIRT resampling emulated through fslpy, labels reproduced the ground truth in
seven T1/QSM geometries:
- neurological, radiological and permuted storage;
- odd and even first dimensions;
- 1 mm and 1.0x0.5x0.5 mm T1;
- 0.8 and 1.0 mm QSM;
- 118 and 144 mm slabs;
- int16 and singleton-4D QSM.

Brain agreement was >= 99.95% and cortex >= 99.9%; the remainder is rounding ties. The
original P8 failed exactly in the three neurological, even-first-dimension cases.

**P1 `flirt -usesqform -dof 6 -cost corratio -nosearch` (`P1:197-203`):**
- A rigid model is appropriate for QSM and T1 from the same session.
- `-usesqform` starts from the scanner world alignment, which is rigid in FSL coordinates
  when header pixdims agree with the sform/qform scaling.
- The matrix is consumed in FLIRT convention by FLIRT itself, so the conventions are
  consistent.
- This cannot be executed here. `tests/test_fsl_real.py` therefore checks, with real FSL:
  - whether `fslmaths -subsamp2` produces the assumed headers;
  - that real `flirt -applyxfm` with P8's matrix reproduces ground truth;
  - that P1's exact call recovers a known rigid misalignment with a rigid matrix.
- P8 now also records `q2t_scale_dev` (deviation of the matrix's singular values from 1)
  and warns above 1%.

**Residual risk:** the `-subsamp2` header model is unconfirmed until those tests run on the
cluster. The fix does not depend on it, because it computes the conversion from the real
headers; the model only affects P8's informational warning.

## Priority 1: methodological issues (no thresholds changed)

### CR-07 P8 QC thresholds (methodological, High risk of bias; decision D1)

Each threshold was evaluated for whether it measures technical failure or biological
variation:
- **Cortical volume 350-700 mL (`P8 915e78b:398, 464`):**
  - It is measured on the label resampled to the QSM grid (`cortex_ml = ribbon voxels x
    QSM voxel volume`), so it excludes cortex outside the QSM slab.
  - 250 of 711 worklist scans have less than 125 mm of coverage, almost all from the
    iso1.0mm protocol. Those scans can fall below 350 mL through truncation alone.
  - Severe atrophy (AD, FTD) also lowers cortical volume.
  - The check therefore mixes protocol, field of view and disease severity.
- **SynthSeg ICV 1200-1900 mL (`:462`):** ICV is largely fixed in adulthood and smaller in
  women and smaller adults, so a 1200 mL floor can exclude valid scans, preferentially from
  women.
- **Registration Dice >= 0.60 (`:456`; also P1's gate):** this is the overlap between the
  QSM nonzero mask and the BET mask of the whole T1. Short slabs and eroded reconstruction
  masks lower Dice without any misregistration.
- **Cortical ribbon coverage >= 0.50 (`:466`):** this is computed within the QSM field of
  view, so it is blind to slab truncation. It mainly reflects mask erosion at the cortex.
- **GM-WM contrast > 0 ppb (`:468`):** unit-free and lenient. Gross failures fail it, but
  millimetre-scale misregistration usually still passes.

**Changes made (reporting only):** P9 now reports `t1_cortex_ml` (full-brain cortical
volume from the T1 segmentation, not truncated by the field of view) next to
`P8_cortex_vol_ml` and `P8_synthseg_icv_ml`, so the thresholds can be re-evaluated against
protocol, sex and diagnosis on the pilot.

### CR-08 ROI coverage (methodological + reporting gap; decision D2)

- `coverage_fraction` (valid / labelled voxels on the QSM grid) is <= 1 but blind to slab
  truncation.
- `coverage_vs_t1` (valid mm3 / T1-space ROI mm3) and `fov_fraction` (QSM-grid label mm3 /
  T1 mm3) include truncation. They can exceed 1 by sampling error, because nearest-neighbour
  resampling of 1 mm labels onto 0.8 mm voxels changes boundary voxel counts. The maximum
  observed on synthetic data was 1.016.
- Larger excesses indicate a label/segmentation mismatch. Zero or non-finite T1 volume is
  now a hard failure (fixed in the previous revision), and P9 now reports
  `max_fov_fraction` per subject.
- `--min-voxels 20` is resolution-dependent. `--min-mm3 20` harmonises the physical
  threshold, so at 0.8 mm (0.512 mm3 voxels) the volume gate dominates.

No threshold was changed.

### CR-09 Support masks (methodological; decision D3)

A mask estimated from QSM intensity can only separate "zero" from "nonzero". Measured
behaviour:
- **Float QSM:** exact interior zeros essentially never occur, so the estimate equals the
  nonzero proxy.
- **Integer-quantized QSM:** many genuine zeros occur. Filling enclosed holes of up to
  10 mm3 recovers them in the interior, but zeros touching the brain edge connect to the
  background. With 1 ppb steps and SD 3 ppb, 5.6% of the outermost 3 mm (where the cortex
  lies) stayed unsupported, versus 0.45% deeper.

Sensitivity to the hole limit is in `METHODOLOGICAL_DECISIONS.md`. Records now state
`support_type`, and P9 reports `support_type`, `support_method`, `support_max_hole_mm3` and
`qsm_nonzero_outside_support`.

### CR-10 Units, scaling and integer data (non-issue for rescaling; reporting added)

- No code rescales QSM. `possible_ppm_units` is a flag only.
- NIfTI `scl_slope`/`scl_inter` are applied by nibabel in P9 and by FSL in P8.
- P8 QC thresholds are unit-free.
- New report columns: `qsm_dtype`, `qsm_scl_slope`, `qsm_scl_inter` and
  `qsm_integer_valued` (P9), and `qsm_scl_slope`/`qsm_scl_inter` (P8 QC).
- Between-subject reference offsets remain a modelling decision (D4).

### CR-11 CSF referencing and reliability wording (confirmed bug; fixed)

**Where:** `P9 915e78b:773-774` added "CSF reference unavailable" to the primary
`analysis_pass`, so `--csf-ref` silently shrank the unreferenced analysis set. NOCSF passed
without `--csf-ref` and failed with it.

**Fix:**
- `analysis_pass` no longer depends on CSF.
- `csfref_analysis_pass` (= `analysis_pass` and a finite CSF reference) selects the rows of
  `*_csfref_analysis.csv`.
- A test confirms the primary `_analysis.csv` is identical with and without `--csf-ref`.

**Wording:** `reliability()` (`:445-473`) printed "Reliability check ... Spearman-Brown
coefficient" for left/right homologous correlations. These measure anatomical concordance
inflated by shared global variance, not test-retest reliability. The function is renamed
`lr_concordance` and its message says so; the computation is unchanged. `test_retest()`,
which uses repeat scans, is the actual reliability estimate.

## Priority 2: engineering robustness

### CR-12 Stale outputs, atomicity, concurrency and reuse provenance (fixed except D7)

- **Non-atomic P8 QC (`P8 915e78b:481`):** an interrupted write could leave a truncated QC
  file. Because `ERROR:` lines are written last, a truncated file could read as passing. QC
  is now written to a temporary file and renamed, its last line is `qc_complete=1`, and
  both P8 (reuse) and P9 (eligibility) require that line.
- **P1 matrix record (`:174-191, 211`):** it held only input hashes, so a matrix, T1_brain
  or P1 QC replaced later would still verify, and the Dice read from P1's QC might belong
  to another run. The record now also holds hashes of the matrix, T1_brain and P1 QC. It is
  deleted before re-registration and written atomically after success.
- **Segmentation record (`:254`):** it lacked the segmentation's own hash and survived an
  interrupted SynthSeg run. It now records `t1seg_sha1`, is deleted before SynthSeg runs,
  and is written atomically.
- **Concurrent runs:** there was no guard. Six duplicated IDs in the worklist make
  concurrent P8 jobs for one ID realistic, and P1 also shares a work directory per ID. P8
  now takes a per-ID `flock` lock and exits 16 if another run holds it.
- **P9 and mask outputs:** all CSVs, masks and records are written atomically.
- **Legacy P1 matrices:** matrices registered before P8 kept records are still reused when
  P1's QC names the same paths and the inputs are not newer than the matrix. This is not
  content-verified, so it is a potential risk. It is recorded per subject as
  `p1_matrix_provenance=legacy_path_mtime` (decision D7: accept, or force re-registration).

### CR-13 ROI labels (non-issue)

- The 68 label IDs (1001-1003, 1005-1035 and their 2000-series counterparts) and names
  match the independent `FreeSurferColorLUT.txt`, and all 68 are unique.
- 1004/2004 (corpus callosum) and 1000/2000 (unknown) are correctly excluded.
- P8's `CORTEX_L` equals P9's `FS_L`.
- P9 also checks CSV names against the LUT at run time; an alphabetically ordered CSV is
  rejected.

### CR-14 Exclusion reasons and reproducibility (confirmed bug, Low; fixed)

- Rows with `runnable != yes` were dropped without a record (`P9 915e78b:597-598`). They
  are now listed in `*_excluded.csv`.
- The support-mask check ran before the QSM-exists check, so a missing QSM was reported as
  "support mask missing". The first missing input in pipeline order is now reported.
- With `--ids-file` or `--sample-root`, worklist-wide exclusions of unselected subjects
  were listed. They are now restricted to the selection.
- `*_run_info.json` records the command, all options, the P9 script hash, the hashes of
  the worklist, label CSV, LUT and failure list, package versions and counts.
  `analysis_pass` is a deterministic function of these inputs.

### CR-15 Language and compatibility (verified)

- No non-ASCII text exists in any code file, and all messages, comments and CSV fields are
  English.
- P1 is unchanged.
- P8 keeps its interface and adds environment variables only.
- P9 adds the `--p1-qc-dir` option, removes `--allow-legacy-labels`, and requires
  `t1_file` (present in `Clinical_cohort/worklist.csv`).
- All earlier output columns are kept. New columns are listed in `CHANGELOG.md`.
- `qsm_matches_P8` changed meaning from path to content.

## Additional findings

- **CR-16 (potential risk):** inconsistent header geometry was not detected. Subject
  00413156 has a reported 1.83 mm slice thickness. If pixdim and the sform/qform scaling
  disagree, FSL coordinates (pixdim) and world coordinates (sform) describe different
  geometries. P8 now reports `qsm_pixdim_sform_rel_diff` and `t1_pixdim_sform_rel_diff` and
  warns above 1%, without failing (decision D10).
- **CR-17 (methodological):** re-scanned patients (13 hospital_ids with two scans) enter
  the analysis set twice. They are flagged by `repeat_patient`; which scan to analyse is
  decision D6.

## External review of revision 3 (round 4)

Five comments from an external static review were checked against the code (commit
`271576d`) and the worklist. Each verdict says whether the comment is technically correct,
and what changed. No threshold, inclusion rule, support definition or ROI statistic was
changed. Every change reports, records or adds an option, and the defaults are unchanged.

| ID | Comment | Verdict | Change |
|---|---|---|---|
| CR-18 | P8 QC thresholds and P1's early Dice stop decide inclusion | Correct, with two corrections (below) | Per-rule QC record; P1 early stop removed (same exclusions); exclusion report by stratum |
| CR-19 | Estimated support is not the reconstruction domain | Correct (already CR-09/D3); the hole-filling half matters most for float QSM | Per-ROI zero counts, `qsm_quantized`, float-hole warning, run comparison tool |
| CR-20 | `--min-reportable-rois 1` does not give complete cases | Correct (already D2) | ROI missingness table by protocol stratum; complete-case count |
| CR-21 | P9 ignores `include_main_model` | Correct | `--main-model-only`; note when the analysis set holds other rows |
| CR-22 | Registration uses a downsampled T1 | Claims correct; a related provenance bug found | Registration settings checked before reuse and recorded; comparison workflow |

### CR-18 QC thresholds and the P1 Dice stop (correct with corrections; reporting changes)

**What is correct:**
- **Cortex volume (350-700 mL):** it is measured on the QSM grid, so it falls with slab
  coverage. 250 of 711 runnable scans (172 of 479 main-model rows) have less than 125 mm.
  This is the clearest case of a rule that tests field of view rather than segmentation.
- **ICV (1200-1900 mL):** it depends on head size and therefore on sex.
- **`reg_mask_dice` (P1, gate 0.60):** it compares the QSM non-zero mask with the BET mask of
  the whole T1. Suppose registration is perfect, the slab holds a fraction c of the BET
  brain, and a fraction e of that is non-zero QSM. Then Dice = 2ec / (1 + ec). With P1's own
  estimate e = 0.73, Dice falls below 0.60 once c < 0.59. A short slab can fail this gate
  with correct registration.
- **P1's early stop:** this was worse than the comment says. In rigid mode, P1 exited before
  keeping the matrix, so P8 stopped with exit 5. The scan had no label and no QC metrics,
  could not be inspected, and could not be counted under any other rule. Also,
  `QSM_ROI_DICE_MIN` could change native-branch inclusion through P1, while P8's own gate
  was fixed at 0.60.

**What needs correcting:**
- `seg_qsm_dice` and ribbon coverage are computed on the QSM grid. Cortex outside the slab
  is not on that grid, so slab truncation does not lower them. Only reconstruction erosion
  and genuine zeros do.
- A ribbon with less than half of its voxels non-zero also means unreliable cortical QSM.
  That rule is not only a registration check.

**Changes:**
- P8 runs P1 with `QSM_ROI_DICE_MIN=0` in rigid mode. Step [5] still applies
  `reg_mask_dice < 0.60`, so a low-Dice scan is excluded exactly as before. It now exits 15
  with a label and full QC instead of 5 with nothing.
- P8 records:
  - `qc_failed_rules`, the machine-readable names of the failed rules;
  - `synthseg_cortex_ml`, the cortex volume of the whole-brain T1 segmentation;
  - `qc_report_only_failed_rules=cortex_vol_t1`, which applies the same 350-700 mL range to
    `synthseg_cortex_ml`. This rule is evaluated but not applied.
- New `qc_exclusion_report.py` counts, for each rule, the scans evaluated, failed, and
  failed by that rule alone. It reports these overall and by protocol, coverage, diagnosis
  and any covariate (sex via `--covariates`). It also compares the QSM-grid and T1 cortex
  rules.
- Thresholds are unchanged (decision D1).

### CR-19 Support masks and zero-valued voxels (correct; reporting changes)

**The comment restates CR-09/D3:**
- An estimated support cannot tell genuine zeros touching the brain edge from padding.
- It fills enclosed zero holes up to 10 mm3.

**The second point is the more important one for non-quantized QSM:**
- An exact 0.0 is essentially never a reconstructed value there.
- So a filled hole is more likely tissue the reconstruction masked out than a genuine zero.

**Reporting gap found:** `qsm_integer_valued` missed quantized data stored as integers with
a non-unit `scl_slope` (e.g. int16 at 0.001 ppm), since the scaled values are not integers.

**Changes:**
- `native_roi_coverage*.csv` gains `valid_zero_voxels` (zeros counted as values) and
  `unsupported_zero_voxels` (zeros left out as padding).
- `*_qc.csv` gains `cortex_valid_zero_fraction` and `cortex_unsupported_zero_voxels`.
- P9 and the support records report `qsm_quantized`, meaning stored as integers or
  integer-valued.
- `make_support_masks.py` warns when holes were filled in non-quantized QSM.
- New `compare_p9_runs.py` compares a run at another `--max-hole-mm3` with the primary run,
  ROI by ROI and by stratum.
- The support definition is unchanged (decision D3).

### CR-20 Completeness of ROI values (correct; reporting change)

- With `--min-reportable-rois 1`, `*_analysis.csv` holds subjects with missing ROIs.
- A partially covered ROI (coverage of at least 0.5) is summarised over its covered part
  only.
- Missingness depends on slab coverage, so values are not missing at random across
  protocols.

**Changes:**
- New `native_roi_missingness*.csv`: per ROI, overall and per protocol stratum, the number
  of analysis-set subjects, the number with a value, the missing fraction and the median
  `coverage_vs_t1`.
- P9 prints complete-case counts per stratum.
- The default is unchanged (decision D2).

### CR-21 Main-model selection (correct; option added)

- P9 and `make_support_masks.py` select `runnable=yes` only.
- 232 of 711 runnable rows have `include_main_model` other than `yes`. These include 119
  rows without `model_label`, plus MCI, other neurological and vascular groups.
- `*_analysis.csv` carries `include_main_model`, but nothing stopped it from being used as
  the main-model cohort.

**Changes:**
- New P9 `--main-model-only` option. Rows outside the main model are listed in
  `*_excluded.csv`.
- Without the option, P9 prints how many analysis-set subjects are outside the main model.
- Extra support masks are harmless, because P9 selects the cohort, so
  `make_support_masks.py` is unchanged.
- Making the option the default changes N (decision D11).

### CR-22 Registration on a downsampled T1 (claims correct; provenance bug fixed)

**The comment is right on both counts:**
- The final extraction stays on the native QSM grid.
- Whole-brain Dice cannot show local ribbon misalignment.

**Provenance bug found:** reuse checks ignored registration settings. A comparison run with
`QSM_ROI_T1_SUBSAMP=0` (or another `QSM_ROI_SEARCH_MODE`) in the same native root would have
silently reused the downsampled-T1 matrix and its label.

**Changes:**
- P8 compares the requested T1 resolution and search mode with those P1 recorded in its
  rigid QC, and re-registers when they differ.
- It records `p1_t1_resolution` and `p1_search_mode` in its QC, and the label is rebuilt
  when they change.
- The README describes the full-resolution comparison: a second native root, then
  `compare_p9_runs.py`.
- The default resolution is unchanged (decision D12).
- `QSM_ROI_SEARCH_DEG` is not recorded by P1. After changing it, use `P8_FORCE_P1=1`.
