# Methodological decisions requiring approval

No scientific threshold, inclusion criterion, referencing scheme or summary statistic has
been changed without approval. Each item below gives:
- the current behaviour;
- the evaluation;
- the options, with a recommendation.

The pilot (README, "Representative pilot") produces the data needed to decide most items.
Each decision lists the columns to examine.

---

## D1. P8 QC thresholds (`P8_native_aparc_qsm_ALL.sh`, QC step)

| Threshold | Measures | Technical failure it catches | Valid variation it can exclude |
|---|---|---|---|
| SynthSeg ICV 1200-1900 mL | Intracranial volume (T1) | Failed segmentation, wrong image | Small heads, which are more common in women and smaller adults. ICV does not shrink with atrophy, so the floor selects by sex and body size, not disease. |
| Cortical volume 350-700 mL | Cortex **within the QSM slab**, on the QSM grid | Mislabelled or missing cortex | Short slabs (250 of 711 scans have < 125 mm coverage, nearly all iso1.0mm) and severe atrophy (AD, FTD). Protocol, field of view and disease severity are confounded. |
| Registration Dice >= 0.60 | QSM nonzero mask vs whole-T1 BET mask | Gross misregistration, failed BET | Short slabs and eroded reconstruction masks lower Dice with correct registration. |
| Cortical ribbon coverage >= 0.50 | Ribbon voxels with nonzero QSM, within the slab | Severe mask erosion, misplacement | Mostly mask erosion at the cortex. Blind to slab truncation, so it does not protect against it. |
| GM-WM contrast > 0 ppb | Median cortex minus median WM QSM | Gross label/QSM mismatch | Unlikely to exclude valid scans (GM is paramagnetic relative to diamagnetic WM). Weak: millimetre-scale misregistration usually passes. |

**Options.**
- (a) Keep all thresholds as they are.
- (b) Compute the cortical-volume check on the T1 segmentation instead of the
  slab-truncated QSM grid. P9 now reports `t1_cortex_ml` for this.
- (c) Replace the ICV range with a sex-specific or cohort-percentile range, or report ICV
  only.
- (d) Report Dice, ribbon coverage and the volume checks, and gate only on gross failures,
  with visual QC for borderline cases.

**Registration Dice and slab coverage.** Suppose registration is perfect, the slab holds a
fraction c of the T1 BET brain, and a fraction e of that is non-zero QSM. Then
`reg_mask_dice` = 2ec / (1 + ec). With e = 0.73 (P1's estimate), Dice < 0.60 once c < 0.59,
so the shortest slabs can fail the gate with correct registration. `seg_qsm_dice` and ribbon
coverage are computed on the QSM grid and do not depend on slab coverage.

**Now recorded (revision 4).**
- P8 writes `qc_failed_rules`, the rule names: `reg_mask_dice`, `synthseg_icv`,
  `cortex_vol_qsm_grid`, `ribbon_coverage`, `gm_wm_contrast`, `seg_qsm_dice`.
- P8 writes `synthseg_cortex_ml`, the whole-brain T1 cortex volume.
- P8 writes `qc_report_only_failed_rules=cortex_vol_t1`: the same 350-700 mL range applied
  to that volume, evaluated and not applied.
- P1 no longer stops at its Dice gate inside P8. A low-Dice scan keeps its label and full
  QC, and still fails the same 0.60 rule.

**Recommendation.** Before fixing thresholds, run `qc_exclusion_report.py` on the pilot,
then the full cohort. Pass sex with `--covariates`, because the worklist has no sex column.
- `qc_exclusions_by_rule.csv`: how many scans each rule excludes, and how many it excludes
  alone, by protocol, coverage, diagnosis and sex.
- `qc_cortex_rule_comparison.csv`: how many scans option (b) would rescue or newly exclude.

Columns: `P8_native_qc_errors`, `P8_cortex_vol_ml`, `t1_cortex_ml`, `P8_synthseg_icv_ml`,
`P8_reg_mask_dice`, `P8_ribbon_in_qsm_frac`, `P8_gm_wm_contrast_ppb`.

I recommend (b) for the cortical-volume check at least. As implemented, it penalises the
short-coverage protocol, which is enriched for AD and DLB.

## D2. ROI coverage criteria (`P9 --min-voxels 20 --min-mm3 20 --min-coverage 0.5`)

**Current rule.** An ROI value is reported when all of these hold:
- at least 20 valid voxels and at least 20 mm3 of valid volume;
- valid volume at least 50% of the ROI's T1-space volume. This includes cortex outside the
  slab: `coverage_vs_t1`.

**Evaluation.**
- **Thresholds:** the physical-volume gate makes the minimum comparable across 0.8 mm and
  1.0 mm protocols. The voxel-count gate is redundant for voxels of 1 mm3 or smaller.
- **Coverage above 1:** `fov_fraction` and `coverage_vs_t1` can exceed 1 slightly through
  nearest-neighbour sampling (up to 1.016 observed on synthetic data). Larger values
  indicate a label/segmentation mismatch, and `max_fov_fraction` is reported for that.
- **Missing values are not random:** ROIs at the top (vertex) or bottom (temporal pole,
  orbitofrontal) of short slabs fail the coverage gate more often. Missing ROI values
  therefore depend on protocol.
- **Unrestricted parts:** the gate does not restrict which part of an ROI is sampled.
  Partially covered ROIs are biased toward the retained part.

**Options.**
- (a) Keep 0.5.
- (b) Choose another coverage threshold after inspecting the coverage distribution by ROI
  and protocol on the pilot (`native_roi_coverage*.csv`).
- (c) Add `max_fov_fraction > 1.1` as a subject-level QC failure (currently report-only).
- (d) Include protocol and coverage as covariates in downstream models.

**Recommendation.** (b) and (d). Decide (c) after inspecting the observed
`max_fov_fraction` values.

**Downstream completeness (`--min-reportable-rois`, default 1).**
- Per-ROI mass-univariate models can use subjects with some missing ROIs.
- Multivariate or complete-case models need 68, or an explicit imputation plan.

Please confirm which model the extraction feeds.

**Missingness by protocol (revision 4).** `native_roi_missingness*.csv` gives, per ROI,
overall and per protocol stratum (QSM voxel size | short/full coverage):
- the number of analysis-set subjects and the number with a value;
- the missing fraction and the median `coverage_vs_t1`.

P9 also prints complete-case counts per stratum. Use these to choose the completeness rule
or imputation model, and to check whether missingness differs by protocol.

## D3. QSM support masks

**Current behaviour.** Formal extraction requires `--support-dir` with masks that carry a
verified provenance record. Two kinds exist:
- **`support_type=reconstruction`:** the reconstruction pipeline's own mask (registered with
  `make_support_masks.py --recon-mask-template`). This is the scientifically correct
  definition of valid QSM support.
- **`support_type=estimated_from_qsm`:** recovered from the QSM image as the nonzero region
  plus enclosed zero holes of at most 10 mm3. This is an estimate, not independent
  evidence of the reconstruction support.

**Measured behaviour of the estimate:**

| QSM storage | Interior zeros | Unsupported, outer 3 mm (cortex) | Unsupported, deeper |
|---|---|---|---|
| Float32 | essentially none | = nonzero proxy | = nonzero proxy |
| Integer, 1 ppb steps, SD 30 ppb | 1.3% | 0.41% | 0.000% |
| Integer, 1 ppb steps, SD 10 ppb | 4.1% | 1.40% | 0.000% |
| Integer, 1 ppb steps, SD 3 ppb | 13.1% | 5.61% | 0.45% |

The integer rows use the 10 mm3 limit. The nonzero proxy (limit 0) leaves all zeros
unsupported, for example 13.1% at SD 3 ppb.

**Sensitivity to the hole limit (SD 3 ppb):**

| Limit | Outer 3 mm unsupported | Deeper unsupported |
|---|---|---|
| 0 mm3 | 13.09% | 13.1% |
| 10 mm3 | 5.61% | 0.45% |
| 50 mm3 | 5.46% | 0.011% |

Above about 10 mm3 the limit barely matters for the cortex, but it starts to admit masked
regions (e.g. CSF pockets) as "valid zeros".

**Limitation.** Zero-valued voxels touching the brain edge cannot be told apart from
background. For quantized QSM, an estimated mask therefore undercounts valid cortex voxels
whose value is exactly 0. This shifts ROI statistics away from 0 for ROIs with values near
0.

**Holes in non-quantized QSM.** An exact 0.0 is essentially never a reconstructed value in
float QSM. A filled hole there is therefore more likely tissue the reconstruction masked out
than a genuine zero. `make_support_masks.py` now warns when this happens; check
`filled_voxels` for those maps.

**Now reported (revision 4).**
- `qsm_quantized`: stored as integers (any `scl_slope`) or integer-valued.
  `qsm_integer_valued` missed int16 data with a slope such as 0.001.
- Per ROI: `valid_zero_voxels` and `unsupported_zero_voxels`.
- Per subject: `cortex_valid_zero_fraction` and `cortex_unsupported_zero_voxels`.

**Sensitivity analysis.**
1. Rebuild the masks with `--max-hole-mm3 0` (and e.g. 50) into another directory.
2. Run P9 on them with `--suffix`.
3. Compare the result with the primary run using `compare_p9_runs.py` (README).

**Options.**
- (a) Obtain the reconstruction masks from the scanner or reconstruction pipeline. This is
  preferred.
- (b) Use estimated masks for float QSM, where they equal the nonzero proxy and the
  limitation does not arise.
- (c) For integer QSM without reconstruction masks, treat results as a sensitivity
  analysis, or report `qsm_quantized` and model it.
- (d) Use a hole limit of 0 for non-quantized QSM, where filled holes cannot be genuine
  zeros, and keep 10 mm3 for quantized QSM.

Columns: `support_type`, `qsm_dtype`, `qsm_quantized`, `support_masks_summary.csv`
(`filled_voxels`), `qsm_nonzero_outside_support`, `valid_zero_voxels`,
`unsupported_zero_voxels`.

## D4. QSM units, scaling and reference

**Current behaviour.**
- No rescaling. NIfTI scaling is applied by nibabel and FSL.
- `possible_ppm_units` flags subjects whose cortical interquartile range is below 1
  (consistent with ppm), and nothing else.
- Values are unreferenced in the primary output. `--csf-ref` adds a lateral-ventricle-CSF
  referenced table with its own eligibility flag.

**Evaluation.**
- **Two protocols:** iso1.0mm and 0.8mm may use different reconstruction versions or
  reference conventions (e.g. whole-brain mean vs none). A between-protocol offset would
  appear as a protocol effect.
- **CSF referencing:**
  - The CSF reference uses SynthSeg ventricles eroded by one voxel. SynthSeg has no choroid
    plexus label, so choroid plexus is inside the ventricle label.
  - Ventricles may be masked out by the reconstruction (P1 notes that reconstruction zeros
    CSF). In that case no CSF reference is available.

**Options.**
- (a) Keep unreferenced values as primary, with protocol as a covariate.
- (b) Use CSF-referenced values as primary. This needs a check of how often CSF is
  unavailable (`csf_voxels`), and a decision about choroid plexus (e.g. restrict to frontal
  horns).
- (c) Use a whole-brain or deep-WM reference.

Before any of these: confirm units with the reconstruction documentation, and inspect
`possible_ppm_units` and value ranges per protocol on the pilot.

## D5. Stage-specific use of the P1 failure list (implemented as instructed; confirm two points)

**Implemented behaviour.** Each failure-list entry is classified from P1's QC file:
- **`native`** (input or rigid-registration failure): excluded.
- **`mni_only`** (FNIRT, T1->MNI Dice, MNI coverage or output dimensions): stays in the
  native analysis.
- **`unknown`** (no P1 QC file, no recognised warning, or a P1 QC showing no failure):
  excluded.

**Please confirm:**
1. that `mni_only` entries should enter the native analysis. This follows your
   instruction, and it increases the sample size relative to the audited code;
2. how to treat `unknown` entries. The current choice excludes them, conservatively. The
   alternative is to rely on P8's own rigid-registration QC, which every native subject
   passes anyway.

Columns: `P1_failure_stage`, `P1_failure_evidence`, and `*_excluded.csv`.

## D6. Repeated scans of the same patient

**Current behaviour.** 13 hospital_ids have two image directories. Both scans enter the
analysis set when they pass, flagged by `repeat_patient`. `test_retest()` uses them to
estimate reliability (scan intervals are unknown to the script).

**Options.**
- (a) Keep both and model the dependence (random effect per patient).
- (b) Keep one scan per patient by a pre-specified rule: first scan, or the scan with the
  better protocol or coverage.
- (c) Keep both for the reliability analysis only.

## D7. Legacy P1 registration matrices

**Current behaviour.** Matrices registered before P8 kept content records are reused if
P1's QC names the same QSM/T1 paths and the files are not newer than the matrix.
`p1_matrix_provenance=legacy_path_mtime` marks them. This evidence is weaker than content
hashes.

**Options.**
- (a) Accept them, as now.
- (b) Run the first full P8 pass with `P8_FORCE_P1=1`. FLIRT re-registration is
  deterministic for the same inputs, so every matrix then has content provenance
  (`sha1_record`), at the cost of one rigid registration per subject.

**Recommendation:** (b), for a fully auditable cohort.

## D8. Primary ROI statistic

The median and the mean are both written. Nothing has changed, and the downstream choice
is yours. The median is robust to veins and calcifications in the ribbon. For integer QSM
the median is discrete, and the mean is smoother.

## D9. Legacy labels

Legacy labels are now prohibited in formal analysis (CR-04), as one of your two specified
options. Revised P8 rebuilds them automatically.

Please confirm whether historical results computed from legacy labels should be:
- re-derived for all subjects, or
- re-derived only for subjects that `check_grid_conversion.py` flags as `label_offset`.

## D10. Header-geometry and rigidity warnings

P8 now reports these, and warns above 1%:
- `qsm_pixdim_sform_rel_diff` and `t1_pixdim_sform_rel_diff`: disagreement between pixdim
  and sform/qform voxel sizes;
- `q2t_scale_dev`: deviation of the registration matrix's singular values from 1.

These are not QC failures. Decide after the pilot whether either should become one. Subject
00413156 (reported 1.83 mm slices) is the known candidate.

## D11. Main-model cohort (`include_main_model`)

**Current behaviour.**
- P9 processes every `runnable=yes` row and carries `include_main_model` into its tables.
- 232 of 711 runnable rows are not in the main model. They include 119 without
  `model_label`, plus MCI, other neurological and vascular groups.
- `--main-model-only` (revision 4) restricts P9 to `include_main_model=yes` and lists the
  other rows in `*_excluded.csv`.
- Without it, P9 prints how many analysis-set subjects are outside the main model.

**Options.**
- (a) Keep the default. Run the main analysis with `--main-model-only` (with `--suffix`),
  and use the full run for secondary analyses.
- (b) Make `--main-model-only` the default.

**Recommendation:** (a), unless the full runnable set is never used. Either way, pass the
same option to `qc_exclusion_report.py`.

## D12. T1 resolution of the QSM->T1 registration

**Current behaviour.**
- P1 registers the QSM to the T1 downsampled by `fslmaths -subsamp2` (about 2 mm).
- P8 composes the matrix with the full-resolution grid conversion. Labels come from the
  full-resolution SynthSeg, and extraction stays on the native QSM grid.
- The downsampled reference limits registration precision. A millimetre of misalignment
  matters for a 2-3 mm cortical ribbon, and whole-brain Dice cannot detect it.

**Comparison (pilot).**
1. Run P8 with `QSM_ROI_T1_SUBSAMP=0` into a second `QSM_ROI_NATIVE_ROOT`. P8 now records
   the registration resolution and never reuses a matrix or label made with other settings.
2. Run `make_support_masks.py` and P9 on that root.
3. Compare the two P9 runs with `compare_p9_runs.py`.
4. Inspect label overlays for both resolutions in each stratum.

**Options.**
- (a) Keep the downsampled registration if ROI differences are small relative to
  between-subject variation.
- (b) Register at full resolution (`QSM_ROI_T1_SUBSAMP=0`) for the whole cohort. This
  costs a slower rigid registration per subject.
