# Changelog: Clinical_cohort_revised

Each change is marked for whether it can affect:
- **Values:** QSM values in output tables;
- **ROIs:** which ROI values are reported;
- **QC:** P8/P9 QC status;
- **N:** which subjects enter the analysis set.

"yes" means the change can alter that quantity for some subjects; "no" means it cannot.

## Revision 4 (this revision): response to the external review (CODE_REVIEW CR-18 to CR-22)

No threshold, inclusion rule, support definition or ROI statistic changed. With identical
inputs and options, `analysis_pass` and ROI values equal those of revision 3. Exit codes
are the one exception: a scan below P8's Dice gate now exits 15 with a label, not 5.

### P8 (`P8_native_aparc_qsm_ALL.sh`)

| Change | Values | ROIs | QC | N |
|---|---|---|---|---|
| P1 runs with its Dice gate off (`QSM_ROI_DICE_MIN=0`) in rigid mode. P8's own 0.60 gate still fails such scans, which now get a label and full QC (exit 15, not 5) (CR-18) | no | no | no (same rule) | no |
| `qc_failed_rules` and `qc_report_only_failed_rules` recorded; `synthseg_cortex_ml` and the report-only `cortex_vol_t1` rule added (CR-18) | no | no | no | no |
| Registration settings (T1 resolution, search mode) compared with P1's rigid QC before a matrix or label is reused; `p1_t1_resolution`, `p1_search_mode` recorded. Labels without them are rebuilt once, reusing verified matrices (CR-22) | no | no | no | no |

### P9 (`P9_extract_native_qsm_ALL.py`)

| Change | Values | ROIs | QC | N |
|---|---|---|---|---|
| `--main-model-only` option (default off); note when the analysis set includes other rows (CR-21) | no | no | no | only with the option |
| Coverage table: `valid_zero_voxels`, `unsupported_zero_voxels`. QC table: `cortex_valid_zero_fraction`, `cortex_unsupported_zero_voxels`, `qsm_quantized`, `qsm_protocol`, `qsm_fov_si_mm`, `qsm_coverage_group`, and P8's `qc_failed_rules`, `qc_report_only_failed_rules`, `synthseg_cortex_ml`, `p1_t1_resolution`, `p1_search_mode` (CR-18, CR-19) | no | no | no | no |
| New `native_roi_missingness*.csv` and complete-case counts per protocol stratum (CR-20) | no | no | no | no |

### Other

- **`make_support_masks.py`:** records `qsm_quantized`, and warns when holes were filled in
  non-quantized QSM (CR-19). Masks are unchanged.
- **New `qc_exclusion_report.py`:** exclusions per rule by protocol, coverage, diagnosis
  and covariates (CR-18).
- **New `compare_p9_runs.py`:** paired comparison of two P9 runs (CR-19, CR-22).
- **New `strata.py`:** protocol strata shared with `select_pilot.py` (whose output is
  unchanged).
- **Tests:** stub P1 now applies the Dice gate and records `t1_pixdim` and `search_mode`
  like the real P1; tests added for every change (`TEST_REPORT.md`).
- **Documents:** CODE_REVIEW (CR-18 to CR-22), METHODOLOGICAL_DECISIONS (D1-D3 updated, D11
  and D12 added), README (sensitivity analyses).

## Revision 3 (`f32283c`): independent audit and provenance hardening

### P9 (`P9_extract_native_qsm_ALL.py`)

| Change | Values | ROIs | QC | N |
|---|---|---|---|---|
| QSM verified by SHA-1 against P8's `qsm_sha1`; a replaced file at the same path is no longer accepted (CR-01). `qsm_matches_P8` now means content; `qsm_path_matches_P8` added | no | no | yes | yes (fewer) |
| `t1_file` required; worklist T1 verified by SHA-1 against P8's `t1_sha1`; missing T1 file is an explicit exclusion (CR-03). `t1_matches_P8`, `t1_path_matches_P8` added | no | no | yes | yes (fewer) |
| `t1_file` included in duplicate-ID conflict detection (CR-03) | no | no | no | yes (fewer) |
| Support-mask provenance record verified (current QSM hash, mask hash, known `support_type`); masks without a valid record rejected (CR-02). `support_type`, `support_method`, `support_max_hole_mm3`, `support_provenance_ok`, `qsm_nonzero_outside_support` added | no | no | yes | yes (fewer) |
| `--allow-legacy-labels` removed; legacy labels never eligible (CR-04) | no | no | yes | yes (fewer, only if the option was used) |
| P8 QC must be complete (`qc_complete=1`) | no | no | yes | yes (fewer; P8 rebuilds such QC) |
| P1 failure list classified by stage from P1 QC (`--p1-qc-dir`); `mni_only` failures no longer exclude; `native` and `unknown` still do (CR-05). `P1_failure_stage`, `P1_failure_evidence` added | no | no | no | **yes (more)**; please confirm (D5) |
| `--csf-ref` no longer changes `analysis_pass`; new `csfref_analysis_pass` selects `*_csfref_analysis.csv` (CR-11) | no | no | no | yes (primary set restored when `--csf-ref` is given) |
| Exclusion records: non-runnable rows listed; first missing input in pipeline order reported; exclusions restricted to the `--ids-file`/`--sample-root` selection (CR-14) | no | no | no | no |
| `*_run_info.json` (command, options, script hash, input hashes, versions, counts) | no | no | no | no |
| All CSV outputs written atomically | no | no | no | no |
| Report-only columns: `qsm_dtype`, `qsm_scl_slope`, `qsm_scl_inter`, `qsm_integer_valued`, `t1_cortex_ml`, `max_fov_fraction`, `P8_qc_complete`, `P8_cortex_vol_ml`, `P8_synthseg_icv_ml`, `P8_p1_matrix_provenance`, `P8_q2t_scale_dev`, `P8_qsm_pixdim_sform_rel_diff`, `P8_t1_pixdim_sform_rel_diff` | no | no | no | no |
| Left/right homologous correlation renamed `lr_concordance` and labelled "not test-retest reliability" (log text only) | no | no | no | no |

ROI statistics (`extract_one`) and all thresholds are unchanged. With identical inputs and
gates, the ROI values equal those of revision 2.

### P8 (`P8_native_aparc_qsm_ALL.sh`)

| Change | Values | ROIs | QC | N |
|---|---|---|---|---|
| QC written atomically, ending with `qc_complete=1`; reuse requires it (CR-12) | no | no | yes (an interrupted QC can no longer read as passing) | indirectly |
| P1 matrix record also holds hashes of the matrix, T1_brain and P1 QC; deleted before re-registration, written atomically; `p1_matrix_provenance` recorded (CR-12) | no | no | no | no |
| Segmentation record holds the segmentation's hash; deleted before SynthSeg, written atomically; reuse and skip require it (CR-12) | no | no | no | no |
| Per-subject `flock` lock; a concurrent run exits 16 (CR-12) | no | no | no | no |
| Report-only QC fields and warnings: `qsm_pixdim_sform_rel_diff`, `t1_pixdim_sform_rel_diff`, `q2t_scale_dev`, `qsm_scl_slope`, `qsm_scl_inter` (CR-06, CR-10, CR-16). Warnings do not change pass/fail | no | no | no | no |
| Unused `json` import in the QC step replaced by `os` | no | no | no | no |

### Support masks (`make_support_masks.py`)

| Change | Values | ROIs | QC | N |
|---|---|---|---|---|
| Provenance record `<IID>_support.json` (`support_type`, method, parameters, QSM and mask SHA-1); stale masks rebuilt automatically; atomic writes (CR-02) | yes (stale masks replaced) | yes | no | yes |
| `--recon-mask-template` registers genuine reconstruction masks (grid-checked) | yes (if used) | yes | no | yes |
| Summary reports QSM data type, integer values, filled voxels, nonzero QSM outside support | no | no | no | no |

### Other

- **`check_grid_conversion.py`:** adds `p8_qc_grid_conversion` per subject.
- **New `tests/`:** pytest suite with stub tools and a synthetic generator (87 tests; 6 need
  real FSL).
- **New documents:** `CODE_REVIEW.md`, `TEST_REPORT.md`, `METHODOLOGICAL_DECISIONS.md`.
- **`README.md`:** rewritten around the pilot and the extraction workflow.

## Revision 2 (`915e78b`): verified output reuse and stricter eligibility

| Change | Values | ROIs | QC | N |
|---|---|---|---|---|
| P8 skips a label only after verifying QC, inputs and hashes; force flags honoured | no | no | yes | yes |
| SynthSeg reused only with a verified source record | no | no | no | no |
| P9: ROI requires a finite, positive T1 volume; label/segmentation disagreement fails the subject | no | yes | yes | yes (fewer) |
| P9: at least `--min-reportable-rois` reported values required (default 1) | no | no | yes | yes (fewer) |
| P9: proxy outputs carry `_proxy`; mode recorded | no | no | no | no |
| `make_support_masks.py`, `select_pilot.py` added | yes | yes | no | yes |

## Revision 1 (`2f998e8`): first revised pipeline

| Change | Values | ROIs | QC | N |
|---|---|---|---|---|
| P8 grid conversion corrected (FLIRT x-flip for neurologically stored T1s). Affected labels (neurological storage, even first dimension) were offset by one T1 voxel | **yes** | yes | yes | yes |
| P9: label names checked against FreeSurfer; conflicting duplicate IDs excluded; corrected `analysis_pass`; coverage relative to T1 volume (`--min-mm3 20`, `--min-coverage 0.5`); per-subject error handling; CSF reference for 4D QSM | yes (4D CSF) | **yes** | yes | **yes** |

Revisions 1 and 2 introduced the coverage gates and analysis-set rules. Under
METHODOLOGICAL_DECISIONS they still await your confirmation, together with the items
listed there.
