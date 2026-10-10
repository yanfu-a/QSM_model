# Test report

## 1. Environment and scope

**Software:** Python 3.13.16, numpy 2.5.3, nibabel 5.4.2, pandas 3.0.6, scipy 1.18.1,
pytest 9.1.1, fslpy 3.29.1, GNU bash 5.2.

**Not available here:**
- FSL (`flirt`, `fslmaths`, `bet`, `fslinfo`);
- FreeSurfer/SynthSeg;
- the module system;
- real MRI data.

**How tools were replaced** (`tests/conftest.py`):
- **`flirt -applyxfm` (nearest neighbour):** reimplemented with fslpy, FSL's own coordinate
  library. It is independent of `native_grid.py`.
- **`fslinfo`:** reads the header with nibabel.
- **`mri_synthseg`:** copies a prepared segmentation.
- **P1 (rigid mode):** a fake with the real P1's skip rule and Dice gate (`QSM_ROI_DICE_MIN`).
  It writes the true registration matrix and the same QC fields as the real P1
  (`reg_mask_dice`, `t1_pixdim`, `search_mode`).

**Synthetic subjects** (`tests/synth.py`) are ellipsoidal brains with the 68 Desikan parcels,
WM and ventricles. Each has a ground-truth label on the QSM grid, obtained by sampling the
T1 labels straight through world space, without any FLIRT convention.

All results below are on synthetic data and say nothing about clinical validity.

**Run:**

```bash
python -m pytest Clinical_cohort_revised/tests -v                 # everything feasible here
python -m pytest Clinical_cohort_revised/tests -m "not slow" -v   # unit tests only (seconds)
module load fsl/6.0.4 && python -m pytest Clinical_cohort_revised/tests -m fsl -v    # on the cluster
QSM_TEST_REAL_FSL=1 python -m pytest Clinical_cohort_revised/tests -v                # P8 tests with real flirt
```

## 2. Results

| Suite | Tests | Passed | Skipped | Failed |
|---|---|---|---|---|
| `test_native_grid.py` (FLIRT conventions) | 9 | 9 | 0 | 0 |
| `test_p9_unit.py` (extraction, gates, parsing, zero counts, quantization, strata) | 29 | 29 | 0 | 0 |
| `test_support_masks.py` | 6 | 6 | 0 | 0 |
| `test_p8.py` (P8 end to end, stubs) | 25 | 25 | 0 | 0 |
| `test_geometry.py` (7 geometries, revised vs original P8) | 15 | 15 | 0 | 0 |
| `test_p9_integration.py` (24-subject synthetic cohort) | 11 | 11 | 0 | 0 |
| `test_reports.py` (QC exclusion report, run comparison) | 4 | 4 | 0 | 0 |
| `test_fsl_real.py` (needs FSL) | 6 | 0 | 6 | 0 |
| **Total** | **105** | **99** | **6** | **0** |

The full run took 246 s (revision 6). Static checks:
- `bash -n` on P8: clean.
- `shellcheck -S warning`: one pre-existing, intentional SC2163 (exporting variables by
  name).
- `pyflakes`: clean on all Python files and on P8's embedded Python blocks.

## 3. Required tests: expected and observed

| Requirement | Test(s) | Expected | Observed |
|---|---|---|---|
| QSM replaced at the same path | `test_p9_integration` (QREP); `test_p8::test_input_content_change...[qsm]` | P9: not eligible ("QSM content differs"); P8: label not reused, P1 re-run | As expected. `qsm_path_matches_P8=True`, `qsm_matches_P8=False` |
| Stale support masks | `test_support_masks::test_provenance_reuse_and_invalidation`; integration STALESUP, NOREC | Rebuilt when QSM content, parameters or mask change or the record is lost; P9 rejects stale or unrecorded masks | As expected |
| T1 provenance mismatch | integration T1SAME (same path, new content), T1MIS (different file); `test_p8...[t1]` | Not eligible; P8 re-runs P1 and SynthSeg | As expected |
| Duplicate worklist IDs with conflicting data | `test_p9_unit::test_duplicate_ids_conflicting_t1_are_excluded`; integration DUP | Excluded, with the conflicting `t1_file` values named | As expected |
| Affected and unaffected legacy labels | integration LEG_AFF, LEG_UNAFF; `test_geometry::test_grid_audit_classifies_legacy_labels` | Both rejected by P9; audit flags `label_offset` only for neurological even-first-dimension T1s | As expected |
| FLIRT transforms across geometries | `test_native_grid` (300 random oblique pairs vs fslpy; 8 subsamp2 geometries); `test_geometry` (7 end-to-end geometries) | Matches fslpy; labels reproduce ground truth | Max error vs fslpy < 1e-4 mm; label agreement >= 99.95% brain and >= 99.9% cortex in all 7 |
| 3D and singleton-4D QSM | `test_p9_unit::test_singleton_4d...`; integration QSM4D; geometry `radio_aniso`, `radio_permuted` | Identical handling, CSF reference available | As expected |
| Corrupt or missing images | integration CORRUPT (truncated label), NOQSM, NOT1, NOLABEL; `test_p9_unit::test_grid_mismatch_raises` | Per-subject error status or explicit exclusion; the run continues | As expected |
| Missing or stale QC files | integration QCMISS, QCSTALE (label edited after QC), QCINC (no `qc_complete`); `test_p8` invalidation cases | Not eligible; P8 rebuilds | As expected |
| Zero supported cortical voxels | `test_p9_unit::test_zero_supported_voxels...`; integration ZEROSUP (all-zero reconstruction mask) | "No valid supported cortical voxels"; label counts kept | As expected |
| ROI coverage and physical-volume thresholds | `test_p9_unit::test_voxel_volume_and_coverage_thresholds_are_inclusive`, `test_roi_without_finite_positive_t1_volume...`, `test_coverage_gate_needs_t1...` | Inclusive thresholds; zero, negative or infinite T1 volume fails the subject | As expected |
| Subject-level eligibility | `test_p9_integration::test_every_subject_has_the_expected_outcome` (24 IDs) | Each ID has its designed outcome; `_analysis.csv` rows equal `analysis_pass` | As expected |
| No QC-passing subjects | `test_run_with_no_passing_subject`; `test_run_with_no_extractable_subject_still_records_why` | Exit 0 with an empty analysis table; or exit 1 with exclusions and run info written | As expected |
| Repeated scans of one patient | integration OK1 and REP; `test_p9_unit::test_concordance_and_test_retest_run` | Both flagged `repeat_patient`; test-retest computed | As expected |
| ROI means, medians, counts vs an independent reference | `test_p9_unit::test_roi_statistics_match_independent_reference` (pandas groupby, 3 random volumes); `test_p9_integration::test_roi_values_match_independent_calculation` (68 ROIs recomputed from the files) | Equal | Counts identical; means within 1e-9 (unit) and 1e-6 (integration); medians within 1e-6 |
| CSF referencing does not change the primary set | `test_csf_reference_does_not_change_the_primary_set` | `_analysis.csv` identical with and without `--csf-ref` | Identical (`assert_frame_equal`) |
| Stage-specific P1 failures | `test_p9_unit::test_p1_failure_stage_classification` (9 cases); integration FAILMNI, FAILRIG, FAILUNK | `mni_only` stays; `native` and `unknown` excluded with reasons | As expected |
| Proxy outputs separate from formal | `test_proxy_outputs_never_touch_formal_outputs` | Formal files byte-identical after a proxy run | As expected |
| P8 verified reuse and invalidation | `test_p8` (8 invalidation cases, verified skip, legacy matrix) | Rebuild only when needed; P1 and SynthSeg reused only when verified | As expected |
| Interrupted and failed runs | `test_failure_after_cleanup...`, `test_failed_registration_removes_matrix_record`, `test_missing_registration_dice_fails_qc` | No stale QC, label or matrix record; a missing Dice fails QC and is not skipped later | As expected |
| Concurrent runs | `test_concurrent_run_is_refused` | Exit 16 while another run holds the lock | As expected |
| Low registration Dice (revision 4, CR-18) | `test_p8::test_low_registration_dice_keeps_label_and_fails_qc` | Label and complete QC are written; QC fails `reg_mask_dice`, as before (exit 15); not reused later | As expected |
| Per-rule QC record (CR-18) | `test_p8::test_fresh_run...` (`qc_failed_rules` empty, `synthseg_cortex_ml`); `test_reports::test_report_reads_p8_qc_files` | Rule names recorded, read by the report | As expected |
| Exclusions per rule and stratum (CR-18) | `test_reports::test_exclusions_are_counted_per_rule_and_stratum` (7 scans, every stage, sex covariate) | Hand-computed counts per rule, stratum, "alone" and cortex-rule comparison | As expected |
| Zero-valued voxels inside and outside support (CR-19) | `test_p9_unit::test_zero_valued_voxels...`; integration `test_report_columns_and_missingness_table` | Exact counts; statistics unchanged | As expected |
| Quantized QSM with `scl_slope` (CR-19) | `test_p9_unit::test_scaled_integer_qsm_is_reported_as_quantized`; `test_support_masks::test_quantization_is_recorded...` | int16 x 0.001: `qsm_integer_valued=False`, `qsm_quantized=True`; filled holes in float QSM flagged | As expected |
| ROI missingness by stratum (CR-20) | integration `test_report_columns_and_missingness_table` | 68 rows per stratum; missing fraction = 1 - reported / analysis set | As expected |
| Main-model selection (CR-21) | integration `test_main_model_only` | Without the option, a note; with it, other rows excluded with reason | As expected |
| Registration settings provenance (CR-22) | `test_p8::test_registration_settings_are_recorded_and_enforced`, `test_legacy_p1_qc_without_settings_is_re_registered`, invalidation case "registration settings" | A changed T1 resolution or search mode re-registers and rebuilds; unchanged settings skip; unrecorded settings re-register | As expected |
| Missing T1 cortex volume is not a pass (revision 6, CR-23) | `test_reports::test_exclusions_are_counted_per_rule_and_stratum` (GRIDNOT1, NOVOL, P1STOP) | Counted as not measured or not evaluated; "alone" only when every other applied rule passed | As expected. Before the fix, `of_which_pass_t1_rule` counted the unmeasured scan (2 instead of 1) |
| Volume table part of segmentation provenance (CR-24) | `test_p8::test_segmentation_reuse_requires_its_volume_table[delete, edit]`, `test_missing_volume_measures_fail_or_stop` | Deleted or edited table: label not reused, SynthSeg rerun, ICV restored; no table: exit 13; no ICV column: ICV rule fails | As expected. Before the fix, the deleted table led to reuse and a QC pass with `synthseg_icv_ml=nan` |
| Paired run comparison (CR-19, CR-22) | `test_reports::test_compare_p9_runs` | Hand-computed differences, r, gained/lost values, strata | As expected |

## 4. Bugs reproduced on the audited code (915e78b) and resolved

The audited and fixed P9 were run on one synthetic cohort.

| Subject (issue) | Audited 915e78b | Fixed |
|---|---|---|
| QREP: QSM replaced at same path (CR-01) | PASS | fail: QSM content differs from the file P8 used |
| T1MIS: worklist names a different T1 (CR-03) | PASS | fail: T1 content differs from the file P8 segmented |
| DUP: rows differ in `t1_file` (CR-03) | PASS | excluded: worklist rows disagree: `t1_file=...` |
| LEG_AFF: affected legacy label, `--allow-legacy-labels` (CR-04) | PASS | fail: P8 label predates grid-conversion fix |
| FAILMNI: FNIRT-only P1 failure (CR-05) | excluded: on P1 failure list | PASS |
| NOCSF: no CSF reference, `--csf-ref` given (CR-11) | fail: CSF reference unavailable | PASS (CSF set only: `csfref_analysis_pass=False`) |
| OK1: control | PASS | PASS |

**Stale support mask (CR-02).** The QSM was replaced at the same path after the mask was
made:
- The audited `make_support_masks.py` kept the old mask (8,000 voxels) with status "ok".
- The fixed script reported "stale (provenance changed)" and rebuilt it (4,000 voxels,
  equal to the current nonzero region).

**Label geometry (CR-06 context).** Two experiments measured this:
- **Full-size geometry runs:**
  - The original P8 reproduced 85.4%, 85.9% and 93.6% of cortical voxels for the three
    neurologically stored, even-first-dimension T1s (1 mm, 1.0x0.5x0.5 mm, and permuted
    0.5x0.5x1 mm).
  - It reproduced >= 99.9999% for the other four geometries.
  - The revised P8 reproduced >= 99.9997% for all seven.
- **Half-size regression tests (`test_geometry.py`):** these assert cortex agreement below
  99% for exactly the affected geometries, and >= 99.9% for the revised P8 in all
  geometries. `test_p8.py` asserts below 95% for the full-size affected subject.

## 5. Findings made by the tests

- **Quantized QSM and estimated masks:** with integer QSM at 1 ppb steps and SD 3 ppb,
  13.1% of interior voxels are exactly zero. The estimated support recovers deep zeros but
  leaves 5.6% of the outermost 3 mm unsupported. Zero voxels touching the brain edge
  cannot be told apart from background. This is documented in METHODOLOGICAL_DECISIONS D3,
  and the test asserts this measured behaviour.
- **Misleading exclusion reasons:** a missing QSM file was reported as "support mask
  missing", because the mask was looked up first. This is fixed: the first missing input in
  pipeline order is reported.
- **Exclusions outside the selection:** with `--ids-file`, worklist-wide exclusions of
  unselected subjects were listed. This is fixed.

## 6. Validation gaps (not completed here)

| Gap | Needs | How to close it |
|---|---|---|
| `fslmaths -subsamp2` header behaviour (the basis of P8's "expected identity" message; not of the fix) | FSL | `pytest -m fsl` (`test_fslmaths_subsamp2_grid_conversion_is_identity`) |
| Real `flirt -applyxfm` with P8's composed matrix | FSL | `pytest -m fsl` (`test_flirt_applyxfm_with_p8_matrix_reproduces_truth`); `QSM_TEST_REAL_FSL=1 pytest` reruns all P8 tests with real flirt |
| P1 `-usesqform -dof 6 -cost corratio` recovers a rigid misalignment and returns a rigid matrix | FSL | `pytest -m fsl` (`test_rigid_registration_recovers_known_misalignment`); on real data, inspect `P8_q2t_scale_dev` |
| Real P1 (BET, BET retry, Dice) on real scans | FSL + MRI data | Pilot (README) |
| SynthSeg `--parc --keepgeom` labels on real T1 (label set, insula, geometry) | FreeSurfer 8.2 + MRI data | Pilot. P9 checks `synthseg_parcellation_labels.npy` and the LUT at run time |
| Label placement on real anatomy, especially lateral cortex in both hemispheres | Real data | Pilot visual QC (overlays) per geometry stratum |
| Contents of the real `failed_ID.txt` and P1 QC files (stage classification) | Cluster files | Pilot: review `P1_failure_stage` and `P1_failure_evidence` |
| Real QSM storage (float vs integer) and units; whether reconstruction masks exist | Real data / scanner documentation | `support_masks_summary.csv` (`qsm_dtype`, `qsm_quantized`); D3, D4 |
| Exclusions per QC rule by protocol, diagnosis and sex | P8 on the cohort + demographics | `qc_exclusion_report.py` with `--covariates` (D1) |
| Effect of downsampled-T1 registration on ROI values and ribbon placement | FSL + MRI data | Full-resolution comparison run and overlays (README "Sensitivity analyses"; D12) |
| Effect of the support-mask hole limit on ROI values | Real data | `--max-hole-mm3 0` run and `compare_p9_runs.py` (D3) |
| `flock` behaviour on the cluster's network file system | Cluster | Run two P8 jobs for one ID; the second should exit 16 |
| Memory and time at cohort scale (705 subjects, 16 threads) | Cluster | Pilot timing; README memory note |
| Clinical validity of QC thresholds and coverage rules | Real cohort | METHODOLOGICAL_DECISIONS D1-D2 |

No full-cohort run was performed, and no historical result was overwritten.
