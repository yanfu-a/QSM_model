# Clinical_cohort_revised

Revised native-space cortical QSM pipeline (P8 labels, P9 extraction). The original
scripts in `../Clinical_cohort/` are unchanged.

| File | Status |
|---|---|
| `P8_native_aparc_qsm_ALL.sh` | revised |
| `P9_extract_native_qsm_ALL.py` | revised |
| `native_grid.py` | new; FLIRT-coordinate grid conversion shared by P8 and the audit |
| `check_grid_conversion.py` | new; reports which existing P8 labels were built with the faulty conversion |
| `P1_qsm_to_mni_ALL.sh` | unchanged copy, so P8 finds P1 in this folder |

## P8 changes

1. **Grid conversion fix.** The T1 → T1_brain conversion now uses FLIRT's coordinate
   convention, which flips the first voxel axis for neurologically stored images (positive
   determinant). The original formula left the flip out. For such T1s with an even first
   dimension it produced the "−1 mm" shift the old comments called typical. That offset
   every label by one T1 voxel along the first axis (left–right for the sagittal T1s). The
   correct conversion for `fslmaths -subsamp2` is the identity.
2. **Grid-conversion tag.** The QC file now records `grid_conversion=flirt_fsl_v2`, the T1
   orientation, and both the corrected and legacy shifts. P8 skips a subject only if its
   label carries this tag. Untagged (legacy) labels are regenerated automatically.
3. **Cheap regeneration.** The P1 `.mat` is reused only when P1's QC shows it came from the
   same QSM and T1 files. Otherwise P1 is forced to register again: P1 reads `FORCE_RERUN`,
   not `P8_FORCE_RERUN`, and skips silently whenever a `.mat` exists. SynthSeg output is
   reused when it lies on the current T1 grid and its recorded source T1 matches.
   Regenerating a legacy label therefore repeats only steps [3]–[5].
4. **No stale outputs.** The old QC and label are deleted before regeneration, so a crash
   cannot leave a passing QC file next to a new or missing label.
5. **Stricter QC.** A missing P1 registration Dice is now a QC error. NaN QSM background no
   longer counts as coverage.
6. **P1 script lookup.** P8 runs `P1_qsm_to_mni.sh`, or `P1_qsm_to_mni_ALL.sh` if the first
   is absent. The original called only `P1_qsm_to_mni.sh`, which is not in the repository.

New environment variables:
- `P8_PYTHON`, `P8_FREESURFER_HOME`, `P8_P1_SCRIPT`: path overrides.
- `P8_FORCE_P1=1`: re-register with P1. An exported `FORCE_RERUN=1` still forces P1 as
  before.
- `P8_FORCE_SYNTHSEG=1`: rerun SynthSeg.

## P9 changes

- **Label names:** each CSV region name is checked against the FreeSurfer LUT name assigned
  to it by ID, and output names must be unique. Mismatches stop the run unless
  `--allow-label-name-mismatch` is given. The mapping is written to
  `native_label_mapping*.csv`.
- **Duplicate worklist IDs:** an ID whose repeated rows disagree (hospital_id, diagnosis,
  include flag, qsm_file) is excluded instead of taking the first row. Identical repeats are
  merged.
- **Analysis set:** a subject is in the analysis set only if all of these hold:
  - extraction succeeded;
  - it is not on the P1 failure list;
  - P8 QC passed with a numeric Dice;
  - the label carries the grid-conversion tag (`--allow-legacy-labels` overrides this);
  - P8 used the same QSM file (compared by path, inode, or content hash);
  - with `--csf-ref`, the CSF reference is finite.

  `--include-failed` subjects are extracted but never pass. The result is saved as
  `analysis_pass` and `fail_reasons` in `*_qc.csv`, and as `*_analysis.csv` (plus
  `*_csfref_analysis.csv`).
- **Coverage:** valid QSM volume is compared with the ROI's T1-space volume (from P8's
  `t1seg`), so cortex cut off by the QSM slab is counted. Gates: `--min-voxels` (20),
  `--min-mm3` (20) and `--min-coverage` (0.5). New coverage columns: `t1_label_mm3`,
  `qsm_grid_label_mm3`, `valid_qsm_mm3`, `fov_fraction`, `coverage_vs_t1`,
  `roi_value_reported`.
- **Robustness:**
  - Each subject is loaded once; any error becomes that subject's status and the cohort
    run continues.
  - The CSF reference now works with 4D `(x,y,z,1)` QSM.
  - A subject with zero valid voxels keeps its label counts.
- **Bookkeeping:**
  - `*_excluded.csv` lists every subject dropped and why.
  - New QC columns: `repeat_patient`, QSM voxel size, worklist coverage, the P8 QC fields
    (Dice, seg Dice, GM–WM contrast, grid tag, T1 orientation), and a ppm-versus-ppb
    heuristic (`possible_ppm_units`).
- **Reliability:** prints the mean non-homologous left/right correlation as a baseline,
  runs on CSF-referenced values when available, and adds a test–retest correlation for
  patients scanned twice (at least 5 needed).
- **Small fixes:**
  - `--stat mean` now writes only mean columns.
  - `--which` without `--ids-file` is an error, and an unknown column is reported.
  - A missing failure list or SynthSeg label file produces a warning.
  - No crash when no subject passes.
  - `--freesurfer-home` sets the FreeSurfer path.

All original output columns are kept. With the coverage and mm³ gates off
(`--min-coverage 0 --min-mm3 0`), ROI medians are identical to the original and means
differ by less than 1e-5 ppb (float64 accumulation).

**Defaults that change results:** stricter analysis set, `--min-coverage 0.5`,
`--min-mm3 20`, exclusion of conflicting worklist IDs, and rejection of legacy labels.

## Running on the cluster

1. Put these files in the working directory used for P8/P9. `native_grid.py` must sit next
   to P8.
2. Optional: `python check_grid_conversion.py --worklist worklist.csv` writes
   `output/native_grid_audit.csv`. It lists how many existing labels (and results built on
   them) were offset.
3. Rerun the usual P8 batch over all subjects. Legacy labels are regenerated, reusing the
   P1 `.mat` files and SynthSeg output.
4. Fix the conflicting rows in `worklist.csv`. Currently: P191936, P194010, P195951,
   P196693, P196837, P201485, each listed under two hospital_ids, some with different
   diagnoses.
5. Run P9. Its defaults (`--worklist`, `--failed-list`, `--out-dir`) are relative to the
   script's folder, so pass them explicitly if that folder is not your working directory.
   P9 now also reads each subject's T1-space segmentation (about 200 MB for a
   176×512×512 T1), so peak memory is roughly `--jobs` × 0.4 GB. Lower `--jobs` on small
   nodes.

## Testing and limits

FSL and FreeSurfer were not available where these changes were made. They were tested with
stub tools:
- a fake `flirt -applyxfm` built on fslpy's FLIRT conventions;
- fake `fslinfo`, `mri_synthseg` and rigid-mode P1;
- synthetic brains with a neurologically stored 1 mm T1 and a radiological 0.8 mm QSM slab.

Results:
- **Labels:** the revised P8 reproduced ground-truth labels for 100% of brain voxels; the
  original reproduced 85.9% of cortical voxels.
- **Conversion:** `native_grid.flirt_grid_conversion` matched fslpy to within 3e-5 mm on
  300 random oblique grid pairs.
- **Paths:** every P8 skip, reuse, provenance and failure path was exercised, as were P9's
  exclusion and gating rules.

One modelling assumption: the synthetic T1_brain header follows `fslmaths -subsamp2`
(centred subsampling in FSL's radiological voxel order). The fix itself does not depend on
it, because it computes the conversion from the real headers.

Unchanged on purpose: P1 registration, the `QSM_ROI_T1_SUBSAMP` default, the worklist,
and the choice of CSF reference region.
