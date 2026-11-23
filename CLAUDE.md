# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

Research code for a T1-MRI "disease likelihood" model: a 7-class 3-D ResNet10 that
scores a brain T1 scan for REF (reference/healthy), AD, MSA, PSP, PD, DLB and FTD.
The model is developed in a hospital clinical cohort and applied, locked, to ADNI,
PPMI and UK Biobank. The study design (cohort roles, estimands, analysis plan) is in
`09_Design/study_design.md`; read it before changing anything that affects model inputs,
scores or evaluation.

Branches differ in content: `new` holds the T1 model code described here; `main`
holds the clinical-cohort QSM preprocessing pipeline (`Clinical_cohort/`).

There is no build system, test suite, linter or dependency file. Scripts are run
directly on a GPU server with FSL; nothing here can be run end-to-end without the
(non-public) data.

## Running the code

Every script imports sibling modules by bare name (`from resnet import ...`,
`from model_utils import ...`), so run each one from inside its own directory.
Data paths are anonymised placeholders (`/path_to/`) that must be filled in, either
in the `config` dict / argparse defaults or, for `model_utils.py`, the module-level
`mask = nib.load('/path_to/')` that runs at import time.

```bash
# Session 1: fine-tune the diffusion super-resolution model (DDP, NCCL)
cd 02_Code/Session1_T1_Reconstruction && torchrun --nproc_per_node=<N> finetune_T1_sr.py
# Session 1: sample reconstructions (DDIM, CFG, averages --average_num samples)
cd 02_Code/Session1_T1_Reconstruction && python evaluation_T1_sr.py --ckpt <ckpt>
# Session 2: train/evaluate the classifier over add_ratio 0.2..1.0 (paths in the config dict in __main__)
cd 02_Code/Session2_Disease_Classification && python train_and_inference_add_0_6.py
# Session 3: UKB inference; --limit N is the smoke test; output CSV is append-only and resumable
cd 02_Code/Session3_UKB_Inference && python inference_ukb_0_6_v2.py --limit 10
# Session 3: PPMI registration to MNI 2 mm with FSL, batched by file index
bash 02_Code/Session3_PPMI/PPMI_T1_MNI.sh <START_INDEX> <END_INDEX> [JOBS]
# Session 3: PPMI scan-time diagnosis labels (output dir must be outside the checkout)
bash 02_Code/Session3_PPMI/run_scan_diagnosis_20261123.sh "" <output_dir>
# Session 3: supplementary nearest-diagnosis match (may use records after the scan)
cd 02_Code/Session3_PPMI && python3 build_nearest_diagnosis.py --scan-list ../../01_Data/PPMI/PPMI_selected_scans.csv \
    --research-dx ../../01_Data/PPMI/Primary_Research_Diagnosis_23Nov2026.csv \
    --clinical-dx ../../01_Data/PPMI/Clinical_Diagnosis_23Nov2026.csv --out-dir <output_dir>
```

The PPMI matching scripts are the only code with tests (stdlib `unittest`):

```bash
cd 02_Code/Session3_PPMI && python3 -m unittest test_scan_diagnosis test_nearest_diagnosis
# a single test:
cd 02_Code/Session3_PPMI && python3 -m unittest test_nearest_diagnosis.NearestDiagnosisTests.test_tie_different_codes_is_unassigned
```

PPMI diagnosis dates are month/year only while scan dates are exact days, so both
matching scripts compare at calendar-month resolution. `build_scan_diagnosis.py`
gives the time-of-scan label (no look-ahead); `build_nearest_diagnosis.py` imports
its parsing and code mapping and may use later records. Only the `PRIMDXPD` page of
the Primary Research Diagnosis export counts by default (`--research-page` to add
others). The labeling design is Appendix A of `09_Design/study_design.md`.

Python dependencies (from imports): torch, torchvision, nibabel, numpy, pandas,
scipy, scikit-learn, tqdm; `Optional_Models/model_multimodal.py` also needs
`rotary_embedding_torch` and `einops`. FSL (`bet`, `flirt`, `fnirt`, `applywarp`) is
needed for PPMI preprocessing.

## Pipeline architecture

1. **Session 1, reconstruction** (`Session1_T1_Reconstruction`): a 1-D Diffusion
   Transformer (`models.py`, `DiT_UKB`, built on the guided-diffusion code in
   `diffusion/`) operates on the vector of brain voxels inside the MNI152 2 mm mask
   (228,453 voxels), not on a 3-D grid. `PatchEmbed_1d` pads and patches this vector.
   `DiT_UKB_finetune` replaces the patch embedder with a 2-channel one so that the
   clinical (low-resolution, single-plane) T1 is concatenated with the noisy target,
   and conditions on age, sex and modality embeddings. Sampling uses classifier-free
   guidance; null conditions are age 0, sex 2 and modality `--modality_num`. Outputs
   are saved per subject as `<id>_<plane>/T1W_brain_nonlinear_recon_average_<n>.npy` (Session 2 reads
   `average_5`).
2. **Session 2, classification** (`Session2_Disease_Classification`): ResNet10 on
   one-channel 3-D volumes. `Multi_disease_dataset_with_clinical` concatenates
   high-resolution training images with a fraction (`add_ratio`) of Session 1
   reconstructions. The model is selected on validation AUROC and written to
   `add_ratio_<r>/best_model.pth`.
3. **Session 3, external inference**: the locked classifier is applied to UKB
   (`inference_ukb_0_6_v2.py`) and, after FSL registration, to PPMI.
   `Session3_UKB_Inference/resnet.py` is an identical copy of the Session 2 file;
   keep them in sync.

## Invariants that must match across sessions

All classifier inputs, in training and in every external cohort, go through the
same transform. A mismatch produces wrong scores without raising an error.
- Grid: MNI152 2 mm, 91×109×91, with the template brain mask (`MNI152_T1_2mm_brain`,
  voxels `!= 0`). Images are brain-only (skull-stripped), not whole-head.
- Intensities: z-score over the mask voxels with population std (`ddof=0`,
  numpy), scattered back into the mask volume, then cropped
  `[9:81, 10:100, 1:77]` to shape (72, 90, 76). Note that Session 1 z-scores with
  `torch.std` (unbiased), so reconstruction and classifier normalisation differ slightly.
- Voxel order: `mask.nonzero()` / boolean mask indexing (C order). Cached `.npy`
  files store only the masked voxel vector.
- Class order comes from `label_id` 0–6 with `disease_label` names in
  `classification_train.csv`. The UKB script checks this and refuses to run on a mismatch.
- UKB 1 mm images are reduced to 2 mm with `[::2]` only after checking that the affine
  matches the template affine. Other grids are resampled with `affine_transform`.
- The UKB output has a `.meta.json` with the checkpoint SHA-256 and class names. A
  rerun against a different checkpoint is refused, so keep this provenance in any new
  inference script.

The dataset classes silently create `.npy` caches in `image_dir_add` when a file
is missing. Account for this when changing paths or preprocessing.

## Data handling

`01_Data/` contains participant-level ADNI and PPMI exports, and
`03_Result/Session3_cohort/` contains participant-level matching results
(`*_PRIVATE.csv`); these were committed by the repository owner. They are covered
by the ADNI and PPMI data-use agreements: do not copy them elsewhere, and do not
add further participant-level data, image paths or per-participant counts unless the
owner asks. Hard-coded server paths have been replaced with `/path_to/` in earlier
commits; keep that convention in code (`02_Code/Session3_PPMI/PPMI_T1_MNI.sh`,
`01_Data/ADNI/ADNI_T1_registration_QC.xlsx` and `01_Data/PPMI/PPMI_selected_scans.csv`
still contain real paths).
