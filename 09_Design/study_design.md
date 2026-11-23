# Study design: T1-MRI disease likelihood model across four cohorts

Status: draft for discussion.

## 1. Aim and overall structure

The model is a 7-class 3-D ResNet10 that outputs softmax probabilities for
REF, AD, MSA, PSP, PD, DLB and FTD from a T1-weighted brain MRI in MNI 2 mm
space. It is **developed only in the clinical cohort**, then **locked** and
**evaluated without retraining** in three external cohorts, each answering a
different question:

| Cohort | Role | Question it answers |
|---|---|---|
| Clinical cohort | Development and internal validation | Can the model separate the 7 classes in the clinical population it was built for? |
| ADNI | External validation (disease cohort) | Does the AD output transfer to an independent research cohort, and does it follow AD biology? |
| PPMI | External validation (disease cohort) | Does the PD output transfer to an independent research cohort, and does it follow dopaminergic biology? |
| UK Biobank | External evaluation (population cohort) | How does the model behave in a mostly healthy population: specificity, false positives, and association with prevalent and future disease? |

Design principle: **one locked model, one preprocessing definition, one label
map, cohort-specific estimands.** ADNI, PPMI and UKB are never used for
training, model selection, threshold choice or calibration in the primary
analysis.

Reporting follows TRIPOD+AI (prediction model development and external
validation), with STARD elements for the diagnostic-accuracy analyses and CLAIM
for the imaging AI parts.

## 2. Cohort characteristics that shape the design

| | Clinical cohort | ADNI | PPMI | UK Biobank |
|---|---|---|---|---|
| Setting | Hospital, routine care | Multi-site research study | Multi-site research study | Population volunteer study |
| T1 acquisition | Heterogeneous clinical protocols, often thick-slice or single-plane, need diffusion reconstruction | Protocolised 3D T1 across sites; vendor and field strength vary by phase (to be verified per phase) | Protocolised 3D T1 across sites; vendors and field strengths vary (to be verified) | Uniform 3D T1 on 3T Siemens Skyra, about 1 mm isotropic |
| Image route to model | Reconstruction (Session 1) → MNI 2 mm | FSL pipeline → MNI 2 mm (brain-only) | FSL pipeline (`PPMI_T1_MNI.sh`) → MNI 2 mm (brain-only) | Released `T1_brain_to_MNI` (1 mm) → decimated to 2 mm |
| Classes present | All 7 | REF (CN), AD; MCI has no matching class | REF (HC), PD; Prodromal and SWEDD have no matching class | Mostly REF-like; few cases of each target disease, with MSA, PSP, FTD and DLB expected to be very rare (counts to be confirmed) |
| Label source | Clinical diagnosis (criteria to be documented) | Research diagnosis + amyloid/tau PET, CSF, APOE | Research diagnosis + DAT-SPECT, MDS-UPDRS, alpha-synuclein assay | ICD-10 hospital records, primary care, death registry, self-report |
| Label quality | Moderate; misdiagnosis across parkinsonian and dementia classes likely | High for AD continuum | High for PD; SWEDD are known mislabels | Low to moderate; incomplete case finding |
| Disease prevalence | Enriched (case mix of a referral centre) | Enriched for AD | Enriched for PD | Very low |
| Age range | To be documented | Older (mostly 55–90) | Broad (PD and HC mostly 40–80) | 45–82 at imaging |
| Follow-up | To be documented | Dense longitudinal imaging and cognition | Dense longitudinal imaging and motor scores | Registry-based, several years since imaging |
| Main biases | Referral and spectrum bias; site/protocol | Selected volunteers; strict exclusion of comorbidity | Selected volunteers; de novo PD | Healthy-volunteer bias; few cases |

Consequences for the design:

- **Spectrum and prevalence differ by cohort.** Discrimination (AUROC) can be
  compared across cohorts with care, but predictive values and calibration
  cannot, and must be reported per cohort.
- **Absent classes.** Only the clinical cohort contains all 7 classes. External
  performance is defined only for the classes each cohort contains (section 6).
- **Image domain differs.** Clinical images are reconstructed, while ADNI, PPMI and
  UKB images are native research-grade 3D T1. External results therefore
  test transfer from the reconstructed-clinical domain to the research domain;
  section 4 makes this explicit.
- **Label quality differs.** Biomarker-defined subgroups in ADNI and PPMI serve
  as a higher-quality reference standard for sensitivity analyses.

## 3. Development in the clinical cohort

### 3.1 Data and splits
- Unit of analysis: participant. One primary scan per participant (the earliest
  diagnostic scan); further scans only as secondary data.
- Split: stratified by class, site and age group into training, validation
  and internal test sets at participant level. If the cohort spans several
  years or sites, also hold out the most recent period or one site as a
  temporal or geographic internal test.
- Leakage checks before any training:
  - Participant overlap between training, validation and test sets, including the
    `clinical_csv` subjects that receive reconstructed images.
  - Overlap between the high-resolution training images (`train_csv`) and
    ADNI, PPMI or UKB participants. If the high-resolution training set comes
    from any of these cohorts, those participants must be excluded from that
    cohort's external evaluation, and the paper must say so.

### 3.2 Model training (existing Session 1–2 code, after the fixes noted in review)
- Reconstruction model (Session 1) trained and frozen before classifier
  training. Reconstruction quality reported on paired data (PSNR, SSIM, MAE,
  regional volumes).
- Classifier (Session 2): the `add_ratio` experiment (0.2–1.0 of
  reconstructed clinical images) with nested subsets, 5 seeds per ratio.
- Model selection on the validation set only. Pre-specify the selection
  metric (macro AUROC) and the rule for choosing the final `add_ratio`:
  the ratio with the highest mean validation macro AUROC across the 5 seeds.
- The locked model is the **5-seed ensemble** at the chosen ratio, not the
  best single seed. Ensemble rule: average the pre-softmax logits of the 5
  networks, then apply softmax. This avoids seed-selection bias and usually
  improves calibration.
- Address class imbalance (class-weighted loss or balanced sampling) as a
  pre-specified option, not after looking at external results.

### 3.3 Internal validation
- Metrics on the internal test set:
  - 7-class: macro AUROC (one-vs-rest), balanced accuracy, confusion matrix.
  - Clinical sub-tasks: patients vs REF; PD vs MSA vs PSP; AD vs DLB vs FTD.
  - Calibration: reliability curves, calibration slope and intercept per
    class, Brier score.
- Uncertainty: 1,000 bootstrap resamples at participant level for CIs.
- Internal test set should contain reconstructed clinical images. The current
  code evaluates on high-resolution images only, which does not test the
  intended use.
- Report internal test performance **stratified by image type**:
  reconstructed clinical scans and native high-resolution scans. Where the
  same participant has both (paired scans), also report the within-participant
  difference in scores. This separates reconstruction-to-native domain shift
  from cohort and spectrum shift when external performance differs.

### 3.4 Model lock
Before any external cohort is scored, freeze and record:
- The 5 checkpoint files and their SHA-256 hashes (as
  `inference_ukb_0_6_v2.py` already does for one), and the ensemble rule.
- Class order (`label_id` 0–6 → names).
- Preprocessing definition: template, mask file and hash, z-score rule,
  crop `[9:81, 10:100, 1:77]`.
- Operating thresholds, chosen on the clinical validation set, **each tied to
  the score it applies to** (see section 6 for the score definitions):
  - one-vs-rest threshold on the raw softmax `P(class)` for each class;
  - pairwise threshold on `P(D) / (P(D) + P(REF))` for each binary
    comparison used externally (AD vs REF, PD vs REF), chosen using only the
    REF and D participants of the validation set.

  The two scores are on different scales (the pairwise score is always at
  least as large as `P(D)`), so a threshold for one is never applied to the
  other.
- The statistical analysis plan in this document, with the date.

## 4. Preprocessing harmonisation across cohorts

All cohorts must reach the same model input: a brain-only T1 on the MNI152
2 mm grid (91×109×91), masked by the training template mask, z-scored over the
mask voxels (`ddof=0`) and cropped to 72×90×76.

| Cohort | Route | Required checks |
|---|---|---|
| Clinical | Reconstruction → MNI 2 mm | Same as training by definition |
| ADNI | Raw 3D T1 → FSL pipeline → brain-only MNI 2 mm | Bias-field correction, brain mask QC, registration QC |
| PPMI | Raw 3D T1 → `PPMI_T1_MNI.sh` → brain-only MNI 2 mm | As ADNI; the current script outputs whole-head T1 and must also write a brain-only image |
| UKB | Released `T1_brain_to_MNI` → `[::2]` | Already implemented |

**UKB 1 mm → 2 mm step.** The inference script takes every second voxel
(`[::2]`) only after checking that the image affine times 2 equals the
template affine. For the FSL MNI152 grids this makes each 2 mm voxel centre
coincide exactly with a 1 mm voxel centre, so there is no sub-voxel shift. The
step does not smooth, however, so it keeps more high-frequency noise than
resampling would. The reference for which route is "correct" is how the
training high-resolution images (`T1w_brain_MNI_nonlinear.npy`) were made,
which needs to be documented. The bridging study therefore compares three
routes: decimation (`[::2]`), resampling with anti-aliasing (Gaussian
smoothing, then trilinear or spline), and the native FSL route straight to
2 mm. Pick the route closest to the training preprocessing.

**Bridging study (required before external scoring):** process about 50 UKB
participants from raw T1 through the ADNI/PPMI route and compare with the
released-image route. Acceptance: high agreement of class probabilities (for
example intraclass correlation above 0.9 for each class) and no systematic
shift. If this fails, fix the pipeline before scoring ADNI and PPMI.

QC applies to every scan: header checks, slice thickness and voxel size,
brain-mask overlay, registration similarity to the template, and visual review
of flagged scans. QC is done blind to diagnosis, and failure rates are reported
per diagnostic group.

## 5. Label mapping and reference standards

A fixed file (`label_mapping.yaml`) maps each cohort diagnosis to a model class
or to "no matching class". The table shows the mapping and what each cohort
uses as its reference standard.

| Cohort | Mapped to a class | Not mapped (secondary analyses) | Reference standard |
|---|---|---|---|
| ADNI | CN → REF; AD dementia → AD. ADNI has no clinical diagnoses of other dementias, so no participant is mapped to DLB, FTD, PSP or MSA | MCI (prognostic analysis); SMC | Clinical diagnosis (primary); amyloid and tau status (biomarker evidence); CSF α-synuclein SAA and NfL (biomarker evidence); neuropathology (pathological evidence). Biomarker and pathology groups are never relabelled as clinical diagnoses (section 7.2) |
| PPMI | Confirmed HC → REF; PD → PD; clinical MSA, PSP, DLB, AD, FTD diagnoses → matching class (descriptive) | Prodromal; SWEDD; CBS; other neurological disorders; unaffected carriers | Research diagnosis at the scan date (Appendix A); DaTscan and SAA as biomarker evidence only |
| UKB | No target disease → REF-like; PD (G20) → PD; AD (G30/F00) → AD | MSA, PSP, FTD, DLB (descriptive only) | ICD-10 / primary care / death registry with dates |

Baseline visit is the primary scan in ADNI and PPMI. A participant's diagnosis
at the scan date (and, for sensitivity analysis, the diagnosis confirmed at
last follow-up) defines the label.

## 6. External evaluation: how scores are used when classes are absent

The model's output is a 7-class softmax, but each external cohort contains only
two or three of the classes. Pre-specify the score used for each comparison:

- **Primary score for a binary comparison (for example AD vs REF):**
  the pairwise score `P(AD) / (P(AD) + P(REF))`. It uses the two classes that
  exist in the cohort and is not distorted by probability assigned to classes
  that cannot occur there.
- **Secondary score:** the raw one-vs-rest probability `P(AD)`, which is what a
  user would see in practice.
- **Specificity of the signal:** report all 7 probabilities by diagnostic group
  (a score-by-diagnosis matrix) in every cohort.
- **Score for regression and survival models:** the log-odds
  `log(P(D) / P(REF))`, which equals the logit difference `z_D − z_REF` of
  the network output. It is unbounded and close to normally distributed,
  while raw probabilities in a healthy population pile up near zero.
  Hazard ratios and odds ratios are reported per SD of this log-odds.
  The log-odds is also the logit of the pairwise score, so the two give the
  same ranking and the same AUROC.

Thresholds come from the clinical validation set only, and each is applied only
to the score it was set on (section 3.4). External cohorts report
threshold-free discrimination first, and threshold-based metrics second.

To compute the log-odds without numerical problems, the inference script must
save the 7 pre-softmax logits (ensemble-averaged) in addition to the
probabilities. Very small probabilities stored after softmax can underflow to
zero, and their logarithm is then undefined.

## 7. Cohort-specific analyses

### 7.1 Clinical cohort (internal validation)
As section 3.3. Primary endpoint: macro AUROC on the internal test set.

### 7.2 ADNI

**What ADNI can and cannot provide.** ADNI enrolled people along the AD
continuum and excluded other dementias at entry. Its clinical diagnoses are
cognitively normal (CN), MCI and AD dementia. It gives **no clinical
diagnoses of DLB, FTD, PSP, MSA or vascular dementia**. ADNI therefore
cannot test the model's accuracy for these diagnoses. What it can do is test
whether the non-AD class scores respond to independent **biomarker** or
**pathological** evidence of non-AD disease. These results are associations,
not diagnostic validation, and no subtype is assumed to be validatable in
advance.

**Three kinds of evidence, kept separate in every table and figure:**

| Evidence type | What it shows | Examples in ADNI | Can it be called a diagnosis? |
|---|---|---|---|
| Clinical diagnosis | The clinical syndrome at a visit | CN, MCI, AD dementia | Yes, for these three categories only |
| In vivo biomarker | Presence of a disease process in a living participant | Amyloid (PET or CSF), tau (CSF p-tau, tau PET), CSF α-synuclein SAA, NfL | No. It shows biology, not a clinical syndrome |
| Neuropathology | Post-mortem presence and stage of each pathology | AD neuropathologic change, Lewy body disease, LATE-NC, FTLD subtypes | No. It shows pathology at death, which may differ from the brain at the time of the scan |

**General rules for ADNI**
- The same locked model scores **every T1 scan that passes QC**, at every
  visit. All 7 raw class probabilities and the 7 ensemble logits are kept for
  every scan.
- Each analysis selects scans by a rule fixed in advance: the baseline scan
  for the primary analysis, the scan closest to the fluid sample for
  biomarker analyses, and the last scan before death for pathology analyses.
  Longitudinal analyses use all scans.
- **Biomarkers and pathology are used only to define comparison groups or to
  validate scores.** They are never model inputs, and they are never used to
  choose thresholds, recalibrate or update the model.
- Clinicians may have seen the MRI when making the ADNI diagnoses, so some
  incorporation bias is possible in the primary analysis and is noted as a
  limitation.

**Score notation.** With ensemble logits `z_k`, all scores below are exact
functions of the logits:
- Pairwise score for classes A and B: `P(A) / (P(A) + P(B))`.
- Log-ratio of A to B: `log[P(A) / P(B)] = z_A − z_B`.
- Decomposition of a ratio against REF, for example
  `z_DLB − z_AD = (z_DLB − z_REF) − (z_AD − z_REF)`. Reporting both terms
  shows whether a change in the ratio comes from the numerator, the
  denominator, or both.

#### 7.2.0 Data availability and linkage audit (required before analysis)
Which markers and pathology classifications can be used depends on what the
current ADNI release contains and how many results link to a QC-passed T1. I
could not inspect the ADNI tables (they require approved access, and
`01_Data/ADNI/` in this repository is empty), so the list below is what to
check, not a statement of what exists. Each item is confirmed against the
ADNI data dictionary and tables before the analysis plan is locked.

| Data | Expected content (to be confirmed) | Linkage rule to T1 | Count to report |
|---|---|---|---|
| Diagnostic summary | CN, MCI, AD dementia by visit | Diagnosis at the scan visit | Participants and scans per diagnosis |
| CSF core AD markers | Aβ42, Aβ42/40 where available, p-tau181, t-tau; assay platform and cut-offs | Scan within ±12 months of the lumbar puncture | Participants per amyloid and tau status |
| Amyloid and tau PET | Amyloid status (centiloid or SUVR cut-off); tau PET (flortaucipir) regional SUVR | Scan within ±12 months of PET | Participants per status |
| CSF α-synuclein SAA | Positive, negative, indeterminate; assay and laboratory | Scan within ±12 months of the CSF sample | Participants per SAA result, by amyloid status and cognitive stage |
| CSF or plasma NfL | Continuous concentration | Scan within ±12 months of the sample | Participants with both |
| Other fluid markers (for example GFAP, plasma p-tau217, proteomic panels) | Record only if present | As above | Participants with both |
| Neuropathology | Per-case staging of AD neuropathologic change (ABC), Lewy body disease (type and distribution), LATE-NC stage, hippocampal sclerosis, FTLD (TDP-43 type; tau subtype such as Pick's disease, CBD, PSP or argyrophilic grain disease), and any MSA | Last QC-passed T1 before death; record the MRI-to-death interval | Cases per pathology category, and per interval band (≤2, 2–5, >5 years) |

A non-AD marker that does not exist in ADNI (for example, a validated CSF
assay for TDP-43 or 4-repeat tau) is not replaced by a proxy. NfL indicates
non-specific neurodegeneration and is never used as evidence of a particular
non-AD disease.

**Minimum numbers, fixed in advance.** A group comparison is reported only if
each group has at least 10 participants (proposed; to be agreed before
lock). Below that, individual scores are shown and no group-level claim is
made. The audit result for each analysis below is reported, including
analyses that could not be run.

#### 7.2.1 Primary external validation (unchanged)
- **Comparison:** AD dementia vs CN at baseline, by clinical diagnosis.
- **Score:** pairwise AD score `P(AD) / (P(AD) + P(REF))`.
- **Outcome measures:** AUROC with 95% bootstrap CI, tested against the
  minimum useful AUROC (section 9); sensitivity and specificity at the locked
  AD-vs-REF threshold; calibration of the pairwise score.
- Remains one of the four primary external tests in section 9.

#### 7.2.2 AD biology (biomarker evidence)
- **Biomarker-defined AD vs controls:** amyloid-positive and tau-positive
  impaired participants (MCI or dementia) vs amyloid-negative CN, using the
  pairwise AD score.
- **Within the AD continuum:** association of `z_AD − z_REF` with amyloid
  burden (centiloid), CSF p-tau and tau PET, adjusted for age, sex and
  cognitive stage.
- **Clinical AD dementia without amyloid:** amyloid-negative participants
  with a clinical diagnosis of AD dementia are reported separately. Their
  scores show how the model treats a clinical AD syndrome without AD biology.
- Reported as estimates with CIs.

#### 7.2.3 Lead secondary analysis: Lewy body co-pathology (biomarker evidence)
- **Question:** in people with AD biology, is CSF α-synuclein SAA positivity
  associated with a shift of the model from the AD score toward the DLB
  score?
- **Population:** amyloid-positive participants (PET or CSF, cut-offs fixed in
  advance) with a valid SAA result and a QC-passed T1 within ±12 months of
  the CSF sample.
- **Groups:** SAA-positive vs SAA-negative. Indeterminate results are excluded
  and counted.
- **Comparable cognitive stage:** comparisons are made within cognitive stage
  (MCI or dementia, or CDR global 0.5 vs ≥1), and the model adjusts for
  CDR-SB, because SAA-positive participants are often more impaired.
- **Outcome:** `z_DLB − z_AD`.
- **Decomposition:** the same model is fitted separately for `z_DLB − z_REF`
  and `z_AD − z_REF`, and the raw `P(DLB)` and `P(AD)` are shown by group.
  Reading the result:
  - higher DLB score with unchanged AD score: the model picks up an
    additional Lewy-body-related pattern;
  - lower AD score with unchanged DLB score: SAA positivity is associated
    with a less typical AD pattern, not with a DLB signal;
  - both: report both.
- **Model:** linear regression of each outcome on SAA status, adjusted for
  age, sex, CDR-SB, and field strength or ADNI phase. Sensitivity analyses:
  additional adjustment for APOE ε4 and for tau burden; restriction to the
  dementia stage.
- **Reporting:** group difference with 95% CI, standardized effect size, and
  the AUROC of the log-ratio for SAA-positive vs SAA-negative as a
  descriptive measure.
- **Interpretation:** an association with **Lewy body co-pathology** only.
  SAA positivity does not establish clinical DLB, and the analysis does not
  validate the DLB class as a diagnostic output.

#### 7.2.4 Conditional: α-synuclein without AD biology (biomarker evidence)
- **Population:** cognitively impaired participants (MCI or dementia) who are
  amyloid-negative, split by SAA result.
- **Comparison:** SAA-positive vs SAA-negative on `z_DLB − z_REF` and
  `z_AD − z_REF`, adjusted as in 7.2.3. This group is the closest ADNI comes
  to Lewy body disease without AD pathology, but it is not a clinical DLB
  group.
- **Condition:** run only if the audit finds at least 10 participants in each
  group; otherwise individual scores are shown.

#### 7.2.5 Specificity check: non-specific neurodegeneration (biomarker evidence)
- **Question:** do the FTD, PSP, MSA and DLB scores simply rise with any
  neurodegeneration, rather than with a specific disease?
- **Analysis:** association of each non-AD log-odds (`z_k − z_REF`) with NfL,
  adjusted for age, sex, cognitive stage, amyloid and tau status.
- **Interpretation:** a strong association with NfL after this adjustment
  means the score reflects general neurodegeneration and has limited
  disease specificity. No association does not validate the score.

#### 7.2.6 Exploratory: neuropathology (pathological evidence)
- **Population:** autopsy cases with a QC-passed T1 before death. The last
  scan before death is used. The primary set is cases with an MRI-to-death
  interval of 5 years or less; all cases are shown in a sensitivity analysis.
- **Pathology categories, defined before analysis from the neuropathology
  staging:**
  - AD neuropathologic change (intermediate or high) without other major
    pathology: the reference group;
  - AD neuropathologic change with Lewy body disease (limbic or neocortical);
  - AD neuropathologic change with LATE-NC (stage 2 or higher);
  - FTLD-TDP;
  - FTLD-tau, with subtype recorded (for example Pick's disease, CBD, PSP);
  - MSA, if any case exists;
  - little or no neurodegenerative pathology.
- **Analyses, by available numbers:**
  - Lewy body co-pathology in AD: if both groups reach the minimum, compare
    `z_DLB − z_AD` and its two components between AD with and without Lewy
    body disease. This is the pathological counterpart of 7.2.3. If CSF SAA
    and pathology are both available, report their agreement.
  - LATE-NC in AD: if both groups reach the minimum, describe how
    `z_AD − z_REF` and `z_FTD − z_REF` differ between AD with and without
    LATE-NC. The model has no LATE class, so this is descriptive.
  - FTLD-TDP, FTLD-tau, PSP, CBD and MSA: individual-case plots of all 7
    scores with pathology and MRI-to-death interval. The expected pattern is
    stated in advance (FTD score raised in FTLD; PSP score raised in PSP; MSA
    score raised in MSA), but no subtype is described as validated.
- **Caveats:** autopsy cases are a selected subgroup; long MRI-to-death
  intervals weaken the link between scan and pathology; mixed pathology is
  common, so single-pathology groups will be small.

#### 7.2.7 Other ADNI secondary analyses (unchanged)
- MCI: time to conversion to AD dementia, Cox model with the baseline AD
  log-odds; C-index and time-dependent AUROC.
- Longitudinal: change in the AD score over visits vs change in CDR-SB.
- Score-by-diagnosis matrix of all 7 scores for CN, MCI and AD dementia, and,
  separately, by biomarker group and by pathology category.

### 7.3 PPMI
Scan-level labels follow Appendix A: each scan is labelled from the clinical
records dated on or before its acquisition, with flags for transitional,
revised, conflicting or uncertain labels.

- **Primary:** PD vs HC discrimination of the pairwise PD score at baseline.
- **Secondary:**
  - DAT-deficit-confirmed PD vs HC.
  - Association of the PD score with DAT-SPECT striatal binding and
    MDS-UPDRS part III.
  - SWEDD as a separate group; prodromal participants: conversion to PD.
  - Off-target scores: MSA and PSP probabilities in PD (expected to rise,
    reported as specificity information).
- Note: PPMI PD participants are de novo and early-stage, so their structural
  changes are small. Lower discrimination than in the clinical cohort is
  expected and does not necessarily indicate failure.

### 7.4 UK Biobank
Exclusions: QC failures; participants with a prevalent diagnosis of any of
the 7 target diseases are excluded from the reference group.

- **Primary 1 (PD, cross-sectional):** prevalent PD vs participants
  without any target disease; AUROC of the pairwise PD score, and
  covariate-adjusted Firth logistic regression on the PD log-odds.
- **Primary 2 (AD, prospective):** incident AD in participants free
  of dementia at imaging; cause-specific Cox model with age as the time scale,
  adjusted for sex and site; hazard ratio per SD of the AD log-odds; C-index.
  Check linearity (restricted cubic spline) and proportional hazards
  (Schoenfeld residuals).
- **Competing risks:** death and non-AD dementia compete with AD. The
  cause-specific hazard ratio answers the association question. Any absolute
  risk or cumulative incidence is estimated with Aalen–Johansen curves, and a
  Fine–Gray model is reported as a sensitivity analysis.
- **Secondary:**
  - Incident PD only if enough events occur after the scan; prevalent cases must
    be excluded first. Confirming how prevalent and incident PD cases are
    defined is a prerequisite: if no PD diagnoses fall after the scan, no PD
    survival analysis is planned.
  - False-positive rate per class at the locked threshold among participants
    without target disease, and the positive predictive value implied at UKB
    prevalence.
  - REF score as a continuous brain-health index: association with age,
    cognition and incident dementia.
  - Rare classes (MSA, PSP, FTD, DLB): descriptive only, as percentiles of each
    case's score in the cohort.
- Sensitivity: exclude events in the first 1–2 years; adjust additionally for
  head size, motion and Townsend index; APOE as stratification.

## 8. Confounding, robustness and harmonisation

- Covariates (pre-specified, kept small): age, sex, site or scanner;
  for UKB also head size and motion.
- Age imbalance between cases and controls in ADNI and PPMI: a 3-D CNN on T1
  is partly a brain-age estimator and will exploit any age gap. Report, for
  each external AUROC:
  - the unadjusted AUROC in the full set;
  - a covariate-adjusted AUROC (age and sex), computed by placing each case's
    score within the score distribution of controls with the same age and sex
    (the control-standardised AROC of Janes and Pepe), which keeps the full
    sample;
  - the AUROC in an age- and sex-matched subset, as a check.
- Site effects: mixed models with site as a random effect; ComBat applied to
  logit scores and fitted on controls only, as a sensitivity analysis.
- Field strength (ADNI, PPMI): stratified results for 1.5T and 3T.
- Test-retest: UKB repeat-imaging subset and stable ADNI/PPMI controls with
  repeat scans; ICC of each score.
- Recalibration (secondary only): logistic recalibration of intercept and
  slope, or prior-shift adjustment for the change in class prevalence, fitted
  in a random half of each external cohort and tested in the other half.
  Primary results always use the unmodified locked model.

## 9. Statistics, precision and multiplicity

Primary descriptive endpoint (estimate with 95% bootstrap CI, no hypothesis
test):
- Clinical internal test: macro AUROC, and the clinical sub-task AUROCs.

Primary external hypothesis tests (one family of 4, Holm correction):
1. ADNI: AD vs CN AUROC.
2. PPMI: PD vs HC AUROC.
3. UKB: prevalent PD vs reference AUROC.
4. UKB: incident AD hazard ratio per SD of the AD log-odds (null HR = 1).

For the three AUROC tests, the null is a pre-specified minimum useful value
rather than 0.5, because almost any brain-MRI model beats chance in enriched
disease cohorts. Proposed null: AUROC = 0.70 (to be agreed before lock).
Results are reported mainly as estimates with CIs; the tests decide only
whether the claim of external validity is supported.

All other analyses are secondary or exploratory, reported with FDR
correction where many tests are run. In ADNI, the Lewy body co-pathology
analysis (section 7.2.3) is the lead secondary analysis: its outcome
(`z_DLB − z_AD`) is pre-specified and reported with its own p-value. The
other biomarker and neuropathology analyses (sections 7.2.2 and 7.2.4–7.2.6)
are reported as estimates with CIs or as individual cases, without claims of
diagnostic validation.

Precision in UKB is computed once the case and event counts are confirmed:

- AUROC: expected 95% CI half-width from the Hanley–McNeil variance, for an
  assumed true AUROC (for example 0.75) and the confirmed number of cases.
  With a very large reference group, precision depends almost entirely on the
  number of cases.
- Cox model: detectable hazard ratio per SD from Schoenfeld's formula,
  `HR = exp((z_(1−α/2) + z_power) / √events)`, at 80% power, both unadjusted
  (α = 0.05) and with the Holm threshold for the smallest p-value
  (α = 0.0125).
- A PD Cox analysis is added only if the incident PD count, after excluding
  prevalent cases, supports it.

ADNI and PPMI sample sizes depend on the number of participants with a usable
baseline 3D T1; compute these after the imaging manifest is built.

Missing data: imaging QC failures are not imputed but counted by diagnosis.
Missing covariates use multiple imputation in the regression analyses.

## 10. Governance and reproducibility

- ADNI, PPMI and UKB data stay under their data-use agreements. Only code and
  aggregate results go into the repository; manifests, images and
  participant-level scores stay on approved storage.
- Hard-coded server paths are replaced by configuration files.
- The analysis plan (this document), label mapping, model hash and thresholds
  are dated and archived before external scoring (for example on OSF).
- Each cohort's scores are written with the checkpoint hash and pipeline
  version in their metadata.

## 11. Implementation outline

```
02_Code/
  Session2_Disease_Classification/   # development, internal validation, model lock
  Session3_External_Inference/
    manifests/                       # per-cohort scan lists (not committed)
    preprocess/                      # PPMI_T1_MNI.sh, ADNI equivalent, QC
    infer.py                         # generalised from inference_ukb_0_6_v2.py;
                                     # 5-model ensemble, saves logits + probabilities
  Session4_External_Evaluation/
    label_mapping.yaml
    01_bridging_ukb.py
    02_adni.py                     # primary AD vs CN, 7.2.7
    02a_adni_audit.py              # 7.2.0 availability and linkage counts
    02b_adni_specificity.py        # biomarker and neuropathology analyses
    03_ppmi.py
    04_ukb_cross_sectional.py
    05_ukb_survival.py
    06_robustness.py
    07_figures_tables.py
```

Order of work:
1. Fix development code issues; leakage checks; lock the model.
2. Build ADNI and PPMI manifests; fix and validate preprocessing; bridging study.
3. Score all external cohorts with the locked model.
4. Run the pre-specified analyses; then secondary and exploratory analyses.

## 12. Open questions

1. Clinical cohort: size per class, sites, scan years, diagnostic criteria,
   and whether paired high-resolution scans exist.
2. Source of the high-resolution training images in Session 2 (risk of overlap
   with ADNI, PPMI or UKB), and the exact pipeline that produced
   `T1w_brain_MNI_nonlinear.npy` (needed to choose the UKB 2 mm route).
3. UKB: definition of `*_status` vs `*_status_r`, event dates relative to the
   scan, and DLB extraction.
4. ADNI and PPMI: number of participants with a usable baseline 3D T1, per
   diagnosis and field strength.
5. ADNI audit (section 7.2.0): which CSF, PET and plasma markers are in the
   current release, with assay platforms and cut-offs; number of
   participants with a valid SAA result and a T1 within the window, by
   amyloid status and cognitive stage; which neuropathology fields and
   staging systems are recorded; and the number of autopsy cases with a
   pre-death T1, per pathology category and MRI-to-death interval.
6. Minimum group size for any group comparison (proposed 10 per group), the
   scan-to-sample window (proposed ±12 months), and the maximum MRI-to-death
   interval for the primary pathology set (proposed 5 years).
7. Whether additional cohorts (for example NIFD/FTLDNI or 4RTNI for FTD and
   PSP) can be obtained; otherwise MSA, PSP, DLB and FTD are validated only
   internally.

## Appendix A. PPMI scan-level disease labeling

This appendix defines how each PPMI T1 scan receives a disease label from the
clinical records dated around its acquisition. It is a design only: no
labels have been computed. Field names and codes below come from the PPMI
harmonized data dictionary and code list in `01_Data/PPMI/` (11 Nov 2026).

### A.1 Source files

**Already in `01_Data/PPMI/`:**

| File | Fields used | Role |
|---|---|---|
| `PPMI_selected_scans.csv` | `subject_id` (= `PATNO`), `exam_date`, `exam_key`, `selected_mni2mm_nifti`, `acceptable_scan_count`, `pending_scan_count`, `flagged_scan_count` | One row per MRI exam; defines the scans to label and their QC status |
| `Participant_Status_11Nov2026.csv` | `COHORT`, `ENROLL_DATE`, `ENROLL_STATUS`, `STATUS_DATE`, enrollment subgroup flags (`ENRLSRDC`, `ENRLLRRK2`, `ENRLGBA`, `ENRLSNCA`, `ENRLPRKN`, `ENRLPINK1`, `ENRLOTHGV`, `ENRLHPSM`, `ENRLRBD`, `ENRLNORM`), `DATELIG` | Enrollment cohort and date; fallback label at enrollment |
| Data dictionary and code list | — | Code definitions for all fields below |

**Required from PPMI but not yet uploaded:**

| PPMI table (module) | Fields used | Role |
|---|---|---|
| Primary Research Diagnosis (`PRIMDXPD`) | `PATNO`, `EVENT_ID`, `INFODT`, `PRIMDIAG`, `DXLVL`, `BIORV`, `DATSCAN`, `ASYNSAA`, `NEWDIAGEXP`, `PAG_NAME` | Site investigator's most likely primary diagnosis and its confidence, by visit |
| Clinical Diagnosis (`NEWCLINDX`) | `PATNO`, `EVENT_ID`, `INFODT`, `NEWDIAG`, `DIAGDT`, `DIAGSINLV`, `DIAGCHANGEXP` | Diagnosis changes during follow-up, with the date of the new diagnosis |
| PD Diagnosis History (`PDDXHIST`) | `PATNO`, `PDDXDT`, `SXDT` | Date of PD diagnosis (including for prodromal participants who convert) |
| Visit Status (`VISSTAT`) | `PATNO`, `EVENT_ID`, `EVENT_START`, `INFODT`, `VISCMPLT` | Visit dates, to check the dates of diagnosis records |
| MRI acquisition metadata (`XI_MRI_ACQUISITION`) | `PATNO`, `EVENT_ID`, `MRI_SCAN_DATE`, `MRI_SEQ_T1_WEIGHTED` | Links each scan to a visit ID |
| Conclusion of Study Participation (`CONCL`) | `WDDT`, `WDDEATH` | End of follow-up |
| DaTscan visual read and SBR (`XI_DATSCAN_VI_RESULTS`, `XI_DATSCAN_SBR`) | `DATSCAN_DATE`, `DATSCAN_VISINTRP`, SBR values | Biomarker evidence only (A.6) |
| α-synuclein SAA (`SAA`) | `CLINICAL_EVENT`, `TYPE`, `SAA_Status`, `RUNDATE` | Biomarker evidence only (A.6) |
| Pathology Core (`PATHCORE`) and any released neuropathology | `PATHDOD`, `AUTDT`, `PCSTATUS` | Death date; pathology, if released, as separate evidence |

### A.2 Principles

1. **The label describes clinical status at the scan date.** It is taken from
   the most recent diagnostic record dated on or before the scan, without
   using information from later visits. A separate *final diagnosis* (A.5) is
   kept for sensitivity analyses that do use later information.
2. **Research diagnosis takes precedence over enrollment cohort.** `COHORT`
   records why a participant was enrolled, not their current diagnosis. It is
   used only when no diagnostic record exists on or before the scan.
3. **Labels are clinical.** DaTscan, SAA and pathology never change a label.
   They are stored next to it as biomarker or pathological evidence, as in
   the ADNI design.
4. **Dates are compared at month resolution.** PPMI clinical dates are
   recorded as month and year (`MM/YYYY`), while scan dates have a day. All
   comparisons use calendar months. A record in the same month as the scan
   counts as concurrent (A.4, rule 4).
5. **Every scan gets a label, a source and flags.** Scans are never silently
   dropped; analysis sets are defined later by eligibility flags (A.8).

### A.3 Diagnosis events per participant

For each participant, build a dated list of diagnosis events:

| Event source | Diagnosis | Effective date |
|---|---|---|
| Enrollment | `COHORT` mapped as in A.7 | `ENROLL_DATE` |
| `PRIMDXPD` | `PRIMDIAG` | `INFODT` |
| `NEWCLINDX` | `NEWDIAG` | `DIAGDT` if recorded; otherwise see rule 3 |
| `PDDXHIST` | Idiopathic PD | `PDDXDT` |

When two sources disagree in the same month, the order of precedence is
`NEWCLINDX` > `PRIMDXPD` > `PDDXHIST` > enrollment `COHORT`, and the
disagreement is flagged (rule 8).

### A.4 Rules for assigning the label at the scan date

1. **Current label.** The label is the diagnosis of the most recent event with
   an effective month on or before the scan month.
2. **Enrollment fallback.** If the only event before the scan is enrollment,
   the label comes from `COHORT`, with `label_source = cohort`. This is
   accepted for scans within 6 months of `ENROLL_DATE` (the baseline
   window). Later scans with only a cohort label are flagged
   `unconfirmed_cohort_only`.
3. **Diagnosis changes between visits.** `NEWCLINDX` records a change at a
   visit, and `DIAGDT` gives the date of the new diagnosis.
   - If `DIAGDT` is recorded, the new label applies from that month.
   - If `DIAGDT` is missing and `DIAGSINLV = 1` (the change happened between
     the previous visit and this one), the change month is unknown. Scans
     between the two visits are labelled with the earlier diagnosis and
     flagged `transitional`.
   - If `DIAGDT` is missing and `DIAGSINLV = 0`, the change is dated to the
     visit's `INFODT`.
4. **Same month as a change.** If the scan falls in the same month as a
   diagnosis change, the day order cannot be known. The scan is flagged
   `transitional` and keeps the earlier label.
5. **Prodromal conversion.** For participants enrolled as prodromal (or
   unaffected genetic carriers) who later receive a PD or other diagnosis,
   scans before the conversion month are labelled prodromal and those from
   that month onward with the new diagnosis. Each scan stores the months to
   conversion (negative before, positive after). Conversion is dated by
   `NEWCLINDX.DIAGDT`, then `PDDXHIST.PDDXDT`, then the `INFODT` of the first
   record of the new diagnosis.
6. **Revised diagnoses.** If a later record replaces a diagnosis (for
   example PD revised to MSA or PSP, explained in `NEWDIAGEXP` or
   `DIAGCHANGEXP`), earlier scans keep the label that was current at the
   time. They are also flagged `later_revised`, with the revised diagnosis
   stored. These scans matter for PD vs atypical parkinsonism, because they
   capture early misdiagnosis.
7. **Unstable diagnoses.** A participant whose diagnosis changes more than
   once, or returns to an earlier diagnosis, is flagged `unstable`. Their
   scans keep the time-of-scan label and are removed in a sensitivity
   analysis.
8. **Conflicts.** If the enrollment cohort and the first research diagnosis
   disagree (for example `COHORT` = Healthy Control but `PRIMDIAG` =
   Idiopathic PD, or `COHORT` = PD but `PRIMDIAG` = 17 "No PD nor other
   neurological disorder"), the research diagnosis is used and the scan is
   flagged `conflict`. Conflicting scans are excluded from primary analyses
   and listed for review.
9. **Diagnostic confidence.** `DXLVL` from `PRIMDXPD` is carried with each
   label. Primary analyses require confidence above 50% (`DXLVL` 1 or 2).
   Labels with `DXLVL` 3 or 4 are flagged `low_confidence`. A sensitivity
   analysis uses `DXLVL` = 1 only.
10. **Stale labels.** If the most recent diagnostic record is more than 24
    months before the scan, the scan is flagged `stale`. It keeps its label
    but is excluded from primary analyses.
11. **End of follow-up.** A scan dated after withdrawal (`WDDT`) or death is
    flagged as a data error and excluded.
12. **Missing information.** A scan with no diagnostic record and no
    enrollment cohort, or a code that cannot be mapped, is labelled
    `unknown` and excluded. Code 97 ("Other neurological disorder(s)") is
    labelled `other_neurological` and the free text kept.

### A.5 Final diagnosis (for sensitivity analyses only)

Each participant also gets a final diagnosis: the most recent diagnosis at
last contact, with its date and the follow-up time after each scan. This is
used to:
- define **confirmed healthy controls** (A.8);
- identify scans taken before a later revision or conversion;
- check whether results change when labels use all available information.

The final diagnosis is never used as the primary label.

### A.6 Biomarker and pathology evidence attached to each scan

Stored next to the label, never used to change it:
- **DaTscan:** result closest to the scan (visual read and SBR), with the
  interval in months. Enrollment eligibility (`DATELIG`) is stored as well.
- **α-synuclein SAA:** CSF result closest to the scan (`SAA_Status`:
  positive, negative or inconclusive), with the interval. Inconclusive
  results are kept as inconclusive.
- **Pathology:** death date and, where released, neuropathology diagnosis.

These fields define biomarker-supported subgroups (for example DaT-deficit
PD, SAA-positive prodromal) for secondary analyses, as described in section
7.3.

### A.7 Mapping to labels and to model classes

| PPMI source | Code | Scan label | Model class |
|---|---|---|---|
| `PRIMDIAG` / `NEWDIAG` | 1 Idiopathic PD | `PD` | PD |
| | 17 No PD nor other neurological disorder | `no_disorder` | REF only if a confirmed healthy control (A.8) |
| | 11 Multiple system atrophy | `MSA` | MSA |
| | 14 Progressive supranuclear palsy | `PSP` | PSP |
| | 5 Dementia with Lewy bodies | `DLB` | DLB |
| | 2 Alzheimer's disease | `AD` | AD |
| | 3 Frontotemporal dementia | `FTD` | FTD |
| | 4 Corticobasal syndrome | `CBS` | No matching class |
| | 23, 24, 25 Prodromal (non-motor, motor, synucleinopathy) | `prodromal` | No matching class |
| | 6–10, 12, 13, 15, 16, 18 (other movement or neurological disorders) | `other_neurological` | No matching class |
| | 97 Other neurological disorder(s) | `other_neurological` | No matching class |
| `COHORT` (fallback only) | 1 Parkinson's Disease | `PD` | PD |
| | 2 Healthy Control | `no_disorder` | REF only if confirmed (A.8) |
| | 3 SWEDD | `SWEDD` | No matching class |
| | 4 Prodromal | `prodromal` | No matching class |
| | 7 Genetic Registry - PD | `PD`, flagged `genetic_registry` | PD (sensitivity analyses only) |
| | 8 Genetic Registry - Unaffected | `unaffected_carrier` | No matching class |

Additional attributes stored per scan: enrollment subgroup flags (sporadic,
LRRK2, GBA, SNCA, PRKN, PINK1, other variant, hyposmia, RBD, normosmic), and
the months since PD diagnosis for PD scans.

Notes:
- `SWEDD` stays `SWEDD` until a later diagnosis replaces it, which is common.
- The genetic registry cohorts were recruited differently from the de novo
  PD and healthy-control cohorts, so they are analysed separately.
- MSA, PSP, DLB, AD, FTD and CBS labels in PPMI are clinical diagnoses made
  during follow-up and are expected to be rare. They are reported
  descriptively unless the minimum group size (10) is reached.

### A.8 Eligibility flags and analysis sets

| Analysis set | Scans included |
|---|---|
| **Primary PD vs HC (section 7.3)** | One scan per participant: the earliest QC-passed scan within the baseline window. PD = label `PD` from `COHORT` 1 (not genetic registry), `DXLVL` 1–2, no `conflict`, `transitional` or `stale` flag. HC = `COHORT` 2 with label `no_disorder` at the scan **and** no later neurological diagnosis at last contact (confirmed healthy control), with at least 12 months of follow-up after the scan |
| DaT-deficit PD vs HC | Primary set, restricted to PD with an abnormal DaTscan |
| Atypical parkinsonism (descriptive) | Scans labelled MSA, PSP, CBS or DLB, plus `later_revised` PD scans whose final diagnosis is one of these |
| SWEDD | Scans labelled SWEDD, with their final diagnosis |
| Prodromal and conversion | All scans labelled `prodromal`, with months to conversion and SAA status |
| Longitudinal | All QC-passed scans with a time-of-scan label; participants flagged `unstable` removed in a sensitivity analysis |
| Sensitivity: final diagnosis | Primary set relabelled with the final diagnosis |
| Sensitivity: high confidence | Primary set with `DXLVL` = 1 only |
| Sensitivity: genetic | Add `genetic_registry` PD and genetic-variant PD |

Scans with `selected_mni2mm_nifti` missing, or with only pending or flagged
QC results, receive a label but are not eligible for any analysis until QC
passes.

### A.9 Output: one row per scan

`PATNO`, `exam_key`, scan date, scan month, visit ID (if linked), enrollment
cohort and subgroup flags, `label_at_scan`, `label_source`, source record
date, months since source record, `DXLVL`, `final_diagnosis` and its date,
months of follow-up after the scan, months to conversion, flags
(`transitional`, `later_revised`, `unstable`, `conflict`, `low_confidence`,
`stale`, `unconfirmed_cohort_only`, `genetic_registry`), model class, DaTscan
result and interval, SAA status and interval, QC status, analysis-set
eligibility, and the reason for any exclusion. The table holds
participant-level data and stays on approved storage, not in this
repository.

### A.10 Values to confirm before lock

- Baseline window: within 6 months of `ENROLL_DATE`.
- Maximum age of a diagnostic record before the scan: 24 months.
- Minimum follow-up for a confirmed healthy control: 12 months.
- Confidence threshold: `DXLVL` 1–2 in primary analyses.
- Same-month rule: scans in the month of a diagnosis change are transitional.
- The diagnosis tables listed in A.1 need to be downloaded from PPMI. Until
  then, only the enrollment fallback (rule 2) can be applied.
