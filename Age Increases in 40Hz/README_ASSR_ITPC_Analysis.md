# README: ASSR_ITPC_Analysis.m

**Statistical analysis pipeline for:**
> "Age-Related Increases in 40Hz Neural Synchrony Are Specific to Typical Development: A Cross-Sectional Study of Autism Spectrum Disorder"
> Beker et al., submitted to *Autism Research*

---

## Requirements

- MATLAB R2019b or later
- Statistics and Machine Learning Toolbox (for `fitlm`, `kruskalwallis`, `lillietest`, `ranksum`, `corr`)
- Data file: `DATA_SB.xlsx` (must be in the MATLAB working directory)

---

## Input Data

**File:** `DATA_SB.xlsx`

| Column | Description |
|--------|-------------|
| `DX` | Diagnosis: `'TD'`, `'ASD'`, or `'PMS'` |
| `Age` | Age in years (continuous) |
| `Gender` | `'M'` or `'F'` |
| `ITPC` | 40Hz inter-trial phase coherence at electrode Cz |
| `Segments` | Number of retained EEG segments (trials) |
| `ID` | Intellectual disability status: `'noID'`, `'ID'`, or `'PMSID'` |
| `IQ` | Full-scale IQ score (numeric; missing for TD and some clinical participants) |

---

## Participant Exclusions (applied at load time)

| Exclusion | Reason |
|-----------|--------|
| PMS participants aged ≥ 18 (n=2) | Insufficient adult PMS representation for trajectory analysis |
| ASD w/ID participants aged ≥ 18 (n=1) | Same rationale |
| One PMS participant (ITPC=0.388, age 13) | Identified as influential observation: Cook's D=0.117, studentized residual=2.36, within-group z=2.88 SD above PMS mean |

**Final analytic sample: N=127** (TD=43, ASD=61, PMS=23)

> **Note:** One TD participant has an empty Age cell in the Excel file. This is handled automatically by `fitlm` (excluded due to missing predictor), but should be explicitly removed if running descriptive statistics manually.

---

## Analytical Approach

Age was modeled as a **continuous log-transformed predictor** (`logAge = log(Age)`) using **ordinary least squares (OLS) linear regression** (`fitlm`). The log transformation linearizes the expected non-linear (decelerating) relationship between age and ITPC.

All regression models include **Sex** as a covariate (coded: F=1, M=0).

Non-parametric **Kruskal-Wallis tests** were used for group comparisons of ITPC collapsed across age. Post-hoc pairwise **Mann-Whitney U** (Wilcoxon rank-sum, `ranksum`) tests with **Bonferroni correction** were applied only for the trial count analysis (the sole significant omnibus KW test with >2 groups).

---

## Models

### Model 0A — KW: All groups, ages collapsed
`kruskalwallis(ITPC, DX)` across TD, ASD (all), PMS (all)

Replicates the approach of prior studies that did not account for age. Included as context to motivate the regression-based approach.

### Model 0B — OLS: All groups, all ages
`ITPC ~ logAge * DX + Sex` (TD as reference group)

Provides a formal test of whether age-related trajectories differ across groups in the full sample. **Note:** Models 0A and 0B are retained in the script for completeness but were not included in the final published manuscript, which focuses on Models 1, 2, and 3S.

### Model 1 — OLS: TD vs ASD w/o ID (all ages)
`ITPC ~ logAge * DX + Sex` (N=80: TD=43, ASD w/o ID=37)

**Primary model.** Tests whether the age-related increase in ITPC differs between TD and iASD without intellectual disability across the full age range.

### Model 2 — OLS: TD vs ASD w/ID vs PMS (under 18 only)
`ITPC ~ logAge * DX + Sex` (N=67: TD=20, ASD w/ID=24, PMS=23)

Tests for group differences in ITPC and age-related trajectories among children and adolescents under 18. Groups are matched on age range. PMS is included as a genetically defined model for E/I-related neural synchrony disruptions in ASD.

### Model 3S — OLS: IQ sensitivity (children under 18 only)
`ITPC ~ logAge * DX + Sex + IQ` (N=39: ASD w/ID=17, PMS=22)

Sensitivity analysis testing whether IQ accounts for ITPC variance beyond age and diagnosis. Restricted to children under 18 with comparable age and cognitive profiles (both groups have severe intellectual disability, mean IQ ≈ 40). Participants with missing IQ excluded.

### Sex Effects Post-hoc
`ITPC ~ logAge * Sex` across all participants (N=131, full dataset including excluded subgroups)

Post-hoc exploratory analysis. Uses the full dataset (including PMS and ASD w/ID adults) since the question is about sex differences broadly. KW tests then compare M vs F within each diagnostic group separately.

### Trial Count Analysis
Kruskal-Wallis test across 5 trial-count subgroups, followed by Bonferroni-corrected pairwise Mann-Whitney U tests. Pearson correlation between segment count and ITPC. Partial correlation controlling for logAge.

---

## Output

The script prints all statistics to the MATLAB console and generates the following figures:

| Figure | Content |
|--------|---------|
| Model 0A | Boxplot + jittered points: ITPC by group, ages collapsed |
| Model 0B | Scatter + regression lines: all groups vs age; boxplot by age group |
| Model 1 | Scatter + regression: TD vs ASD w/o ID; boxplots adults/children |
| Model 2 | Scatter + regression: TD/ASD w/ID/PMS under 18; boxplot |
| Model 3S | Scatter + regression: ASD w/ID vs PMS (IQ controlled); boxplot |
| Sex Effects | Scatter by sex; TD M vs F boxplot; ASD/PMS M vs F boxplot |

---

## Key Results (for reference)

| Model | F / χ² | df | R² | p |
|-------|---------|----|----|---|
| Model 1 | F=6.592 | 4,75 | 0.260 | <0.001 |
| Model 2 | F=1.252 | 6,60 | 0.111 | 0.293 |
| Model 3S | F=0.342 | 5,33 | 0.049 | 0.884 |
| Sex (full model) | F=6.617 | 3,127 | 0.135 | <0.001 |
| Trial count KW | χ²=27.3 | 4 | — | <0.001 |

**Within-group slopes (Model 1):**
- TD: β=0.080, SE=0.029, p=0.008 ✓
- ASD w/o ID: β=0.024, SE=0.020, p=0.241 (ns)

---

## Color Scheme

| Group | Color (RGB) | MATLAB variable |
|-------|-------------|-----------------|
| TD | Blue [0.13 0.47 0.71] | `col_TD` |
| ASD | Red [0.84 0.15 0.16] | `col_ASD` |
| PMS | Green [0.17 0.63 0.17] | `col_PMS` |

---

## Notes

- The script uses `fitlm` for OLS regression (Statistics Toolbox). Output includes coefficient estimates, SE, t-statistics, and p-values.
- `kruskalwallis` is called with `'off'` to suppress the interactive figure.
- Pairwise Mann-Whitney tests use `ranksum` (equivalent to Mann-Whitney U for independent samples).
- The sex model reloads the raw data file to include the two excluded PMS adults and one excluded ASD w/ID adult, giving N=131.
- Models 0A and 0B in the script use the older 3-group structure (TD + all ASD + PMS) which was removed from the final manuscript; they are retained here for reproducibility.
