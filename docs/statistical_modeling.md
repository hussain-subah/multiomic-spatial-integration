# Statistical Modeling

This repository implements a structured multi-contrast statistical framework to analyze spatially inferred cell-type proportions derived from Bayesian deconvolution (Cell2Location / SpaceJam).

---

## Overview

Cell-type abundances inferred from spatial transcriptomics are:

- Continuous
- Bounded between 0 and 1
- Heteroscedastic

To appropriately model these properties, we use:

```r
glmmTMB::beta_family(link = "logit")
```

Formal statistical testing is performed using beta mixed-effects models.

### Why Beta Models?

Cell-type proportions:
- Continuous
- Bounded between 0 and 1
→ Beta distribution is appropriate

---

## Model
All models are fit using beta mixed-effects regression, accounting for repeated measurements at the ROI/sample level.

### General Model

```r
rel_abundance ~ contrast_variable + (1 | Scan_ID)
```
Where:
  - rel_abundance: inferred cell-type proportion
  - contrast_variable: depends on the biological question
  - (1 | Scan_ID): random intercept for ROI/sample

---
## Data Preprocessing

Before modeling:
  - Proportions are clipped to avoid boundary issues:

```r
rel_abundance <- pmin(pmax(rel_abundance, 1e-6), 1 - 1e-6)
```
- Variables are encoded as factors:
    * disease_status: Control, AD-CAA
    * pathology: AmyloidFree, Amyloid
    * region: spatial region
    * Scan_ID: ROI identifier
      
# Multi-Contrast Framework

To disentangle disease and pathology effects, we implement four complementary contrasts.

## 1. Amyloid Effect (within AD)

 * Goal: isolate the effect of amyloid deposition within diseased tissue. *

Subset:
  - AD-CAA samples only

Model:
``` r
rel_abundance ~ pathology + (1 | Scan_ID)
```
Contrast:
```r
Amyloid - AmyloidFree
```

Interpretation:
- Captures amyloid-associated changes independent of disease baseline

## 2. Disease Effect (amyloid-free regions)

 * Goal: measure disease-associated changes independent of amyloid pathology. *

Subset:
- Only AmyloidFree ROIs

Model:
```r
rel_abundance ~ disease_status + (1 | Scan_ID)
```
Contrast:
```r
AD-CAA - Control
```
Interpretation:
- Reflects disease-driven changes not confounded by amyloid
  
## 3. Maximal Pathology Effect

 * Goal: capture the strongest possible pathological contrast. *

Groups:
- Control_AmyloidFree
- AD_Amyloid

Model:
```r
rel_abundance ~ group2 + (1 | Scan_ID)
```
where : 
```r
group2 = c("Control_AmyloidFree", "AD_Amyloid")
```
Contrast:
```r
AD_Amyloid - Control_AmyloidFree
```
Interpretation:
- Represents extreme phenotype differences
- Useful for identifying high-confidence shifts


## 4. Overall Disease Effect (weighted model)

 * Goal: estimate the true disease effect accounting for heterogeneous amyloid burden. *

We define:
- Control_AmyloidFree
- AD_AmyloidFree
- AD_Amyloid

Model:
```r
rel_abundance ~ group + (1 | Scan_ID)
```
Weighted Contrast
We compute a region-specific weight:
```r
w = proportion of AD ROIs that are Amyloid
```
Then define:
```r
AD_overall = (1 − w) × AD_AmyloidFree + w × AD_Amyloid
```
Final contrast:
```r
AD_overall - Control_AmyloidFree
```
Interpretation:
- Integrates spatial heterogeneity of pathology
- Reflects realistic disease composition
- Avoids bias from uneven amyloid distribution

---

## Effect Size

Effect sizes are reported as:
```r
log2_OR = log2(odds ratio)
```
Derived from model estimates:
- odds.ratio (if available)
- or transformed log-scale coefficients
  
## Post-hoc Testing

Post-hoc contrasts are computed using:
```r
emmeans::contrast()
```
Providing:
- Estimated differences
- Confidence intervals
- p-values
  
## Multiple Testing Correction

P-values are corrected using:
```r
p.adjust(method = "BH")
```
Applied within:
- each region
- each contrast type

---
  
# Output Structure

Each contrast produces:
- means → estimated marginal means
- contrasts → statistical comparisons
- weights → (for overall model only)

A unified summary table is generated:

```text
combined_spatial_contrast_summary.csv
```
Containing:
- celltype
- region
- contrast
- contrast_type
- log2_OR
- p-value
- adjusted p-value

---
### Exploratory vs Formal Statistics

Exploratory analyses (Notebook 05):
- Mann–Whitney tests
- visualization

Formal inference:

- performed exclusively using beta mixed-effects models in R

--- 

# Summary

This framework enables:
- separation of disease and pathology effects
- modeling of heterogeneous amyloid burden
- robust inference across spatial compartments
- biologically interpretable contrasts
  
### Key Strengths

Unlike standard approaches, this pipeline:

* avoids confounding between disease and pathology
* incorporates spatial heterogeneity
* provides multiple complementary biological contrasts

This allows a more accurate representation of disease biology in complex multi-microenvironment systems such as AD/CAA

