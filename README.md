---

editor_options: 
  markdown: 
    wrap: 72
---

# SFinxSDFA: Spatial Fingerprint Analytics & Spatial Functional Data Analysis

[![R-CMD-check](https://img.shields.io/badge/R--package-SFinxSDFA-blue.svg)](https://github.com/mortezahaji/SFinxSDFA) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![bioRxiv](https://img.shields.io/badge/bioRxiv-10.64898%2F2026.07.15.738836-red.svg)](https://doi.org/10.64898/2026.07.15.738836)

**`SFinxSDFA`** is an R package for reconstructing continuous molecular landscapes, evaluating spatial pathway operators, and fitting spatially varying regression models on spatial transcriptomics (ST) data.

By combining **Spatial Fingerprint Analytics (SFinx)** and **Spatial Functional Data Analysis (SDFA)**, this platform moves beyond gene- or domain-level clustering to model continuous 2D/3D pathway activity and pathway–phenotype crosstalk across complex tissue sections.

------------------------------------------------------------------------

\## Key Capabilities

- **Spatial Fingerprints (SFinx):** Quantitative, programmable representations of molecular signatures across 2D tissue spaces.
- **Spatial Operators:**
  - **Unary Operators:** Spatial autocorrelation (Moran's $I$), activation scaling ($|Z| \ge \text{threshold}$), and 2D kernel density estimations.
  - **Binary Operators:** Union, Intersection (Jaccard Index), Difference, and Maximum (Dominance) maps between pathway signatures.
- **Spatial Functional Data Analysis (SDFA):** Continuous 2D surface modeling and Generalized Additive Models (GAMs) using tensor product B-splines (`mgcv::bam`).
- **Multi-Slice 3D Integration:** Evaluates spatial relationships across consecutive tissue sections while accounting for slice-specific random effects along the z-axis.

------------------------------------------------------------------------

## Workflow Overview

------------------------------------------------------------------------

## Installation

Install the latest version directly from GitHub:

``` r
# Install remotes if not already available
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install SFinxSDFA
remotes::install_github("mortezahaji/SFinxSDFA")
```

## Quick Start Example

### Evaluate Single Pathway Fingerprints (Unary Operator)

Evaluate spatial autocorrelation (Moran's $I$) and 2D density landscapes of active spots ($|Z| \ge 1.0$) over tissue histology:
``` r
library(SFinxSDFA)

# Run Unary Analysis

u_plot <- SFinx_UO( Data = normalized_scores, Finpre = "GOBP_REGULATION_OF_NEUTROPHIL_ACTIVATION", Thresh = 1.0, Pos = spatial_coords, Image = histology_image, Title = "Neutrophil Activation Surface" )

# Render combined plot

print(u_plot)
```
### Compare Pathway Crosstalk (Binary Operators)
Quantify co-localization, asymmetric activation, and dominance between two pathway fingerprints:

``` r
# Run Binary Operators (Union, Intersection, Difference, Dominance)
binary_res <- SFinx_BO(
  Data    = normalized_scores, 
  Finpre1 = "GOBP_REGULATION_OF_NEUTROPHIL_ACTIVATION", 
  Finpre2 = "HP_LUPUS_NEPHRITIS", 
  Thresh  = 1.0, 
  Pos     = spatial_coords, 
  Image   = histology_image
)

# Display arranged multi-panel plot
print(binary_res$combined)
```

### Spatial Function-on-Function Regression (SDFA)
Fit a spatially varying coefficient GAM with tensor product splines ($k=20$) to model continuous pathway interactions across coordinates:

``` r
# Fit spatial functional regression model
sdfa_fit <- FOF_fit(
  Data      = combined_df, 
  Response  = "HP_LUPUS_NEPHRITIS", 
  Predictor = "GOBP_REGULATION_OF_NEUTROPHIL_ACTIVATION", 
  k         = 20, 
  grid_n    = 50
)

# Extract spatial effect surface and significance mask
beta_matrix <- sdfa_fit$beta_surface
sig_mask    <- sdfa_fit$significance_mask
```