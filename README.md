
# SFinxSDFA: Spatial Fingerprint Analytics & Spatial Functional Data Analysis

[![R-CMD-check](https://img.shields.io/badge/R--package-SFinxSDFA-blue.svg)](https://github.com/mortezahaji/SFinxSDFA)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![bioRxiv](https://img.shields.io/badge/bioRxiv-10.64898%2F2026.07.15.738836-red.svg)](https://doi.org/10.64898/2026.07.15.738836)

**`SFinxSDFA`** is an R package for reconstructing continuous molecular landscapes, evaluating spatial pathway operators, and fitting spatially varying regression models on spatial transcriptomics (ST) data.

By combining **Spatial Fingerprint Analytics (SFinx)** and **Spatial Functional Data Analysis (SDFA)**, this platform moves beyond gene- or domain-level clustering to model continuous 2D/3D pathway activity and pathway–phenotype crosstalk across complex tissue sections.

---
  
  ## Key Capabilities
  
  * **Spatial Fingerprints (SFinx):** Quantitative, programmable representations of molecular signatures across 2D tissue spaces.
* **Spatial Operators:**
  * **Unary Operators:** Spatial autocorrelation (Moran's $I$), activation scaling ($|Z| \ge \text{threshold}$), and 2D kernel density estimations.
  * **Binary Operators:** Union, Intersection (Jaccard Index), Difference, and Maximum (Dominance) maps between pathway signatures.
* **Spatial Functional Data Analysis (SDFA):** Continuous 2D surface modeling and Generalized Additive Models (GAMs) using tensor product B-splines (`mgcv::bam`)[cite: 3].
* **Multi-Slice 3D Integration:** Evaluates spatial relationships across consecutive tissue sections while accounting for slice-specific random effects along the z-axis[cite: 3].

---

## Workflow Overview