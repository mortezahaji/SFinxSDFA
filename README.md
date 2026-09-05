# SFinxSDFA: Spatial Fingerprint Analytics & Spatial Functional Data Analysis

[![R-CMD-check](https://img.shields.io/badge/R--package-SFinxSDFA-blue.svg)](https://github.com/mortezahaji/SFinxSDFA) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![bioRxiv](https://img.shields.io/badge/bioRxiv-10.64898%2F2026.07.15.738836-red.svg)](https://doi.org/10.64898/2026.07.15.738836)

**`SFinxSDFA`** is an R package for reconstructing continuous molecular landscapes, evaluating spatial pathway operators, and fitting spatially varying regression models on spatial transcriptomics (ST) data.

By combining **Spatial Fingerprint Analytics (SFinx)** and **Spatial Functional Data Analysis (SDFA)**, this platform moves beyond gene- or domain-level clustering to model continuous 2D/3D pathway activity and pathway–phenotype crosstalk across complex tissue sections.

------------------------------------------------------------------------

## Key Capabilities

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

# Complete End-to-End Example

The following script runs a synthetic spatial transcriptomics workflow covering:

- **Compatibility patching (ggpubr & ggplot2 bindings).**

- **Synthetic coordinate (Longitude, Latitude) and pathway score generation.**

- **SFinx Unary Operator (SFinx_UO).**

- **SFinx Binary Operators (SFinx_BO).**

- **SDFA Spatial Functional Regression (FOF_fit).**

- **3D Interactive & Static Coefficient Surface Plotting.**

### 1. Load Libraries & Apply Patches

``` r
library(SFinxSDFA)
library(ggplot2)
library(ggpubr)
library(plotly)
library(mgcv)

set.seed(42)

# Patch hardcoded namespace calls for local session compatibility
body_UO <- deparse(body(SFinx_UO))
body_UO_fixed <- gsub("egg::background_image", "ggpubr::background_image", body_UO)
body_UO_fixed <- gsub("viridis::scale_fill_viridis_c", "ggplot2::scale_fill_viridis_c", body_UO_fixed)
body(SFinx_UO) <- parse(text = body_UO_fixed)

body_BO <- deparse(body(SFinx_BO))
body_BO_fixed <- gsub("egg::background_image", "ggpubr::background_image", body_BO)
body_BO_fixed <- gsub("viridis::scale_fill_viridis_c", "ggplot2::scale_fill_viridis_c", body_BO_fixed)
body(SFinx_BO) <- parse(text = body_BO_fixed)

# Create 1x1 blank image matrix for runs without histology images
blank_image <- matrix("#FFFFFF", nrow = 1, ncol = 1)
```

### 2. Generate Synthetic Dataset

``` r
n_spots <- 150

# Spatial coordinates MUST be named Longitude and Latitude for FOF_fit
spatial_coords <- data.frame(
  Longitude = runif(n_spots, min = 0, max = 100),
  Latitude  = runif(n_spots, min = 0, max = 100)
)

# Pathway activity scores
normalized_scores <- data.frame(
  GOBP_REGULATION_OF_NEUTROPHIL_ACTIVATION = rnorm(n_spots, mean = 0, sd = 1.5),
  HP_LUPUS_NEPHRITIS                      = rnorm(n_spots, mean = 0, sd = 1.2)
)

# Merged data frame for SDFA regression modeling
combined_df <- cbind(spatial_coords, normalized_scores)

# Add x and y aliases to spatial_coords
spatial_coords$x <- spatial_coords$Longitude
spatial_coords$y <- spatial_coords$Latitude
```

### 3. Test Unary Operator (SFinx_UO)

``` r
cat("--- Running Unary Operator (SFinx_UO) ---\n")

u_plot <- SFinx_UO(
  Data   = normalized_scores, 
  Finpre = "GOBP_REGULATION_OF_NEUTROPHIL_ACTIVATION", 
  Thresh = 1.0, 
  Pos    = spatial_coords, 
  Image  = blank_image, 
  Title  = "Neutrophil Activation Surface"
)

# Render plot output
print(u_plot)
```

![Unary Operator Surface](https://github.com/mortezahaji/SFinxSDFA/blob/main/images/Unary_Operator_Plot.png)

### 4. Test Binary Operators (SFinx_BO)

``` r
cat("--- Running Binary Operators (SFinx_BO) ---\n")

binary_res <- SFinx_BO(
  Data    = normalized_scores, 
  Finpre1 = "GOBP_REGULATION_OF_NEUTROPHIL_ACTIVATION", 
  Finpre2 = "HP_LUPUS_NEPHRITIS", 
  Thresh  = 1.0, 
  Pos     = spatial_coords, 
  Image   = blank_image
)

# Render multi-panel output
print(binary_res$combined)
```
![Binary Operators Crosstalk](https://github.com/mortezahaji/SFinxSDFA/blob/main/images/Binary_Operators_Plot.png)

### 5. Fit Spatial Functional Regression (FOF_fit)

``` r
cat("--- Fitting SDFA Spatial Model ---\n")

grid_resolution <- 30 # Resolution of spatial prediction surface

sdfa_fit <- FOF_fit(
  Data      = combined_df, 
  Response  = "HP_LUPUS_NEPHRITIS", 
  Predictor = "GOBP_REGULATION_OF_NEUTROPHIL_ACTIVATION", 
  k         = 12,                  # Spline basis dimension
  grid_n    = grid_resolution      # Grid dimension (30x30)
)
```

<a href="https://mortezahaji.github.io/SFinxSDFA/images/3d_surface_plotly.html" target="_blank">
  <img src="images/3d_coefficient_surface.png" width="100%" alt="Click for 3D Interactive Surface" />
</a>

<p align="center">
  👉 <b><a href="https://mortezahaji.github.io/SFinxSDFA/images/3d_surface_plotly.html">Launch Interactive 3D Surface Model in Browser</a></b>
</p>

### 6. Extract Surface Coefficients & Generate 3D Plots

``` r
# Extract the estimated spatio-temporal coefficient matrix beta(s_1, s_2)
beta_matrix <- sdfa_fit$beta_surface 

# Create grid sequence vectors based on grid_n
grid_x <- seq(min(combined_df$Longitude), max(combined_df$Longitude), length.out = grid_resolution)
grid_y <- seq(min(combined_df$Latitude),  max(combined_df$Latitude),  length.out = grid_resolution)

# Interactive 3D Plot (Plotly)
p_3d_interactive <- plot_ly(
  x = ~grid_x, 
  y = ~grid_y, 
  z = ~beta_matrix
) %>% 
  add_surface(
    colorscale = "Viridis",
    contours = list(
      z = list(show = TRUE, usecolormap = TRUE, highlightcolor = "#ff0000", project = list(z = TRUE))
    )
  ) %>% 
  layout(
    title = list(text = "3D Spatially Varying Coefficient Surface β(s)"),
    scene = list(
      xaxis = list(title = "Longitude (s1)"),
      yaxis = list(title = "Latitude (s2)"),
      zaxis = list(title = "Beta Coefficient Value")
    )
  )

print(p_3d_interactive)

# Static 3D Perspective Plot (persp)
ncz <- ncol(beta_matrix)
nrz <- nrow(beta_matrix)
zfacet <- beta_matrix[-1, -1] + beta_matrix[-1, -ncz] + beta_matrix[-nrz, -1] + beta_matrix[-nrz, -ncz]
facetcol <- cut(zfacet, 100)

persp(
  x = grid_x, 
  y = grid_y, 
  z = beta_matrix,
  theta = 35,
  phi = 25,
  expand = 0.6,
  col = terrain.colors(100)[facetcol],
  shade = 0.4,
  ticktype = "detailed",
  xlab = "Longitude (s1)",
  ylab = "Latitude (s2)",
  zlab = "Beta Coefficient",
  main = "3D Functional Coefficient Surface β(s1, s2)"
)
```

### Exporting & Saving Generated Figures

Depending on whether your output is a ggplot2 object, a static 3D graphics device object, or an interactive plotly surface, use the following methods to export high-resolution publication figures:

Saving 2D Plots (SFinx_UO & SFinx_BO) with ggsave

``` r
# Save Unary Operator output
ggplot2::ggsave(
  filename = "Unary_Operator_Plot.png", 
  plot     = u_plot, 
  width    = 10, 
  height   = 6, 
  dpi      = 300
)

# Save Binary Operators multi-panel output
ggplot2::ggsave(
  filename = "Binary_Operators_Plot.png", 
  plot     = binary_res$combined, 
  width    = 12, 
  height   = 10
)
```

Saving Static 3D Perspective Plots

``` r
# Save static 3D perspective plot to PNG
png("3D_Coefficient_Surface.png", width = 2400, height = 2000, res = 300)

persp(
  x = grid_x, 
  y = grid_y, 
  z = beta_matrix,
  theta = 35, phi = 25, expand = 0.6,
  col = terrain.colors(100)[facetcol],
  shade = 0.4, ticktype = "detailed",
  xlab = "Longitude (s1)", ylab = "Latitude (s2)", zlab = "Beta Coefficient",
  main = "3D Functional Coefficient Surface β(s1, s2)"
)

dev.off() # Close graphics device to write file
```

Saving Interactive 3D Plots

``` r
# Option A: Export as standalone interactive HTML page
htmlwidgets::saveWidget(
  widget = p_3d_interactive, 
  file   = "Interactive_3D_Surface.html"
)

# Option B: Export as static image (PNG/PDF/SVG)
# Requires 'kaleido' or 'orca' installed
plotly::save_image(
  p = p_3d_interactive, 
  file = "3D_Surface_Plotly.png", 
  width = 1000, 
  height = 800
)
```

# Citation

### If you use SFinxSDFA in your research, please cite our preprint

Hajihosseini, M., Patino-Martinez, E., Ghosal, R., Kaplan, M. J., & Pyne, S. (2026). An automated platform for spatial functional modeling and fingerprint analysis of tissue molecular landscapes. bioRxiv 2026.07.15.738836; doi: <https://doi.org/10.64898/2026.07.15.738836>
