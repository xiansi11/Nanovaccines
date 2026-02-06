# Nanovaccines (private)

This repository contains R code and data used to generate figures and statistical results for the Nanovaccines manuscript.

## Repository structure

- `R/figures/`  
  Figure-generation scripts (e.g., `Fig2.R`, `Fig7.R`, `FigS4.R`).
- `R/functions/`  
  Helper functions sourced by figure scripts (e.g., `plot_pca.R`, `plot_combined.R`, `data_prep.R`, `my_fviz_cluster.R`).
- `data/`  
  Input data files (CSV).
- `res/`  
  Generated figures/tables/results data (created locally when running; not required to commit).

> If your files are currently in the repo root, you can move them into these folders (GitHub: edit file name → add path like `R/figures/Fig7.R`).

## Requirements

- R (>= 4.1 recommended)
- R packages used in this project include (non-exhaustive):
  `ggplot2`, `dplyr`, `tidyr`, `lme4`, `lmerTest`, `ggrepel`, `psych`, `ggh4x`, `ComplexHeatmap`, `circlize`, `factoextra`, `devtools`

Install packages (example):
```r
install.packages(c(
  "ggplot2","dplyr","tidyr","lme4","lmerTest","ggrepel","psych","ggh4x",
  "ComplexHeatmap","circlize","factoextra","devtools"
))
