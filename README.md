# Nanovaccines

R code to reproduce key statistical analyses and figures for the Nanovaccines manuscript.

> **For peer review / quick validation**, we provide a fully runnable **Fig.7 minimal example** under `example/`.
> The GitHub Actions workflow (`run-example`) runs this example end-to-end and uploads the outputs as an artifact.

## Contents
- Overview
- Repository contents
- System requirements
- Installation guide
- Demo: reproduce Fig.7 (minimal example)
  - A) Run locally
  - B) Run via GitHub Actions (no local setup)
- Expected outputs
- Notes on data availability
- License
- Citation

## Overview
This repository contains analysis scripts and helper functions used to generate figures and statistical results reported in the manuscript.
A minimal reproducible example is provided for Fig.7 to facilitate editor/reviewer evaluation.

## Repository contents
- `example/`  
  Minimal, runnable Fig.7 example (script + required helper functions + minimal dataset).
- `.github/workflows/`  
  CI workflow to run the Fig.7 example and upload outputs as an artifact.
- `R/`  
  Full analysis scripts and helper functions used in the project.
- `data/`  
  Project data (full raw datasets may not be publicly hosted here; see “Notes on data availability” below).

## System requirements
### OS requirements
Tested on:
- GitHub Actions: `ubuntu-latest` (see Actions logs)

### Software requirements
- R (>= 4.1 recommended; GitHub Actions uses the default R from `r-lib/actions/setup-r`)

### Hardware requirements
The Fig.7 minimal example should run on a standard laptop/desktop.

## Installation guide
### Install R
Install R from CRAN: https://cran.r-project.org/

### Install required R packages (Fig.7 minimal example)
The Fig.7 example uses the following R packages (minimum set):
- ggplot2, dplyr, tidyr, ggrepel
- emmeans, boot
- cowplot, gridExtra, knitr, stringr, tidyverse
- bayesplot, rstanarm

Install them in R:
```r
install.packages(c(
  "ggplot2","dplyr","tidyr","ggrepel",
  "emmeans","boot",
  "cowplot","gridExtra","knitr","stringr","tidyverse",
  "bayesplot","rstanarm"
))
```r
## Demo: reproduce Fig.7 (minimal example)
###  A) Run locally

- Clone the repository
     - git clone https://github.com/<YOUR-ORG-OR-USERNAME>/Nanovaccines.git cd Nanovaccines
- Run the Fig.7 example script
     - Rscript example/Fig7.R
- Outputs will be written to:
    - example/res/ (figures + tables)

### B) Run via GitHub Actions (recommended for reviewers)

- Go to the repository → Actions tab
- Select workflow run-example
- Click Run workflow (workflow_dispatch)
- After it finishes, download the artifact:
       - example-output (contains example/res/*)

### Expected outputs
After a successful run, the following files (or similarly named outputs) should appear in `example/res/`:
- PCA plots:
    Fig7_PCA_Aged.png
    Fig7_PCA_Young.png
- Forest plots:
    Fig7_forest_overall.png
    Fig7_forest_Aged.png
    Fig7_forest_Young.png
- Objective plot:
    Fig7_E_objective_plot.png
- Summary tables (CSV):
   Fig7_contrasts_overall.csv
   Fig7_contrasts_by_age.csv
   PCA_combine_Aged.csv, PCA_combine_Young.csv
   PCA_combine_detail_Aged.csv
   PCA_combine_detail_Young.csv



## License

See LICENSE.

## Citation

If you use this code, please cite the Nanovaccines manuscript (and/or add DOI here once available).
