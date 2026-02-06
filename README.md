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
