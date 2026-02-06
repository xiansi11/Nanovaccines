# Nanovaccines

This repository contains R code (and limited example data) used to generate figures and statistical results for the Nanovaccines manuscript.

> **Note on data availability:** Full raw/complete datasets may not be included in this repository.  
> We provide either (i) a minimal reproducible subset, or (ii) example data sufficient to run representative analyses and reproduce example outputs.  
> If additional data are required for full reproduction, they will be made available via an appropriate data-sharing mechanism (e.g., controlled access, institutional repository, or upon reasonable request), consistent with ethical and privacy constraints.

---

## Repository structure (current)

- `.github/workflows/`  
  GitHub Actions workflows for automated runs (e.g., running the example pipeline).

- `R/`  
  Main analysis scripts and helper functions (project-wide).

- `data/`  
  Project input data files (CSV).  
  **Tip:** Large raw data should generally not be committed to GitHub. Use a minimal subset or a link to an external archive if needed.

- `example/`  
  A self-contained reproducible example (currently includes `Fig7.R`, helper scripts, and a small input CSV).

- `example/res/`  
  Output folder used by the example run (generated during execution).  
  This folder can exist in the repo with a `.gitkeep`, but typically outputs do **not** need to be committed.

---

## Quick start (run the Fig.7 example locally)

### 1) Clone the repository
```bash
git clone https://github.com/<YOUR-ORG-OR-USERNAME>/<YOUR-REPO>.git
cd <YOUR-REPO>
