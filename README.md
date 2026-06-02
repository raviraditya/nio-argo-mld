# Mixed Layer Depth Definition Sensitivity and Barrier Layer Structure in the North Indian Ocean

MATLAB code accompanying the manuscript submitted to *Ocean Science* (Copernicus).
This repository processes Argo profiles over the North Indian Ocean, computes
mixed layer depth (MLD) under multiple density and temperature criteria, derives
isothermal layer depth (ILD) and barrier layer thickness (BLT), and produces the
summary tables and figures used in the paper.

## Authors
- Raviraditya Singh — Department of Civil Engineering, Indian Institute of Technology Bombay, Powai, Mumbai, Maharashtra 400076, India
- Manasa Ranjan Behera — Department of Civil Engineering and Centre for Climate Studies, Indian Institute of Technology Bombay, Powai, Mumbai, Maharashtra 400076, India

Correspondence: raviraditya.singh@gmail.com

## Overview
The workflow has two stages:

1. **`code/process_argo_mld.m`** — walks every Argo `*_prof.nc` file, applies
   quality control, computes per-profile MLD/ILD/BLT, and writes one flat table
   plus the summary outputs in `data/`.
2. **`code/make_figures_nio_mld.m`** — reads the per-profile table and produces
   the manuscript figures (profile-density map, example profiles, bias/RMSE bars,
   bias-vs-BLT scatter, seasonal panels) and the error-metrics table, with a
   sidecar `.txt` of every plotted value for full reproducibility.

## Method summary
- Reference depth: 10 dbar.
- Density-threshold MLD: delta-sigma = 0.01, 0.02, 0.03, 0.04, 0.05, 0.125 kg m-3.
- Temperature-threshold MLD: delta-T = 0.2, 0.5, 0.8 deg C.
- ILD: delta-T = 0.5 deg C; BLT = ILD minus density MLD at 0.03 kg m-3.
- Reference criterion for bias/RMSE: density MLD at 0.03 kg m-3.
- Density computed with TEOS-10 (Gibbs SeaWater toolbox).
- Region: 0-30 N, 30-100 E; Arabian Sea west of 78 E, Bay of Bengal east of 78 E;
  north/south split at 12 N.
- Physical caps before summaries: MLD and ILD <= 300 dbar, |BLT| <= 100 dbar.

Criteria and definitions follow de Boyer Montegut et al. (2004), Sprintall and
Tomczak (1992), and related literature cited in the manuscript.

## Requirements
- MATLAB (R2019b or later recommended; uses `histcounts2`, `exportgraphics`,
  `datetime`, `readtable`/`writetable`).
- Statistics and Machine Learning Toolbox is not required.
- Gibbs SeaWater (GSW) Oceanographic Toolbox for MATLAB (TEOS-10):
  https://www.teos-10.org/software.htm

## Input data
The raw input is the publicly available Argo profile dataset (NetCDF `*_prof.nc`
files) for the North Indian Ocean. Argo data are freely available from the Argo
Global Data Assembly Centres; see https://argo.ucsd.edu and
https://www.seanoe.org/data/00311/42182/. The raw Argo archive is large and is
NOT redistributed here; only the processed summary outputs needed to reproduce
the figures are included in `data/`.

## How to run
1. Install the GSW toolbox and note its folder path.
2. Open `code/process_argo_mld.m` and edit the paths at the top:
   - `ROOT` — parent folder containing the Argo `YYYY/MM/*_prof.nc` files.
   - `BASE` — output folder for tables and summaries.
   - `GSW_PATH` — folder of the installed GSW toolbox.
3. Run `process_argo_mld.m` (press Run / F5; do not use Run Section, the local
   function needs the whole file). This creates `nio_mld_profiles.mat`,
   `nio_mld_profiles.xlsx`, and the summary files.
   To rebuild only the summaries from a saved run, set `RUN_STAGE1 = false`.
4. Open `code/make_figures_nio_mld.m`, set `CSV` to the path of the per-profile
   table (`nio_mld_profiles.xlsx`) and `OUTD` to the desired figure folder, then
   run it. Figures are written as vector/raster PDF and 300-dpi PNG, each with a
   sidecar `.txt` of plotted values.

## Outputs
- `nio_mld_profiles.mat` / `.xlsx` — per-profile table (not redistributed; large).
- `data/nio_summary_overall.*` — overall min/mean/max/std by ALL / AS / BoB.
- `data/nio_summary_region_season.*` — region x 4-season summary.
- `data/nio_summary_monsoon.*` — pre-/post-monsoon vs other summary.
- `data/nio_summary_bias.txt` — temperature-criterion bias and BLT correlation.
- Figure PDFs/PNGs and per-figure value `.txt` files (created on running the
  figures script).

## License
Code is released under the MIT License (see `LICENSE`).

## How to cite
If you use this code, please cite the accompanying paper and this repository.
See `CITATION.cff`. A Zenodo DOI will be added here on archival.
