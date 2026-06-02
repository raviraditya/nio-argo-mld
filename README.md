# Mixed Layer Depth Definition Sensitivity and Barrier Layer Structure in the North Indian Ocean

MATLAB code and processed data accompanying the manuscript submitted to *Ocean Science* (Copernicus).
This repository processes Argo profiles over the North Indian Ocean, computes mixed layer
depth (MLD) under multiple density and temperature criteria, derives isothermal layer depth
(ILD) and barrier layer thickness (BLT), and produces the tables and figures used in the paper.

## Authors
- Raviraditya Singh — Department of Civil Engineering, Indian Institute of Technology Bombay, Powai, Mumbai, Maharashtra 400076, India
- Manasa Ranjan Behera — Department of Civil Engineering and Centre for Climate Studies, Indian Institute of Technology Bombay, Powai, Mumbai, Maharashtra 400076, India

Correspondence: raviraditya.singh@gmail.com

## Files in this repository
- `process_argo_mld.m` — reads the Argo `*_prof.nc` files, applies quality control, computes
  per-profile MLD/ILD/BLT, and writes the per-profile table plus summary outputs.
- `make_figures_nio_mld.m` — reads the per-profile table and produces the manuscript figures
  and the error-metrics table, with a sidecar `.txt` of every plotted value.
- `nio_mld_profiles.xlsx` — the processed per-profile table that underlies every figure and
  number in the paper. (All summary tables are regenerated from this file by the code.)

## Method summary
- Reference depth: 10 dbar.
- Density-threshold MLD: delta-sigma = 0.01, 0.02, 0.03, 0.04, 0.05, 0.125 kg m-3.
- Temperature-threshold MLD: delta-T = 0.2, 0.5, 0.8 deg C.
- ILD: delta-T = 0.5 deg C; BLT = ILD minus density MLD at 0.03 kg m-3.
- Reference criterion for bias/RMSE: density MLD at 0.03 kg m-3.
- Density computed with TEOS-10 (Gibbs SeaWater toolbox).
- Region: 0-30 N, 30-100 E; Arabian Sea west of 78 E, Bay of Bengal east of 78 E; N/S split at 12 N.
- Physical caps before summaries: MLD and ILD <= 300 dbar, |BLT| <= 100 dbar.

Criteria follow de Boyer Montegut et al. (2004), Sprintall and Tomczak (1992), and the
literature cited in the manuscript.

## Requirements
- MATLAB (R2019b or later recommended).
- Gibbs SeaWater (GSW) Oceanographic Toolbox for MATLAB (TEOS-10):
  https://www.teos-10.org/software.htm

## Input data
The raw input is the publicly available Argo profile dataset (NetCDF `*_prof.nc` files) for
the North Indian Ocean, available from the Argo Global Data Assembly Centres
(https://argo.ucsd.edu; https://www.seanoe.org/data/00311/42182/). The raw Argo archive is
large and is not redistributed here; only the processed per-profile table is included.

## How to run
1. Install the GSW toolbox and note its folder path.
2. In `process_argo_mld.m`, set the paths at the top: `ARGO_ROOT` (folder of the Argo
   `*_prof.nc` files), `OUT_DIR` (output folder), and `GSW_PATH`. Run it to create
   `nio_mld_profiles.xlsx` and the summary files.
3. In `make_figures_nio_mld.m`, set `IN_TABLE` to `nio_mld_profiles.xlsx` and `OUT_DIR` to a
   figure folder, then run it. Figures are written as PDF and 300-dpi PNG, each with a sidecar
   `.txt` of plotted values.

If you only want to reproduce the figures, you can skip step 2 and point `make_figures_nio_mld.m`
directly at the included `nio_mld_profiles.xlsx`.

## License
MIT License (see `LICENSE`).
