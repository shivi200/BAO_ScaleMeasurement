# BAO Scale Measurement

Measuring the Baryon Acoustic Oscillation (BAO) scale from SDSS eBOSS DR16 galaxy clustering data, validated against a mock catalog pipeline built with `nbodykit`.

## Overview

This project implements a full pipeline to detect and measure the BAO feature — a subtle imprint from early-universe sound waves that shows up as a preferred separation scale (~150 Mpc) in how galaxies are clustered today.

The pipeline was first validated on a synthetic (mock) galaxy catalog with a known, built-in BAO signal, then applied to real observational data from the SDSS eBOSS DR16 Luminous Red Galaxy (LRG) sample.

## Contents

- `mock_pipeline.py` — Generates a mock galaxy catalog (`LogNormalCatalog`) in a periodic simulation box, computes the power spectrum (`FFTPower`) and two-point correlation function (`SimulationBox2PCF`), and isolates the BAO wiggle by comparing against a smooth no-wiggle reference spectrum.
- `real_pipeline.py` — Loads real eBOSS DR16 LRG clustering data and randoms, computes systematic weights, converts redshifts to comoving distances (Planck15 cosmology), and measures the real two-point correlation function via the Landy-Szalay estimator (`SurveyData2PCF`), comparing it against smooth no-wiggle theory.
- `verify_packages.py` — Environment/dependency check script.
- `r_values.npy`, `xi_values.npy` — Saved output arrays (separation bins and correlation function values) from the real-data run, so results can be replotted without rerunning the ~30 minute pair-counting computation.
- `bao_mock_diagnostics.png` — Mock catalog validation figure (power spectrum, correlation function, and wiggle isolation).
- `bao_comparison_plot.png` — Real-data correlation function vs. no-wiggle theory, with candidate BAO feature regions highlighted.

## Data

Real data used: SDSS DR16 eBOSS LRG NGC clustering catalogs (`eBOSS_LRG_clustering_data-NGC-vDR16.fits` and `eBOSS_LRG_clustering_random-NGC-vDR16.fits`), publicly available from the SDSS DR16 archive. These files are not included in this repository due to size; they can be re-downloaded directly from the SDSS DR16 public archive.

## Environment

- Python (conda environment `bao`)
- `nbodykit`, `astropy`, `numpy` (pinned <1.24 for Corrfunc compatibility), `scipy`, `matplotlib`

## Results Summary

- **Mock validation**: BAO bump successfully recovered at ~105 Mpc/h in both the power spectrum and correlation function, consistent with theoretical expectations.
- **Real data**: A localized feature appears around r≈80-110 Mpc, with a smaller secondary feature near the theoretically expected ~150 Mpc BAO scale. Given the modest sample size (107,500 galaxies), subsampled random catalog, and absence of formal error bars, this result is treated as a preliminary/illustrative finding rather than a confirmed detection.

## Notes

This project was completed as part of a two-week research internship (Project 45).
