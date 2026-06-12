# fluxcapacitoR

<!-- badges: start -->
[![R-CMD-check](https://github.com/nickmckay/fluxcapacitoR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/nickmckay/fluxcapacitoR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Flux-focused age modeling for layered sediments. fluxcapacitoR builds
varve-informed age-ensemble priors that combine low-frequency monotone
spline fits through calibrated date distributions with high-frequency
simulated annual-layer (varve) thickness variability, and evaluates the
resulting age ensembles against the dated distributions.

**Status: experimental.** The API is in flux (sorry).

## Installation

```r
# install.packages("remotes")
remotes::install_github("nickmckay/fluxcapacitoR")
```

## Key functions

- `createVarveAgePriors()` — build an ensemble of varve-informed age-depth
  priors from a LiPD distribution table
- `simulateVarves()` — simulate gamma-distributed varve thickness series
  with AR(1) or long-memory (Hurst) persistence
- `gammify()` — transform data to a gamma distribution via the inverse
  Rosenblatt transform
- `sampleAge()` / `ageProbsDT()` — sample from and evaluate likelihoods
  against calibrated age distributions
- `removeOutliers()` — heuristic age-reversal removal for sampled
  age-depth sequences

## Package family

| Package | Role |
|---|---|
| [ens](https://github.com/nickmckay/ens) | Ensemble methods for time-uncertain data |
| [lipdViz](https://github.com/nickmckay/lipdViz) | Visualization for LiPD data and ensemble results |
| [geoChronR](https://github.com/nickmckay/GeoChronR) | Age modeling and age-uncertain analysis |
| [actR](https://github.com/LinkedEarth/actR) | Abrupt change detection |
| [compositeR](https://github.com/nickmckay/compositeR) | Ensemble compositing |
| fluxcapacitoR | Flux-focused age modeling (this package) |
