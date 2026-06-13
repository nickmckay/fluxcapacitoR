# CLAUDE.md — fluxcapacitoR

`fluxcapacitoR` does **flux-focused age modeling** for layered sediments: it models
sedimentation/accumulation rate (and its uncertainty) rather than age directly, then
integrates to ages. Extracted June 2026 from the `sedflux` branch of the monolithic
geoChronR — the novel work lived entirely in `R/sedFluxR.R` + `R/syntheticVarves.R`. It is
the **sixth, loosely-coupled** member of the paleogeoscience package family.

Repo: nickmckay/fluxcapacitoR. Branch: `main`.

## Package family (DAG: ens ← lipdViz ← geoChronR; actR, compositeR & fluxcapacitoR on top)

| Repo (`~/GitHub/...`) | GitHub | Branch | Role |
|---|---|---|---|
| ens | nickmckay/ens | main | Ensemble methods + UQ engine; `computeSpectraEns` |
| lipdViz | nickmckay/lipdViz | main | Plotting (`plotModelDistributions`, `plotTimeseriesEns*`, `plotSpectraEns`) |
| geoChronR-chronOnly | nickmckay/geoChronR-chronOnly | main | Age modeling |
| actR | LinkedEarth/actR | refactor | Abrupt-change detection |
| compositeR | nickmckay/compositeR | refactor | Record compositing |
| **fluxcapacitoR** (this repo) | nickmckay/fluxcapacitoR | main | Flux-focused varve age modeling |

**Loosely coupled:** the core has NO hard dependency on ens/geoChronR. It *Suggests* ens +
lipdViz + lipdR for the testing vignette only.

## What lives here

- **Varve age priors** (`R/varve-age-priors.R`):
  - `createVarveAgePriors(DT, ...)` — top-level convenience: low-frequency `scam` monotone
    spline through sampled calibrated dates + high-frequency simulated varve variability.
  - `addVarves(ages, model.depths, yrPerDepth, totalDepth, varveMean, ar1|H, ...)` — **the
    core call**: adds high-frequency varve variability to an age-depth ensemble, scales to
    total depth + total age, scores against the dated distributions. Called repeatedly in
    the fitting loop.
  - `ageProbsDT(DT, ageEstimates, depths)` — log-likelihood of age ensembles vs the LiPD
    `distributionTable`. `sampleAge()`, `removeOutliers()` helpers.
- **Varve simulation** (`R/simulate-varves.R`):
  - `simulateVarves(n, n.ens, ar1=|H=, ...)` — gamma-distributed varve thickness series.
  - `fgnEnsemble(n, n.ens, H)` — exact fractional Gaussian noise via Davies-Harte / Wood-Chan
    circulant embedding (power-of-2 padded; eigenvalues computed once, whole ensemble in one
    `mvfft`).
  - `gammify(X, shape, mean)` — inverse-Rosenblatt transform to a gamma marginal.

## The workflow (the vignette is the spec)

Real usage (from `~/Dropbox/fluxCapacitoR.R`, data `~/Dropbox/RAW_lakes/*.lpd`) is a
Metropolis-Hastings MCMC tuning `ar1` + `varveScalingFactor`: each iteration calls
`addVarves()` and scores with `ageProbsDT()` vs the LiPD `distributionTable`, keeping the
nBest ensemble members. `vignettes/fluxcapacitoR.Rmd` runs a reduced version on a bundled
lake (`inst/extdata/Imandra.Holtzman.2024.lpd`).

## Performance (June 2026)

Per-iteration cost is ~all inside `addVarves()` (the MH loop is negligible orchestration).
On Imandra (12.3k yr, 1053 depths, 100 ens):
- **AR(1): 2.2s → 0.67s (3.3x)** — `qgamma` memoized over the fixed plotting positions +
  one `matrixStats::colRanks` (gammify), `rowsum` binning (replaced dplyr pivot/group), and
  a single recursive `filter` for the whole ensemble (replaced per-member `arima.sim`).
- **Long memory: ~33x** — `fracdiff.sim` (~12.6s/ensemble, non-viable) replaced by exact fGn
  circulant embedding (~0.39s); H-path `addVarves` ~0.9s, on par with AR(1).
- **C/Rcpp is NOT worth it**: remaining hotspots (colRanks/rowsum/rnorm/filter/ageProb) are
  all already C-backed; `qgamma` (the old hotspot) is already C. Next R-level target if
  needed: vectorize the per-date `approx()` likelihood in `ageProbsDT`.

## Gotchas

- **Markdown roxygen** (`Roxygen: list(markdown = TRUE)`).
- **Do NOT depend on the geoChronR monolith** — rbacon/Bchron/oxcAAR/JAGS won't resolve on
  CI (pak: "Can't find package called geoChronR"). Use lipdViz (plots) + ens
  (`computeSpectraEns`) directly. The vignette uses `computeSpectraEns(method = "mtm")`
  (default) so `nuspectral` (an ens Suggest) isn't needed.
- `addVarves()` clamps extrapolated `yrPerDepth >= 0` — `Hmisc::approxExtrap(method =
  "linear")` can go negative at the grid ends, causing age reversals.
- Tests in `tests/testthat/` (synthetic + a real-Imandra `addVarves` monotonicity regression
  test guarded by `skip_if_not_installed("lipdR")`). CI: R CMD check Windows + macOS with the
  vignette built against ens + lipdViz.

## Refactor TODO

- Wrap the MH driver into a package function (`fluxCapacitor()`/`runFluxModel()`); move
  `getSedRatesFromChronEnsemble()` + the `ar1`/`varveScaling` priors from the vignette into
  the package.
- `gammify()` duplicates ens's gaussianize-style util — candidate to share if/when wired onto
  the ens stack.

## Dev

`devtools::load_all()` · `devtools::document()` · `devtools::test()` · `devtools::check()`.
Commit work when complete. Co-author trailer: `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`.
