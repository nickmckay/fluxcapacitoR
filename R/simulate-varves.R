#' Simulate synthetic varve thickness series
#'
#' Simulates gamma-distributed annual-layer (varve) thickness series with
#' temporal persistence from either an AR(1) process or long-memory
#' fractional differencing (fBm-like, via the Hurst parameter).
#'
#' Both paths are fully vectorized across the ensemble: AR(1) via a single
#' recursive filter over a matrix of innovations (with a burn-in to reach
#' stationarity), and long-memory via exact circulant-embedding simulation of
#' fractional Gaussian noise (see [fgnEnsemble()]).
#'
#' @param n Number of years (layers) to simulate.
#' @param n.ens Number of ensemble members (columns) to simulate.
#' @param ar1 AR(1) coefficient. Specify `ar1` or `H`, not both.
#' @param H Hurst parameter for long-memory (fractional Gaussian noise)
#'   simulation, `0.5 < H < 1`.
#' @param shape Shape parameter of the gamma distribution of varve
#'   thicknesses.
#' @param mean Mean varve thickness.
#' @param length.out If not `NA`, pad the output matrix with `NA` rows to this
#'   length.
#'
#' @return An `n` (or `length.out`) x `n.ens` matrix of simulated varve
#'   thicknesses.
#' @export
simulateVarves <- function(n, n.ens = 100, ar1 = NULL, H = NULL, shape = 2, mean = 1, length.out = NA) {
  if (!is.null(ar1) && !is.null(H)) {
    stop("Specify either AR(1) coefficient 'ar1' or Hurst parameter 'H', not both.")
  }

  if (!is.null(ar1)) {
    # AR(1) process, vectorized across the whole ensemble in one recursive
    # filter. A burn-in is discarded so each column starts at stationarity
    # (mirroring arima.sim's start-up innovations).
    burn <- 100L
    innov <- matrix(stats::rnorm((n + burn) * n.ens), nrow = n + burn, ncol = n.ens)
    sim <- stats::filter(innov, filter = ar1, method = "recursive")
    sim <- matrix(as.numeric(sim), nrow = n + burn, ncol = n.ens)[-seq_len(burn), , drop = FALSE]
    gamma_series <- gammify(sim, shape = shape, mean = mean)
  } else if (!is.null(H)) {
    # Long memory: exact fractional Gaussian noise (the increments of
    # fractional Brownian motion) via circulant embedding -- vectorized across
    # the whole ensemble, no per-member simulation.
    gamma_series <- gammify(fgnEnsemble(n, n.ens, H), shape = shape, mean = mean)
  } else {
    stop("You must specify either 'ar1' for AR(1) or 'H' for fBM.")
  }

  if (!is.na(length.out)) {
    out <- matrix(NA, nrow = length.out, ncol = n.ens)
    out[seq_len(nrow(gamma_series)), ] <- gamma_series
    return(out)
  } else {
    return(gamma_series)
  }
}

#' Transform data to a gamma distribution
#'
#' Transforms each column of a data matrix to a gamma distribution using the
#' inverse Rosenblatt transform, preserving rank order.
#'
#' Inspired by gaussianize.R, and split.m in normal.m by Van Albada, S.J.,
#' Robinson P.A. (2006), Transformation of arbitrary distributions to the
#' normal distribution with application to EEG test-retest reliability
#' (J Neurosci Meth, doi:10.1016/j.jneumeth.2006.11.004). Modified from
#' matlab code written 26/06/2015 by Julien Emile-Geay (USC); translated to R
#' and jitter option added 7/06/2017 by Nick McKay (NAU).
#'
#' @param X A vector or matrix of data; each column is transformed
#'   independently.
#' @param shape Shape parameter of the target gamma distribution.
#' @param mean Mean of the target gamma distribution.
#' @param jitter Add tiny random noise to break ties?
#'
#' @return A matrix the same size as `X`, gamma-distributed by column.
#' @export
gammify <- function (X, shape = 1.5, mean = 1, jitter = FALSE){
  if(!is.matrix(X)){
    X = as.matrix(X)
  }
  p = NCOL(X)
  n = NROW(X)

  if(jitter){
    #add tiny random numbers to avoid ties
    X = array(stats::rnorm(p*n, mean = 0, sd = stats::sd(as.vector(X))/1e6), c(n,p)) + X
  }

  # The inverse-Rosenblatt target values are qgamma() evaluated at the fixed
  # plotting positions (seq_len(n) - 0.5)/n. Those positions are identical for
  # every column -- only their assignment to data points (the ranks) differs --
  # so qgamma is computed once here and gathered by rank, instead of once per
  # column. All column ranks come from a single C call. Tie ranks (which yield
  # fractional plotting positions) fall back to an explicit qgamma; continuous
  # data essentially never ties.
  qg = stats::qgamma((seq_len(n) - 0.5)/n, shape = shape, rate = shape/mean)
  R = matrixStats::colRanks(X, ties.method = "average", preserveShape = TRUE)

  Xn = matrix(qg[round(R)], n, p)
  ties = R != floor(R)
  if (any(ties)){
    Xn[ties] = stats::qgamma(R[ties]/n - 1/(2*n), shape = shape, rate = shape/mean)
  }

  return(Xn)
}

#' Simulate fractional Gaussian noise by circulant embedding
#'
#' Exact simulation of stationary fractional Gaussian noise (fGn, the
#' increments of fractional Brownian motion) with Hurst parameter `H`, using
#' the Davies-Harte / Wood-Chan circulant embedding. The embedding eigenvalues
#' are computed once and reused for every ensemble member -- all members are
#' drawn with a single column-wise FFT -- and the circulant is padded to a
#' power-of-2 length so the FFT is fast. For `0.5 < H < 1` the embedding is
#' non-negative definite; tiny negative eigenvalues from round-off are clamped
#' to zero as a guard.
#'
#' @param n Number of samples (years) per member.
#' @param n.ens Number of ensemble members (columns).
#' @param H Hurst parameter, `0.5 < H < 1`.
#'
#' @return An `n` x `n.ens` matrix of unit-variance fGn.
#' @export
fgnEnsemble <- function(n, n.ens, H){
  if (n < 3) {                       # too short for embedding; H is irrelevant
    return(matrix(stats::rnorm(n * n.ens), nrow = n, ncol = n.ens))
  }
  m <- as.integer(2^ceiling(log2(2 * (n - 1))))   # power-of-2 embedding length
  half <- m %/% 2L
  k <- 0:half
  rk <- 0.5 * (abs(k + 1)^(2*H) - 2*abs(k)^(2*H) + abs(k - 1)^(2*H))  # fGn autocovariance
  firstRow <- c(rk, rk[half:2])      # first row of the symmetric circulant
  lambda <- Re(stats::fft(firstRow)) # embedding eigenvalues (>= 0 for fGn)
  lambda[lambda < 0] <- 0
  sl <- sqrt(lambda / m)
  Z <- matrix(stats::rnorm(m * n.ens), m, n.ens) +
    1i * matrix(stats::rnorm(m * n.ens), m, n.ens)
  X <- Re(stats::mvfft(sl * Z))
  X[seq_len(n), , drop = FALSE]
}
