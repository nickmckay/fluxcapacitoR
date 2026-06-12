#' Simulate synthetic varve thickness series
#'
#' Simulates gamma-distributed annual-layer (varve) thickness series with
#' temporal persistence from either an AR(1) process or long-memory
#' fractional differencing (fBm-like, via the Hurst parameter).
#'
#' @param n Number of years (layers) to simulate.
#' @param n.ens Number of ensemble members (columns) to simulate.
#' @param ar1 AR(1) coefficient. Specify `ar1` or `H`, not both.
#' @param H Hurst parameter for long-memory simulation (much slower than
#'   AR(1)).
#' @param shape Shape parameter of the gamma distribution of varve
#'   thicknesses.
#' @param mean Mean varve thickness.
#' @param length.out If not `NA`, pad the output matrix with `NA` rows to this
#'   length.
#'
#' @return An `n` (or `length.out`) x `n.ens` matrix of simulated varve
#'   thicknesses.
#' @export
simulateVarves <- function(n, n.ens = 100, ar1 = NULL, H = NULL, shape = 2, mean = 1,length.out = NA) {
  if (!is.null(ar1) && !is.null(H)) {
    stop("Specify either AR(1) coefficient 'ar1' or Hurst parameter 'H', not both.")
  }

  if (!is.null(ar1)) {
    # AR(1) Process with given ar1
    gamma_series <-  purrr::map(rep(n,times = n.ens),\(x) stats::arima.sim(n = x, list(ar = ar1,ma = 0), innov = stats::rnorm(n))) |>
      purrr::list_c()  |>
      matrix(nrow = n, ncol = n.ens,byrow = FALSE) |>
      gammify(shape = shape, mean = mean)  # Convert to Gamma
  } else if (!is.null(H)) {
    # Fractal Brownian Motion (fBM) using fractional differencing
    # This is MUCH slower...
    gamma_series <- purrr::map(rep(n,times = n.ens),\(x) fracdiff::fracdiff.sim(n = x,d = H - 0.5)$series) |>
      purrr::list_c()  |>
      matrix(nrow = n, ncol = n.ens,byrow = FALSE) |>
      gammify(shape = shape, mean = mean)  # Convert to Gamma
  } else {
    stop("You must specify either 'ar1' for AR(1) or 'H' for fBM.")
  }

  if(!is.na(length.out)){
    out <- matrix(NA, nrow = length.out,ncol = n.ens)
    out[seq_len(nrow(gamma_series)), ] <- gamma_series
    return(out)
  }else{
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
gammify <- function (X,shape = 1.5, mean = 1,jitter=FALSE){
  if(!is.matrix(X)){
    X=as.matrix(X)
  }
  p=NCOL(X)
  n=NROW(X)

  if(jitter){
    #add tiny random numbers to avoid ties
    X=array(stats::rnorm(p*n,mean=0,sd=stats::sd(as.vector(X))/1e6),c(n,p))+X
  }

  Xn    = matrix(0,n,p);
  for (j in 1:p){
    # Sort the data in ascending order and retain permutation indices
    R=rank(X[,j])
    # The cumulative distribution function
    CDF = R/n - 1/(2*n);
    # Apply the inverse Rosenblatt transformation
    Xn[,j] = stats::qgamma(CDF,shape = shape,rate = shape/mean)  # Xn is now gamma distributed
  }

  return(Xn)
}
