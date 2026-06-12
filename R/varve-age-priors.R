#' Sample one age from a LiPD distribution table object
#'
#' Draws a single age from the calibrated probability distribution stored in
#' one entry of a LiPD `distributionTable` (e.g.
#' `L$chronData[[1]]$model[[1]]$distributionTable[[i]]`), by inverting the
#' empirical CDF.
#'
#' @param dto A single LiPD distribution table object: a list with
#'   `age$values` and `probabilityDensity$values` components.
#'
#' @return A single sampled age.
#' @export
sampleAge <- function(dto){
  #calculate the cumulative probabilites, after adding in a tiny slope to avoid ties, and normalizing to sum to 1.
  cdf <- (cumsum(dto$probabilityDensity$values) + seq(0,0.00000001,length.out = length(dto$probabilityDensity$values)))/sum(dto$probabilityDensity$values)
  sampledAge <- stats::approx(x = cdf, y = dto$age$values, xout = stats::runif(1))$y
  return(sampledAge)
}

#' @noRd
fitMonotonicSpline <- function(age,depth,model.depths,kFrac){
  k <- min(max(round(kFrac * sum(is.finite(age))),5),sum(is.finite(age)))

  scam_fit <- scam::scam(age ~ s(depth, bs = "mpi",k = k), sp = 1) |>
    stats::predict(newdata = data.frame(depth = model.depths))

  return(scam_fit)
}

#' Remove age reversals from a sampled age-depth sequence
#'
#' Heuristically removes ages until the sequence is monotonically increasing
#' with depth, then tries to re-add removed ages (in random order) wherever
#' they don't reintroduce a reversal. Removed ages are returned as `NA`.
#'
#' @param age A vector of sampled ages.
#' @param depth The corresponding depths.
#'
#' @return The age vector, in the original order, with reversal-causing ages
#'   set to `NA`.
#' @export
removeOutliers <- function(age,depth) {
  #first sort by depth
  di <- sort(depth,index.return = TRUE)$ix

  age <- age[di]

  ageOrig <- age

  removeWeights <- rep(1,times = length(age))

  while(any(diff(stats::na.omit(age)) <= 0,na.rm = TRUE)){
    tr <- sample(seq_along(age),prob = removeWeights,size = 1)
    age[tr] <- NA
    removeWeights[tr] <- 0
  }

  #now let's try to add as many back in as we can.
  removed <- which(is.na(age))
  shuffled <- sample(removed)

  for(i in seq_along(shuffled)){
    age[shuffled[i]] <- ageOrig[shuffled[i]]
    if(any(diff(stats::na.omit(age)) <= 0,na.rm = TRUE)){
      #reversal! take it back out
      age[shuffled[i]] <- NA
    }
  }

  return(age[order(di)])
}

#' @noRd
ageProb <- function(dto,ageToTest){
  pdf <- dto$probabilityDensity$values
  prob <- stats::approx(y = pdf, x = dto$age$values, xout = ageToTest)$y
  min.prob <- min(pdf[pdf > 0])

  if(is.na(prob)){#then its off the scale, interpolate between 0 and the lowest point

    if(ageToTest < min(dto$age$values)){#too young
      prob <- stats::approx(x = c(0,dto$age$values[which.min(dto$age$values)]), y = c(0,min.prob), xout = ageToTest)$y
    }else{
      prob <- stats::approx(x = c(10* dto$age$values[which.min(dto$age$values)],dto$age$values[which.min(dto$age$values)]), y = c(0,min.prob), xout = ageToTest)$y
    }
  }

  logProb <- log(prob)
  if(!is.finite(logProb)){
    logProb <- log(min.prob)
  }
  return(logProb)
}

#' Create varve-informed age-ensemble priors
#'
#' Builds an ensemble of age-depth priors from a LiPD distribution table by
#' combining (1) low-frequency structure, from monotone spline fits through
#' ages sampled from the calibrated date distributions, with (2)
#' high-frequency structure, from simulated annual-layer (varve) thickness
#' variability with AR(1) or long-memory persistence.
#'
#' @param DT A LiPD distribution table: a list of distribution table objects,
#'   each with `age`, `probabilityDensity`, and `depth` components.
#' @param model.depths Depths at which to evaluate the age model. If `NA`
#'   (default), a regular sequence from 0 to the deepest dated sample in steps
#'   of `model.depth.step`.
#' @param model.depth.step Depth increment used to construct `model.depths`
#'   when they're not supplied.
#' @param varveScalingFactor Prior on the ratio of mean varve thickness to the
#'   mean spline-derived sedimentation rate.
#' @param H Prior Hurst parameter for long-memory varve simulation. Specify
#'   `H` or `ar1`, not both.
#' @param ar1 Prior AR(1) coefficient for varve simulation.
#' @param kFrac Fraction of the number of dates to use as the spline basis
#'   dimension `k` (bounded between 5 and the number of dates).
#' @param n.ms.ens Number of monotone-spline ensemble members.
#' @param n.varve.ens Number of varve ensemble members (currently must equal
#'   `n.ms.ens`).
#' @param heuristicOutlierRemoval Apply [removeOutliers()] to sampled ages?
#' @param outlierRemovedFraction Fraction of ensemble members that get the
#'   outlier-removal treatment.
#' @param progress Show progress bars and messages?
#'
#' @return A list with `agePriors` (matrix of ensemble age priors),
#'   `ageDepths` (the depths they're evaluated at), and `varvedPriorLogObj`
#'   (log-likelihood of each ensemble member given the dated distributions).
#' @export
createVarveAgePriors <- function(DT,
                                 model.depths  = NA,
                                 model.depth.step = 1,
                                 varveScalingFactor = 30,#prior
                                 H = NULL, #prior
                                 ar1 = NULL, #prior
                                 kFrac = 1/3, #prior
                                 n.ms.ens = 100,
                                 n.varve.ens = 100,
                                 heuristicOutlierRemoval = TRUE,
                                 outlierRemovedFraction = .95,
                                 progress = TRUE){


  sampledAgeEns <- purrr::map(1:n.ms.ens,\(x) purrr::map_dbl(DT,sampleAge)) |>
    purrr::list_c() |>
    matrix(nrow = length(DT),ncol = n.ms.ens)


  sampleDepths <- purrr::map_dbl(DT,"depth")

  if(heuristicOutlierRemoval){
    filteredSampledAgeEns <- apply(sampledAgeEns,2, removeOutliers, sampleDepths)
    wc <- sample(seq_len(ncol(sampledAgeEns)),size = round(ncol(sampledAgeEns) * outlierRemovedFraction))
    sampledAgeEns[,wc] <- filteredSampledAgeEns[,wc]
  }

  if(all(is.na(model.depths))){
    model.depths <- seq(0,ceiling(max(sampleDepths)),by = model.depth.step)
  }

  if(progress){message(crayon::blue("Estimating low-frequency variability in the model..."))}

  if(progress){
    aeig <- pbapply::pbapply(X = sampledAgeEns,
                             MARGIN = 2,
                             FUN = fitMonotonicSpline,
                             depth = sampleDepths,
                             model.depth = model.depths,
                             kFrac = kFrac)
  }else{
    aeig <- apply(X = sampledAgeEns,
                  MARGIN = 2,
                  FUN = fitMonotonicSpline,
                  depth = sampleDepths,
                  model.depth = model.depths,
                  kFrac = kFrac)
  }
  #calculate accum rates
  mgDepthSteps <- matrix(rep(diff(model.depths),times = n.ms.ens),ncol = n.ms.ens,nrow = length(model.depths) - 1)
  mgAgeSteps <- apply(X = aeig, MARGIN = 2,FUN = diff)
  yrPerDepth <- mgAgeSteps/mgDepthSteps
  meanScamSedrate <- 1/mean(yrPerDepth)
  varveMean <- meanScamSedrate * varveScalingFactor


  #add the high frequency bit
  if(progress){message(crayon::blue("Estimating high-frequency variability in the model..."))}
  varvedPrior <- addVarves(ages = aeig,
                           model.depths = model.depths,
                           yrPerDepth = yrPerDepth,
                           totalDepth = diff(range(model.depths)),
                           varveMean = varveMean,
                           H = H,
                           ar1 = ar1,
                           n.varve.ens = n.ms.ens,
                           DT = DT,
                           progress = progress)

  return(varvedPrior)
}

#' Add high-frequency varve variability to age-depth ensembles
#'
#' The core of the flux model: takes an ensemble of low-frequency age-depth
#' estimates and a per-member sedimentation-rate (years-per-depth) estimate,
#' simulates annual-layer (varve) thickness variability with the requested
#' persistence, scales it to be consistent with total depth and total age, and
#' returns an updated ensemble of age priors plus their log-likelihood against
#' the dated distributions. Called repeatedly inside the parameter-fitting loop
#' (see the package vignette).
#'
#' @param ages Matrix of low-frequency age estimates (depth x ensemble member).
#' @param model.depths Depths corresponding to the rows of `ages` /
#'   `yrPerDepth`.
#' @param yrPerDepth Matrix of years-per-depth (inverse sedimentation rate),
#'   same shape as `ages`.
#' @param totalDepth Total depth spanned by the record, used to rescale
#'   simulated varve depths.
#' @param varveMean Mean simulated varve thickness.
#' @param H Hurst parameter for long-memory varve simulation. Specify `H` or
#'   `ar1`, not both.
#' @param ar1 AR(1) coefficient for varve simulation.
#' @param n.varve.ens Number of varve ensemble members (must match the number
#'   of columns of `ages`).
#' @param DT A LiPD distribution table (list of distribution table objects)
#'   used to evaluate the fit.
#' @param progress Show progress bars?
#'
#' @return A list with `agePriors` (matrix of updated age priors), `ageDepths`
#'   (depths they're evaluated at), and `varvedPriorLogObj` (log-likelihood of
#'   each ensemble member).
#' @export
addVarves <- function(ages, model.depths,  yrPerDepth, totalDepth, varveMean, H, ar1, n.varve.ens, DT, progress = TRUE){
  nYears <- apply(ages, 2, \(x) ceiling(max(x)) - floor(min(x)))
  maxN <- max(nYears)

  # simulate the whole ensemble at the maximum length in one call, then trim
  # each member back to its own number of years (NA-padding the tail, as the
  # per-member length.out did before)
  v1 <- simulateVarves(n = maxN, n.ens = length(nYears), ar1 = ar1, H = H,
                       mean = varveMean)
  for (k in which(nYears < maxN)) v1[(nYears[k] + 1L):maxN, k] <- NA

  vInv <- 1/v1

  depths <- apply(v1,2,cumsum)
  totalDepths <- apply(depths,2,max,na.rm = TRUE)

  #scale depth
  adjustFactor <- totalDepth/totalDepths
  adjustedDepths <- matrix(adjustFactor,ncol = ncol(depths),nrow = nrow(depths), byrow = TRUE) * depths
  depthStep <- stats::median(abs(diff(adjustedDepths)),na.rm = TRUE)
  modelDepthStep <- stats::median(diff(model.depths))

  regularGrid <- stats::sd(diff(model.depths))/modelDepthStep < 0.1
  if(regularGrid){#typical use case, standard gaps much faster
    nearestDepths <- roundAny(adjustedDepths,accuracy = modelDepthStep)
  }else{#irregular depths, much slower
    nearestDepths <- apply(adjustedDepths,2,
                           \(col) model.depths[purrr::map_dbl(col, \(d) ifelse(is.na(d),NA,which.min(abs(d - model.depths))))])
  }

  if(depthStep > modelDepthStep){#then interpolate
    stop("Using depth steps that are finer than annual is not yet supported.")
  }else{#then bin: mean accumulation rate within each (depth, ensemble) cell
    depthVec <- as.vector(nearestDepths)
    accVec <- as.vector(vInv)
    ensVec <- rep(seq_len(ncol(vInv)), each = nrow(vInv))
    keep <- !is.na(depthVec) & !is.na(accVec)
    depthVec <- depthVec[keep]; accVec <- accVec[keep]; ensVec <- ensVec[keep]
    nEns <- ncol(vInv)

    # Map each observation to a depth-bin index `di` and the sorted unique
    # depths. On a regular grid the depths are exact multiples of the step, so
    # integer indexing + tabulate beats a unique()/match() over ~1e6 values.
    if(regularGrid){
      k <- as.integer(round(depthVec / modelDepthStep))
      kmin <- min(k)
      counts <- tabulate(k - kmin + 1L)
      occ <- which(counts > 0L)                       # bins that actually occur
      new.depths.to.model <- (kmin + occ - 1L) * modelDepthStep
      dense2sparse <- integer(length(counts))
      dense2sparse[occ] <- seq_along(occ)
      di <- dense2sparse[k - kmin + 1L]
    }else{
      new.depths.to.model <- sort(unique(depthVec))
      di <- match(depthVec, new.depths.to.model)
    }
    nDepth <- length(new.depths.to.model)

    # one integer cell id per (depth, ensemble); bin sums and counts in a single
    # C-backed rowsum over a two-column matrix (so the grouping is done once),
    # then place means into a (depth x ensemble) matrix. Cells with no
    # observations stay 0 (matching the previous NA -> 0 fill).
    cell <- (ensVec - 1L) * nDepth + di
    agg <- rowsum(cbind(accVec, 1), cell)
    bvInv <- numeric(nDepth * nEns)
    bvInv[as.integer(rownames(agg))] <- agg[, 1] / agg[, 2]
    bvInv <- matrix(bvInv, nrow = nDepth, ncol = nEns)
  }

  depthStep <- stats::median(diff(new.depths.to.model))

  #interpolate yrPerDepth to match new model.depth
  nypd <- apply(yrPerDepth,
                MARGIN = 2,
                \(sr) Hmisc::approxExtrap(x = model.depths, y = sr,xout = new.depths.to.model,method = "linear")$y)

  # years-per-depth is physically non-negative, but linear extrapolation can dip
  # below zero at the grid ends; clamp so the prior stays >= 0 and the resulting
  # age model is monotone (no reversals).
  nypd[!is.finite(nypd) | nypd < 0] <- 0

  if(all(dim(bvInv) == dim(nypd))){
    prior <- bvInv * nypd
  }else{
    stop("Feature is not set up yet - n.varve.ens and n.ms.ens need to be the same for now")
  }

  totalAge <- colSums(prior * depthStep)

  #scale to total age
  adjustFactor <- nYears/totalAge
  adjustedForTotalAgePrior <- matrix(adjustFactor,ncol = ncol(prior),nrow = nrow(prior), byrow = TRUE) * prior


  if(!all(round(colSums(adjustedForTotalAgePrior * depthStep)) == nYears)){
    stop("Mismatch")
  }

  #sum up the ages
  cumAges <- apply(adjustedForTotalAgePrior * depthStep,MARGIN = 2,cumsum)

  #adjust for starting age
  minAge <- apply(X = ages,MARGIN = 2,min)
  agePriors <- matrix(minAge,ncol = ncol(prior),nrow = nrow(prior),byrow = TRUE) + cumAges

  #assess fit
  #ages and depths may not match?
  varvedPriorLogObj <- ageProbsDT(DT,ageEstimates = agePriors,depths = new.depths.to.model)


  return(list(agePriors = agePriors,
              ageDepths = new.depths.to.model,
              varvedPriorLogObj = varvedPriorLogObj))
}

#' @noRd
getAgeLikelihoodFromEns <- function(i,DT,sampleAE){
  oneEns <- sampleAE[,i]
  logObj <- sum(purrr::map_dbl(seq_along(DT),\(x) ageProb(DT[[x]],ageToTest = oneEns[x])))
  return(logObj)
}

#' Log-likelihood of age ensembles given dated distributions
#'
#' Evaluates each column of an age-estimate ensemble against the calibrated
#' age distributions in a LiPD distribution table, summing log probability
#' density across dated horizons (matched to the nearest model depth).
#'
#' @param DT A LiPD distribution table (list of distribution table objects).
#' @param ageEstimates A matrix of age estimates (depth x ensemble member).
#' @param depths The depths corresponding to the rows of `ageEstimates`.
#'
#' @return A vector of log-likelihoods, one per ensemble member.
#' @export
ageProbsDT <- function(DT, ageEstimates,depths){
  sampleDepths <- purrr::map_dbl(DT,"depth")
  whichSample <- purrr::map_dbl(sampleDepths,\(x) which.min(abs(x-depths)))

  sampleAE <- ageEstimates[whichSample,]

  ensObj <- purrr::map_dbl(seq_len(ncol(ageEstimates)), getAgeLikelihoodFromEns,DT,sampleAE)

  return(ensObj)
}

# In-progress alternative likelihood path: these expect ATdepths, depthBins,
# and cdfFuns to exist in the calling environment and aren't wired into the
# rest of the package yet.

#' @noRd
ageProbs <- function(ATdepth, ageEstimates, cdfFuns, ATdepths){
  whichFun <- which(dplyr::near(ATdepth,ATdepths,tol = 0.5))
  if(length(whichFun) == 0){
    stop("Depth is not within a cm of an age ensemble table depth")
  }
  if(length(whichFun) > 1){
    whichFun <- whichFun[1]
  }
  ap <- cdfFuns[[whichFun]](ageEstimates)
  if(ap <= 0){ap <- 1e-16}

  return(log(ap))
}

#' @noRd
modelProbs <- function(newModel, depthBins, ATdepths, cdfFuns){
  agesEstimatesAtATDepths <- Hmisc::approxExtrap(x = depthBins,y = newModel,xout = ATdepths)$y

  logObj <- sum(purrr::map2_dbl(ATdepths,agesEstimatesAtATDepths,ageProbs,cdfFuns,ATdepths = ATdepths))
  return(logObj)
}
