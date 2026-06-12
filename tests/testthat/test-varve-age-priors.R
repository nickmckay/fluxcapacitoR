# Build a small synthetic LiPD-style distribution table: gaussian-ish date
# distributions at a few depths.
makeDto <- function(meanAge, sdAge, depth){
  ages <- seq(meanAge - 4 * sdAge, meanAge + 4 * sdAge, length.out = 200)
  list(age = list(values = ages, units = "BP", variableName = "age"),
       probabilityDensity = list(values = dnorm(ages, meanAge, sdAge),
                                 variableName = "probabilityDensity"),
       depth = depth)
}

makeDT <- function(){
  list(makeDto(100, 20, 10),
       makeDto(500, 30, 50),
       makeDto(1000, 40, 100),
       makeDto(1500, 50, 150))
}

test_that("sampleAge draws from the date distribution", {
  set.seed(7)
  dto <- makeDto(1000, 40, 100)
  draws <- replicate(500, sampleAge(dto))
  expect_true(all(draws >= min(dto$age$values) & draws <= max(dto$age$values)))
  expect_equal(mean(draws), 1000, tolerance = 10)
  expect_equal(sd(draws), 40, tolerance = 10)
})

test_that("removeOutliers yields monotonically increasing non-NA ages", {
  set.seed(7)
  depth <- 1:20
  age <- sort(runif(20, 0, 2000))
  age[c(5, 12)] <- age[c(5, 12)] - 500 # introduce reversals
  cleaned <- removeOutliers(age, depth)
  expect_true(all(diff(na.omit(cleaned)) > 0))
  expect_true(any(is.na(cleaned)))
  # untouched values are preserved
  expect_equal(cleaned[!is.na(cleaned)], age[!is.na(cleaned)])
})

test_that("ageProbsDT returns one finite log-likelihood per ensemble member", {
  set.seed(7)
  DT <- makeDT()
  depths <- 0:160
  agePerDepth <- 10
  ageEstimates <- cbind(depths * agePerDepth,
                        depths * agePerDepth + 50)
  ll <- ageProbsDT(DT, ageEstimates = ageEstimates, depths = depths)
  expect_length(ll, 2)
  expect_true(all(is.finite(ll)))
})

test_that("createVarveAgePriors runs end to end on synthetic dates", {
  set.seed(7)
  DT <- makeDT()
  out <- createVarveAgePriors(DT,
                              model.depth.step = 1,
                              ar1 = 0.3,
                              n.ms.ens = 3,
                              n.varve.ens = 3,
                              progress = FALSE)
  expect_named(out, c("agePriors", "ageDepths", "varvedPriorLogObj"))
  expect_equal(ncol(out$agePriors), 3)
  expect_equal(nrow(out$agePriors), length(out$ageDepths))
  # age priors should be monotonically increasing with depth
  expect_true(all(apply(out$agePriors, 2, \(x) all(diff(x) >= 0))))
  expect_length(out$varvedPriorLogObj, 3)
})

test_that("addVarves yields monotone age priors on real data (no end reversals)", {
  # Regression test: yrPerDepth used to be linearly extrapolated onto the depth
  # grid and could go negative at the ends, producing small age reversals. It
  # is now clamped to >= 0. The Imandra ensemble triggered this on both paths.
  skip_if_not_installed("lipdR")
  set.seed(1)

  L <- lipdR::readLipd(system.file("extdata", "Imandra.Holtzman.2024.lpd",
                                   package = "fluxcapacitoR"))
  et <- L$chronData[[1]]$model[[1]]$ensembleTable[[1]]
  ae <- et$ageEnsemble$values
  depths <- et$depth$values
  DT <- L$chronData[[1]]$model[[1]]$distributionTable

  srEns <- apply(ae, 2, \(x) diff(x) / diff(depths))
  srDepths <- rowMeans(cbind(depths[-1], depths[-length(depths)]))
  meanSR <- apply(ae, 2, \(a, d) diff(range(d)) / diff(range(a)), depths)
  ds <- 0.5
  dtm <- seq(min(depths) - ds / 2, max(depths) + ds / 2, by = ds)
  dtmsr <- rowMeans(cbind(dtm[-1], dtm[-length(dtm)]))
  iypd <- apply(srEns, 2, \(sr) Hmisc::approxExtrap(srDepths, sr, xout = dtm, method = "linear")$y)
  iypd[!is.finite(iypd) | iypd < 0] <- 0
  sym <- matrix(apply(ae, 2, min, na.rm = TRUE), ncol = ncol(iypd), nrow = nrow(iypd), byrow = TRUE)
  iae <- apply(iypd, 2, \(x) cumsum(x * ds)) + sym

  nE <- 40
  for (params in list(list(ar1 = 0.5, H = NULL), list(ar1 = NULL, H = 0.8))) {
    r <- addVarves(ages = iae[, seq_len(nE)], model.depths = dtmsr,
                   yrPerDepth = iypd[, seq_len(nE)],
                   totalDepth = diff(range(depths)),
                   varveMean = mean(meanSR) * 25,
                   ar1 = params$ar1, H = params$H,
                   n.varve.ens = nE, DT = DT, progress = FALSE)
    expect_true(all(apply(r$agePriors, 2, \(x) all(diff(x) >= 0))))
  }
})
