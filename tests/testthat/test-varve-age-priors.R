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
