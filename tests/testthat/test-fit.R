test_that("maximum likelihood agrees with sn::selm when the estimate is finite", {
  skip_if_not_installed("sn")
  set.seed(10)
  x <- runif(200)
  X <- cbind(1, x)
  y <- -2 + 2 * x + rsnorm(200, 0, 1, 2)
  f <- sn_mle(X, y)
  expect_true(f$finite)
  s <- sn::selm.fit(X, y, family = "SN")
  expect_equal(unname(f$theta), unname(s$param$dp), tolerance = 1e-5)
})

test_that("an all-positive residual pattern gives an infinite shape estimate", {
  set.seed(11)
  x <- runif(30)
  X <- cbind(1, x)
  y <- 1 + x + abs(rnorm(30))   # half-normal errors
  f <- sn_mle(X, y)
  expect_false(f$finite)
  expect_equal(f$theta[4], Inf, ignore_attr = TRUE)
})

test_that("the profile likelihood approaches the half-normal limit", {
  set.seed(12)
  x <- runif(30)
  X <- cbind(1, x)
  y <- 1 + x + abs(rnorm(30))
  L <- snbr:::.loglik_limit(X, y, 1)
  prof <- function(a) {
    o <- optim(c(L$beta - c(1e-3, 0), log(L$sigma)),
               function(p) -sn_loglik(c(p[1:2], exp(p[3]), a), X, y),
               method = "BFGS", control = list(reltol = 1e-14, maxit = 5000))
    -o$value
  }
  pr <- sapply(c(10, 100, 1000, 1e4), prof)
  expect_true(all(diff(pr) > 0))           # increasing in alpha
  expect_true(all(pr < L$loglik))          # below the limit
  expect_lt(L$loglik - pr[4], 0.05)        # and approaching it
})

test_that("bias-reduced estimates are finite when the MLE is not", {
  set.seed(13)
  x <- runif(40)
  X <- cbind(1, x)
  y <- 1 + x + abs(rnorm(40)) + 0.05 * rnorm(40)
  expect_false(sn_mle(X, y)$finite)
  bm <- sn_br(X, y, "mean")
  bd <- sn_br(X, y, "median")
  expect_true(bm$converged && all(is.finite(bm$theta)))
  expect_true(bd$converged && all(is.finite(bd$theta)))
})

test_that("fits are equivariant under rescaling of a covariate", {
  set.seed(14)
  x <- runif(60)
  y <- -1 + 3 * x + rsnorm(60, 0, 1, 4)
  for (type in c("mean", "median")) {
    a <- sn_br(cbind(1, x), y, type)$theta
    b <- sn_br(cbind(1, 1000 * x), y, type)$theta
    expect_equal(unname(b), unname(a * c(1, 1e-3, 1, 1)), tolerance = 1e-6)
  }
})

test_that("the adjusted score vanishes at the mean bias-reduced estimate", {
  set.seed(15)
  x <- runif(80)
  X <- cbind(1, x)
  y <- 2 - x + rsnorm(80, 0, 0.5, 3)
  b <- sn_br(X, y, "mean")
  th <- b$theta
  adj <- sn_score(th, X, y) - drop(sn_info(th, X) %*% sn_bias(th, X))
  expect_lt(max(abs(adj)), 1e-5)
})
