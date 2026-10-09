test_that("Bartlett identities hold for the unit cumulants", {
  for (al in c(-4, 0.5, 3, 12)) {
    cu <- sn_cumulants(al)
    e11 <- -cu$k2
    ## third Bartlett identity
    b3 <- array(0, c(3, 3, 3))
    for (a in 1:3) for (b in 1:3) for (c in 1:3) {
      b3[a, b, c] <- cu$k3[a, b, c] + cu$k21[a, b, c] + cu$k21[a, c, b] +
        cu$k21[b, c, a] + cu$k111[a, b, c]
    }
    expect_lt(max(abs(b3)), 1e-10)
    expect_true(isSymmetric(unname(e11)))
  }
})

test_that("the information is singular at alpha = 0", {
  W <- -sn_cumulants(0)$k2
  expect_equal(W[1, 3]^2, W[1, 1] * W[3, 3], tolerance = 1e-10)
  expect_equal(W[1, 3], sqrt(2 / pi), tolerance = 1e-10)
})

test_that("reflection alpha -> -alpha changes signs as expected", {
  a <- sn_cumulants(2.7)
  b <- sn_cumulants(-2.7)
  s1 <- c(-1, 1, -1)
  s2 <- outer(s1, s1)
  s3 <- outer(s2, s1)
  expect_equal(unname(b$k2), unname(s2 * a$k2), tolerance = 1e-10)
  expect_equal(unname(b$k3), unname(s3 * a$k3), tolerance = 1e-10)
  expect_equal(unname(b$k21), unname(s3 * a$k21), tolerance = 1e-10)
})

test_that("interpolated cumulants agree with quadrature", {
  for (al in c(0.37, -1.9, 7.3, 150)) {
    e <- sn_cumulants(al)
    f <- sn_cumulants_fast(al)
    ev <- c(e$k2, e$k3, e$k21, e$k111)
    fv <- c(f$k2, f$k3, f$k21, f$k111)
    expect_lt(max(abs(ev - fv) / pmax(abs(ev), 1e-3)), 1e-7)
  }
})

test_that("scale enters through powers of 1/sigma", {
  a <- sn_cumulants_fast(3, sigma = 2)
  d <- function(w) snbr:::sn_loglik_derivs(w, 2, 3)
  direct <- snbr:::.sn_integrate(function(w) { x <- d(w); x$es * x$a }, 3)
  expect_equal(unname(a$k21[1, 2, 3]), direct, tolerance = 1e-8)
})
