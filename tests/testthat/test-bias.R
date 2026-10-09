## Closed-form bias of Lopes and Vasconcellos (quantities Q1..Q18 and the
## vector delta), transcribed from the earlier manuscript, used here as an
## independent check of the general Cox-Snell implementation.
bias_closed_form <- function(theta, X) {
  A <- function(n, m, alpha) {
    snbr:::.sn_integrate(function(w) {
      z1 <- exp(stats::dnorm(alpha * w, log = TRUE) - stats::pnorm(alpha * w, log.p = TRUE))
      z1^m * w^n
    }, alpha)
  }
  n <- nrow(X); p <- ncol(X)
  s <- theta[p + 1]; a <- theta[p + 2]; b <- sqrt(2 / pi)
  A02 <- A(0, 2, a); A12 <- A(1, 2, a); A22 <- A(2, 2, a); A32 <- A(3, 2, a); A42 <- A(4, 2, a)
  A03 <- A(0, 3, a); A13 <- A(1, 3, a); A23 <- A(2, 3, a); A33 <- A(3, 3, a)
  r <- 1 + a^2
  Q <- numeric(18)
  Q[1] <- a^3 / (2 * s^3) * (3 * a * A12 + 2 * A03 - b * r^(-3/2))
  Q[2] <- 1 / s^3 * (a^2 * A02 + 1 + 1.5 * a^4 * A22 + a^3 * A13)
  Q[3] <- a / s^2 * (a^2 / 2 * A22 - A02)
  Q[4] <- 1 / (2 * s^3) * (3 * a^4 * A22 - 2 * a^2 * A02 + 2 * a^3 * A13 - 2)
  Q[5] <- a / s^2 * (A02 - a * A13 - 1.5 * a^2 * A22)
  Q[6] <- a / s^3 * (b * a^2 * (a^2 - 0.5) * r^(-5/2) + a^2 * A23 + 1.5 * a^3 * A32)
  Q[7] <- 1 / (2 * s^2) * (a^3 * A32 - 3 * b * (1 + 2 * a^2) * r^(-5/2) - a * A12)
  Q[8] <- 1 / s^2 * (b * (2 * a^2 + 0.5) * r^(-5/2) + a / 2 * A12 - a^2 * A23 - 1.5 * a^3 * A32)
  Q[9] <- 1 / (2 * s) * (3 * b * a * r^(-5/2) - a^2 * A32)
  Q[10] <- n / s^3 * (1.5 * a^4 * A42 + a^3 * A33 - 1 - a^2 * A22)
  Q[11] <- n * a^3 / (2 * s^2) * A42
  Q[12] <- n * a / s^2 * (A22 - 1.5 * a^2 * A42 - a * A33)
  Q[13] <- -n * a^2 / (2 * s) * A42
  Q[14] <- n / s * (1.5 * a^2 * A42 + a * A33 - A22)
  Q[15] <- n * a / 2 * A42
  Q[16] <- a / (2 * s^3) * (3 * a^3 * A32 + 2 * a^2 * A23 - 4 * a * A12 - b * (6 * a^4 + 13 * a^2 + 4) * r^(-5/2))
  Q[17] <- 1 / (2 * s^2) * (3 * a * A12 - 3 * a^3 * A32 - 2 * a^2 * A23 - b * (1 - 2 * a^2) * r^(-5/2))
  Q[18] <- 1 / (2 * s) * (3 * a^2 * A32 + 2 * a * A23 - 2 * A12 - 3 * b * a * r^(-5/2))
  Ki <- solve(sn_info(theta, X, sn_cumulants(a, s)))
  ib <- 1:p; is <- p + 1; ia <- p + 2
  Kbb <- Ki[ib, ib]; Kbs <- Ki[ib, is]; Kba <- Ki[ib, ia]
  Kss <- Ki[is, is]; Ksa <- Ki[is, ia]; Kaa <- Ki[ia, ia]
  dbb <- rowSums((X %*% Kbb) * X); trb <- sum(dbb)
  d11 <- Q[1] * dbb + (Q[2] + Q[4]) * drop(X %*% Kbs) + (Q[3] + Q[5]) * drop(X %*% Kba) +
    ((Q[7] + Q[8]) * Ksa + Q[6] * Kss + Q[9] * Kaa)
  d21 <- Q[4] * trb + (Q[6] + Q[16]) * sum(X %*% Kbs) + (Q[7] + Q[17]) * sum(X %*% Kba) +
    (Q[11] + Q[12]) * Ksa + Q[10] * Kss + Q[13] * Kaa
  d31 <- Q[5] * trb + (Q[8] + Q[17]) * sum(X %*% Kbs) + (Q[9] + Q[18]) * sum(X %*% Kba) +
    (Q[13] + Q[14]) * Ksa + Q[12] * Kss + Q[15] * Kaa
  drop(Ki %*% c(crossprod(X, d11), d21, d31))
}

test_that("general Cox-Snell bias matches the closed form of the paper", {
  set.seed(1)
  x <- runif(40)
  X <- cbind(1, x, x^2)
  for (th in list(c(-2, 2, 1, 1, 3), c(0.5, -1, 0.3, 2, -6))) {
    g <- sn_bias(th, X, cu = sn_cumulants(th[5], th[4]))
    expect_equal(unname(g), unname(bias_closed_form(th, X)), tolerance = 1e-9)
  }
})

test_that("penalized bias reduces to Cox-Snell when lambda = 0", {
  set.seed(2)
  X <- cbind(1, runif(30))
  th <- c(1, 2, 0.7, 4)
  S <- diag(c(0, 1, 0, 0))
  expect_equal(sn_bias(th, X, S, lambda = 0), sn_bias(th, X))
})

test_that("bias is linear-equivariant in beta", {
  set.seed(3)
  x <- runif(50)
  X <- cbind(1, x)
  th <- c(-2, 2, 1, 5)
  X2 <- cbind(1, 100 * x)
  th2 <- c(-2, 0.02, 1, 5)
  b1 <- sn_bias(th, X)
  b2 <- sn_bias(th2, X2)
  expect_equal(unname(b2), unname(b1 * c(1, 0.01, 1, 1)), tolerance = 1e-9)
})
