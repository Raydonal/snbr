## Second-order bias of (penalized) maximum likelihood estimators in the
## skew-normal regression model y_i = x_i' beta + eps_i,
## eps_i ~ SN(0, sigma, alpha), theta = (beta, sigma, alpha).
##
## Each observation depends on theta only through the reduced parameter
## phi_i = (eta_i, sigma, alpha), eta_i = x_i' beta, with Jacobian
## J_i = blockdiag(x_i', 1, 1).  Writing H for the expected negative
## Hessian of the estimating function g = U - lambda S theta and K for the
## Fisher information, the O(1/n) bias is
##   b = H^{-1} sum_i J_i' d_i,
##   d_i[a] = sum_{b,c} k21[a,b,c] (J_i H^{-1} J_i')[b,c]
##            + 1/2 sum_{b,c} k3[a,b,c] (J_i H^{-1} K H^{-1} J_i')[b,c].
## With lambda = 0, H = K and b is the Cox and Snell (1968) bias.

## block entries of J_i A J_i' for all i, as an n x 3 x 3 array
.jaj <- function(X, A) {
  p <- ncol(X)
  ib <- seq_len(p)
  is <- p + 1L
  ia <- p + 2L
  n <- nrow(X)
  out <- array(0, c(n, 3, 3))
  out[, 1, 1] <- rowSums((X %*% A[ib, ib, drop = FALSE]) * X)
  out[, 1, 2] <- out[, 2, 1] <- drop(X %*% A[ib, is])
  out[, 1, 3] <- out[, 3, 1] <- drop(X %*% A[ib, ia])
  out[, 2, 2] <- A[is, is]
  out[, 2, 3] <- out[, 3, 2] <- A[is, ia]
  out[, 3, 3] <- A[ia, ia]
  out
}

## sum_i J_i' v_i for an n x 3 matrix v
.jtv <- function(X, v) {
  c(drop(crossprod(X, v[, 1])), sum(v[, 2]), sum(v[, 3]))
}

#' Expected information of the skew-normal regression model
#'
#' @param theta parameter vector `c(beta, sigma, alpha)`.
#' @param X model matrix of the location (n x p).
#' @param cu optional output of [sn_cumulants()] at `theta`'s shape and
#'   scale, to avoid recomputation.
#' @return The \eqn{(p+2) \times (p+2)} Fisher information matrix.
#' @export
sn_info <- function(theta, X, cu = NULL) {
  p <- ncol(X)
  sigma <- theta[p + 1L]
  alpha <- theta[p + 2L]
  if (is.null(cu)) cu <- sn_cumulants_fast(alpha, sigma)
  W <- -cu$k2
  n <- nrow(X)
  sx <- colSums(X)
  K <- matrix(0, p + 2L, p + 2L)
  ib <- seq_len(p)
  K[ib, ib] <- W[1, 1] * crossprod(X)
  K[ib, p + 1L] <- K[p + 1L, ib] <- W[1, 2] * sx
  K[ib, p + 2L] <- K[p + 2L, ib] <- W[1, 3] * sx
  K[p + 1L, p + 1L] <- n * W[2, 2]
  K[p + 1L, p + 2L] <- K[p + 2L, p + 1L] <- n * W[2, 3]
  K[p + 2L, p + 2L] <- n * W[3, 3]
  nm <- .coef_names(X)
  dimnames(K) <- list(nm, nm)
  K
}

#' Second-order bias of the (penalized) maximum likelihood estimator
#'
#' Evaluates the \eqn{O(n^{-1})} bias of the maximum likelihood estimator
#' of \eqn{\theta = (\beta, \sigma, \alpha)} in the skew-normal regression
#' model (Cox and Snell, 1968). When a quadratic penalty
#' \eqn{\lambda \theta^\top S \theta / 2} is present, the bias is that of
#' the penalized estimator around its own target, obtained from the same
#' expansion with the penalized information \eqn{H = K + \lambda S}.
#'
#' @inheritParams sn_info
#' @param S optional penalty matrix of dimension \eqn{(p+2) \times (p+2)}
#'   (zero rows and columns for unpenalized parameters).
#' @param lambda smoothing parameter (non-negative scalar).
#' @return Numeric vector of length \eqn{p+2} with the estimated biases.
#' @references Cox, D. R. and Snell, E. J. (1968). A general definition of
#'   residuals. *Journal of the Royal Statistical Society B*, 30, 248--275.
#' @export
sn_bias <- function(theta, X, S = NULL, lambda = 0, cu = NULL) {
  p <- ncol(X)
  sigma <- theta[p + 1L]
  alpha <- theta[p + 2L]
  if (is.null(cu)) cu <- sn_cumulants_fast(alpha, sigma)
  K <- sn_info(theta, X, cu)
  H <- K
  if (!is.null(S) && lambda > 0) H <- K + lambda * S
  Hi <- solve(H)
  M1 <- .jaj(X, Hi)
  M2 <- if (identical(H, K)) M1 else .jaj(X, Hi %*% K %*% Hi)
  n <- nrow(X)
  d <- matrix(0, n, 3)
  for (a in 1:3) {
    for (b in 1:3) for (c in 1:3) {
      d[, a] <- d[, a] + cu$k21[a, b, c] * M1[, b, c] + 0.5 * cu$k3[a, b, c] * M2[, b, c]
    }
  }
  out <- drop(Hi %*% .jtv(X, d))
  names(out) <- .coef_names(X)
  out
}
