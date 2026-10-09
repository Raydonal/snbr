## Per-observation joint cumulants of the derivatives of the skew-normal
## log-density, in the reduced parameter phi = (eta, sigma, alpha).
##
## For a single observation y ~ SN(eta, sigma, alpha) let l be its
## log-density and l_a, l_ab, l_abc its derivatives.  The functions below
## return, at sigma = 1,
##   k2[a, b]     = E(l_ab)              (minus the unit Fisher information)
##   k3[a, b, c]  = E(l_abc)
##   k21[a, b, c] = E(l_ab l_c)
##   k111[a, b, c]= E(l_a l_b l_c)
## Location invariance makes these free of eta.  Every derivative with
## respect to eta or sigma carries a factor 1/sigma, so the values at a
## general sigma are obtained by rescaling (see sn_cumulants()).

.sn_index <- c("e", "s", "a")

.sn_key <- function(idx) paste(.sn_index[sort(idx)], collapse = "")

.sn_integrate <- function(f, alpha, rel.tol = 1e-11) {
  ## integrand f(w) * 2 phi(w) Phi(alpha w); the mass sits on w > 0 when
  ## alpha is large and positive, so the real line is split at 0 and at
  ## +/- 1/|alpha| to help the adaptive quadrature.
  dens <- function(w) {
    2 * exp(stats::dnorm(w, log = TRUE) + stats::pnorm(alpha * w, log.p = TRUE))
  }
  g <- function(w) f(w) * dens(w)
  cut <- sort(unique(c(-1 / max(abs(alpha), 1e-8), 0, 1 / max(abs(alpha), 1e-8))))
  cut <- cut[abs(cut) < 40]
  brk <- c(-Inf, cut, Inf)
  tot <- 0
  for (j in seq_len(length(brk) - 1)) {
    ## integrals that vanish exactly (odd moments at alpha = 0) can make the
    ## pure relative criterion fail; retry with an absolute tolerance
    v <- tryCatch(
      stats::integrate(g, brk[j], brk[j + 1], rel.tol = rel.tol,
                       abs.tol = 0, subdivisions = 2000L)$value,
      error = function(e) {
        stats::integrate(g, brk[j], brk[j + 1], rel.tol = rel.tol,
                         abs.tol = 1e-13, subdivisions = 2000L)$value
      })
    tot <- tot + v
  }
  tot
}

#' Exact unit cumulants of the skew-normal log-likelihood derivatives
#'
#' Computes, by adaptive quadrature, the expected values of products of
#' derivatives of the skew-normal log-density with respect to the reduced
#' parameter \eqn{(\eta, \sigma, \alpha)}, where \eqn{\eta} is the location.
#' These are the building blocks of the Fisher information, of the
#' Cox and Snell second-order bias and of the bias-reducing adjusted scores.
#'
#' @param alpha shape parameter (scalar).
#' @param sigma scale parameter (scalar, default 1).
#' @param rel.tol relative tolerance passed to [stats::integrate()].
#' @return A list with components `k2` (3 x 3), `k3`, `k21` and `k111`
#'   (3 x 3 x 3 arrays), indexed by `c("eta", "sigma", "alpha")`. For
#'   `k21` the first two indices refer to the second derivative and the
#'   third to the first derivative.
#' @examples
#' cu <- sn_cumulants(3)
#' -cu$k2   # unit Fisher information at sigma = 1
#' @export
sn_cumulants <- function(alpha, sigma = 1, rel.tol = 1e-11) {
  stopifnot(length(alpha) == 1L, is.finite(alpha), length(sigma) == 1L, sigma > 0)
  dfun <- function(w) sn_loglik_derivs(w, 1, alpha)
  ev <- function(h) .sn_integrate(h, alpha, rel.tol)
  dn <- list(c("eta", "sigma", "alpha"), c("eta", "sigma", "alpha"),
             c("eta", "sigma", "alpha"))
  k2 <- matrix(0, 3, 3, dimnames = dn[1:2])
  k3 <- k21 <- k111 <- array(0, c(3, 3, 3), dimnames = dn)
  for (a in 1:3) for (b in a:3) {
    kab <- .sn_key(c(a, b))
    k2[a, b] <- k2[b, a] <- ev(function(w) dfun(w)[[kab]])
  }
  for (a in 1:3) for (b in 1:3) for (c in 1:3) {
    k3[a, b, c] <- ev(function(w) dfun(w)[[.sn_key(c(a, b, c))]])
    kab <- .sn_key(c(a, b))
    kc <- .sn_index[c]
    k21[a, b, c] <- ev(function(w) { d <- dfun(w); d[[kab]] * d[[kc]] })
    k111[a, b, c] <- ev(function(w) {
      d <- dfun(w); d[[.sn_index[a]]] * d[[.sn_index[b]]] * d[[kc]]
    })
  }
  if (sigma != 1) {
    ## number of (eta, sigma) indices in each entry gives the power of 1/sigma
    p2 <- outer(1:3 < 3, 1:3 < 3, "+")
    p3 <- outer(outer(1:3 < 3, 1:3 < 3, "+"), 1:3 < 3, "+")
    k2 <- k2 / sigma^p2
    k3 <- k3 / sigma^p3
    k21 <- k21 / sigma^p3
    k111 <- k111 / sigma^p3
  }
  list(k2 = k2, k3 = k3, k21 = k21, k111 = k111, alpha = alpha, sigma = sigma)
}
