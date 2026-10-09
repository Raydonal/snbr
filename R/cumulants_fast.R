## Fast evaluation of the unit cumulants by cubic-spline interpolation on
## the grid stored in R/sysdata.rda (see data-raw/cumulant_grid.R).
## The grid is regular in t = asinh(alpha) on [0, asinh(2000)]; beyond
## that range, or when the grid is unavailable, the exact quadrature of
## sn_cumulants() is used.

.snbr_cache <- new.env(parent = emptyenv())

.grid_splines <- function() {
  if (is.null(.snbr_cache$spl)) {
    g <- get(".cumgrid", envir = asNamespace("snbr"), inherits = FALSE)
    .snbr_cache$spl <- lapply(seq_len(ncol(g$vals)), function(j) {
      stats::splinefun(g$t, g$vals[, j], method = "fmm")
    })
    .snbr_cache$tmax <- max(g$t)
    ## reflection signs: (-1)^(number of eta and alpha indices)
    s1 <- c(-1, 1, -1)
    .snbr_cache$sign <- c(as.vector(outer(s1, s1)),
                          rep(as.vector(outer(outer(s1, s1), s1)), 3))
  }
  .snbr_cache
}

#' Unit cumulants by interpolation
#'
#' Same output as [sn_cumulants()], obtained by spline interpolation of a
#' precomputed grid. The interpolation error is below \eqn{10^{-8}} in
#' relative terms over \eqn{|\alpha| \le 2000} (checked in the package
#' tests); outside that range the exact quadrature is used.
#'
#' @inheritParams sn_cumulants
#' @return See [sn_cumulants()].
#' @export
sn_cumulants_fast <- function(alpha, sigma = 1) {
  ok <- tryCatch({ .grid_splines(); TRUE }, error = function(e) FALSE)
  t <- asinh(abs(alpha))
  if (!ok || t > .snbr_cache$tmax) return(sn_cumulants(alpha, sigma))
  v <- vapply(.snbr_cache$spl, function(f) f(t), numeric(1))
  if (alpha < 0) v <- v * .snbr_cache$sign
  dn <- list(c("eta", "sigma", "alpha"), c("eta", "sigma", "alpha"),
             c("eta", "sigma", "alpha"))
  k2 <- matrix(v[1:9], 3, 3, dimnames = dn[1:2])
  k3 <- array(v[10:36], c(3, 3, 3), dimnames = dn)
  k21 <- array(v[37:63], c(3, 3, 3), dimnames = dn)
  k111 <- array(v[64:90], c(3, 3, 3), dimnames = dn)
  if (sigma != 1) {
    p2 <- outer(1:3 < 3, 1:3 < 3, "+")
    p3 <- outer(outer(1:3 < 3, 1:3 < 3, "+"), 1:3 < 3, "+")
    k2 <- k2 / sigma^p2
    k3 <- k3 / sigma^p3
    k21 <- k21 / sigma^p3
    k111 <- k111 / sigma^p3
  }
  list(k2 = k2, k3 = k3, k21 = k21, k111 = k111, alpha = alpha, sigma = sigma)
}
