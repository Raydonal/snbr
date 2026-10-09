## Penalized B-spline (P-spline) representation of the nonparametric
## component of the partially linear skew-normal model
##   y_i = x_i' beta + f(t_i) + eps_i.
## f is written as Z gamma, with Z a cubic B-spline basis reparameterized
## to satisfy sum_i f(t_i) = 0 (so that the intercept stays in x_i), and
## the roughness penalty is lambda gamma' S gamma / 2 with S built from
## second-order differences (Eilers and Marx, 1996).

#' P-spline basis with a sum-to-zero constraint
#'
#' @param t covariate entering nonparametrically.
#' @param k number of B-spline basis functions before the constraint.
#' @param degree spline degree.
#' @param diff order of the difference penalty.
#' @return A list with the constrained basis `Z` (n x (k-1)), the
#'   penalty `S` ((k-1) x (k-1)) and the objects needed by
#'   [predict_pspline()].
#' @references Eilers, P. H. C. and Marx, B. D. (1996). Flexible smoothing
#'   with B-splines and penalties. *Statistical Science*, 11, 89--121.
#' @export
pspline_basis <- function(t, k = 12L, degree = 3L, diff = 2L) {
  rng <- range(t)
  nint <- k - degree
  dx <- diff(rng) / nint
  knots <- seq(rng[1] - degree * dx, rng[2] + degree * dx, by = dx)
  B <- splines::splineDesign(knots, t, ord = degree + 1L, outer.ok = TRUE)
  cvec <- colSums(B)
  qrc <- qr(matrix(cvec, ncol = 1))
  N <- qr.Q(qrc, complete = TRUE)[, -1, drop = FALSE]
  D <- diff(diag(k), differences = diff)
  S <- crossprod(D %*% N)
  Z <- B %*% N
  colnames(Z) <- paste0("s(", seq_len(ncol(Z)), ")")
  list(Z = Z, S = S, knots = knots, degree = degree, N = N, range = rng)
}

#' Evaluate a fitted P-spline component at new points
#'
#' @param basis output of [pspline_basis()].
#' @param gamma coefficients of the constrained basis.
#' @param tnew new values of the covariate.
#' @return Vector of fitted values of \eqn{f} at `tnew`.
#' @export
predict_pspline <- function(basis, gamma, tnew) {
  B <- splines::splineDesign(basis$knots, tnew, ord = basis$degree + 1L, outer.ok = TRUE)
  drop(B %*% basis$N %*% gamma)
}

## full (p+2) x (p+2) penalty for theta = (beta, gamma, sigma, alpha)
.full_penalty <- function(p_lin, S) {
  q <- ncol(S)
  P <- p_lin + q + 2L
  Sf <- matrix(0, P, P)
  idx <- p_lin + seq_len(q)
  Sf[idx, idx] <- S
  Sf
}

#' Partially linear skew-normal regression with bias reduction
#'
#' Fits \eqn{y_i = x_i^\top\beta + f(t_i) + \varepsilon_i},
#' \eqn{\varepsilon_i \sim SN(0, \sigma, \alpha)}, with \eqn{f} a P-spline.
#' For a fixed smoothing parameter the estimating equation is the
#' penalized score, optionally corrected by the mean bias-reducing
#' adjustment computed with the penalized information
#' \eqn{H = K + \lambda S}. The smoothing parameter is chosen on a grid by
#' AIC with effective degrees of freedom \eqn{\mathrm{tr}(H^{-1}K)}.
#'
#' @param X linear part of the model matrix (including the intercept).
#' @param t covariate of the smooth term.
#' @param y response.
#' @param method `"ML"` (penalized likelihood) or `"BR"` (penalized
#'   likelihood with the mean bias-reducing adjustment).
#' @param k basis dimension.
#' @param lambda optional fixed smoothing parameter; if `NULL` it is
#'   selected from `lambda_grid`.
#' @param lambda_grid grid of candidate values.
#' @return A list with the coefficients, the smoothing parameter, the
#'   effective degrees of freedom, the fitted smooth, AIC and the basis.
#' @export
snplm <- function(X, t, y, method = c("BR", "ML"), k = 12L, lambda = NULL,
                  lambda_grid = 10^seq(-3, 4, by = 0.5)) {
  method <- match.arg(method)
  bs <- pspline_basis(t, k)
  Xf <- cbind(X, bs$Z)
  p_lin <- ncol(X)
  Sf <- .full_penalty(p_lin, bs$S)
  P <- ncol(Xf) + 2L
  fit_one <- function(lam, start = NULL) {
    if (method == "BR") {
      f <- sn_br(Xf, y, "mean", start = start, S = Sf, lambda = lam)
      th <- f$theta
      ok <- f$converged
    } else {
      f <- .penalized_mle(Xf, y, Sf, lam, start)
      th <- f$theta
      ok <- f$converged
    }
    if (!ok || any(!is.finite(th))) return(NULL)
    cu <- sn_cumulants_fast(th[P], th[P - 1L])
    K <- sn_info(th, Xf, cu)
    edf <- sum(diag(solve(K + lam * Sf, K)))
    ll <- sn_loglik(th, Xf, y)
    list(theta = th, lambda = lam, edf = edf, loglik = ll, aic = -2 * ll + 2 * edf)
  }
  if (is.null(lambda)) {
    fits <- list()
    start <- NULL
    for (lam in rev(lambda_grid)) {   # from smooth to rough, warm starts
      f <- fit_one(lam, start)
      if (!is.null(f)) { fits[[length(fits) + 1L]] <- f; start <- f$theta }
    }
    if (!length(fits)) stop("no convergent fit on the lambda grid")
    best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "aic"))]]
  } else {
    best <- fit_one(lambda)
    if (is.null(best)) stop("fit did not converge")
  }
  th <- best$theta
  names(th) <- .coef_names(X, bs$Z, extra_prefix = "f")
  gam <- th[p_lin + seq_len(ncol(bs$Z))]
  list(theta = th, beta = th[seq_len(p_lin)], gamma = gam,
       sigma = th[P - 1L], alpha = th[P], lambda = best$lambda, edf = best$edf,
       loglik = best$loglik, aic = best$aic, fitted_smooth = drop(bs$Z %*% gam),
       basis = bs, X = Xf, S = Sf, method = method)
}

## Penalized maximum likelihood for fixed lambda, with the same exact
## treatment of infinite shape estimates as sn_mle(): as alpha -> +-Inf the
## penalized log-likelihood converges to
##   sup_{sigma} [ n log 2 - n log sigma - n log(2 pi)/2
##                 - min_{b: sign*(y - X b) >= 0} { |y - X b|^2/(2 sigma^2) + lambda b'S b/2 } ],
## a quadratic programme for each sigma followed by a one-dimensional search.
.loglik_limit_pen <- function(X, y, S, lambda, sign = 1) {
  n <- nrow(X)
  p <- ncol(X)
  Sb <- S[seq_len(p), seq_len(p), drop = FALSE]
  inner <- function(ls) {
    s2 <- exp(2 * ls)
    D <- crossprod(X) / s2 + lambda * Sb
    D <- D + diag(1e-10 * max(abs(diag(D))), p)
    qp <- tryCatch(quadprog::solve.QP(Dmat = D, dvec = drop(crossprod(X, y)) / s2,
                                      Amat = t(-sign * X), bvec = -sign * y),
                   error = function(e) NULL)
    if (is.null(qp)) return(list(val = -Inf, beta = rep(NA_real_, p)))
    b <- qp$solution
    r <- y - drop(X %*% b)
    val <- n * log(2) - n * ls - n / 2 * log(2 * pi) - sum(r^2) / (2 * s2) -
      lambda / 2 * drop(crossprod(b, Sb %*% b))
    list(val = val, beta = b)
  }
  r0 <- y - stats::lm.fit(X, y)$fitted.values
  ls0 <- log(stats::sd(r0) + 1e-8)
  o <- stats::optimize(function(ls) -inner(ls)$val, c(ls0 - 4, ls0 + 3), tol = 1e-10)
  fin <- inner(o$minimum)
  list(loglik = fin$val, beta = fin$beta, sigma = exp(o$minimum))
}

.penalized_mle <- function(X, y, S, lambda, start = NULL, tol = 1e-8) {
  P <- ncol(X) + 2L
  p <- ncol(X)
  plik <- function(th) sn_loglik(th, X, y) - 0.5 * lambda * drop(crossprod(th, S %*% th))
  obj <- function(par) -plik(c(par[seq_len(p)], exp(par[p + 1L]), par[P]))
  gr <- function(par) {
    th <- c(par[seq_len(p)], exp(par[p + 1L]), par[P])
    g <- sn_score(th, X, y) - lambda * drop(S %*% th)
    g[p + 1L] <- g[p + 1L] * th[p + 1L]
    -g
  }
  st0 <- if (is.null(start)) .sn_start(X, y) else start
  if (!all(is.finite(st0))) st0 <- .sn_start(X, y)
  best <- NULL
  for (a0 in c(st0[P], -3, -1, 1, 3)) {
    s <- st0
    s[P] <- a0
    o <- tryCatch(stats::optim(c(s[seq_len(p)], log(s[p + 1L]), s[P]), obj, gr,
                               method = "BFGS", control = list(maxit = 500, reltol = 1e-12)),
                  error = function(e) NULL)
    if (is.null(o) || !is.finite(o$value)) next
    if (is.null(best) || -o$value > best$val) {
      best <- list(theta = c(o$par[seq_len(p)], exp(o$par[p + 1L]), o$par[P]), val = -o$value)
    }
  }
  lim <- c(.loglik_limit_pen(X, y, S, lambda, -1)$loglik,
           .loglik_limit_pen(X, y, S, lambda, 1)$loglik)
  finite <- !is.null(best) && best$val > max(lim) + tol
  theta <- if (finite) best$theta else {
    sg <- if (lim[2] >= lim[1]) 1 else -1
    L <- .loglik_limit_pen(X, y, S, lambda, sg)
    c(L$beta, L$sigma, sg * Inf)
  }
  list(theta = theta, converged = finite, finite = finite,
       value = if (finite) best$val else max(lim))
}
