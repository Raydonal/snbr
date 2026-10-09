## Likelihood, maximum likelihood with exact detection of infinite shape
## estimates, and the bias-reduced estimators.

#' Random generation from the skew-normal distribution
#'
#' Uses the stochastic representation
#' \eqn{Z = \delta |U| + (1-\delta^2)^{1/2} V}, \eqn{\delta = \alpha/(1+\alpha^2)^{1/2}},
#' with \eqn{U, V} independent standard normal (Henze, 1986).
#'
#' @param n number of observations.
#' @param xi,omega,alpha location, scale and shape (recycled).
#' @return Numeric vector of length `n`.
#' @examples
#' z <- rsnorm(1000, alpha = 5)
#' mean(z)  # close to sqrt(2/pi) * 5 / sqrt(26)
#' @export
rsnorm <- function(n, xi = 0, omega = 1, alpha = 0) {
  delta <- alpha / sqrt(1 + alpha^2)
  u <- abs(stats::rnorm(n))
  v <- stats::rnorm(n)
  xi + omega * (delta * u + sqrt(1 - delta^2) * v)
}

#' Log-likelihood and score of the skew-normal regression model
#'
#' @param theta `c(beta, sigma, alpha)`.
#' @param X model matrix (n x p).
#' @param y response vector.
#' @name sn_loglik
#' @return `sn_loglik()` returns a scalar, `sn_score()` a vector of length
#'   \eqn{p+2}, `sn_hessian()` the observed Hessian matrix.
#' @export
sn_loglik <- function(theta, X, y) {
  p <- ncol(X)
  s <- theta[p + 1L]
  a <- theta[p + 2L]
  if (!is.finite(s) || s <= 0) return(-Inf)
  w <- drop(y - X %*% theta[seq_len(p)]) / s
  sum(log(2) - log(s) + stats::dnorm(w, log = TRUE) + stats::pnorm(a * w, log.p = TRUE))
}

#' @rdname sn_loglik
#' @export
sn_score <- function(theta, X, y) {
  p <- ncol(X)
  s <- theta[p + 1L]
  a <- theta[p + 2L]
  w <- drop(y - X %*% theta[seq_len(p)]) / s
  d <- sn_loglik_derivs(w, s, a)
  c(drop(crossprod(X, d$e)), sum(d$s), sum(d$a))
}

#' @rdname sn_loglik
#' @export
sn_hessian <- function(theta, X, y) {
  p <- ncol(X)
  s <- theta[p + 1L]
  a <- theta[p + 2L]
  w <- drop(y - X %*% theta[seq_len(p)]) / s
  d <- sn_loglik_derivs(w, s, a)
  ib <- seq_len(p)
  H <- matrix(0, p + 2L, p + 2L)
  H[ib, ib] <- crossprod(X * d$ee, X)
  H[ib, p + 1L] <- H[p + 1L, ib] <- drop(crossprod(X, d$es))
  H[ib, p + 2L] <- H[p + 2L, ib] <- drop(crossprod(X, d$ea))
  H[p + 1L, p + 1L] <- sum(d$ss)
  H[p + 1L, p + 2L] <- H[p + 2L, p + 1L] <- sum(d$sa)
  H[p + 2L, p + 2L] <- sum(d$aa)
  H
}

## Supremum of the log-likelihood as alpha -> +Inf (sign = 1) or -Inf
## (sign = -1).  In the limit the model becomes a half-normal regression
## with all residuals of the given sign, so the supremum is attained by
## least squares under the constraints sign * (y - X beta) >= 0.
.loglik_limit <- function(X, y, sign = 1) {
  n <- nrow(X)
  out <- tryCatch({
    qp <- quadprog::solve.QP(Dmat = crossprod(X), dvec = drop(crossprod(X, y)),
                             Amat = t(-sign * X), bvec = -sign * y)
    beta <- qp$solution
    rss <- sum((y - X %*% beta)^2)
    list(loglik = n * log(2) - n / 2 * log(2 * pi * rss / n) - n / 2,
         beta = beta, sigma = sqrt(rss / n))
  }, error = function(e) list(loglik = -Inf, beta = rep(NA_real_, ncol(X)), sigma = NA_real_))
  out
}

## Moment-based starting values (centred-parameter inversion).
.sn_start <- function(X, y) {
  fit <- stats::lm.fit(X, y)
  r <- fit$residuals
  s <- stats::sd(r)
  g1 <- mean((r - mean(r))^3) / mean((r - mean(r))^2)^1.5
  g1 <- sign(g1) * min(abs(g1), 0.95)
  c0 <- sign(g1) * (2 * abs(g1) / (4 - pi))^(1 / 3)
  muz <- c0 / sqrt(1 + c0^2)
  b <- sqrt(2 / pi)
  delta <- max(min(muz / b, 0.99), -0.99)
  alpha <- delta / sqrt(1 - delta^2)
  sigma <- s / sqrt(1 - b^2 * delta^2)
  beta <- fit$coefficients
  beta[is.na(beta)] <- 0
  icpt <- which(apply(X, 2, function(z) all(abs(z - 1) < 1e-12)))
  if (length(icpt)) beta[icpt[1]] <- beta[icpt[1]] - sigma * b * delta
  c(beta, sigma, alpha)
}

## Newton polishing on (beta, log sigma, alpha) with step halving.
.newton_polish <- function(theta, X, y, maxit = 50L, tol = 1e-10) {
  p <- ncol(X)
  ll <- sn_loglik(theta, X, y)
  for (it in seq_len(maxit)) {
    U <- sn_score(theta, X, y)
    Hs <- sn_hessian(theta, X, y)
    step <- tryCatch(-solve(Hs, U), error = function(e) NULL)
    if (is.null(step) || any(!is.finite(step))) break
    if (sum(step * U) <= 0) step <- U / max(1, sqrt(sum(U^2)))  # ascent direction
    h <- 1
    repeat {
      new <- theta + h * step
      lln <- if (new[p + 1L] > 0) sn_loglik(new, X, y) else -Inf
      if (is.finite(lln) && lln >= ll - 1e-12) break
      h <- h / 2
      if (h < 1e-10) break
    }
    if (h < 1e-10) break
    theta <- new
    conv <- abs(lln - ll) < tol && max(abs(h * step)) < 1e-7
    ll <- lln
    if (conv) break
  }
  theta
}

#' Maximum likelihood fit of the skew-normal regression model
#'
#' Maximizes the log-likelihood from several starting points and then
#' decides whether the shape estimate is finite. The decision is exact: as
#' \eqn{\alpha \to \pm\infty} the log-likelihood converges to that of a
#' half-normal regression whose residuals all have the sign of
#' \eqn{\alpha}, and its supremum is computed by constrained least squares.
#' When this limit is not smaller than the best finite value, the maximum
#' likelihood estimate of \eqn{\alpha} is reported as `Inf` or `-Inf`.
#'
#' @param X model matrix (n x p).
#' @param y response.
#' @param start optional starting value `c(beta, sigma, alpha)`.
#' @param tol tolerance (on the log-likelihood scale) used in the
#'   comparison with the limit.
#' @return A list with `theta`, `loglik`, `finite` (logical) and
#'   `loglik_limit` (the two limits, for `alpha -> -Inf, +Inf`).
#' @export
sn_mle <- function(X, y, start = NULL, tol = 1e-8) {
  sc <- .col_scale(X)
  Xs <- sweep(X, 2, sc, "*")
  st <- if (is.null(start)) NULL else .to_scaled(start, sc)
  f <- .sn_mle_core(Xs, y, st, tol)
  f$theta <- .from_scaled(f$theta, sc)
  names(f$theta) <- .coef_names(X)
  f
}

## Column scaling used internally by the fitting functions: non-constant
## columns are divided by their standard deviation, so that the model
## matrix is well conditioned.  A diagonal rescaling of beta leaves the
## maximum likelihood and both bias-reduced estimators equivariant.
.col_scale <- function(X) {
  s <- apply(X, 2, stats::sd)
  s[!is.finite(s) | s < 1e-12] <- 1
  1 / s
}
.to_scaled <- function(theta, sc) {
  p <- length(sc)
  theta[seq_len(p)] <- theta[seq_len(p)] / sc
  theta
}
.from_scaled <- function(theta, sc) {
  p <- length(sc)
  theta[seq_len(p)] <- theta[seq_len(p)] * sc
  theta
}

.sn_mle_core <- function(X, y, start = NULL, tol = 1e-8) {
  p <- ncol(X)
  st0 <- if (is.null(start)) .sn_start(X, y) else start
  starts <- list(st0)
  for (a0 in c(-3, -1, 1, 3)) {
    s1 <- st0
    s1[p + 2L] <- a0
    starts[[length(starts) + 1L]] <- s1
  }
  negll <- function(par) {
    th <- c(par[seq_len(p)], exp(par[p + 1L]), par[p + 2L])
    -sn_loglik(th, X, y)
  }
  neggr <- function(par) {
    th <- c(par[seq_len(p)], exp(par[p + 1L]), par[p + 2L])
    g <- sn_score(th, X, y)
    g[p + 1L] <- g[p + 1L] * th[p + 1L]
    -g
  }
  best <- NULL
  for (s in starts) {
    par0 <- c(s[seq_len(p)], log(s[p + 1L]), s[p + 2L])
    o <- tryCatch(stats::optim(par0, negll, neggr, method = "BFGS",
                               control = list(maxit = 500, reltol = 1e-12)),
                  error = function(e) NULL)
    if (is.null(o) || !is.finite(o$value)) next
    th <- c(o$par[seq_len(p)], exp(o$par[p + 1L]), o$par[p + 2L])
    if (abs(th[p + 2L]) < 500) th <- .newton_polish(th, X, y)
    ll <- sn_loglik(th, X, y)
    if (is.null(best) || ll > best$loglik) best <- list(theta = th, loglik = ll)
  }
  lim <- c(.loglik_limit(X, y, -1)$loglik, .loglik_limit(X, y, 1)$loglik)
  finite <- TRUE
  theta <- best$theta
  if (max(lim) >= best$loglik - tol) {
    finite <- FALSE
    sg <- if (lim[2] >= lim[1]) 1 else -1
    L <- .loglik_limit(X, y, sg)
    theta <- c(L$beta, L$sigma, sg * Inf)
  }
  names(theta) <- .coef_names(X)
  list(theta = theta, loglik = if (finite) best$loglik else max(lim),
       finite = finite, loglik_limit = lim)
}

#' Corrective (Cox and Snell) bias-corrected estimator
#'
#' Returns \eqn{\tilde\theta = \hat\theta - b(\hat\theta)}, with
#' \eqn{b} the second-order bias of [sn_bias()]. Not defined when the
#' maximum likelihood estimate of the shape is infinite.
#'
#' @param theta maximum likelihood estimate.
#' @param X model matrix.
#' @return Numeric vector (all `NA` if `theta` is not finite).
#' @export
sn_correct <- function(theta, X) {
  if (any(!is.finite(theta))) return(theta * NA)
  theta - sn_bias(theta, X)
}

## full-model arrays T[r, s, t] = sum_i sum_{abc} k[a,b,c] J_i[a,r] J_i[b,s] J_i[c,t]
.full3 <- function(X, k) {
  n <- nrow(X)
  p <- ncol(X)
  P <- p + 2L
  J <- array(0, c(n, 3, P))
  J[, 1, seq_len(p)] <- X
  J[, 2, p + 1L] <- 1
  J[, 3, p + 2L] <- 1
  Tt <- array(0, c(P, P, P))
  for (a in 1:3) for (b in 1:3) for (c in 1:3) {
    if (k[a, b, c] == 0) next
    Ja <- J[, a, ]
    Jb <- J[, b, ]
    Jc <- J[, c, ]
    for (r in seq_len(P)) {
      if (all(Ja[, r] == 0)) next
      Tt[r, , ] <- Tt[r, , ] + k[a, b, c] * crossprod(Ja[, r] * Jb, Jc)
    }
  }
  Tt
}

## Median bias-reducing adjustment M1 (Kenne Pagui, Salvan and Sartori,
## 2017, expressions (7) and (10)): the adjusted score is U + K M1.
.median_adjustment <- function(theta, X, cu = NULL) {
  p <- ncol(X)
  P <- p + 2L
  if (is.null(cu)) cu <- sn_cumulants_fast(theta[P], theta[p + 1L])
  K <- sn_info(theta, X, cu)
  T21 <- .full3(X, cu$k21)    # T21[a, b, c] = E(U_ab U_c)
  T111 <- .full3(X, cu$k111)  # E(U_a U_b U_c)
  M1 <- numeric(P)
  for (r in seq_len(P)) {
    lam <- setdiff(seq_len(P), r)
    Kll <- K[lam, lam, drop = FALSE]
    gam <- solve(Kll, K[lam, r])
    v <- numeric(P)
    v[r] <- 1
    v[lam] <- -gam
    nu_inv <- solve(Kll)
    A21 <- matrix(0, P, P)
    A111 <- matrix(0, P, P)
    for (t in seq_len(P)) {
      if (v[t] == 0) next
      A21 <- A21 + v[t] * T21[, , t]
      A111 <- A111 + v[t] * T111[t, , ]
    }
    k1 <- -0.5 * sum(nu_inv * (A21[lam, lam] + A111[lam, lam]))
    k2 <- drop(crossprod(v, K %*% v))
    k3 <- 0
    for (t in seq_len(P)) if (v[t] != 0) k3 <- k3 + v[t] * drop(crossprod(v, T111[t, , ] %*% v))
    M1[r] <- (-k1 + k3 / (6 * k2)) / k2
  }
  M1
}

#' Bias-reducing adjusted score estimators
#'
#' Solves the adjusted score equations of Firth (1993) (`type = "mean"`)
#' or of Kenne Pagui, Salvan and Sartori (2017) (`type = "median"`) for
#' all parameters of the skew-normal regression model, by the quasi
#' Fisher-scoring iteration of Kosmidis and Firth (2010). A quadratic
#' penalty \eqn{\lambda \theta^\top S \theta/2} can be added (mean
#' adjustment only), as used for the semiparametric model.
#'
#' @param X model matrix.
#' @param y response.
#' @param type `"mean"` or `"median"`.
#' @param start optional starting value.
#' @param S,lambda optional penalty.
#' @param maxit,tol iteration controls.
#' @return A list with `theta`, `converged`, `iterations` and the value of
#'   the adjusted score at the solution.
#' @references
#' Firth, D. (1993). Bias reduction of maximum likelihood estimates.
#' *Biometrika*, 80, 27--38.
#'
#' Kosmidis, I. and Firth, D. (2010). A generic algorithm for reducing
#' bias in parametric estimation. *Electronic Journal of Statistics*, 4,
#' 1097--1112.
#'
#' Kenne Pagui, E. C., Salvan, A. and Sartori, N. (2017). Median bias
#' reduction of maximum likelihood estimates. *Biometrika*, 104, 923--938.
#' @export
sn_br <- function(X, y, type = c("mean", "median"), start = NULL,
                  S = NULL, lambda = 0, maxit = 200L, tol = 1e-8) {
  type <- match.arg(type)
  sc <- .col_scale(X)
  Xs <- sweep(X, 2, sc, "*")
  st <- if (is.null(start)) NULL else .to_scaled(start, sc)
  Ss <- if (is.null(S)) NULL else {
    dd <- c(sc, 1, 1)
    S * outer(dd, dd)
  }
  f <- .sn_br_core(Xs, y, type, st, Ss, lambda, maxit, tol)
  f$theta <- .from_scaled(f$theta, sc)
  names(f$theta) <- .coef_names(X)
  f
}

.sn_br_core <- function(X, y, type, start, S, lambda, maxit, tol) {
  p <- ncol(X)
  P <- p + 2L
  if (type == "median" && lambda > 0) stop("median adjustment with a penalty is not implemented")
  pen <- !is.null(S) && lambda > 0
  adj_score <- function(th) {
    tryCatch({
      cu <- sn_cumulants_fast(th[P], th[p + 1L])
      U <- sn_score(th, X, y)
      K <- sn_info(th, X, cu)
      H <- if (pen) K + lambda * S else K
      g <- if (pen) U - lambda * drop(S %*% th) else U
      A <- if (type == "mean") -drop(H %*% sn_bias(th, X, S, lambda, cu))
           else drop(K %*% .median_adjustment(th, X, cu))
      step <- solve(H, g + A)
      if (any(!is.finite(step))) NULL else list(val = g + A, step = step)
    }, error = function(e) NULL)
  }
  merit <- function(as) sum(as$step^2)
  ## Stage 0: plain quasi Fisher scoring (fast when it works)
  run_plain <- function(theta) {
    for (it in seq_len(maxit)) {
      as <- adj_score(theta)
      if (is.null(as)) return(list(theta = theta, conv = FALSE, it = it, val = NA))
      step <- as$step
      h <- 1
      while (theta[p + 1L] + h * step[p + 1L] <= 0) h <- h / 2
      if (abs(h * step[P]) > 5) h <- 5 / abs(step[P])
      theta <- theta + h * step
      if (max(abs(h * step)) < tol) {
        fin <- adj_score(theta)
        ok <- !is.null(fin) && max(abs(fin$val)) < 1e-5
        return(list(theta = theta, conv = ok, it = it, val = if (ok) fin$val else NA))
      }
    }
    list(theta = theta, conv = FALSE, it = maxit, val = NA)
  }
  ## Stage 1: quasi Fisher scoring theta <- theta + h H^{-1}(g + A) (Kosmidis
  ## and Firth, 2010; for the mean adjustment theta + H^{-1} g - b(theta)),
  ## with step halving on the length of the quasi-Fisher step.  Steps keep
  ## sigma > 0 and are capped in the shape.
  run <- function(theta) {
    as <- adj_score(theta)
    if (is.null(as)) return(list(theta = theta, conv = FALSE, it = 0L, val = NA))
    for (it in seq_len(maxit)) {
      step <- as$step
      h <- 1
      while (theta[p + 1L] + h * step[p + 1L] <= 0) h <- h / 2
      if (abs(h * step[P]) > 2) h <- 2 / abs(step[P])
      repeat {
        new <- theta + h * step
        asn <- adj_score(new)
        if (!is.null(asn) && merit(asn) < merit(as)) break
        h <- h / 2
        if (h < 1e-6) break
      }
      if (is.null(asn) || h < 1e-6) return(list(theta = theta, conv = FALSE, it = it, val = NA))
      theta <- new
      as <- asn
      if (max(abs(as$step)) < tol) {
        ok <- max(abs(as$val)) < 1e-5
        return(list(theta = theta, conv = ok, it = it, val = if (ok) as$val else NA))
      }
    }
    list(theta = theta, conv = FALSE, it = maxit, val = NA)
  }
  ## Stage 2 (fallback): profile the adjusted score in the shape.  For fixed
  ## alpha the remaining components are solved by quasi Fisher scoring;
  ## the alpha component of the adjusted score is then scanned on a grid in
  ## asinh(alpha) and a decreasing sign change is located by uniroot().
  rest <- seq_len(P - 1L)
  solve_rest <- function(alpha, th0) {
    th <- th0
    th[P] <- alpha
    for (it in seq_len(100L)) {
      as <- adj_score(th)
      if (is.null(as)) return(NULL)
      Hs <- tryCatch({
        cu <- sn_cumulants_fast(th[P], th[p + 1L])
        K <- sn_info(th, X, cu)
        if (pen) K + lambda * S else K
      }, error = function(e) NULL)
      if (is.null(Hs)) return(NULL)
      st <- tryCatch(solve(Hs[rest, rest], as$val[rest]), error = function(e) NULL)
      if (is.null(st)) return(NULL)
      h <- 1
      while (th[p + 1L] + h * st[p + 1L] <= 0) h <- h / 2
      th[rest] <- th[rest] + h * st
      if (max(abs(h * st)) < 1e-10) break
    }
    as <- adj_score(th)
    if (is.null(as)) return(NULL)
    list(theta = th, h = as$val[P] / sqrt(max(Hs[P, P], 1e-12)))
  }
  profile_root <- function(th0) {
    tg <- seq(-5.125, 5.125, by = 0.25)
    hv <- rep(NA_real_, length(tg))
    sols <- vector("list", length(tg))
    i0 <- which.min(abs(tg - asinh(th0[P])))
    cur <- th0
    for (i in c(i0:length(tg))) {
      s <- solve_rest(sinh(tg[i]), cur)
      if (!is.null(s)) { hv[i] <- s$h; sols[[i]] <- s$theta; cur <- s$theta }
    }
    cur <- if (!is.null(sols[[i0]])) sols[[i0]] else th0
    for (i in rev(seq_len(i0 - 1L))) {
      s <- solve_rest(sinh(tg[i]), cur)
      if (!is.null(s)) { hv[i] <- s$h; sols[[i]] <- s$theta; cur <- s$theta }
    }
    ## decreasing sign changes, excluding the interval that contains alpha = 0:
    ## a sign change across the singular point is not a root
    nt <- length(tg)
    br <- which(!is.na(hv[-1]) & !is.na(hv[-nt]) & hv[-nt] > 0 & hv[-1] < 0 &
                !(tg[-nt] < 0 & tg[-1] > 0))
    if (!length(br)) {
      ## no regular root; record whether the sign change happens only across
      ## the singular point alpha = 0
      neg <- tg < 0 & !is.na(hv)
      pos <- tg > 0 & !is.na(hv)
      ## the only decreasing sign change is across alpha = 0
      sing <- any(neg) && any(pos) && hv[max(which(neg))] > 0 && hv[min(which(pos))] < 0
      diag_reason <<- if (sing) "sign change only across alpha = 0" else "no sign change found"
      return(NULL)
    }
    j <- br[which.min(abs(tg[br] - tg[i0]))]
    cache <- sols[[j]]
    fz <- function(t) {
      s <- solve_rest(sinh(t), cache)
      if (is.null(s)) return(NA_real_)
      cache <<- s$theta
      s$h
    }
    u <- tryCatch(stats::uniroot(fz, c(tg[j], tg[j + 1]), f.lower = hv[j], f.upper = hv[j + 1],
                                 tol = 1e-10), error = function(e) NULL)
    if (is.null(u)) return(NULL)
    solve_rest(sinh(u$root), cache)$theta
  }
  diag_reason <- NA_character_
  st0 <- if (is.null(start)) .sn_start(X, y) else start
  res <- run_plain(st0)
  if (!res$conv) res <- run(st0)
  if (!res$conv) {
    th2 <- profile_root(st0)
    if (!is.null(th2)) res <- run(th2)
    if (!res$conv && !is.null(th2)) {
      fin <- adj_score(th2)
      if (!is.null(fin) && max(abs(fin$val)) < 1e-5) res <- list(theta = th2, conv = TRUE, it = NA, val = fin$val)
    }
  }
  if (!res$conv) res <- list(theta = rep(NA_real_, P), conv = FALSE, it = maxit, val = NA)
  theta <- res$theta
  names(theta) <- .coef_names(X)
  list(theta = theta, converged = res$conv, iterations = res$it,
       adjusted_score = res$val, type = type,
       failure = if (res$conv) NA_character_ else
         if (is.na(diag_reason)) "iteration did not converge" else diag_reason)
}

#' Parametric bootstrap bias correction
#'
#' Generates `B` samples from the model fitted by maximum likelihood,
#' refits, and returns \eqn{2\hat\theta - \bar\theta^*}. Bootstrap
#' replicates whose shape estimate is infinite are counted and, following
#' the policy chosen by `infinite`, either discarded (the policy of the
#' earlier versions of this work) or make the correction undefined.
#'
#' @param theta maximum likelihood estimate (finite).
#' @param X model matrix.
#' @param B number of bootstrap samples.
#' @param infinite `"discard"` or `"na"`.
#' @return A list with the corrected estimate `theta`, the number of
#'   infinite replicates `n_inf` and the bootstrap mean `mean_star`.
#' @export
sn_boot <- function(theta, X, B = 500L, infinite = c("discard", "na")) {
  infinite <- match.arg(infinite)
  p <- ncol(X)
  if (any(!is.finite(theta))) return(list(theta = theta * NA, n_inf = NA, mean_star = NA))
  mu <- drop(X %*% theta[seq_len(p)])
  st <- matrix(NA_real_, B, p + 2L)
  ninf <- 0L
  for (b in seq_len(B)) {
    yb <- mu + rsnorm(length(mu), 0, theta[p + 1L], theta[p + 2L])
    f <- sn_mle(X, yb, start = theta)
    if (f$finite) st[b, ] <- f$theta else ninf <- ninf + 1L
  }
  ms <- colMeans(st, na.rm = TRUE)
  out <- 2 * theta - ms
  if (infinite == "na" && ninf > 0) out <- out * NA
  list(theta = out, n_inf = ninf, mean_star = ms)
}
